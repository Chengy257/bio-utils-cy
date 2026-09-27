#!/usr/bin/env bash
#########################################################################
# File Name: bin/plot_go_heatmap_batch.sh
# Author: ChengYu
# Description: Batch generate heatmaps for GO term gene sets from
#              expression data. For each GO ID in the input list, genes
#              are extracted from a GMT file and used to subset an
#              expression matrix. The resulting per-GO matrices are then
#              passed to the R helper script (plot_heatmap_multi.R) for
#              combined and individual heatmap generation.
# Created Time: 2026
# Changelog:
#   v1.1.0  2026-09-28
#   - FIX: gene/GO matching used regex patterns (grep); gene IDs
#     containing "." (e.g. AT1G01010.1) over-matched unrelated rows
#     (AT1G01010X1). Matching is now exact first-field (no regex).
#   - FIX: subset matching also accepted word-boundary neighbours
#     (gene "bHLH38" matched row "bHLH38-like") and could match the
#     header line. Only exact first-column equality matches now, and
#     the header is never re-matched into the data.
#   - FIX: group labels were emitted in sample-info row order, not in
#     expression-matrix column order; plot_heatmap_multi.R aligns labels
#     to columns positionally, so a different order silently mislabeled
#     every heatmap. Labels are now emitted in matrix column order via a
#     sample->group lookup, and an unknown sample is a hard error.
#   - FIX: GO term descriptions containing "." broke the cluster-name
#     contract (filename field 2 after "." split), and two GO IDs whose
#     descriptions sanitize identically overwrote each other's subsets
#     and PDFs. Filenames now embed the sanitized GO ID plus a fully
#     sanitized description (safe chars only), making names unique.
#   - FIX: CRLF input files: CR is stripped from the GO ID list, sample
#     info and expression header before matching.
#   - CHANGE: all INFO/WARN logs go to stderr; stdout is reserved for
#     pipeline data.
#   - GENE_LIST is deduplicated (sort -u) and empty fields dropped.
#########################################################################

set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf

# ---------------------------------------------------------------------------
# Default values
# ---------------------------------------------------------------------------
GMT=""
R_SCRIPT="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/plot_heatmap_multi.R"
OUT_PREFIX="GO_Heatmap"

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
    cat <<'EOF'
Usage: plot_go_heatmap_batch.sh [OPTIONS] GO_ID_LIST EXPRESSION_MATRIX SAMPLE_INFO

Batch generate heatmaps for GO term gene sets from expression data.

For each GO ID listed in GO_ID_LIST, genes are extracted from the
provided GMT file. The expression matrix is subset to those genes, and
all per-GO subsets are passed to the R heatmap plotting script.

Positional arguments:
  GO_ID_LIST         File with one GO term ID per line (e.g. GO:0008150)
  EXPRESSION_MATRIX  Tab-delimited expression matrix (genes x samples)
  SAMPLE_INFO        CSV file mapping sample column names to groups
                     (header row + rows of: SampleID,GroupName)

Required options:
  -g GMT_FILE        Path to a GMT file containing GO term gene sets

Optional options:
  -r R_SCRIPT        Path to the R heatmap plotting script
                     (default: plot_heatmap_multi.R alongside this script)
  -o OUT_PREFIX      Prefix for output files (default: GO_Heatmap)
  -h                 Show this help message and exit

Notes:
  Every sample column of the expression matrix must have a group in
  SAMPLE_INFO (entries are matched in matrix column order). Gene and
  GO ID matching is exact first-column equality — IDs are not treated
  as regular expressions.

Examples:
  # Basic usage
  plot_go_heatmap_batch.sh -g go_terms.gmt go_ids.txt expr_matrix.tsv sample_info.csv

  # Custom R script and output prefix
  plot_go_heatmap_batch.sh -g go_terms.gmt -r /path/to/plot_heatmap_multi.R \
      -o my_results go_ids.txt expr_matrix.tsv sample_info.csv
EOF
}

# ---------------------------------------------------------------------------
# Parse options
# ---------------------------------------------------------------------------
while getopts ":g:r:o:h" opt; do
    case "${opt}" in
        g) GMT="${OPTARG}" ;;
        r) R_SCRIPT="${OPTARG}" ;;
        o) OUT_PREFIX="${OPTARG}" ;;
        h) usage; exit 0 ;;
        :) echo "ERROR: Option -${OPTARG} requires an argument." >&2; usage; exit 1 ;;
        \?) echo "ERROR: Unknown option -${OPTARG}." >&2; usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

# ---------------------------------------------------------------------------
# Validate positional arguments
# ---------------------------------------------------------------------------
if [[ $# -ne 3 ]]; then
    echo "ERROR: Expected exactly 3 positional arguments (GO_ID_LIST EXPRESSION_MATRIX SAMPLE_INFO), got $#." >&2
    usage
    exit 1
fi

GO_ID_LIST="$1"
EXPRESSION_MATRIX="$2"
SAMPLE_INFO="$3"

# ---------------------------------------------------------------------------
# Validate required flag
# ---------------------------------------------------------------------------
if [[ -z "${GMT}" ]]; then
    echo "ERROR: -g GMT_FILE is required." >&2
    usage
    exit 1
fi

# ---------------------------------------------------------------------------
# Validate all input files exist
# ---------------------------------------------------------------------------
for f in "${GO_ID_LIST}" "${EXPRESSION_MATRIX}" "${SAMPLE_INFO}" "${GMT}" "${R_SCRIPT}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Required file not found: ${f}" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# Create and manage temporary directory
# ---------------------------------------------------------------------------
TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/plot_go_heatmap_batch.XXXXXX")
cleanup() {
    echo "[INFO] Cleaning up temporary directory: ${TMPDIR}" >&2
    rm -rf "${TMPDIR}"
}
trap cleanup EXIT

# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------
log_info() { echo "[INFO] $*" >&2; }
log_warn() { echo "[WARN] $*" >&2; }

log_info "Starting GO heatmap batch generation"
log_info "GO ID list       : ${GO_ID_LIST}"
log_info "Expression matrix : ${EXPRESSION_MATRIX}"
log_info "Sample info       : ${SAMPLE_INFO}"
log_info "GMT file          : ${GMT}"
log_info "R script          : ${R_SCRIPT}"
log_info "Output prefix     : ${OUT_PREFIX}"

# Normalize the GO ID list once (strip CR from CRLF files, drop blanks,
# drop duplicate IDs so a term is never processed twice).
GO_IDS="${TMPDIR}/go_ids.txt"
tr -d '\r' < "${GO_ID_LIST}" | sed '/^[[:space:]]*$/d' | awk '!seen[$0]++' > "${GO_IDS}"

# ---------------------------------------------------------------------------
# Build the sample-group file in EXPRESSION-MATRIX COLUMN ORDER.
# plot_heatmap_multi.R aligns group labels to matrix columns positionally,
# so the label order here must follow the matrix header, not the row
# order of SAMPLE_INFO.
# ---------------------------------------------------------------------------
EXP_SAMPLES="${TMPDIR}/exp_samples.txt"
# matrix header -> one sample per line, skipping the gene-ID column,
# stripping CR and empty (trailing-tab) fields.
head -1 "${EXPRESSION_MATRIX}" | tr '\t' '\n' | tail -n +2 \
    | sed 's/\r$//' | sed '/^[[:space:]]*$/d' > "${EXP_SAMPLES}"

SAMPLE_HEADER="${TMPDIR}/sample_groups.txt"
awk -F',' '
    FNR == NR {
        if (FNR > 1 && NF >= 2) {
            sub(/\r$/, "", $1); gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
            sub(/\r$/, "", $2); gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
            map[$1] = $2
        }
        next
    }
    {
        if (!($0 in map)) {
            printf "ERROR: sample %s has no group in sample info file %s\n", $0, ARGV[1] > "/dev/stderr"
            exit 1
        }
        print map[$0]
    }
' "${SAMPLE_INFO}" "${EXP_SAMPLES}" > "${SAMPLE_HEADER}"

if [[ ! -s "${SAMPLE_HEADER}" ]]; then
    echo "ERROR: No sample columns found in the expression matrix header, or no groups resolved." >&2
    exit 1
fi
log_info "Found $(wc -l < "${SAMPLE_HEADER}") sample group entries (matrix column order)"

# For each GO ID, extract gene set from GMT and subset expression matrix.
FILELIST="${TMPDIR}/exp_file_list.txt"
: > "${FILELIST}"

GO_COUNT=0
SKIP_COUNT=0
while IFS= read -r GO_ID; do
    GO_TAG=$(printf '%s' "${GO_ID}" | sed 's/[^A-Za-z0-9_-]/_/g')
    GMT_HITS="${TMPDIR}/gmt_hits_${GO_TAG}.txt"
    # Exact first-column match on the GMT (no regex, CR-safe).
    awk -F'\t' -v id="${GO_ID}" '{ sub(/\r$/, "", $1); if ($1 == id) print }' \
        "${GMT}" > "${GMT_HITS}"

    GO_NAME=$(cut -f2 "${GMT_HITS}" | sort -u | head -1)
    if [[ -z "${GO_NAME}" ]]; then
        log_warn "GO ID '${GO_ID}' not found in GMT file -- skipping"
        ((SKIP_COUNT++)) || true
        continue
    fi
    # Description must be filename-safe AND free of "." (the R helper
    # derives the cluster name from filename field 2 after a "." split).
    GO_NAME=$(printf '%s' "${GO_NAME}" | sed 's/[^A-Za-z0-9_-]/_/g')

    # Extract gene list for this GO term: columns 3+, literal fields,
    # deduplicated, empty fields (e.g. trailing tabs) dropped.
    GENE_LIST="${TMPDIR}/genes_${GO_TAG}_${GO_NAME}.txt"
    cut -f3- "${GMT_HITS}" | tr '\t' '\n' | sed 's/\r$//' \
        | sed '/^[[:space:]]*$/d' | sort -u > "${GENE_LIST}"

    GENE_COUNT=$(wc -l < "${GENE_LIST}")
    if [[ "${GENE_COUNT}" -eq 0 ]]; then
        log_warn "No genes found for GO ID '${GO_ID}' (${GO_NAME}) -- skipping"
        ((SKIP_COUNT++)) || true
        continue
    fi

    # Subset expression matrix: keep header + rows whose first field
    # equals one of the gene list entries (exact literal match — regex
    # and word-boundary heuristics would pull in wrong genes). The
    # header line itself is never re-matched into the data.
    SUBSET_EXP="${TMPDIR}/tempExp.${GO_TAG}_${GO_NAME}"
    head -1 "${EXPRESSION_MATRIX}" > "${SUBSET_EXP}"
    tail -n +2 "${EXPRESSION_MATRIX}" | awk -F'\t' '
        FNR == NR { sub(/\r$/, ""); if ($0 != "") p[$0] = 1; next }
        { sub(/\r$/, "", $1); if ($1 in p) print }
    ' "${GENE_LIST}" - >> "${SUBSET_EXP}"

    MATCHED=$(($(wc -l < "${SUBSET_EXP}") - 1))
    if [[ "${MATCHED}" -eq 0 ]]; then
        log_warn "No genes from ${GO_ID} (${GO_NAME}) matched in expression matrix -- skipping"
        rm -f "${SUBSET_EXP}"
        ((SKIP_COUNT++)) || true
        continue
    fi

    echo "${SUBSET_EXP}" >> "${FILELIST}"
    log_info "${GO_ID} (${GO_NAME}) -- ${MATCHED}/${GENE_COUNT} genes matched"
    ((GO_COUNT++)) || true

done < "${GO_IDS}"

log_info "Processed ${GO_COUNT} GO terms (${SKIP_COUNT} skipped)"

if [[ "${GO_COUNT}" -eq 0 ]]; then
    echo "ERROR: No GO terms produced valid gene subsets. Check inputs and try again." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Call R heatmap plotting script
# ---------------------------------------------------------------------------
log_info "Calling R heatmap script: ${R_SCRIPT}"
RSCRIPT_BIN="$(buc_resolve_bin BUC_RSCRIPT_BIN Rscript)" || { echo "ERROR: Rscript not found. Set BUC_RSCRIPT_BIN in config/env.local.sh." >&2; exit 1; }
# plot_heatmap_multi.R takes flag options only (-f/-g/-o); a positional
# call makes getopt() reject the file-list path.
"${RSCRIPT_BIN}" "${R_SCRIPT}" \
    -f "${FILELIST}" -g "${SAMPLE_HEADER}" -o "${OUT_PREFIX}"
log_info "R script completed"

log_info "Done. Output prefix: ${OUT_PREFIX}"
