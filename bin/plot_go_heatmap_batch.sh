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
    echo "[INFO] Cleaning up temporary directory: ${TMPDIR}"
    rm -rf "${TMPDIR}"
}
trap cleanup EXIT

# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------
echo "[INFO] Starting GO heatmap batch generation"
echo "[INFO] GO ID list       : ${GO_ID_LIST}"
echo "[INFO] Expression matrix : ${EXPRESSION_MATRIX}"
echo "[INFO] Sample info       : ${SAMPLE_INFO}"
echo "[INFO] GMT file          : ${GMT}"
echo "[INFO] R script          : ${R_SCRIPT}"
echo "[INFO] Output prefix     : ${OUT_PREFIX}"

EXP_HEADER="${TMPDIR}/exp_headers"
head -1 "${EXPRESSION_MATRIX}" | tr '\t' '\n' > "${EXP_HEADER}"

# Build sample group file: extract group labels for columns present in the
# expression matrix header.
SAMPLE_HEADER="${TMPDIR}/sample_groups.txt"
tail -n +2 "${SAMPLE_INFO}" \
    | grep -w -f "${EXP_HEADER}" \
    | cut -d',' -f2 \
    > "${SAMPLE_HEADER}"
echo "[INFO] Found $(wc -l < "${SAMPLE_HEADER}") sample group entries"

# For each GO ID, extract gene set from GMT and subset expression matrix.
FILELIST="${TMPDIR}/exp_file_list.txt"
: > "${FILELIST}"

GO_COUNT=0
SKIP_COUNT=0
while IFS= read -r GO_ID || [[ -n "${GO_ID}" ]]; do
    # Skip blank lines
    [[ -z "${GO_ID}" ]] && continue

    GO_NAME=$(grep -w "${GO_ID}" "${GMT}" | cut -f2 | sort -u | head -1 | sed 's/ /_/g')

    if [[ -z "${GO_NAME}" ]]; then
        echo "[WARN] GO ID '${GO_ID}' not found in GMT file -- skipping"
        ((SKIP_COUNT++)) || true
        continue
    fi

    # Extract gene list for this GO term
    GENE_LIST="${TMPDIR}/genes_${GO_NAME}.txt"
    grep -w "${GO_ID}" "${GMT}" | cut -f3- | tr '\t' '\n' | sort -k1 > "${GENE_LIST}"

    GENE_COUNT=$(wc -l < "${GENE_LIST}")
    if [[ "${GENE_COUNT}" -eq 0 ]]; then
        echo "[WARN] No genes found for GO ID '${GO_ID}' (${GO_NAME}) -- skipping"
        ((SKIP_COUNT++)) || true
        continue
    fi

    # Subset expression matrix: keep header + matching gene rows
    SUBSET_EXP="${TMPDIR}/tempExp.${GO_NAME}"
    head -1 "${EXPRESSION_MATRIX}" > "${SUBSET_EXP}"
    grep -w -f "${GENE_LIST}" "${EXPRESSION_MATRIX}" >> "${SUBSET_EXP}" || true

    MATCHED=$(($(wc -l < "${SUBSET_EXP}") - 1))
    if [[ "${MATCHED}" -eq 0 ]]; then
        echo "[WARN] No genes from GO:${GO_ID} (${GO_NAME}) matched in expression matrix -- skipping"
        rm -f "${SUBSET_EXP}"
        ((SKIP_COUNT++)) || true
        continue
    fi

    echo "${SUBSET_EXP}" >> "${FILELIST}"
    echo "[INFO] GO:${GO_ID} (${GO_NAME}) -- ${MATCHED} genes matched"
    ((GO_COUNT++)) || true

done < "${GO_ID_LIST}"

echo "[INFO] Processed ${GO_COUNT} GO terms (${SKIP_COUNT} skipped)"

if [[ "${GO_COUNT}" -eq 0 ]]; then
    echo "ERROR: No GO terms produced valid gene subsets. Check inputs and try again." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Call R heatmap plotting script
# ---------------------------------------------------------------------------
echo "[INFO] Calling R heatmap script: ${R_SCRIPT}"
RSCRIPT_BIN="$(mys_resolve_bin MYS_RSCRIPT_BIN Rscript)" || { echo "ERROR: Rscript not found. Set MYS_RSCRIPT_BIN in config/env.local.sh." >&2; exit 1; }
"${RSCRIPT_BIN}" "${R_SCRIPT}" "${FILELIST}" "${SAMPLE_HEADER}" "${OUT_PREFIX}"
echo "[INFO] R script completed"

echo "[INFO] Done. Output prefix: ${OUT_PREFIX}"
