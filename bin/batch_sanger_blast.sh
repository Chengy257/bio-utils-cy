#!/usr/bin/env bash
#########################################################################
# File Name: bin/batch_sanger_blast.sh
# Author: ChengYu
# Description: Batch extract Sanger sequencing results from .zip files,
#              assemble sequences into FASTA, create BLAST database, and
#              run BLASTn in two output formats (pairwise + tabular).
# Created Time: 2026
# Changelog:
#   v1.1.0  2026-09-28
#   - FIX: query names came from the first "."-separated token of the
#     basename, so identically named files in different subdirectories
#     produced duplicate FASTA headers and ambiguous BLAST hits. Names
#     are now derived from the relative path (collision-checked), with
#     the .seq suffix stripped and unsafe characters replaced.
#   - FIX: every *seq file was cat'ed verbatim into the query FASTA; a
#     file holding FASTA format polluted the query with '>' header
#     lines, and non-sequence files (e.g. "refseq") silently entered
#     the BLAST. Candidate files are now validated as nucleotide
#     sequence (IUPAC), embedded '>' header lines are dropped, and
#     invalid files are skipped with a warning.
#   - FIX: THREADS was hardcoded; now defaults to $BUC_THREADS and is
#     validated as a positive integer.
#   - CHANGE: dropped the fragile perl/util-linux rename probing; a
#     single portable shell loop sanitizes filenames.
#   - CHANGE: unzip is verified up front (a missing unzip previously
#     surfaced as a misleading "No *seq files found").
#   - CHANGE: INFO logs go to stderr; stdout is reserved for data.
#   - *seq candidates are processed in sorted order for reproducible
#     FASTA/BLAST output.
#########################################################################

set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf

## ---- Defaults ----
DATABASE=""
WORK_DIR=""
THREADS="${BUC_THREADS:-4}"
OUT_PREFIX="sanger_BLAST"

log_info() { echo "[INFO] $*" >&2; }
log_warn() { echo "[WARN] $*" >&2; }

## ---- Help function ----
usage() {
    cat <<EOF
Usage: bash batch_sanger_blast.sh -d <database.fa> -w <work_dir> [options]

Batch extract Sanger sequencing .zip archives, assemble into FASTA,
and run BLASTn against a reference database.

Required:
  -d <path>    BLAST database FASTA file
  -w <path>    Working directory containing .zip files

Options:
  -t <int>     Number of BLAST threads (default: \$BUC_THREADS or 4)
  -o <prefix>  Output file prefix (default: sanger_BLAST)
  -h           Show this help message

Notes:
  Files whose names end in "seq" are treated as Sanger sequence
  exports. Each candidate is validated as nucleotide sequence (IUPAC
  alphabet); embedded FASTA ">" header lines are stripped and invalid
  files are skipped with a warning. Query names are derived from each
  file's relative path, so identically named files in different
  subdirectories stay distinct.

Outputs (written into the working directory):
  {prefix}_result.txt    BLAST pairwise output (format 1)
  {prefix}_result.tsv    BLAST tabular output (format 6)
  allseq.fa              assembled query FASTA
  sanger.blastdb.*       BLAST database files

Example:
  bash batch_sanger_blast.sh -d ref_sequence.fa -w sanger_zips/ -t 8 -o blast_results

EOF
    exit 0
}

## ---- Parse options ----
while getopts ":d:w:t:o:h" opt; do
    case "${opt}" in
        d) DATABASE="${OPTARG}" ;;
        w) WORK_DIR="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        o) OUT_PREFIX="${OPTARG}" ;;
        h) usage ;;
        :)
            echo "Error: Option -${OPTARG} requires an argument." >&2
            exit 1
            ;;
        \?)
            echo "Error: Invalid option -${OPTARG}." >&2
            exit 1
            ;;
    esac
done

## ---- Validate required arguments ----
if [[ -z "${DATABASE}" ]]; then
    echo "Error: -d (BLAST database FASTA) is required. Use -h for help." >&2
    exit 1
fi

if [[ -z "${WORK_DIR}" ]]; then
    echo "Error: -w (work directory) is required. Use -h for help." >&2
    exit 1
fi

if [[ ! "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: -t must be a positive integer, got '${THREADS}'." >&2
    exit 1
fi

if [[ ! -f "${DATABASE}" ]]; then
    echo "Error: Database file not found: ${DATABASE}" >&2
    exit 1
fi
# The script later cd's into the work directory; a relative -d path must
# keep pointing at the file given on the command line (the documented
# example "-d ref.fa -w sanger_zips/" is exactly this case).
DATABASE="$(readlink -f "${DATABASE}")"

if [[ ! -d "${WORK_DIR}" ]]; then
    echo "Error: Work directory not found: ${WORK_DIR}" >&2
    exit 1
fi

## ---- Resolve BLAST+ tools (BUC_*_BIN from config/env.sh > PATH) ----
MAKEBLASTDB_BIN="$(buc_resolve_bin BUC_MAKEBLASTDB_BIN makeblastdb)" || { echo "Error: makeblastdb not found. Install BLAST+ or set BUC_MAKEBLASTDB_BIN." >&2; exit 1; }
BLASTN_BIN="$(buc_resolve_bin BUC_BLASTN_BIN blastn)" || { echo "Error: blastn not found. Install BLAST+ or set BUC_BLASTN_BIN." >&2; exit 1; }
UNZIP_BIN="$(command -v unzip || true)"
if [[ -z "${UNZIP_BIN}" ]]; then
    echo "Error: unzip not found. Install unzip to extract Sanger .zip archives." >&2
    exit 1
fi

## ---- Cleanup trap ----
ORIG_DIR="$(pwd)"
cleanup() {
    cd "${ORIG_DIR}"
    echo "[INFO] Returned to original directory." >&2
}
trap cleanup EXIT

## ---- Enter work directory ----
cd "${WORK_DIR}"
log_info "Working directory: $(pwd)"

## ---- Step 1: Unzip archives ----
log_info "Step 1: Extracting .zip files"
zip_count=0
if ls ./*.zip &>/dev/null; then
    for zf in ./*.zip; do
        log_info "  Unzipping: ${zf}"
        unzip -O GBK -o "${zf}" >/dev/null 2>&1 || unzip -o "${zf}" >/dev/null 2>&1 || {
            log_warn "  Failed to unzip ${zf}, skipping."
            continue
        }
        zip_count=$((zip_count + 1))
    done
else
    log_info "  No .zip files found in work directory."
fi
log_info "  Extracted ${zip_count} archive(s)."

## ---- Step 2: Sanitize filenames (strip non-printable bytes) ----
log_info "Step 2: Sanitizing filenames"
for f in ./*; do
    [[ -f "${f}" ]] || continue
    clean_name=$(printf '%s' "${f}" | LC_ALL=C sed 's/[^[:print:]]/_/g')
    if [[ "${f}" != "${clean_name}" && ! -e "${clean_name}" ]]; then
        mv "${f}" "${clean_name}"
    fi
done
log_info "  Filename sanitization complete."

## ---- Step 3: Collect *seq files into allseq.fa ----
log_info "Step 3: Assembling sequence files into allseq.fa"
ALLSEQ="allseq.fa"
rm -f "${ALLSEQ}"

seq_count=0
skip_count=0
declare -A seen_names=()
# Sorted for reproducible FASTA record order.
while IFS= read -r -d '' seqfile; do
    # Query name: relative path with the .seq suffix stripped and unsafe
    # characters replaced; kept unique across identical basenames.
    name="${seqfile#./}"
    name="${name%.[sS][eE][qQ]}"
    name="${name//[^A-Za-z0-9._-]/_}"
    if [[ -n "${seen_names[${name}]+x}" ]]; then
        n=2
        while [[ -n "${seen_names[${name}_${n}]+x}" ]]; do
            n=$((n + 1))
        done
        name="${name}_${n}"
    fi
    seen_names["${name}"]=1

    # Validate content: drop embedded FASTA headers and whitespace, then
    # require the remainder to be a nucleotide sequence.
    seq=$(LC_ALL=C tr -d '\r' < "${seqfile}" \
          | sed '/^>/d' | tr -d '[:space:]')
    if [[ -z "${seq}" ]] || ! [[ "${seq}" =~ ^[ACGTURYSWKMBDHVNacgturyswkmbdhvn.-]+$ ]]; then
        log_warn "  Skipping ${seqfile}: not recognisable as nucleotide sequence data."
        skip_count=$((skip_count + 1))
        continue
    fi

    printf '>%s\n%s\n' "${name}" "${seq}" >> "${ALLSEQ}"
    seq_count=$((seq_count + 1))
done < <(find . -type f -name "*seq" -print0 2>/dev/null | sort -z)

if [[ "${seq_count}" -eq 0 ]]; then
    echo "Error: No valid *seq files found after extraction. Nothing to BLAST." >&2
    exit 1
fi
log_info "  Assembled ${seq_count} sequence(s) into ${ALLSEQ} (${skip_count} skipped)."

## ---- Step 4: Create BLAST database ----
log_info "Step 4: Creating BLAST database"
"${MAKEBLASTDB_BIN}" -in "${DATABASE}" -dbtype nucl -parse_seqids \
    -out sanger.blastdb -logfile blastdb_log.txt >/dev/null
log_info "  BLAST database created from: ${DATABASE}"

## ---- Step 5: Run BLASTn (format 1 - pairwise) ----
log_info "Step 5: Running BLASTn (pairwise format)"
RESULT_TXT="${OUT_PREFIX}_result.txt"
"${BLASTN_BIN}" -query "${ALLSEQ}" -db sanger.blastdb \
    -out "${RESULT_TXT}" -outfmt 1 \
    -num_threads "${THREADS}"
[[ -s "${RESULT_TXT}" ]] || log_warn "  Pairwise BLAST produced no hits."
log_info "  Pairwise results: ${RESULT_TXT}"

## ---- Step 6: Run BLASTn (format 6 - tabular) ----
log_info "Step 6: Running BLASTn (tabular format)"
RESULT_TSV="${OUT_PREFIX}_result.tsv"
"${BLASTN_BIN}" -query "${ALLSEQ}" -db sanger.blastdb \
    -out "${RESULT_TSV}" -outfmt 6 \
    -num_threads "${THREADS}"
log_info "  Tabular results: ${RESULT_TSV}"

## ---- Summary ----
{
    echo ""
    echo "=== Summary ==="
    echo "  Sequences processed: ${seq_count}"
    echo "  Sequences skipped:   ${skip_count}"
    echo "  BLAST database:      ${DATABASE}"
    echo "  Threads used:        ${THREADS}"
    echo "  Pairwise output:     ${RESULT_TXT}"
    echo "  Tabular output:      ${RESULT_TSV}"
    echo "  BLAST completed successfully."
} >&2
