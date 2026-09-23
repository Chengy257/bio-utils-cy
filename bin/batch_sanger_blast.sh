#!/usr/bin/env bash
#########################################################################
# File Name: bin/batch_sanger_blast.sh
# Author: ChengYu
# Description: Batch extract Sanger sequencing results from .zip files,
#              assemble sequences into FASTA, create BLAST database, and
#              run BLASTn in two output formats (pairwise + tabular).
# Created Time: 2026
#########################################################################

set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf

## ---- Defaults ----
DATABASE=""
WORK_DIR=""
THREADS=4
OUT_PREFIX="sanger_BLAST"

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
  -t <int>     Number of BLAST threads (default: 4)
  -o <prefix>  Output file prefix (default: sanger_BLAST)
  -h           Show this help message

Outputs:
  {prefix}_result.txt    BLAST pairwise output (format 1)
  {prefix}_result.tsv    BLAST tabular output (format 6)
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

if [[ ! -f "${DATABASE}" ]]; then
    echo "Error: Database file not found: ${DATABASE}" >&2
    exit 1
fi

if [[ ! -d "${WORK_DIR}" ]]; then
    echo "Error: Work directory not found: ${WORK_DIR}" >&2
    exit 1
fi

## ---- Check BLAST+ tools are available ----
for tool in makeblastdb blastn; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "Error: ${tool} not found in PATH. Please install BLAST+." >&2
        exit 1
    fi
done

## ---- Cleanup trap ----
ORIG_DIR="$(pwd)"
cleanup() {
    cd "${ORIG_DIR}"
    echo "Cleanup: returned to original directory."
}
trap cleanup EXIT

## ---- Enter work directory ----
cd "${WORK_DIR}"
echo "Working directory: $(pwd)"

## ---- Step 1: Unzip archives ----
echo "=== Step 1: Extracting .zip files ==="
zip_count=0
if ls ./*.zip &>/dev/null; then
    for zf in ./*.zip; do
        echo "  Unzipping: ${zf}"
        unzip -O GBK -o "${zf}" 2>/dev/null || unzip -o "${zf}" 2>/dev/null || {
            echo "  Warning: Failed to unzip ${zf}, skipping."
            continue
        }
        zip_count=$((zip_count + 1))
    done
else
    echo "  No .zip files found in work directory."
fi
echo "  Extracted ${zip_count} archive(s)."

## ---- Step 2: Rename files to remove problematic characters ----
echo "=== Step 2: Sanitizing filenames ==="
rename_count=0
# Use Perl rename if available, otherwise use a shell loop
if command -v rename &>/dev/null; then
    # Try perl-based rename first
    if rename --version 2>/dev/null | grep -qi perl; then
        rename 's/[^[:ascii:]]/_/g' ./* 2>/dev/null || true
        rename_count=1
    else
        # util-linux rename - limited; use shell fallback
        rename_count=0
    fi
fi

# Shell-based fallback: rename files with non-ASCII characters
if [[ "${rename_count}" -eq 0 ]]; then
    for f in ./*; do
        [[ -f "${f}" ]] || continue
        clean_name=$(echo "${f}" | LC_ALL=C sed 's/[^[:print:]]/_/g')
        if [[ "${f}" != "${clean_name}" ]]; then
            mv "${f}" "${clean_name}" 2>/dev/null || true
        fi
    done
fi
echo "  Filename sanitization complete."

## ---- Step 3: Collect *seq files into allseq.fa ----
echo "=== Step 3: Assembling sequence files into allseq.fa ==="
ALLSEQ="allseq.fa"
if [[ -f "${ALLSEQ}" ]]; then
    rm -f "${ALLSEQ}"
fi

seq_count=0
# Find all files matching *seq pattern (including in subdirs from unzip)
while IFS= read -r -d '' seqfile; do
    name=$(basename "${seqfile}" | cut -d"." -f1)
    echo ">${name}" >> "${ALLSEQ}"
    cat "${seqfile}" >> "${ALLSEQ}"
    echo "" >> "${ALLSEQ}"
    seq_count=$((seq_count + 1))
done < <(find . -type f -name "*seq" -print0 2>/dev/null)

if [[ "${seq_count}" -eq 0 ]]; then
    echo "Error: No *seq files found after extraction. Nothing to BLAST." >&2
    exit 1
fi
echo "  Assembled ${seq_count} sequence(s) into ${ALLSEQ}."

## ---- Step 4: Create BLAST database ----
echo "=== Step 4: Creating BLAST database ==="
makeblastdb -in "${DATABASE}" -dbtype nucl -parse_seqids \
    -out sanger.blastdb -logfile blastdb_log.txt
echo "  BLAST database created from: ${DATABASE}"

## ---- Step 5: Run BLASTn (format 1 - pairwise) ----
echo "=== Step 5: Running BLASTn (pairwise format) ==="
RESULT_TXT="${OUT_PREFIX}_result.txt"
blastn -query "${ALLSEQ}" -db sanger.blastdb \
    -out "${RESULT_TXT}" -outfmt 1 \
    -num_threads "${THREADS}"
echo "  Pairwise results: ${RESULT_TXT}"

## ---- Step 6: Run BLASTn (format 6 - tabular) ----
echo "=== Step 6: Running BLASTn (tabular format) ==="
RESULT_TSV="${OUT_PREFIX}_result.tsv"
blastn -query "${ALLSEQ}" -db sanger.blastdb \
    -out "${RESULT_TSV}" -outfmt 6 \
    -num_threads "${THREADS}"
echo "  Tabular results: ${RESULT_TSV}"

## ---- Summary ----
echo ""
echo "=== Summary ==="
echo "  Sequences processed: ${seq_count}"
echo "  BLAST database:      ${DATABASE}"
echo "  Threads used:        ${THREADS}"
echo "  Pairwise output:     ${RESULT_TXT}"
echo "  Tabular output:      ${RESULT_TSV}"
echo "  BLAST completed successfully."
