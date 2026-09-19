#!/usr/bin/env bash
#
# File Name: fetch_bioproject_metadata.sh
# Author: ChengYu
# Description: Batch download ENA filereport metadata for BioProject accessions.
#              Auto-detects axel or wget for parallel downloading.
# Created Time: 2026
#

set -euo pipefail

readonly SCRIPT_NAME="$(basename "$0")"
readonly VERSION="1.0.0"

# Defaults
DEFAULT_FIELDS="run_accession,experiment_accession,sample_accession,study_accession,experiment_title,instrument_model,library_layout,library_strategy,library_source,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes"
DEFAULT_THREADS=4
DEFAULT_TOOL="auto"

ENA_BASE_URL="https://www.ebi.ac.uk/ena/portal/api/filereport"

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [OPTIONS] -i <input_file> -o <outdir>

Fetch ENA filereport metadata for BioProject accessions listed in a file.

Required:
  -i <file>       Input file with one BioProject accession per line
  -o <dir>        Output directory for downloaded TSV files

Options:
  -f <fields>     Comma-separated ENA fields to retrieve
                  (default: ${DEFAULT_FIELDS})
  -t <threads>    Number of parallel threads for axel (default: ${DEFAULT_THREADS})
  -T <tool>       Download tool: auto|axel|wget (default: ${DEFAULT_TOOL})
  --version       Print version and exit
  -h              Show this help message

Examples:
  ${SCRIPT_NAME} -i bioprojects.txt -o ./ena_data
  ${SCRIPT_NAME} -i projects.txt -o out -f "run_accession,sample_accession,fastq_ftp" -t 8
  ${SCRIPT_NAME} -i projects.txt -o out -T wget
EOF
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

INPUT_FILE=""
OUT_DIR=""
FIELDS="${DEFAULT_FIELDS}"
THREADS="${DEFAULT_THREADS}"
TOOL="${DEFAULT_TOOL}"

while getopts ":i:o:f:t:T:h-:" opt; do
    case "${opt}" in
        i) INPUT_FILE="${OPTARG}" ;;
        o) OUT_DIR="${OPTARG}" ;;
        f) FIELDS="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        T) TOOL="${OPTARG}" ;;
        h) usage; exit 0 ;;
        -)
            case "${OPTARG}" in
                version) echo "${SCRIPT_NAME} ${VERSION}"; exit 0 ;;
                *) echo "Error: unknown option --${OPTARG}" >&2; usage; exit 1 ;;
            esac
            ;;
        :) echo "Error: option -${OPTARG} requires an argument." >&2; usage; exit 1 ;;
        \?) echo "Error: unknown option -${OPTARG}." >&2; usage; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${INPUT_FILE}" ]]; then
    echo "Error: -i <input_file> is required." >&2
    usage
    exit 1
fi

if [[ -z "${OUT_DIR}" ]]; then
    echo "Error: -o <outdir> is required." >&2
    usage
    exit 1
fi

if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "Error: input file not found: ${INPUT_FILE}" >&2
    exit 1
fi

# Create output directory
mkdir -p "${OUT_DIR}"

# ---------------------------------------------------------------------------
# Detect download tool
# ---------------------------------------------------------------------------

detect_tool() {
    if [[ "${TOOL}" == "auto" ]]; then
        if command -v axel &>/dev/null; then
            echo "axel"
        elif command -v wget &>/dev/null; then
            echo "wget"
        else
            echo "Error: neither axel nor wget found in PATH." >&2
            exit 1
        fi
    else
        if ! command -v "${TOOL}" &>/dev/null; then
            echo "Error: tool '${TOOL}' not found in PATH." >&2
            exit 1
        fi
        echo "${TOOL}"
    fi
}

DL_TOOL="$(detect_tool)"
echo "[INFO] Using download tool: ${DL_TOOL}"

# ---------------------------------------------------------------------------
# Build URL for a given BioProject accession
# ---------------------------------------------------------------------------

build_url() {
    local accession="$1"
    local fields="$2"
    local base_url="${ENA_BASE_URL}"
    echo "${base_url}?accession=${accession}&result=read_run&fields=${fields}&format=tsv&download=txt"
}

# ---------------------------------------------------------------------------
# Download function
# ---------------------------------------------------------------------------

download_with_axel() {
    local url="$1"
    local output="$2"
    local n_threads="$3"
    axel -n "${n_threads}" -o "${output}" "${url}" 2>&1 || {
        echo "[WARN] axel failed for ${url}" >&2
        return 1
    }
}

download_with_wget() {
    local url="$1"
    local output="$2"
    wget -q -O "${output}" "${url}" 2>&1 || {
        echo "[WARN] wget failed for ${url}" >&2
        return 1
    }
}

download_one() {
    local url="$1"
    local output="$2"
    case "${DL_TOOL}" in
        axel) download_with_axel "${url}" "${output}" "${THREADS}" ;;
        wget) download_with_wget "${url}" "${output}" ;;
        *) echo "[ERROR] Unknown tool: ${DL_TOOL}" >&2; return 1 ;;
    esac
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

SUCCESS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

while IFS= read -r line || [[ -n "${line}" ]]; do
    # Strip whitespace and skip empty/comment lines
    acc="$(echo "${line}" | xargs)"
    [[ -z "${acc}" ]] && continue
    [[ "${acc}" == \#* ]] && continue

    outfile="${OUT_DIR}/${acc}.tsv"

    # Skip if file already exists and is non-empty
    if [[ -f "${outfile}" ]] && [[ -s "${outfile}" ]]; then
        echo "[SKIP] ${acc}: file already exists (${outfile})"
        (( SKIP_COUNT++ )) || true
        continue
    fi

    url="$(build_url "${acc}" "${FIELDS}")"
    echo "[INFO] Fetching ${acc} -> ${outfile}"

    if download_one "${url}" "${outfile}"; then
        # Validate: file should have at least a header row
        if [[ -f "${outfile}" ]] && [[ -s "${outfile}" ]]; then
            echo "[OK]   ${acc}: downloaded successfully"
            (( SUCCESS_COUNT++ )) || true
        else
            echo "[WARN] ${acc}: output file is empty"
            rm -f "${outfile}"
            (( FAIL_COUNT++ )) || true
        fi
    else
        echo "[FAIL] ${acc}: download error"
        rm -f "${outfile}"
        (( FAIL_COUNT++ )) || true
    fi

done < "${INPUT_FILE}"

echo ""
echo "==========================================="
echo " Summary"
echo "==========================================="
echo "  Total:   $(( SUCCESS_COUNT + FAIL_COUNT + SKIP_COUNT ))"
echo "  Success: ${SUCCESS_COUNT}"
echo "  Failed:  ${FAIL_COUNT}"
echo "  Skipped: ${SKIP_COUNT}"
echo "  Output:  ${OUT_DIR}/"
echo "==========================================="

if [[ "${FAIL_COUNT}" -gt 0 ]]; then
    exit 1
fi

exit 0
