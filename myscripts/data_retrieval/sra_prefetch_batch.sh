#!/usr/bin/env bash
# ==============================================================================
# File Name:    sra_prefetch_batch.sh
# Author:       ChengYu
# Description:  Batch download FASTQ files from NCBI SRA using prefetch and
#               fasterq-dump (or fastq-dump for bulk mode). Supports single-cell
#               and bulk RNA-seq modes, parallel execution via ParaFly, and
#               optional FastQC quality control.
# Created Time: 2026
# ==============================================================================
set -euo pipefail

readonly SCRIPT_NAME="$(basename "$0")"
readonly VERSION="2.0.0"

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [OPTIONS]

Download FASTQ files from NCBI SRA in batch.

Required:
  -i <file>       Input file with one SRR accession per line

Options:
  -d <dir>        Output directory (default: ./sra_download)
  -t <int>        Number of threads for parallel execution (default: 4)
  --sc            Single-cell mode: use fasterq-dump with --split-files
  --fastqc        Run FastQC on downloaded FASTQ files after transfer
  --prefetch-bin  Path to prefetch executable (default: auto-detect)
  --dump-bin      Path to fasterq-dump/fastq-dump executable (default: auto-detect)
  --parafly-bin   Path to ParaFly executable (default: auto-detect)
  --fastqc-bin    Path to fastqc executable (default: auto-detect)
  --version       Print version and exit
  -h, --help      Show this help message

Examples:
  # Bulk mode (fastq-dump), 8 threads
  ${SCRIPT_NAME} -i srr_list.txt -d ./fastq_output -t 8

  # Single-cell mode with FastQC
  ${SCRIPT_NAME} -i srr_list.txt -d ./sc_fastq --sc --fastqc -t 12

  # Specify custom tool paths
  ${SCRIPT_NAME} -i srr_list.txt --prefetch-bin /opt/sratoolkit/bin/prefetch \\
      --dump-bin /opt/sratoolkit/bin/fasterq-dump

Notes:
  - SRA Toolkit must be installed and configured (prefetch, fasterq-dump/fastq-dump).
  - ParaFly is used for parallel execution and must be available on PATH
    unless --parafly-bin is provided.
EOF
}

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
log()     { printf "[%(%Y-%m-%d %H:%M:%S)T] [INFO]  %s\n" -1 "$*"; }
warn()    { printf "[%(%Y-%m-%d %H:%M:%S)T] [WARN]  %s\n" -1 "$*" >&2; }
error()   { printf "[%(%Y-%m-%d %H:%M:%S)T] [ERROR] %s\n" -1 "$*" >&2; }
fatal()   { error "$@"; exit 1; }

# ---------------------------------------------------------------------------
# Auto-detect a tool on PATH, or validate a user-supplied path
# ---------------------------------------------------------------------------
resolve_bin() {
    local name="$1" user_path="$2"
    if [[ -n "${user_path}" ]]; then
        if [[ ! -x "${user_path}" ]]; then
            fatal "Provided ${name} path is not executable: ${user_path}"
        fi
        echo "${user_path}"
    else
        local found
        found="$(command -v "${name}" 2>/dev/null)" || true
        if [[ -z "${found}" ]]; then
            fatal "Could not find '${name}' on PATH. Install it or supply --${name}-bin."
        fi
        echo "${found}"
    fi
}

# ---------------------------------------------------------------------------
# Default values
# ---------------------------------------------------------------------------
INPUT_FILE=""
OUTDIR="./sra_download"
THREADS=4
SC_MODE=false
RUN_FASTQC=false
PREFETCH_BIN=""
DUMP_BIN=""
PARAFLY_BIN=""
FASTQC_BIN=""

# ---------------------------------------------------------------------------
# Parse options
# ---------------------------------------------------------------------------
OPTS="hi:d:t:"
LONGOPTS="help,version,sc,fastqc,prefetch-bin:,dump-bin:,parafly-bin:,fastqc-bin:"

PARSED="$(getopt -o "${OPTS}" -l "${LONGOPTS}" -n "${SCRIPT_NAME}" -- "$@")"
eval set -- "${PARSED}"

while true; do
    case "$1" in
        -h|--help)
            usage; exit 0 ;;
        --version)
            echo "${SCRIPT_NAME} ${VERSION}"; exit 0 ;;
        -i)
            INPUT_FILE="$2"; shift 2 ;;
        -d)
            OUTDIR="$2"; shift 2 ;;
        -t)
            THREADS="$2"; shift 2 ;;
        --sc)
            SC_MODE=true; shift ;;
        --fastqc)
            RUN_FASTQC=true; shift ;;
        --prefetch-bin)
            PREFETCH_BIN="$2"; shift 2 ;;
        --dump-bin)
            DUMP_BIN="$2"; shift 2 ;;
        --parafly-bin)
            PARAFLY_BIN="$2"; shift 2 ;;
        --fastqc-bin)
            FASTQC_BIN="$2"; shift 2 ;;
        --)
            shift; break ;;
        *)
            fatal "Unexpected option: $1"
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate required arguments
# ---------------------------------------------------------------------------
if [[ -z "${INPUT_FILE}" ]]; then
    fatal "Input file is required. Use -i <file>."
fi
if [[ ! -f "${INPUT_FILE}" ]]; then
    fatal "Input file does not exist: ${INPUT_FILE}"
fi

# ---------------------------------------------------------------------------
# Resolve tool paths
# ---------------------------------------------------------------------------
log "Resolving tool paths..."

if "${SC_MODE}"; then
    DUMP_BIN="$(resolve_bin fasterq-dump "${DUMP_BIN}")"
else
    DUMP_BIN="$(resolve_bin fastq-dump "${DUMP_BIN}")"
fi

PREFETCH_BIN="$(resolve_bin prefetch "${PREFETCH_BIN}")"
PARAFLY_BIN="$(resolve_bin ParaFly "${PARAFLY_BIN}")"

if "${RUN_FASTQC}"; then
    FASTQC_BIN="$(resolve_bin fastqc "${FASTQC_BIN}")"
fi

log "  prefetch    : ${PREFETCH_BIN}"
log "  dump tool   : ${DUMP_BIN}"
log "  parafly     : ${PARAFLY_BIN}"
if "${RUN_FASTQC}"; then
    log "  fastqc      : ${FASTQC_BIN}"
fi

# ---------------------------------------------------------------------------
# Prepare output directory
# ---------------------------------------------------------------------------
mkdir -p "${OUTDIR}"
SRA_CACHE="${OUTDIR}/sra_cache"
mkdir -p "${SRA_CACHE}"
CMD_LOG="${OUTDIR}/commands.log"
FAILED_LOG="${OUTDIR}/failed_accessions.txt"

# Read accessions, stripping blank lines and comments
mapfile -t ACCESSIONS < <(sed '/^[[:space:]]*$/d; /^#/d' "${INPUT_FILE}")
TOTAL="${#ACCESSIONS[@]}"

if [[ "${TOTAL}" -eq 0 ]]; then
    fatal "No accessions found in ${INPUT_FILE}"
fi

log "Found ${TOTAL} accession(s) to process"
log "Mode: $("${SC_MODE}" && echo "single-cell (fasterq-dump --split-files)" || echo "bulk (fastq-dump)")"
log "Threads: ${THREADS}"

# ---------------------------------------------------------------------------
# Build per-accession commands
# ---------------------------------------------------------------------------
log "Building download commands..."

> "${CMD_LOG}"

for srr in "${ACCESSIONS[@]}"; do
    srr_trim="$(echo "${srr}" | xargs)"  # trim whitespace
    if [[ -z "${srr_trim}" ]]; then
        continue
    fi

    # prefetch step
    prefetch_cmd="${PREFETCH_BIN} ${srr_trim} -O ${SRA_CACHE} --max-size 100G"

    # dump step depends on mode
    sra_file="${SRA_CACHE}/${srr_trim}/${srr_trim}.sra"
    if "${SC_MODE}"; then
        dump_cmd="${DUMP_BIN} ${sra_file} --split-files --threads ${THREADS} --outdir ${OUTDIR}"
    else
        dump_cmd="${DUMP_BIN} ${sra_file} --gzip --split-3 --threads ${THREADS} --outdir ${OUTDIR}"
    fi

    # Combine into one compound command with error checking
    echo "${prefetch_cmd} && ${dump_cmd}" >> "${CMD_LOG}"
done

# ---------------------------------------------------------------------------
# Execute commands in parallel via ParaFly
# ---------------------------------------------------------------------------
log "Starting parallel download (${THREADS} threads)..."

if ! "${PARAFLY_BIN}" -c "${CMD_LOG}" -CPU "${THREADS}" -shfailed "${FAILED_LOG}"; then
    warn "ParaFly reported some failures. Check ${FAILED_LOG} for details."
    if [[ -f "${FAILED_LOG}" && -s "${FAILED_LOG}" ]]; then
        NUM_FAILED="$(wc -l < "${FAILED_LOG}")"
        warn "  Failed accession count: ${NUM_FAILED}"
    fi
fi

log "Download phase complete."

# ---------------------------------------------------------------------------
# Optional FastQC
# ---------------------------------------------------------------------------
if "${RUN_FASTQC}"; then
    log "Running FastQC on downloaded FASTQ files..."
    mapfile -t FASTQ_FILES < <(find "${OUTDIR}" -maxdepth 1 -name '*.fastq*' -type f)
    if [[ ${#FASTQ_FILES[@]} -gt 0 ]]; then
        "${FASTQC_BIN}" -t "${THREADS}" -o "${OUTDIR}/fastqc_output" "${FASTQ_FILES[@]}" || true
        log "FastQC complete. Results in ${OUTDIR}/fastqc_output/"
    else
        warn "No FASTQ files found for FastQC."
    fi
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
DOWNLOADED="$(find "${OUTDIR}" -maxdepth 1 -name '*.fastq*' -type f | wc -l)"
log "=== Summary ==="
log "  Output directory : ${OUTDIR}"
log "  FASTQ files      : ${DOWNLOADED}"
if [[ -f "${FAILED_LOG}" && -s "${FAILED_LOG}" ]]; then
    log "  Failed accessions: see ${FAILED_LOG}"
fi
log "Done."
