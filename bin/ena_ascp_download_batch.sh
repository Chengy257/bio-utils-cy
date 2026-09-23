#!/usr/bin/env bash
# ==============================================================================
# File Name:    ena_ascp_download_batch.sh
# Author:       ChengYu
# Description:  Batch download FASTQ files from ENA using Aspera ascp.
#               Constructs the ENA remote path from accession IDs, supporting
#               both single-end (SE) and paired-end (PE) file naming
#               conventions.
# Created Time: 2026
# ==============================================================================
set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf

readonly SCRIPT_NAME="$(basename "$0")"
readonly VERSION="2.0.0"

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [OPTIONS]

Batch download FASTQ from ENA via Aspera ascp by constructing remote paths
from accession IDs (SRR/ERR/DRR).

Required:
  -i <file>       Input file with one accession per line

Options:
  -o <dir>        Output directory (default: ./ena_download)
  -b <rate>       Bandwidth limit, e.g. 500M or 1G (default: 500M)
  -s              Single-end mode: download _1.fastq.gz only
  --ascp-bin      Path to ascp executable (default: auto-detect)
  --key           Path to Aspera private key (default: auto-detect)
  --host          ENA Aspera host prefix (default: era-fasp@fasp.sra.ebi.ac.uk:)
  -t <int>        Number of parallel ascp processes (default: 4)
  --version       Print version and exit
  -h, --help      Show this help message

ENA path convention:
  For an accession like SRR1234567 the ENA path is built from the accession
  prefix and length:
    SRR      ->  /vol1/fastq/SRR123/SRR1234567/SRR1234567_1.fastq.gz
  For 3-character suffixes (e.g. SRR12345678):
    SRR      ->  /vol1/fastq/SRR123/SRR12345678/SRR12345678_1.fastq.gz
  The script handles SRR/ERR/DRR prefixes automatically.

Examples:
  # Paired-end download
  ${SCRIPT_NAME} -i accessions.txt -o ./fastq -b 1G

  # Single-end download with 8 parallel transfers
  ${SCRIPT_NAME} -i accessions.txt -o ./fastq -s -t 8

  # Custom ascp location
  ${SCRIPT_NAME} -i accessions.txt --ascp-bin ~/aspera/bin/ascp \\
      --key ~/aspera/etc/aspera_id_rsa
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
# Resolve a binary: user path or auto-detect via command -v
# ---------------------------------------------------------------------------
resolve_bin() {
    local name="$1" user_path="$2"
    if [[ -n "${user_path}" ]]; then
        if [[ ! -x "${user_path}" ]]; then
            fatal "Provided ${name} is not executable: ${user_path}"
        fi
        echo "${user_path}"
    else
        local found
        found="$(command -v "${name}" 2>/dev/null)" || true
        if [[ -z "${found}" ]]; then
            fatal "Cannot find '${name}' on PATH. Install it or supply --${name}."
        fi
        echo "${found}"
    fi
}

# ---------------------------------------------------------------------------
# Resolve Aspera key: search near ascp binary or use provided path
# ---------------------------------------------------------------------------
resolve_key() {
    local ascp_bin="$1" user_key="$2"
    if [[ -n "${user_key}" ]]; then
        if [[ ! -f "${user_key}" ]]; then
            fatal "Provided key file does not exist: ${user_key}"
        fi
        echo "${user_key}"
        return
    fi

    local bindir rootdir
    bindir="$(dirname "$(readlink -f "${ascp_bin}")")"
    rootdir="$(dirname "${bindir}")"

    local candidates=(
        "${bindir}/asperaworkshop.pem"
        "${bindir}/aspera_id_rsa"
        "${bindir}/aspera.openssh"
        "${rootdir}/etc/asperaworkshop.pem"
        "${rootdir}/etc/aspera_id_rsa"
        "${rootdir}/etc/aspera.openssh"
    )

    for c in "${candidates[@]}"; do
        if [[ -f "${c}" ]]; then
            echo "${c}"
            return
        fi
    done

    fatal "Cannot find Aspera private key near ${ascp_bin}. Provide --key explicitly."
}

# ---------------------------------------------------------------------------
# Build the ENA remote directory path from an accession
# ---------------------------------------------------------------------------
# ENA uses the following convention:
#   Accession SRR1234567 (length 10, 3-char numeric suffix) ->
#     /vol1/fastq/SRR123/00{last_digit}/SRR1234567/
#   Accession SRR12345678 (length 11, 4-char numeric suffix) ->
#     /vol1/fastq/SRR123/0{first_two_of_suffix}/SRR12345678/
#   Accession SRR123456789 (length 12, 5-char numeric suffix) ->
#     /vol1/fastq/SRR123/{first_three_of_suffix}/SRR123456789/
# General rule:
#   prefix = first 6 characters (e.g. SRR123)
#   subdir = last 3 digits of accession, zero-padded to 3 digits for path
#            Actually ENA uses: take the last 3 digits of accession and place
#            them as "00X" where X is the last digit.
#   Simpler rule that matches ENA:
#     For accessions of length 10 (standard): /vol1/fastq/{prefix}/{acc}/
#     For accessions of length 11:            /vol1/fastq/{prefix}/00{last_digit}/
#     For accessions of length 12+:           /vol1/fastq/{prefix}/0{last_two}/
#
# The most reliable approach used by ENA documentation:
#   1. Take first 6 chars as prefix directory (e.g., SRR123)
#   2. The sub-directory uses the numeric portion:
#      - For 7-digit numeric (length=10): /vol1/fastq/{prefix}/{acc}/
#      - For 8-digit numeric (length=11): /vol1/fastq/{prefix}/00{last_digit}/{acc}/
#      - For 9-digit numeric (length=12): /vol1/fastq/{prefix}/0{last_two}/{acc}/
# ---------------------------------------------------------------------------
build_ena_path() {
    local acc="$1"
    local prefix subdir acc_len

    # Extract the 6-char prefix (e.g., SRR123, ERR456)
    prefix="${acc:0:6}"
    acc_len="${#acc}"

    case "${acc_len}" in
        10)
            # Standard: /vol1/fastq/SRR123/SRR1234567/
            echo "/vol1/fastq/${prefix}/${acc}"
            ;;
        11)
            # 8-digit suffix: /vol1/fastq/SRR123/007/SRR12345678/
            local last_digit="${acc: -1}"
            echo "/vol1/fastq/${prefix}/00${last_digit}/${acc}"
            ;;
        12)
            # 9-digit suffix: /vol1/fastq/SRR123/078/SRR123456789/
            local last_two="${acc: -2}"
            echo "/vol1/fastq/${prefix}/0${last_two}/${acc}"
            ;;
        *)
            # Fallback: assume standard layout
            warn "Unexpected accession length ${acc_len} for ${acc}, using fallback path"
            echo "/vol1/fastq/${prefix}/${acc}"
            ;;
    esac
}

# ---------------------------------------------------------------------------
# Default values
# ---------------------------------------------------------------------------
INPUT_FILE=""
OUTDIR="./ena_download"
BANDWIDTH="500M"
SE_MODE=false
ASCP_BIN="${MYS_ASCP_BIN:-}"
ASCP_KEY=""
HOST="era-fasp@fasp.sra.ebi.ac.uk:"
THREADS=4

# ---------------------------------------------------------------------------
# Parse options
# ---------------------------------------------------------------------------
OPTS="hi:o:b:st:"
LONGOPTS="help,version,ascp-bin:,key:,host:"

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
        -o)
            OUTDIR="$2"; shift 2 ;;
        -b)
            BANDWIDTH="$2"; shift 2 ;;
        -s)
            SE_MODE=true; shift ;;
        -t)
            THREADS="$2"; shift 2 ;;
        --ascp-bin)
            ASCP_BIN="$2"; shift 2 ;;
        --key)
            ASCP_KEY="$2"; shift 2 ;;
        --host)
            HOST="$2"; shift 2 ;;
        --)
            shift; break ;;
        *)
            fatal "Unexpected option: $1"
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
if [[ -z "${INPUT_FILE}" ]]; then
    fatal "Input file is required. Use -i <file>."
fi
if [[ ! -f "${INPUT_FILE}" ]]; then
    fatal "Input file does not exist: ${INPUT_FILE}"
fi

# ---------------------------------------------------------------------------
# Resolve tools
# ---------------------------------------------------------------------------
log "Resolving tool paths..."
ASCP_BIN="$(resolve_bin ascp "${ASCP_BIN}")"
ASCP_KEY="$(resolve_key "${ASCP_BIN}" "${ASCP_KEY}")"

log "  ascp      : ${ASCP_BIN}"
log "  key       : ${ASCP_KEY}"
log "  host      : ${HOST}"
log "  bandwidth : ${BANDWIDTH}"
log "  threads   : ${THREADS}"
log "  mode      : $("${SE_MODE}" && echo "single-end" || echo "paired-end")"

# ---------------------------------------------------------------------------
# Prepare output directory and command list
# ---------------------------------------------------------------------------
mkdir -p "${OUTDIR}"

CMD_LOG="${OUTDIR}/.ascp_commands.txt"
FAILED_LOG="${OUTDIR}/failed_accessions.txt"
> "${CMD_LOG}"
> "${FAILED_LOG}"

# ---------------------------------------------------------------------------
# Build ascp commands for each accession
# ---------------------------------------------------------------------------
mapfile -t ACCESSIONS < <(sed '/^[[:space:]]*$/d; /^#/d' "${INPUT_FILE}")
TOTAL="${#ACCESSIONS[@]}"

if [[ "${TOTAL}" -eq 0 ]]; then
    fatal "No accessions found in ${INPUT_FILE}"
fi

log "Building download commands for ${TOTAL} accession(s)..."

for acc_raw in "${ACCESSIONS[@]}"; do
    acc="$(echo "${acc_raw}" | xargs)"
    [[ -z "${acc}" ]] && continue

    ena_dir="$(build_ena_path "${acc}")"
    remote_base="${HOST}${ena_dir}/${acc}"

    if "${SE_MODE}"; then
        # Single-end: only _1.fastq.gz
        remote_file="${remote_base}_1.fastq.gz"
        echo "${ASCP_BIN} -T -l ${BANDWIDTH} -P 33001 -i ${ASCP_KEY} -Q ${remote_file} ${OUTDIR}"
    else
        # Paired-end: both _1 and _2
        remote_r1="${remote_base}_1.fastq.gz"
        remote_r2="${remote_base}_2.fastq.gz"
        echo "${ASCP_BIN} -T -l ${BANDWIDTH} -P 33001 -i ${ASCP_KEY} -Q ${remote_r1} ${OUTDIR}"
        echo "${ASCP_BIN} -T -l ${BANDWIDTH} -P 33001 -i ${ASCP_KEY} -Q ${remote_r2} ${OUTDIR}"
    fi
done >> "${CMD_LOG}"

NUM_CMDS="$(wc -l < "${CMD_LOG}")"
log "Total ascp commands: ${NUM_CMDS}"

# ---------------------------------------------------------------------------
# Execute in parallel using xargs
# ---------------------------------------------------------------------------
log "Starting parallel downloads..."

FAILED=0
while IFS= read -r cmd; do
    # Extract the filename for logging
    remote_ref="$(echo "${cmd}" | grep -oP '(?<=-Q )\S+')"
    fname="$(basename "${remote_ref}" 2>/dev/null || echo "unknown")"

    log "Downloading ${fname} ..."
    if eval "${cmd}" >> "${OUTDIR}/ascp.log" 2>&1; then
        log "  OK: ${fname}"
    else
        warn "  FAILED: ${fname}"
        echo "${acc}" >> "${FAILED_LOG}"
        FAILED=$((FAILED + 1))
    fi
done < <(cat "${CMD_LOG}") &
PID_MAIN=$!

# For true parallelism, we use xargs as alternative
# Kill the serial loop above and use xargs instead
kill "${PID_MAIN}" 2>/dev/null || true
wait "${PID_MAIN}" 2>/dev/null || true

log "Launching parallel ascp via xargs (-P ${THREADS})..."

if ! xargs -P "${THREADS}" -I {} bash -c '
    cmd="{}";
    remote_ref=$(echo "${cmd}" | grep -oP "(?<=-Q )\S+");
    fname=$(basename "${remote_ref}" 2>/dev/null || echo "unknown");
    printf "[%(%Y-%m-%d %H:%M:%S)T] [INFO]  Downloading %s ...\n" -1 "${fname}";
    if eval "${cmd}" >> "'"${OUTDIR}"'/ascp.log" 2>&1; then
        printf "[%(%Y-%m-%d %H:%M:%S)T] [INFO]  OK: %s\n" -1 "${fname}";
    else
        printf "[%(%Y-%m-%d %H:%M:%S)T] [WARN]  FAILED: %s\n" -1 "${fname}" >&2;
        # Extract accession from filename
        acc=$(echo "${fname}" | sed "s/_[12]\.fastq\.gz$//");
        echo "${acc}" >> "'"${FAILED_LOG}"'";
    fi
' < "${CMD_LOG}"; then
    warn "Some downloads failed. See ${FAILED_LOG}"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
DOWNLOADED="$(find "${OUTDIR}" -maxdepth 1 -name '*.fastq.gz' -type f | wc -l)"
FAILED_COUNT="$(sort -u "${FAILED_LOG}" | sed '/^$/d' | wc -l)"

log "=== Summary ==="
log "  Output directory : ${OUTDIR}"
log "  FASTQ files      : ${DOWNLOADED}"
log "  Failed accessions: ${FAILED_COUNT}"
if [[ "${FAILED_COUNT}" -gt 0 ]]; then
    log "  Failed list      : ${FAILED_LOG}"
fi
log "Done."
