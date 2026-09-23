#!/usr/bin/env bash
# ==============================================================================
# File Name:    sra_prefetch_batch.sh
# Author:       ChengYu
# Description:  Batch download FASTQ files from NCBI SRA using prefetch and
#               fastq-dump/fasterq-dump. Parallel execution via GNU parallel.
#               Supports single-cell and bulk RNA-seq modes, resume from
#               partial runs, and optional FastQC quality control.
#
# Created Time: 2026
#
# Changelog (vs v2.0.0):
#   v2.1.0  2026-09-19
#   - FIX: bulk mode no longer passes --threads to fastq-dump. fastq-dump has
#     no --threads option (it is a fasterq-dump option); the old script made
#     every conversion fail with "param incorrect ... --threads" and then a
#     bogus accession lookup of the thread count (404). Dump flags are now
#     chosen from the actual dump binary flavor (fastq-dump vs fasterq-dump).
#   - Parallel execution now uses GNU parallel (replaces ParaFly and the
#     --parafly-bin option; a --parallel-bin override is provided instead).
#   - FIX: resolve_bin requires an executable FILE (directories rejected) and
#     error hints name the correct option (--dump-bin, not --fastq-dump-bin).
#   - Resume support: accessions whose FASTQ already exists are skipped; a
#     cached .sra is converted directly (no re-download) with an automatic
#     prefetch fallback in case the cache is incomplete.
#   - Atomic outputs: conversion writes to <outdir>/.staging/<SRR>/ and only
#     moves finished FASTQ into <outdir>, so a crashed dump can never leave a
#     partial file that looks "already done".
#   - Accessions are deduplicated (order preserved); CRLF and inline
#     whitespace tolerated.
#   - New options: --max-size, --dump-threads, --dry-run.
#   - Per-accession tool logs under <outdir>/worker_logs/<SRR>.log and
#     per-accession exit codes in <outdir>/worker_logs/<SRR>.rc, plus a GNU
#     parallel joblog at <outdir>/parallel_joblog.tsv.
#   - Exit code 2 if any accession failed.
# ==============================================================================
set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf
shopt -s nullglob

readonly SCRIPT_NAME="$(basename "$0")"
readonly VERSION="2.1.0"

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [OPTIONS]

Download FASTQ files from NCBI SRA in batch (prefetch + fastq-dump/fasterq-dump,
parallel via GNU parallel, resume-aware).

Required:
  -i <file>       Input file with one SRR accession per line
                  (blank lines and #-comments ignored, duplicates removed)

Options:
  -d <dir>        Output directory (default: ./sra_download)
  -t <int>        Number of parallel download workers (default: 4)
  --sc            Single-cell mode: use fasterq-dump with --split-files
  --fastqc        Run FastQC on downloaded FASTQ files after transfer
  --prefetch-bin  Path to prefetch executable (default: auto-detect on PATH)
  --dump-bin      Path to fasterq-dump/fastq-dump executable (default: auto-detect)
  --parallel-bin  Path to GNU parallel executable (default: auto-detect)
  --fastqc-bin    Path to fastqc executable (default: auto-detect)
  --max-size <s>  prefetch --max-size value (default: 100G)
  --dump-threads <int>  Only for fasterq-dump: --threads per worker
                  (default: unset = tool default). fastq-dump ignores this.
  --dry-run       Print the planned actions per accession and exit
  --version       Print version and exit
  -h, --help      Show this help message

Examples:
  # Bulk mode (fastq-dump), 12 workers
  ${SCRIPT_NAME} -i srr_list.txt -d ./fastq_output -t 12 \\
      --prefetch-bin /opt/anaconda3/envs/download/bin/prefetch \\
      --dump-bin /opt/anaconda3/envs/download/bin/fastq-dump

  # Single-cell mode with FastQC
  ${SCRIPT_NAME} -i srr_list.txt -d ./sc_fastq --sc --fastqc -t 12

  # Preview what would be done (downloads nothing)
  ${SCRIPT_NAME} -i srr_list.txt -d ./fastq_output --dry-run

Notes:
  - Parallel execution uses GNU parallel. A joblog with per-job exit codes is
    written to <outdir>/parallel_joblog.tsv.
  - Re-running the same command is safe: finished accessions are skipped and
    .sra files already in <outdir>/sra_cache are converted without
    re-downloading. Do NOT run two instances on the same output directory at
    the same time (they would fight over the same cache).
  - fastq-dump is single-threaded; parallelism comes from -t workers, not
    from --threads. Use --dump-threads only in fasterq-dump (--sc) workflows.
  - Tool output per accession is logged to <outdir>/worker_logs/<SRR>.log.
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
    local name="$1" user_path="$2" opt_hint="$3"
    if [[ -n "${user_path}" ]]; then
        if [[ -d "${user_path}" ]]; then
            fatal "Provided ${name} path is a DIRECTORY, not an executable file: ${user_path} (check --${opt_hint})"
        fi
        if [[ ! -f "${user_path}" || ! -x "${user_path}" ]]; then
            fatal "Provided ${name} path is not an executable file: ${user_path} (check --${opt_hint})"
        fi
        echo "${user_path}"
    else
        local found
        found="$(command -v "${name}" 2>/dev/null)" || true
        if [[ -z "${found}" ]]; then
            fatal "Could not find '${name}' on PATH. Install it or supply --${opt_hint}."
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
PARALLEL_BIN=""
FASTQC_BIN=""
MAX_SIZE="100G"
DUMP_THREADS=""
DRY_RUN=false

# ---------------------------------------------------------------------------
# Parse options
# ---------------------------------------------------------------------------
OPTS="hi:d:t:"
LONGOPTS="help,version,sc,fastqc,dry-run,prefetch-bin:,dump-bin:,parallel-bin:,fastqc-bin:,max-size:,dump-threads:"

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
        --parallel-bin)
            PARALLEL_BIN="$2"; shift 2 ;;
        --fastqc-bin)
            FASTQC_BIN="$2"; shift 2 ;;
        --max-size)
            MAX_SIZE="$2"; shift 2 ;;
        --dump-threads)
            DUMP_THREADS="$2"; shift 2 ;;
        --dry-run)
            DRY_RUN=true; shift ;;
        --)
            shift; break ;;
        *)
            fatal "Unexpected option: $1"
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate arguments
# ---------------------------------------------------------------------------
if [[ -z "${INPUT_FILE}" ]]; then
    fatal "Input file is required. Use -i <file>."
fi
if [[ ! -f "${INPUT_FILE}" ]]; then
    fatal "Input file does not exist: ${INPUT_FILE}"
fi
if ! [[ "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    fatal "Invalid thread count: ${THREADS} (expected a positive integer)"
fi
if [[ -n "${DUMP_THREADS}" ]] && ! [[ "${DUMP_THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    fatal "Invalid --dump-threads value: ${DUMP_THREADS} (expected a positive integer)"
fi
if [[ -z "${MAX_SIZE}" ]]; then
    fatal "--max-size must not be empty (e.g. 100G)"
fi

# ---------------------------------------------------------------------------
# Resolve tool paths and pick dump flags by the actual binary flavor
# ---------------------------------------------------------------------------
log "Resolving tool paths..."

if "${SC_MODE}"; then
    DUMP_NAME="fasterq-dump"
else
    DUMP_NAME="fastq-dump"
fi
DUMP_BIN="$(resolve_bin "${DUMP_NAME}" "${DUMP_BIN}" "dump-bin")"
PREFETCH_BIN="$(resolve_bin prefetch "${PREFETCH_BIN}" "prefetch-bin")"
PARALLEL_BIN="$(resolve_bin parallel "${PARALLEL_BIN}" "parallel-bin")"

if "${RUN_FASTQC}"; then
    FASTQC_BIN="$(resolve_bin fastqc "${FASTQC_BIN}" "fastqc-bin")"
fi

DUMP_FLAVOR="fastq-dump"
case "$(basename "${DUMP_BIN}")" in
    *fasterq*) DUMP_FLAVOR="fasterq-dump" ;;
    *fastq*)   DUMP_FLAVOR="fastq-dump" ;;
    *) warn "Cannot recognize dump tool flavor from '${DUMP_BIN}'; assuming fastq-dump flags." ;;
esac

DUMP_FLAGS=()
if [[ "${DUMP_FLAVOR}" == "fastq-dump" ]]; then
    if "${SC_MODE}"; then
        DUMP_FLAGS=(--split-files)
        warn "SC mode with fastq-dump is very slow; fasterq-dump (--sc without --dump-bin override) is recommended."
    else
        DUMP_FLAGS=(--gzip --split-3)
    fi
    if [[ -n "${DUMP_THREADS}" ]]; then
        warn "--dump-threads ignored: fastq-dump does not support --threads (parallelism comes from -t workers)."
    fi
else
    if "${SC_MODE}"; then
        DUMP_FLAGS=(--split-files)
    else
        DUMP_FLAGS=(--split-files --gzip)
        warn "Bulk mode with fasterq-dump: --split-3 is unavailable, using --split-files instead."
    fi
    if [[ -n "${DUMP_THREADS}" ]]; then
        DUMP_FLAGS+=(--threads "${DUMP_THREADS}")
    fi
fi

log "  prefetch    : ${PREFETCH_BIN}"
log "  dump tool   : ${DUMP_BIN} [flavor: ${DUMP_FLAVOR}]"
log "  dump flags  : ${DUMP_FLAGS[*]}"
log "  parallel    : ${PARALLEL_BIN}"
if "${RUN_FASTQC}"; then
    log "  fastqc      : ${FASTQC_BIN}"
fi

# ---------------------------------------------------------------------------
# Prepare paths (nothing is created in --dry-run mode)
# ---------------------------------------------------------------------------
SRA_CACHE="${OUTDIR}/sra_cache"
CMD_LOG="${OUTDIR}/commands.log"
FAILED_LOG="${OUTDIR}/failed_accessions.txt"
STAGING_DIR="${OUTDIR}/.staging"
WORKER_LOG_DIR="${OUTDIR}/worker_logs"
JOBS_FILE="${OUTDIR}/jobs_to_run.txt"
JOBLOG="${OUTDIR}/parallel_joblog.tsv"

# ---------------------------------------------------------------------------
# Read accessions: strip CR/blank lines/comments, take first token, dedup
# ---------------------------------------------------------------------------
mapfile -t CLEANED < <(sed 's/\r$//' "${INPUT_FILE}" \
    | grep -v '^[[:space:]]*$' \
    | grep -v '^[[:space:]]*#' \
    | awk '{print $1}')
mapfile -t ACCESSIONS < <(printf '%s\n' ${CLEANED[@]+"${CLEANED[@]}"} | awk '!seen[$0]++')

TOTAL="${#ACCESSIONS[@]}"
if [[ "${TOTAL}" -eq 0 ]]; then
    fatal "No accessions found in ${INPUT_FILE}"
fi

# ---------------------------------------------------------------------------
# Classify: skip finished accessions, remember which have a cached .sra
# ---------------------------------------------------------------------------
is_done() {
    # A finished accession has <SRR>.fastq* or <SRR>_*.fastq* in OUTDIR
    local srr="$1" f
    for f in "${OUTDIR}/${srr}.fastq"* "${OUTDIR}/${srr}_"*".fastq"*; do
        [[ -f "${f}" ]] && return 0
    done
    return 1
}

JOB_SRRS=()
JOB_CACHED=()
SKIP_DONE=0
for srr in "${ACCESSIONS[@]}"; do
    if is_done "${srr}"; then
        SKIP_DONE=$((SKIP_DONE + 1))
        log "  [skip] ${srr}: FASTQ already present in ${OUTDIR}"
    else
        sra_file="${SRA_CACHE}/${srr}/${srr}.sra"
        if [[ -s "${sra_file}" ]]; then
            JOB_SRRS+=("${srr}")
            JOB_CACHED+=("yes")
        else
            JOB_SRRS+=("${srr}")
            JOB_CACHED+=("no")
        fi
    fi
done

N_JOBS="${#JOB_SRRS[@]}"
if "${SC_MODE}"; then MODE_DESC="single-cell"; else MODE_DESC="bulk"; fi

log "Found ${TOTAL} unique accession(s) in ${INPUT_FILE}"
log "Mode: ${MODE_DESC} | Workers: ${THREADS} | prefetch --max-size ${MAX_SIZE}"
log "Already done: ${SKIP_DONE} | To process this run: ${N_JOBS}"

# ---------------------------------------------------------------------------
# Dry run: print planned actions and exit
# ---------------------------------------------------------------------------
if "${DRY_RUN}"; then
    log "=== DRY RUN plan (${N_JOBS} action(s), ${SKIP_DONE} skipped) ==="
    for idx in "${!JOB_SRRS[@]}"; do
        srr="${JOB_SRRS[$idx]}"
        sra_file="${SRA_CACHE}/${srr}/${srr}.sra"
        dump_cmd="${DUMP_BIN} ${sra_file} ${DUMP_FLAGS[*]} --outdir ${STAGING_DIR}/${srr} && mv ${STAGING_DIR}/${srr}/*.fastq* ${OUTDIR}/"
        if [[ "${JOB_CACHED[$idx]}" == "yes" ]]; then
            printf '  [convert cached .sra] %s\n    %s\n    (fallback on failure: %s %s -O %s --max-size %s, then retry)\n' \
                "${srr}" "${dump_cmd}" "${PREFETCH_BIN}" "${srr}" "${SRA_CACHE}" "${MAX_SIZE}"
        else
            printf '  [prefetch + dump] %s\n    %s %s -O %s --max-size %s &&\n    %s\n' \
                "${srr}" "${PREFETCH_BIN}" "${srr}" "${SRA_CACHE}" "${MAX_SIZE}" "${dump_cmd}"
        fi
    done
    log "Dry run complete. No directories were created, nothing was downloaded."
    exit 0
fi

# ---------------------------------------------------------------------------
# Real run: create directories and record the planned commands
# ---------------------------------------------------------------------------
mkdir -p "${OUTDIR}" "${SRA_CACHE}" "${STAGING_DIR}" "${WORKER_LOG_DIR}"

> "${CMD_LOG}"
for idx in "${!JOB_SRRS[@]}"; do
    srr="${JOB_SRRS[$idx]}"
    sra_file="${SRA_CACHE}/${srr}/${srr}.sra"
    dump_cmd="${DUMP_BIN} ${sra_file} ${DUMP_FLAGS[*]} --outdir ${STAGING_DIR}/${srr} && mv ${STAGING_DIR}/${srr}/*.fastq* ${OUTDIR}/"
    if [[ "${JOB_CACHED[$idx]}" == "yes" ]]; then
        echo "${dump_cmd}" >> "${CMD_LOG}"
    else
        echo "${PREFETCH_BIN} ${srr} -O ${SRA_CACHE} --max-size ${MAX_SIZE} && ${dump_cmd}" >> "${CMD_LOG}"
    fi
done

# ---------------------------------------------------------------------------
# Worker: process one accession (invoked by GNU parallel in child shells)
# ---------------------------------------------------------------------------
process_accession() {
    local srr="$1"
    shopt -s nullglob

    # DUMP_FLAGS crosses into this child as a plain string (arrays cannot be
    # exported); split it back. Flags themselves never contain spaces.
    local -a FLAGS=()
    read -ra FLAGS <<< "${DUMP_FLAGS_STR:-}"

    local sra_file="${SRA_CACHE}/${srr}/${srr}.sra"
    local staging="${STAGING_DIR}/${srr}"
    local wlog="${WORKER_LOG_DIR}/${srr}.log"

    dump_once() {
        # Clean staging so anything left in it belongs to THIS attempt
        rm -rf "${staging}"
        mkdir -p "${staging}"
        if ! "${DUMP_BIN}" "${sra_file}" "${FLAGS[@]}" --outdir "${staging}" >>"${wlog}" 2>&1; then
            rm -rf "${staging}"
            return 1
        fi
        local produced=( "${staging}"/*.fastq* )
        if [[ ${#produced[@]} -eq 0 ]]; then
            error "[${srr}] dump exited 0 but produced no FASTQ files; see ${wlog}"
            rm -rf "${staging}"
            return 1
        fi
        if ! mv "${staging}"/*.fastq* "${OUTDIR}/" >>"${wlog}" 2>&1; then
            rm -rf "${staging}"
            return 1
        fi
        rmdir "${staging}" 2>/dev/null || true
        return 0
    }

    : > "${wlog}"
    log "[RUN]  ${srr}"

    # Fast path: .sra already cached -> convert without re-downloading
    if [[ -s "${sra_file}" ]]; then
        if dump_once; then
            log "[OK]   ${srr} (converted from cached .sra)"
            return 0
        fi
        error "[${srr}] conversion from cached .sra failed; cache may be incomplete, re-running prefetch to repair"
    fi

    if ! "${PREFETCH_BIN}" "${srr}" -O "${SRA_CACHE}" --max-size "${MAX_SIZE}" >>"${wlog}" 2>&1; then
        error "[${srr}] prefetch failed; see ${wlog}"
        return 1
    fi
    if dump_once; then
        log "[OK]   ${srr}"
        return 0
    fi
    error "[${srr}] fastq conversion failed; see ${wlog}"
    return 1
}

worker_main() {
    local srr="$1" rc=0
    if [[ -z "${DUMP_FLAGS_STR:-}" ]]; then
        error "[${srr}] internal error: DUMP_FLAGS_STR not transferred to worker"
        rc=9
    else
        process_accession "${srr}" || rc=$?
    fi
    printf '%s' "${rc}" > "${WORKER_LOG_DIR}/${srr}.rc"
    return "${rc}"
}

# Arrays cannot cross into GNU parallel workers; ship the flags as a plain
# string and re-split them inside the worker.
DUMP_FLAGS_STR="${DUMP_FLAGS[*]}"
export OUTDIR SRA_CACHE STAGING_DIR WORKER_LOG_DIR PREFETCH_BIN DUMP_BIN MAX_SIZE DUMP_FLAGS_STR
export -f process_accession worker_main log warn error

# ---------------------------------------------------------------------------
# Execute in parallel via GNU parallel
# ---------------------------------------------------------------------------
printf '%s\n' ${JOB_SRRS[@]+"${JOB_SRRS[@]}"} > "${JOBS_FILE}"
rm -f "${WORKER_LOG_DIR}"/*.rc

PARALLEL_PID=""
on_interrupt() {
    error "Interrupted; stopping dispatcher and running workers..."
    if [[ -n "${PARALLEL_PID}" ]]; then
        # setsid makes parallel a process-group leader: negative pid kills the
        # whole tree (parallel + all prefetch/fastq-dump children). GNU
        # parallel itself also propagates INT/TERM to running jobs.
        kill -TERM -- -"${PARALLEL_PID}" 2>/dev/null || kill -TERM "${PARALLEL_PID}" 2>/dev/null || true
    fi
    exit 130
}
trap on_interrupt INT TERM

if [[ "${N_JOBS}" -gt 0 ]]; then
    log "Starting parallel download (${THREADS} workers, ${N_JOBS} accession(s)) via GNU parallel..."
    rm -f "${JOBLOG}"

    PARALLEL_RUN=(
        "${PARALLEL_BIN}" --will-cite
        -j "${THREADS}"
        --halt never
        --joblog "${JOBLOG}"
        --env OUTDIR --env SRA_CACHE --env STAGING_DIR --env WORKER_LOG_DIR
        --env PREFETCH_BIN --env DUMP_BIN --env MAX_SIZE --env DUMP_FLAGS_STR
        --env worker_main --env log --env warn --env error
        worker_main {}
        :::: "${JOBS_FILE}"
    )
    if command -v setsid >/dev/null 2>&1; then
        setsid "${PARALLEL_RUN[@]}" &
    else
        "${PARALLEL_RUN[@]}" &
    fi
    PARALLEL_PID=$!

    rc=0
    wait "${PARALLEL_PID}" || rc=$?
    PARALLEL_PID=""
    if [[ "${rc}" -ne 0 ]]; then
        warn "GNU parallel exited with status ${rc} (non-zero simply means some jobs failed)."
    fi
else
    log "Nothing to process: all ${TOTAL} accession(s) already have FASTQ output."
    FAILED_SRRS=()
fi

# ---------------------------------------------------------------------------
# Collect per-accession results from worker status files
# ---------------------------------------------------------------------------
FAILED_SRRS=()
while IFS= read -r srr; do
    [[ -z "${srr}" ]] && continue
    rc_file="${WORKER_LOG_DIR}/${srr}.rc"
    if [[ ! -f "${rc_file}" ]]; then
        error "[FAIL] ${srr}: no exit status recorded (worker killed?)"
        FAILED_SRRS+=("${srr}")
        continue
    fi
    wrc="$(cat "${rc_file}")"
    if [[ "${wrc}" != "0" ]]; then
        error "[FAIL] ${srr} (worker exit code ${wrc})"
        FAILED_SRRS+=("${srr}")
    fi
done < "${JOBS_FILE}"

N_FAILED="${#FAILED_SRRS[@]}"
log "Download phase complete."

# ---------------------------------------------------------------------------
# Record failures
# ---------------------------------------------------------------------------
if [[ "${N_FAILED}" -gt 0 ]]; then
    printf '%s\n' ${FAILED_SRRS[@]+"${FAILED_SRRS[@]}"} > "${FAILED_LOG}"
    warn "Failed accession count: ${N_FAILED} (written to ${FAILED_LOG})"
fi

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
N_OK=$(( N_JOBS - N_FAILED ))
DOWNLOADED="$(find "${OUTDIR}" -maxdepth 1 -name '*.fastq*' -type f | wc -l)"
log "=== Summary ==="
log "  Output directory       : ${OUTDIR}"
log "  Unique accessions      : ${TOTAL}"
log "  Skipped (already done) : ${SKIP_DONE}"
log "  Processed this run     : ${N_JOBS} (ok: ${N_OK}, failed: ${N_FAILED})"
log "  FASTQ files in outdir  : ${DOWNLOADED}"
if [[ "${N_FAILED}" -gt 0 ]]; then
    log "  Failed accessions      : see ${FAILED_LOG} and ${WORKER_LOG_DIR}/<SRR>.log"
    exit 2
fi
log "Done."
exit 0
