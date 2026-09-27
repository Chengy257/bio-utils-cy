#!/usr/bin/env bash
#########################################################################
# File Name: batch_dssp.sh
# Author: ChengYu
# Description: Batch run DSSP secondary structure prediction on PDB files.
# Created Time: 2026
# Changelog:
#   v1.1.0  2026-09-27
#   - FIX: the input PDB files were PERMANENTLY MODIFIED: sed -i inserted
#     HEADER/CRYST1 lines into the user's original files (and in the wrong
#     order, CRYST1 before HEADER). A temporary staging copy with the
#     header lines in PDB-spec order is now used; input files are never
#     touched.
#   - FIX: batch failures were silent -- partial/total failures exited 0;
#     the script now verifies every queued output and exits non-zero when
#     any are missing or empty.
#   - FIX: batch command lines are shell-quoted (%q) so paths with spaces
#     survive all three runners; xargs runs one line per invocation
#     (xargs -d '\n') instead of unsafe bash -c re-quoting.
#   - FIX: the 'xssp' candidate (old DSSP 2.x CLI, incompatible args) is
#     no longer auto-detected; resolution goes through buc_resolve_bin
#     (BUC_MKDSSP_BIN) with mkdssp/dssp v4-compatible names only.
#   - FIX: dropped --calculate-accessibility (a no-op for classic dssp
#     output that older mkdssp 2.x builds reject with a hard error).
#   - FIX: all logs go to stderr; -t validated as a positive integer;
#     uppercase *.PDB files are picked up too.
#########################################################################

set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf

VERSION="1.1.0"

log()  { echo "[INFO] $*" >&2; }
warn() { echo "[WARN] $*" >&2; }
err()  { echo "[ERROR] $*" >&2; }

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] -i <PDB_DIR> -o <OUT_DIR>

Batch run DSSP secondary structure prediction on PDB files.
Input files are never modified; temporary staging copies carry the
HEADER/CRYST1 lines required by some DSSP builds.

Options:
  -i <dir>       Input directory containing PDB files (required)
  -o <dir>       Output directory for DSSP files (required)
  -d <path>      Path to DSSP executable (default: BUC_MKDSSP_BIN, then PATH)
  -t <int>       Number of parallel threads (default: 4, or \$BUC_THREADS)
  -p <path>      Path to parallel runner (parallel/ParaFly, default: auto)
  -v             Print version and exit
  -h             Show this help message

Examples:
  $(basename "$0") -i ./pdb_files -o ./dssp_output
  $(basename "$0") -i ./pdbs -o ./out -d /usr/bin/mkdssp -t 8
  $(basename "$0") -i ./pdbs -o ./out -t 20 -p /usr/bin/parallel
EOF
}

# -----------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------
PDB_DIR=""
OUT_DIR=""
DSSP_BIN=""
NUM_CPU="${BUC_THREADS:-4}"
PARALLEL_RUNNER=""

# -----------------------------------------------------------------------
# Option parsing
# -----------------------------------------------------------------------
while getopts ":i:o:d:t:p:vh" opt; do
    case "${opt}" in
        i) PDB_DIR="${OPTARG}" ;;
        o) OUT_DIR="${OPTARG}" ;;
        d) DSSP_BIN="${OPTARG}" ;;
        t) NUM_CPU="${OPTARG}" ;;
        p) PARALLEL_RUNNER="${OPTARG}" ;;
        v) echo "$(basename "$0") version ${VERSION}"; exit 0 ;;
        h) usage; exit 0 ;;
        :) err "Option -${OPTARG} requires an argument."; usage; exit 1 ;;
        \?) err "Invalid option -${OPTARG}."; usage; exit 1 ;;
    esac
done

# -----------------------------------------------------------------------
# Validate required arguments
# -----------------------------------------------------------------------
if [ -z "${PDB_DIR}" ] || [ -z "${OUT_DIR}" ]; then
    err "Both -i (input directory) and -o (output directory) are required."
    usage
    exit 1
fi

if [ ! -d "${PDB_DIR}" ]; then
    err "Input directory '${PDB_DIR}' does not exist."
    exit 1
fi

if ! [[ "${NUM_CPU}" =~ ^[1-9][0-9]*$ ]]; then
    err "Thread count (-t) must be a positive integer, got '${NUM_CPU}'."
    exit 1
fi

mkdir -p "${OUT_DIR}"

# -----------------------------------------------------------------------
# Resolve tools
# -----------------------------------------------------------------------
if [ -n "${DSSP_BIN}" ]; then
    if [ ! -x "${DSSP_BIN}" ]; then
        err "Provided DSSP path is not executable: ${DSSP_BIN}"
        exit 1
    fi
else
    # v4 CLI names only; 'xssp' (DSSP 2.x) uses incompatible arguments.
    if DSSP_BIN="$(buc_resolve_bin BUC_MKDSSP_BIN mkdssp)"; then
        :
    elif DSSP_BIN="$(buc_resolve_bin BUC_MKDSSP_BIN dssp)"; then
        :
    else
        err "DSSP binary not found. Install mkdssp or provide a path with -d."
        exit 1
    fi
fi
log "Using DSSP: ${DSSP_BIN}"

if [ -z "${PARALLEL_RUNNER}" ]; then
    if [ -n "${BUC_PARALLEL_BIN:-}" ] && [ -x "${BUC_PARALLEL_BIN}" ]; then
        PARALLEL_RUNNER="${BUC_PARALLEL_BIN}"
    elif command -v parallel &>/dev/null; then
        PARALLEL_RUNNER="parallel"
    elif command -v ParaFly &>/dev/null; then
        PARALLEL_RUNNER="ParaFly"
    else
        PARALLEL_RUNNER="xargs"
    fi
fi
log "Using parallel runner: ${PARALLEL_RUNNER}"

# -----------------------------------------------------------------------
# Build job list (staging copies; user inputs are never modified)
# -----------------------------------------------------------------------
BATCH_FILE="$(mktemp "${TMPDIR:-/tmp}/batch_dssp_jobs_XXXXXX.txt")"
TMP_PDBS=()
OUT_FILES=()
cleanup() {
    rm -f "${BATCH_FILE}"
    if [ ${#TMP_PDBS[@]} -gt 0 ]; then
        rm -f "${TMP_PDBS[@]+${TMP_PDBS[@]}}"
    fi
}
trap cleanup EXIT

HEADER_LINE='HEADER    ALPHAFOLD PREDICTION                      01-JAN-25   PRED'
CRYST1_LINE='CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1'

count=0
skipped=0

for PDB in "${PDB_DIR}"/*.pdb "${PDB_DIR}"/*.PDB; do
    [ -e "${PDB}" ] || continue

    BASE_NAME="$(basename "${PDB}")"
    OUT_FILE="${OUT_DIR}/${BASE_NAME}.dssp"

    if [ -s "${OUT_FILE}" ]; then
        log "Skipping existing: ${OUT_FILE}"
        skipped=$((skipped + 1))
        continue
    fi

    # Staging copy with minimal PDB header lines in spec order (HEADER
    # first, then CRYST1), inserted only when missing.
    WORK_PDB="$(mktemp "${TMPDIR:-/tmp}/batch_dssp_XXXXXX.pdb")"
    if grep -q '^CRYST1' "${PDB}" 2>/dev/null; then
        { grep -q '^HEADER' "${PDB}" 2>/dev/null || printf '%s\n' "${HEADER_LINE}"; cat -- "${PDB}"; } > "${WORK_PDB}"
    elif grep -q '^HEADER' "${PDB}" 2>/dev/null; then
        awk -v cryst="${CRYST1_LINE}" '
            { print }
            /^HEADER/ && !done { print cryst; done = 1 }
        ' "${PDB}" > "${WORK_PDB}"
    else
        { printf '%s\n' "${HEADER_LINE}" "${CRYST1_LINE}"; cat -- "${PDB}"; } > "${WORK_PDB}"
    fi
    TMP_PDBS+=("${WORK_PDB}")

    printf -v CMDLINE '%q --output-format dssp %q %q' \
        "${DSSP_BIN}" "${WORK_PDB}" "${OUT_FILE}"
    echo "${CMDLINE}" >> "${BATCH_FILE}"
    OUT_FILES+=("${OUT_FILE}")
    count=$((count + 1))
done

if [ "${count}" -eq 0 ]; then
    log "No new PDB files to process (skipped ${skipped} existing)."
    exit 0
fi

log "Queued ${count} PDB files for DSSP analysis (skipped ${skipped} existing)."

# -----------------------------------------------------------------------
# Execute in parallel
# -----------------------------------------------------------------------
set +e
case "${PARALLEL_RUNNER}" in
    parallel)
        parallel -j "${NUM_CPU}" --bar < "${BATCH_FILE}"
        runner_rc=$?
        ;;
    ParaFly)
        ParaFly -c "${BATCH_FILE}" -CPU "${NUM_CPU}"
        runner_rc=$?
        ;;
    xargs)
        xargs -P "${NUM_CPU}" -d '\n' -I{} bash -c '{}' < "${BATCH_FILE}"
        runner_rc=$?
        ;;
    *)
        err "Unknown parallel runner '${PARALLEL_RUNNER}'."
        exit 1
        ;;
esac
set -e

# -----------------------------------------------------------------------
# Verify every queued output
# -----------------------------------------------------------------------
failed=0
for f in ${OUT_FILES[@]+"${OUT_FILES[@]}"}; do
    if [ ! -s "${f}" ]; then
        err "DSSP output missing or empty: ${f}"
        failed=$((failed + 1))
    fi
done

if [ "${failed}" -gt 0 ]; then
    err "${failed} of ${count} DSSP analyses failed (runner exit code: ${runner_rc})."
    exit 1
fi
if [ "${runner_rc}" -ne 0 ]; then
    warn "Runner exited with status ${runner_rc}, but all ${count} outputs were produced."
fi

log "DSSP batch analysis complete. Results in: ${OUT_DIR}"
