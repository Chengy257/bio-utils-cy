#!/usr/bin/env bash
#########################################################################
# File Name: batch_dssp.sh
# Author: ChengYu
# Description: Batch run DSSP secondary structure prediction on PDB files.
# Created Time: 2026
#########################################################################

set -euo pipefail

VERSION="1.0.0"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] -i <PDB_DIR> -o <OUT_DIR>

Batch run DSSP secondary structure prediction on PDB files.

Options:
  -i <dir>       Input directory containing PDB files (required)
  -o <dir>       Output directory for DSSP files (required)
  -d <path>      Path to DSSP executable (default: auto-detect)
  -t <int>       Number of parallel threads (default: 4)
  -p <path>      Path to parallel runner (default: auto-detect)
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
NUM_CPU=4
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
        :) echo "Error: Option -${OPTARG} requires an argument." >&2; usage; exit 1 ;;
        \?) echo "Error: Invalid option -${OPTARG}." >&2; usage; exit 1 ;;
    esac
done

# -----------------------------------------------------------------------
# Validate required arguments
# -----------------------------------------------------------------------
if [ -z "${PDB_DIR}" ] || [ -z "${OUT_DIR}" ]; then
    echo "Error: Both -i (input directory) and -o (output directory) are required." >&2
    usage
    exit 1
fi

if [ ! -d "${PDB_DIR}" ]; then
    echo "Error: Input directory '${PDB_DIR}' does not exist." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

# -----------------------------------------------------------------------
# Auto-detect DSSP binary
# -----------------------------------------------------------------------
if [ -z "${DSSP_BIN}" ]; then
    for candidate in mkdssp dssp xssp; do
        if command -v "${candidate}" &>/dev/null; then
            DSSP_BIN="$(command -v "${candidate}")"
            echo "[INFO] Auto-detected DSSP: ${DSSP_BIN}"
            break
        fi
    done
fi

if [ -z "${DSSP_BIN}" ] || [ ! -x "${DSSP_BIN}" ]; then
    echo "Error: DSSP binary not found. Install DSSP or provide path with -d." >&2
    exit 1
fi

# -----------------------------------------------------------------------
# Auto-detect parallel runner
# -----------------------------------------------------------------------
if [ -z "${PARALLEL_RUNNER}" ]; then
    if command -v parallel &>/dev/null; then
        PARALLEL_RUNNER="parallel"
    elif command -v ParaFly &>/dev/null; then
        PARALLEL_RUNNER="ParaFly"
    else
        PARALLEL_RUNNER="xargs"
    fi
    echo "[INFO] Using parallel runner: ${PARALLEL_RUNNER}"
fi

# -----------------------------------------------------------------------
# Build job list
# -----------------------------------------------------------------------
BATCH_FILE="$(mktemp)"
trap 'rm -f "${BATCH_FILE}"' EXIT

count=0
skipped=0

for PDB in "${PDB_DIR}"/*.pdb; do
    [ -e "${PDB}" ] || continue

    BASE_NAME="$(basename "${PDB}")"
    OUT_FILE="${OUT_DIR}/${BASE_NAME}.dssp"

    # Skip existing output
    if [ -s "${OUT_FILE}" ]; then
        echo "[INFO] Skipping existing: ${OUT_FILE}"
        skipped=$((skipped + 1))
        continue
    fi

    # Insert minimal PDB header if missing (required by some DSSP versions)
    if ! grep -q "^HEADER" "${PDB}" 2>/dev/null; then
        sed -i '1i\HEADER    ALPHAFOLD PREDICTION                      01-JAN-25   PRED' "${PDB}"
    fi
    if ! grep -q "^CRYST1" "${PDB}" 2>/dev/null; then
        sed -i '1i\CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1' "${PDB}"
    fi

    echo "${DSSP_BIN} --calculate-accessibility --output-format dssp ${PDB} ${OUT_FILE}" >> "${BATCH_FILE}"
    count=$((count + 1))
done

if [ "${count}" -eq 0 ]; then
    echo "[INFO] No new PDB files to process (skipped ${skipped} existing)."
    exit 0
fi

echo "[INFO] Queued ${count} PDB files for DSSP analysis (skipped ${skipped} existing)."

# -----------------------------------------------------------------------
# Execute in parallel
# -----------------------------------------------------------------------
case "${PARALLEL_RUNNER}" in
    parallel)
        parallel -j "${NUM_CPU}" --bar < "${BATCH_FILE}"
        ;;
    ParaFly)
        ParaFly -c "${BATCH_FILE}" -CPU "${NUM_CPU}"
        ;;
    xargs)
        xargs -P "${NUM_CPU}" -I{} bash -c '{}' < "${BATCH_FILE}"
        ;;
    *)
        echo "Error: Unknown parallel runner '${PARALLEL_RUNNER}'." >&2
        exit 1
        ;;
esac

echo "[INFO] DSSP batch analysis complete. Results in: ${OUT_DIR}"
