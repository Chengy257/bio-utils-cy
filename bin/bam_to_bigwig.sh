#!/bin/bash
#########################################################################
# File Name: bam_to_bigwig.sh
# Author: ChengYu
# Description: Convert BAM files to bigWig with configurable
#              normalization methods.
# Created Time: 2026
#########################################################################
set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf

VERSION="1.0.0"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] -i <BAM_LIST>

Convert BAM files to bigWig using bamCoverage (deepTools).

Options:
  -i, --input       <FILE>  File with BAM paths (one per line), or single BAM
  -o, --outdir      <DIR>   Output directory (default: same as input BAM)
  -n, --norm        <METHOD> Normalization method (default: RPKM)
                            Options: RPKM, CPM, BPM, RPGC, None
  -g, --genome-size <INT>   Effective genome size for RPGC normalization
                            (default: auto-detect if RPGC, otherwise not needed)
  -b, --bin-size    <INT>   Bin size in bp (default: 1)
  -t, --threads     <INT>   Threads per bamCoverage run (default: 8)
  -p, --parallel    <INT>   Number of parallel BAM files (default: 1)
      --bamcoverage <PATH>  Path to bamCoverage (default: auto-detect)
      --samtools    <PATH>  Path to samtools (default: auto-detect)
  -f, --force               Overwrite existing output files
  -h, --help                Show this help
  -v, --version             Show version

Normalization methods:
  RPKM  - Reads Per Kilobase per Million mapped reads
  CPM   - Counts Per Million mapped reads
  BPM   - Bins Per Million (similar to CPM but bin-based)
  RPGC  - Reads Per Genomic Content (requires --genome-size)
  None  - No normalization (raw counts)

Examples:
  # Basic RPKM normalization
  $(basename "$0") -i bam_list.txt

  # CPM normalization with custom bin size
  $(basename "$0") -i bam_list.txt -n CPM -b 10

  # Single BAM file, RPGC normalization
  $(basename "$0") -i sample.bam -n RPGC -g 373128865
EOF
}

log_info()  { echo "[INFO]  $1"; }
log_warn()  { echo "[WARN]  $1" >&2; }
log_error() { echo "[ERROR] $1" >&2; exit 1; }

# Defaults
INPUT=""
OUTDIR=""
NORM="RPKM"
GENOME_SIZE=""
BIN_SIZE=1
THREADS=8
PARALLEL=1
BAMCOV=""
SAMTOOLS=""
FORCE=false

# Parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)       INPUT="$2"; shift 2 ;;
        -o|--outdir)      OUTDIR="$2"; shift 2 ;;
        -n|--norm)        NORM="$2"; shift 2 ;;
        -g|--genome-size) GENOME_SIZE="$2"; shift 2 ;;
        -b|--bin-size)    BIN_SIZE="$2"; shift 2 ;;
        -t|--threads)     THREADS="$2"; shift 2 ;;
        -p|--parallel)    PARALLEL="$2"; shift 2 ;;
        --bamcoverage)    BAMCOV="$2"; shift 2 ;;
        --samtools)       SAMTOOLS="$2"; shift 2 ;;
        -f|--force)       FORCE=true; shift ;;
        -h|--help)        usage; exit 0 ;;
        -v|--version)     echo "$(basename "$0") ${VERSION}"; exit 0 ;;
        *)                log_error "Unknown option: $1" ;;
    esac
done

[[ -z "${INPUT}" ]] && log_error "Missing required option: -i/--input"

# Find tools
BAMCOV=$(command -v ${BAMCOV:-bamCoverage} 2>/dev/null) || log_error "bamCoverage not found. Install deepTools."
SAMTOOLS=$(command -v ${SAMTOOLS:-samtools} 2>/dev/null) || log_error "samtools not found."

# Validate normalization
VALID_NORMS="RPKM CPM BPM RPGC None"
echo "${VALID_NORMS}" | grep -qw "${NORM}" || log_error "Invalid normalization: ${NORM}. Valid: ${VALID_NORMS}"

# Build normalization args
NORM_ARGS=()
case "${NORM}" in
    RPKM|CPM|BPM)
        NORM_ARGS=(--normalizeUsing "${NORM}")
        ;;
    RPGC)
        [[ -z "${GENOME_SIZE}" ]] && log_error "RPGC normalization requires --genome-size"
        NORM_ARGS=(--normalizeUsing RPGC --effectiveGenomeSize "${GENOME_SIZE}")
        ;;
    None)
        NORM_ARGS=()
        ;;
esac

# Resolve BAM list
if [[ -f "${INPUT}" && "$(head -1 "${INPUT}" | grep -c '\.bam')" -gt 0 ]]; then
    BAM_LIST=$(cat "${INPUT}")
else
    BAM_LIST="${INPUT}"
fi

# Process each BAM
total=0
success=0
for bam in ${BAM_LIST}; do
    total=$((total + 1))
    [[ ! -f "${bam}" ]] && log_warn "BAM not found: ${bam}, skipping." && continue

    # Determine output path
    BASENAME=$(basename "${bam}" .bam)
    if [[ -n "${OUTDIR}" ]]; then
        mkdir -p "${OUTDIR}"
        OUTPUT="${OUTDIR}/${BASENAME}.bw"
    else
        OUTPUT="${bam%.*}.bw"
    fi

    # Check overwrite
    if [[ -f "${OUTPUT}" && "${FORCE}" != "true" ]]; then
        log_warn "Output exists (skipping): ${OUTPUT}. Use -f to overwrite."
        continue
    fi

    # Index if needed
    if [[ ! -f "${bam}.bai" ]]; then
        log_info "Indexing: ${bam}"
        "${SAMTOOLS}" index -@ "${THREADS}" "${bam}"
    fi

    # Run bamCoverage
    log_info "Converting: ${bam} -> ${OUTPUT} (norm=${NORM}, binSize=${BIN_SIZE})"
    "${BAMCOV}" \
        -b "${bam}" \
        -o "${OUTPUT}" \
        -of bigwig \
        -p "${THREADS}" \
        -bs "${BIN_SIZE}" \
        "${NORM_ARGS[@]}"

    success=$((success + 1))
done

log_info "Completed: ${success}/${total} files converted."
