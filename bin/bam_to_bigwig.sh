#!/bin/bash
#########################################################################
# File Name: bam_to_bigwig.sh
# Author: ChengYu
# Description: Convert BAM files to bigWig with configurable
#              normalization methods.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-26
#   - FIX: bash 4.2 + `set -u` crashed with "unbound variable" when
#     --norm None left NORM_ARGS empty ("${NORM_ARGS[@]}" is now
#     expanded with the ${arr[@]+${arr[@]}} idiom).
#   - FIX: -p/--parallel is now implemented (previously a documented
#     no-op): up to N bamCoverage jobs run concurrently; per-job exit
#     codes are collected.
#   - FIX: a list file whose first line is a comment/blank was fed to
#     bamCoverage as if it were a BAM. Input is now classified by gzip
#     magic bytes (BAM/bgzf) vs. list file, and list lines are read
#     with IFS-safe parsing (paths with spaces no longer split).
#   - FIX: a failing conversion no longer aborts the whole batch under
#     `set -e`; failures are counted, partial outputs removed, and the
#     script exits 1 at the end.
#   - FIX: .bai indexes are rebuilt when older than the BAM (stale
#     indexes previously produced misaligned coverage).
#   - FIX: -n None passed no normalization flag to bamCoverage, which
#     then silently applied its default (RPKM); --normalizeUsing None
#     is now passed explicitly.
#   - log_info now goes to stderr; THREADS defaults to BUC_THREADS
#     from config/env.sh; numeric options and RPGC genome size are
#     validated; input validation happens before tool resolution.
#########################################################################
set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf

VERSION="1.1.0"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] -i <BAM_LIST>

Convert BAM files to bigWig using bamCoverage (deepTools).

Options:
  -i, --input       <FILE>  File with BAM paths (one per line; '#' comments
                            and blank lines ignored), or a single BAM
  -o, --outdir      <DIR>   Output directory (default: same as input BAM)
  -n, --norm        <METHOD> Normalization method (default: RPKM)
                            Options: RPKM, CPM, BPM, RPGC, None
  -g, --genome-size <INT>   Effective genome size for RPGC normalization
  -b, --bin-size    <INT>   Bin size in bp (default: 1)
  -t, --threads     <INT>   Threads per bamCoverage run (default: BUC_THREADS or 8)
  -p, --parallel    <INT>   Number of BAM files converted concurrently (default: 1)
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

  # Convert 4 BAMs concurrently, 2 threads each
  $(basename "$0") -i bam_list.txt -p 4 -t 2
EOF
}

log_info()  { echo "[INFO]  $1" >&2; }
log_warn()  { echo "[WARN]  $1" >&2; }
log_error() { echo "[ERROR] $1" >&2; exit 1; }

# Defaults
INPUT=""
OUTDIR=""
NORM="RPKM"
GENOME_SIZE=""
BIN_SIZE=1
THREADS="${BUC_THREADS:-8}"
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
for opt_name in BIN_SIZE THREADS PARALLEL; do
    opt_val="${!opt_name}"
    [[ "${opt_val}" =~ ^[0-9]+$ && "${opt_val}" -ge 1 ]] \
        || log_error "--${opt_name,,} must be a positive integer (got: ${opt_val})"
done

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
        [[ "${GENOME_SIZE}" =~ ^[0-9]+$ ]] || log_error "--genome-size must be a positive integer (got: ${GENOME_SIZE})"
        NORM_ARGS=(--normalizeUsing RPGC --effectiveGenomeSize "${GENOME_SIZE}")
        ;;
    None)
        # Must be explicit: bamCoverage defaults to RPKM when the option
        # is omitted entirely.
        NORM_ARGS=(--normalizeUsing None)
        ;;
esac

# Resolve input: BAM file (bgzf/gzip magic) vs. list of BAM paths.
declare -a BAMS=()
if [[ -f "${INPUT}" ]]; then
    magic="$(head -c 2 "${INPUT}" | od -An -tx1 | tr -d ' \n')"
    if [[ "${magic}" == "1f8b" ]]; then
        BAMS=("${INPUT}")
    else
        while IFS= read -r line || [[ -n "${line}" ]]; do
            line="${line%$'\r'}"
            [[ -z "${line}" || "${line}" == \#* ]] && continue
            BAMS+=("${line}")
        done < "${INPUT}"
    fi
elif [[ "${INPUT}" == *.bam ]]; then
    BAMS=("${INPUT}")
else
    log_error "Input file not found: ${INPUT}"
fi
[[ ${#BAMS[@]} -gt 0 ]] || log_error "No BAM paths found in: ${INPUT}"

# Prepare the job list (bam<TAB>output); validation before tool resolution.
declare -a JOBS=()
total=0
skipped=0
for bam in "${BAMS[@]}"; do
    if [[ ! -f "${bam}" ]]; then
        log_warn "BAM not found: ${bam}, skipping."
        continue
    fi
    total=$((total + 1))
    BASENAME=$(basename "${bam}" .bam)
    if [[ -n "${OUTDIR}" ]]; then
        OUTPUT="${OUTDIR}/${BASENAME}.bw"
    else
        OUTPUT="${bam%.*}.bw"
    fi
    if [[ -f "${OUTPUT}" && "${FORCE}" != "true" ]]; then
        log_warn "Output exists (skipping): ${OUTPUT}. Use -f to overwrite."
        skipped=$((skipped + 1))
        continue
    fi
    JOBS+=("${bam}"$'\t'"${OUTPUT}")
done

if [[ "${total}" -eq 0 ]]; then
    log_error "No valid BAM files to convert."
fi
if [[ -n "${OUTDIR}" ]]; then
    mkdir -p "${OUTDIR}"
fi

# Per-file conversion; returns non-zero on failure (never aborts the batch).
run_one() {
    local bam="$1" output="$2"
    # Index if missing or older than the BAM (stale index misaligns coverage).
    if [[ ! -f "${bam}.bai" || "${bam}" -nt "${bam}.bai" ]]; then
        log_info "Indexing: ${bam}"
        if ! "${SAMTOOLS}" index -@ "${THREADS}" "${bam}"; then
            log_warn "samtools index failed: ${bam}"
            return 1
        fi
    fi
    log_info "Converting: ${bam} -> ${output} (norm=${NORM}, binSize=${BIN_SIZE})"
    if ! "${BAMCOV}" \
            -b "${bam}" \
            -o "${output}" \
            -of bigwig \
            -p "${THREADS}" \
            -bs "${BIN_SIZE}" \
            "${NORM_ARGS[@]+${NORM_ARGS[@]}}"; then
        rm -f "${output}"  # remove partial output
        log_warn "bamCoverage failed: ${bam}"
        return 1
    fi
    return 0
}

# Resolve tools (only needed when there is something to convert).
if [[ ${#JOBS[@]} -gt 0 ]]; then
    if [[ -z "${BAMCOV}" ]]; then
        BAMCOV="$(buc_resolve_bin BUC_BAMCOVERAGE_BIN bamCoverage)" \
            || log_error "bamCoverage not found. Install deepTools or set BUC_BAMCOVERAGE_BIN in config/env.local.sh."
    fi
    SAMTOOLS="$(buc_resolve_bin BUC_SAMTOOLS_BIN samtools)" \
        || log_error "samtools not found. Install samtools or set BUC_SAMTOOLS_BIN in config/env.local.sh."
fi

success=0
failed=0
if [[ "${PARALLEL}" -le 1 ]]; then
    for job in "${JOBS[@]+${JOBS[@]}}"; do
        bam="${job%%$'\t'*}"
        output="${job#*$'\t'}"
        if run_one "${bam}" "${output}"; then
            success=$((success + 1))
        else
            failed=$((failed + 1))
        fi
    done
else
    log_info "Parallel mode: ${PARALLEL} concurrent conversions."
    declare -a PIDS=()
    for job in "${JOBS[@]+${JOBS[@]}}"; do
        bam="${job%%$'\t'*}"
        output="${job#*$'\t'}"
        while [[ "$(jobs -rp | wc -l)" -ge "${PARALLEL}" ]]; do
            sleep 0.2
        done
        run_one "${bam}" "${output}" &
        PIDS+=("$!")
    done
    for pid in ${PIDS[@]+${PIDS[@]}}; do
        if wait "${pid}"; then
            success=$((success + 1))
        else
            failed=$((failed + 1))
        fi
    done
fi

log_info "Completed: ${success}/${total} files converted (failed: ${failed}, skipped: ${skipped})."
[[ "${failed}" -eq 0 ]] || exit 1
