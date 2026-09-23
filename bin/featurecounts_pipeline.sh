#!/bin/bash
#########################################################################
# File Name: featurecounts_pipeline.sh
# Author: ChengYu
# Description: Automated featureCounts pipeline with strandedness
#              detection and PE/SE auto-detection.
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
Usage: $(basename "$0") [OPTIONS] -b <BAM_DIR> -g <GTF> -o <OUT_DIR>

Automated featureCounts pipeline with strandedness and PE/SE detection.

Options:
  -b <DIR>    BAM directory (required)
  -o <DIR>    Output directory (required)
  -g <FILE>   GTF annotation file (required)
  -r <FILE>   BED file for strandedness inference (required)
  -t <INT>    Threads for featureCounts (required)
      --samtools <PATH>     Path to samtools (default: auto-detect)
      --infer-exp <PATH>    Path to infer_experiment.py (default: auto-detect)
      --fc-script <PATH>    Path to run-featurecounts.R (default: auto-detect)
      --bam-pattern <PAT>   BAM filename pattern (default: *_Aligned.sortedByCoord.out.bam)
  -h          Show help
  -v          Show version

Examples:
  $(basename "$0") -b bam/ -g annotation.gtf -r genes.bed -o results/ -t 8
EOF
}

log_info()  { echo "[INFO]  $1"; }
log_warn()  { echo "[WARN]  $1" >&2; }
log_error() { echo "[ERROR] $1" >&2; exit 1; }

# Defaults
BAM_DIR=""
OUT_DIR=""
GTF=""
BED=""
THREADS=""
SAMTOOLS=""
INFER_EXP=""
FC_SCRIPT=""
BAM_PATTERN="*_Aligned.sortedByCoord.out.bam"

while getopts "b:o:g:r:t:hv-:" opt; do
    case $opt in
        b) BAM_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        g) GTF="$OPTARG" ;;
        r) BED="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage; exit 0 ;;
        v) echo "$(basename "$0") ${VERSION}"; exit 0 ;;
        -)
            case "${OPTARG}" in
                samtools=*)    SAMTOOLS="${OPTARG#*=}" ;;
                infer-exp=*)   INFER_EXP="${OPTARG#*=}" ;;
                fc-script=*)   FC_SCRIPT="${OPTARG#*=}" ;;
                bam-pattern=*) BAM_PATTERN="${OPTARG#*=}" ;;
                *)             log_error "Unknown option: --${OPTARG}" ;;
            esac ;;
        *) log_error "Unknown option: -$opt" ;;
    esac
done

[[ -z "${BAM_DIR}" ]] && log_error "Missing required: -b (BAM directory)"
[[ -z "${OUT_DIR}" ]] && log_error "Missing required: -o (output directory)"
[[ -z "${GTF}" ]]    && log_error "Missing required: -g (GTF file)"
[[ -z "${BED}" ]]    && log_error "Missing required: -r (BED file for strandedness)"
[[ -z "${THREADS}" ]] && log_error "Missing required: -t (threads)"
[[ ! -d "${BAM_DIR}" ]] && log_error "BAM directory not found: ${BAM_DIR}"
[[ ! -f "${GTF}" ]]     && log_error "GTF file not found: ${GTF}"
[[ ! -f "${BED}" ]]     && log_error "BED file not found: ${BED}"

# Find tools (MYS_*_BIN from config/env.sh > PATH)
SAMTOOLS="$(mys_resolve_bin MYS_SAMTOOLS_BIN samtools)" || log_error "samtools not found. Install samtools or set MYS_SAMTOOLS_BIN."
INFER_EXP="$(mys_resolve_bin MYS_INFER_EXP_BIN infer_experiment.py)" || log_error "infer_experiment.py not found. Install RSeQC or set MYS_INFER_EXP_BIN."
FEATURECOUNTS="$(mys_resolve_bin MYS_FEATURECOUNTS_BIN featureCounts)" || log_error "featureCounts not found. Install subread or set MYS_FEATURECOUNTS_BIN."
RSCRIPT_BIN="$(mys_resolve_bin MYS_RSCRIPT_BIN Rscript)" || log_error "Rscript not found. Set MYS_RSCRIPT_BIN in config/env.local.sh."

# Create output directory
EXPR_DIR="${OUT_DIR}/expression"
mkdir -p "${EXPR_DIR}"

# Find BAM files
shopt -s nullglob
BAM_FILES=("${BAM_DIR}"/${BAM_PATTERN})
shopt -u nullglob

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    log_error "No BAM files matching '${BAM_PATTERN}' in ${BAM_DIR}"
fi

log_info "Found ${#BAM_FILES[@]} BAM files."

for BAM in "${BAM_FILES[@]}"; do
    SAMPLE=$(basename "${BAM}" | sed 's/_Aligned.sortedByCoord.out.bam//')
    log_info "Processing: ${SAMPLE}"

    # Index if needed
    if [[ ! -f "${BAM}.bai" ]]; then
        log_info "  Indexing..."
        "${SAMTOOLS}" index -@ "${THREADS}" "${BAM}"
    fi

    # Detect strandedness
    log_info "  Detecting strandedness..."
    STRAND_RESULT=$("${INFER_EXP}" -r "${BED}" -i "${BAM}" 2>/dev/null || echo "failed")
    STRAND_CODE=0  # default: unstranded
    if echo "${STRAND_RESULT}" | grep -q "1\.strand-specific"; then
        FWD=$(echo "${STRAND_RESULT}" | grep "1\.strand-specific" | grep -oP '[\d.]+')
        REV=$(echo "${STRAND_RESULT}" | grep "2\.strand-specific" | grep -oP '[\d.]+')
        FWD=${FWD:-0}; REV=${REV:-0}
        if (( $(echo "${FWD} > 0.3" | bc -l) )); then
            STRAND_CODE=2  # firststrand in featureCounts
            log_info "  Strandedness: firststrand (fwd=${FWD}, rev=${REV})"
        elif (( $(echo "${REV} > 0.3" | bc -l) )); then
            STRAND_CODE=1  # secondstrand in featureCounts
            log_info "  Strandedness: secondstrand (fwd=${FWD}, rev=${REV})"
        else
            log_info "  Strandedness: unstranded (fwd=${FWD}, rev=${REV})"
        fi
    else
        log_warn "  Strandedness detection failed, defaulting to unstranded."
    fi

    # Detect PE/SE
    SAMVIEW=$("${SAMTOOLS}" view -c -f 0x1 "${BAM}" 2>/dev/null || echo "0")
    if [[ "${SAMVIEW}" -gt 0 ]]; then
        PE_FLAG="-p"
        log_info "  Paired-end detected."
    else
        PE_FLAG=""
        log_info "  Single-end detected."
    fi

    # Run featureCounts
    log_info "  Running featureCounts..."
    if [[ -n "${FC_SCRIPT}" && -f "${FC_SCRIPT}" ]]; then
        "${RSCRIPT_BIN}" "${FC_SCRIPT}" "${BAM}" "${GTF}" "${STRAND_CODE}" "${THREADS}" "${EXPR_DIR}/${SAMPLE}"
    else
        "${FEATURECOUNTS}" ${PE_FLAG} -s "${STRAND_CODE}" -t exon -g gene_id \
            -a "${GTF}" -o "${EXPR_DIR}/${SAMPLE}.fc.tsv" \
            -T "${THREADS}" "${BAM}"
    fi
    log_info "  Done: ${SAMPLE}"
done

log_info "All samples processed. Results in: ${EXPR_DIR}"
