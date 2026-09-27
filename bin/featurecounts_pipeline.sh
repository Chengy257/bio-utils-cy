#!/bin/bash
#########################################################################
# File Name: featurecounts_pipeline.sh
# Author: ChengYu
# Description: Automated featureCounts pipeline with strandedness
#              detection and PE/SE auto-detection.
# Created Time: 2026
#########################################################################
# Changelog:
#   v1.1.0  2026-09-27
#   - FIX: strandedness detection never fired — the grep pattern
#     "1\.strand-specific" does not exist in real infer_experiment.py
#     output, so every sample silently quantified as unstranded. Parsing
#     now matches the actual output lines:
#       PE: Fraction of reads explained by "1++,1--,2+-,2-+": X   (fwd)
#           Fraction of reads explained by "1+-,1-+,2++,2--": Y   (rev)
#       SE: Fraction of reads explained by "++,--": X              (fwd)
#           Fraction of reads explained by "+-,-+": Y              (rev)
#   - FIX: strandedness mapping was inverted. FWD majority (reads/read1 on
#     the gene's sense strand = secondstrand protocol) now maps to
#     featureCounts "-s 1"; REV majority (firststrand/dUTP protocol) maps
#     to "-s 2" (verified empirically: with read1 on the sense strand of a
#     "+" gene, -s 1 counts the fragment and -s 2 does not; with read1
#     antisense the reverse).
#   - FIX: --samtools/--infer-exp CLI values were silently overwritten by
#     auto-detection; CLI > BUC_*_BIN > PATH is now honored.
#   - FIX: sample name is now derived from --bam-pattern instead of a
#     hardcoded "_Aligned.sortedByCoord.out.bam" suffix.
#   - FIX: PE detection excludes unmapped/secondary/supplementary reads
#     (samtools view -c -f 0x1 -F 0x904).
#   - FIX: infer_experiment.py stderr is captured and reported on failure
#     instead of being discarded (2>/dev/null).
#   - CHANGE: bc replaced with awk for threshold comparison (no extra
#     dependency); strandedness is called when one orientation leads the
#     other by a margin > 0.3 (the old elif order mis-called mixed
#     fractions such as 0.45/0.48 as stranded); --fc-script help text
#     corrected (there is no auto-detect — empty means call featureCounts
#     directly); THREADS validated as numeric; post-run check that the
#     counts file exists.
#   - DOC: version banner on -v.
set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf

VERSION="1.1.0"

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
      --samtools <PATH>     Path to samtools (default: BUC_SAMTOOLS_BIN > PATH)
      --infer-exp <PATH>    Path to infer_experiment.py (default: BUC_INFER_EXP_BIN > PATH)
      --fc-script <PATH>    Post-processing R script run as
                            "Rscript <script> <bam> <gtf> <strand> <threads> <out_prefix>"
                            (empty: call featureCounts directly — there is no
                            auto-detection for this option)
      --bam-pattern <PAT>   BAM filename pattern (default: *_Aligned.sortedByCoord.out.bam)
  -h          Show help
  -v          Show version

Examples:
  $(basename "$0") -b bam/ -g annotation.gtf -r genes.bed -o results/ -t 8
EOF
}

log_info()  { echo "[INFO]  $1" >&2; }
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
            # Accept both "--opt=value" and "--opt value" forms (help shows the space form)
            _lopt="${OPTARG%%=*}"
            if [[ "${OPTARG}" == *=* ]]; then
                _val="${OPTARG#*=}"
            else
                case "${_lopt}" in
                    samtools|infer-exp|fc-script|bam-pattern)
                        _val="${!OPTIND:-}"
                        if [[ -z "${_val}" ]]; then
                            log_error "Option --${_lopt} requires a value"
                        fi
                        OPTIND=$((OPTIND + 1))
                        ;;
                    *) _val="" ;;
                esac
            fi
            case "${_lopt}" in
                samtools)    SAMTOOLS="${_val}" ;;
                infer-exp)   INFER_EXP="${_val}" ;;
                fc-script)   FC_SCRIPT="${_val}" ;;
                bam-pattern) BAM_PATTERN="${_val}" ;;
                *)           log_error "Unknown option: --${_lopt}" ;;
            esac ;;
        *) log_error "Unknown option: -$opt" ;;
    esac
done

[[ -z "${BAM_DIR}" ]] && log_error "Missing required: -b (BAM directory)"
[[ -z "${OUT_DIR}" ]] && log_error "Missing required: -o (output directory)"
[[ -z "${GTF}" ]]    && log_error "Missing required: -g (GTF file)"
[[ -z "${BED}" ]]    && log_error "Missing required: -r (BED file for strandedness)"
[[ -z "${THREADS}" ]] && log_error "Missing required: -t (threads)"
[[ ! "${THREADS}" =~ ^[0-9]+$ ]] && log_error "Threads must be a positive integer: ${THREADS}"
[[ "${THREADS}" -lt 1 ]] && log_error "Threads must be >= 1: ${THREADS}"
[[ ! -d "${BAM_DIR}" ]] && log_error "BAM directory not found: ${BAM_DIR}"
[[ ! -f "${GTF}" ]]     && log_error "GTF file not found: ${GTF}"
[[ ! -f "${BED}" ]]     && log_error "BED file not found: ${BED}"

# Find tools: CLI option > BUC_*_BIN (config/env.sh) > PATH
if [[ -n "${SAMTOOLS}" ]]; then
    [[ -x "${SAMTOOLS}" ]] || log_error "--samtools is not executable: ${SAMTOOLS}"
else
    SAMTOOLS="$(buc_resolve_bin BUC_SAMTOOLS_BIN samtools)" || log_error "samtools not found. Install samtools or set BUC_SAMTOOLS_BIN."
fi
if [[ -n "${INFER_EXP}" ]]; then
    [[ -x "${INFER_EXP}" ]] || log_error "--infer-exp is not executable: ${INFER_EXP}"
else
    INFER_EXP="$(buc_resolve_bin BUC_INFER_EXP_BIN infer_experiment.py)" || log_error "infer_experiment.py not found. Install RSeQC or set BUC_INFER_EXP_BIN."
fi
FEATURECOUNTS="$(buc_resolve_bin BUC_FEATURECOUNTS_BIN featureCounts)" || log_error "featureCounts not found. Install subread or set BUC_FEATURECOUNTS_BIN."
RSCRIPT_BIN="$(buc_resolve_bin BUC_RSCRIPT_BIN Rscript)" || log_error "Rscript not found. Set BUC_RSCRIPT_BIN in config/env.local.sh."
if [[ -n "${FC_SCRIPT}" && ! -f "${FC_SCRIPT}" ]]; then
    log_error "--fc-script not found: ${FC_SCRIPT}"
fi

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

# Sample name = BAM basename minus the literal tail of the glob pattern
# (e.g. pattern "*_Aligned.sortedByCoord.out.bam" strips that suffix;
# "*_star.bam" strips "_star.bam"; pattern without "*" keeps the basename).
SAMPLE_SUFFIX="${BAM_PATTERN##*\*}"

for BAM in "${BAM_FILES[@]}"; do
    SAMPLE="$(basename "${BAM}")"
    if [[ -n "${SAMPLE_SUFFIX}" && "${SAMPLE}" == *"${SAMPLE_SUFFIX}" ]]; then
        SAMPLE="${SAMPLE%"${SAMPLE_SUFFIX}"}"
    fi
    if [[ -z "${SAMPLE}" ]]; then
        SAMPLE="$(basename "${BAM}")"
    fi
    log_info "Processing: ${SAMPLE}"

    # Index if needed
    if [[ ! -f "${BAM}.bai" ]]; then
        log_info "  Indexing..."
        "${SAMTOOLS}" index -@ "${THREADS}" "${BAM}"
    fi

    # Detect strandedness.
    # infer_experiment.py output (RSeQC >= 3), fraction after the colon:
    #   PE: Fraction of reads explained by "1++,1--,2+-,2-+": <fwd>
    #       Fraction of reads explained by "1+-,1-+,2++,2--": <rev>
    #   SE: Fraction of reads explained by "++,--": <fwd>
    #       Fraction of reads explained by "+-,-+": <rev>
    # fwd = reads (SE) / read1 (PE) on the gene's sense strand
    #     -> secondstrand protocol -> featureCounts -s 1
    # rev = reads/read1 antisense -> firststrand (dUTP) -> -s 2
    # The two fractions sum to ~1; the library is called stranded when one
    # orientation leads the other by a margin > 0.3 (i.e. dominant fraction
    # > ~0.65), otherwise unstranded.
    # (mapping verified empirically against featureCounts -s 1/2 on a
    #  synthetic PE BAM; see Changelog v1.1.0)
    log_info "  Detecting strandedness..."
    STRAND_ERR="$(mktemp)"
    STRAND_RESULT="$("${INFER_EXP}" -r "${BED}" -i "${BAM}" 2>"${STRAND_ERR}")" \
        || STRAND_RESULT=""
    STRAND_CODE=0  # unstranded
    if [[ -n "${STRAND_RESULT}" ]]; then
        FWD=$(printf '%s\n' "${STRAND_RESULT}" | awk -F': *' '/1\+\+,1--,2\+-,2-\+|\+\+,--/ {print $2; exit}')
        REV=$(printf '%s\n' "${STRAND_RESULT}" | awk -F': *' '/1\+-,1-\+,2\+\+,2--|\+-,-\+/ {print $2; exit}')
        FWD="${FWD:-0}"; REV="${REV:-0}"
        if awk -v f="${FWD}" -v r="${REV}" 'BEGIN{f+=0; r+=0; exit !((f - r) > 0.3)}'; then
            STRAND_CODE=1  # secondstrand in featureCounts
            log_info "  Strandedness: secondstrand (fwd=${FWD}, rev=${REV})"
        elif awk -v f="${FWD}" -v r="${REV}" 'BEGIN{f+=0; r+=0; exit !((r - f) > 0.3)}'; then
            STRAND_CODE=2  # firststrand in featureCounts
            log_info "  Strandedness: firststrand (fwd=${FWD}, rev=${REV})"
        else
            log_info "  Strandedness: unstranded (fwd=${FWD}, rev=${REV})"
        fi
    else
        log_warn "  Strandedness detection failed (see below), defaulting to unstranded."
        sed 's/^/    infer_experiment.py: /' "${STRAND_ERR}" | head -n 5 >&2 || true
    fi
    rm -f "${STRAND_ERR}"

    # Detect PE/SE (paired primary mapped reads; exclude unmapped 0x4,
    # secondary 0x100, supplementary 0x800)
    SAMVIEW=$("${SAMTOOLS}" view -c -f 0x1 -F 0x904 "${BAM}" 2>/dev/null || echo "0")
    if [[ "${SAMVIEW}" -gt 0 ]]; then
        PE_FLAG="-p"
        log_info "  Paired-end detected."
    else
        PE_FLAG=""
        log_info "  Single-end detected."
    fi

    # Run featureCounts
    log_info "  Running featureCounts..."
    if [[ -n "${FC_SCRIPT}" ]]; then
        "${RSCRIPT_BIN}" "${FC_SCRIPT}" "${BAM}" "${GTF}" "${STRAND_CODE}" "${THREADS}" "${EXPR_DIR}/${SAMPLE}"
        log_info "  Done: ${SAMPLE}"
        continue
    fi
    "${FEATURECOUNTS}" ${PE_FLAG} -s "${STRAND_CODE}" -t exon -g gene_id \
        -a "${GTF}" -o "${EXPR_DIR}/${SAMPLE}.fc.tsv" \
        -T "${THREADS}" "${BAM}"
    if [[ ! -s "${EXPR_DIR}/${SAMPLE}.fc.tsv" ]]; then
        log_error "featureCounts produced no output for ${SAMPLE}: ${EXPR_DIR}/${SAMPLE}.fc.tsv"
    fi
    log_info "  Done: ${SAMPLE}"
done

log_info "All samples processed. Results in: ${EXPR_DIR}"
