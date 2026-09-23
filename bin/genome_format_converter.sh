#!/bin/bash
#########################################################################
# File Name: genome_format_converter.sh
# Author: ChengYu
# Description: Convert between genome annotation formats (GFF3/GTF/BED)
#              using UCSC tools.
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
Usage: $(basename "$0") [OPTIONS] -i <INPUT> -m <MODE>

Convert between genome annotation formats using UCSC tools.

Modes:
  gff2gtf    Convert GFF3 to GTF
  gtf2bed    Convert GTF to BED
  gff2bed    Convert GFF3 to BED (via GTF)
  gtf2gp     Convert GTF to genePred
  gp2bed     Convert genePred to BED

Options:
  -i, --input   <FILE>    Input file (required)
  -o, --output  <FILE>    Output file (default: auto-generated from input name)
  -m, --mode    <MODE>    Conversion mode (required)
  -t, --tmpdir  <DIR>     Temporary directory (default: /tmp)
      --gff3togenepred  <PATH>  Path to gff3ToGenePred (default: auto-detect)
      --genepredtogtf   <PATH>  Path to genePredToGtf (default: auto-detect)
      --genepredtobed   <PATH>  Path to genePredToBed (default: auto-detect)
      --gtftogenepred   <PATH>  Path to gtfToGenePred (default: auto-detect)
  -h, --help              Show this help message
  -v, --version           Show version

Examples:
  $(basename "$0") -i annotation.gff3 -m gff2gtf
  $(basename "$0") -i annotation.gtf -m gtf2bed -o annotation.bed
  $(basename "$0") -i annotation.gff3 -m gff2bed
EOF
}

# Logging
log_info()  { echo "[INFO]  $1"; }
log_warn()  { echo "[WARN]  $1" >&2; }
log_error() { echo "[ERROR] $1" >&2; exit 1; }

# Find UCSC tool
find_tool() {
    local tool_name="$1"
    local explicit="$2"
    if [[ -n "${explicit}" && -x "${explicit}" ]]; then
        echo "${explicit}"
        return
    fi
    # Config-provided path (BUC_<TOOL>_BIN from config/env.sh) beats PATH
    local slot_var="BUC_$(printf '%s' "${tool_name}" | tr 'a-z' 'A-Z' | tr '-' '_')_BIN"
    local slot_val="${!slot_var:-}"
    if [[ -n "${slot_val}" && -x "${slot_val}" ]]; then
        echo "${slot_val}"
        return
    fi
    local found
    found=$(command -v "${tool_name}" 2>/dev/null) || true
    if [[ -n "${found}" ]]; then
        echo "${found}"
        return
    fi
    log_error "Tool not found: ${tool_name}. Install UCSC tools or specify path."
}

# Variables
INPUT=""
OUTPUT=""
MODE=""
TMPDIR="/tmp"
TOOL_GFF3TOGP=""
TOOL_GP2GTF=""
TOOL_GP2BED=""
TOOL_GTF2GP=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)   INPUT="$2"; shift 2 ;;
        -o|--output)  OUTPUT="$2"; shift 2 ;;
        -m|--mode)    MODE="$2"; shift 2 ;;
        -t|--tmpdir)  TMPDIR="$2"; shift 2 ;;
        --gff3togenepred) TOOL_GFF3TOGP="$2"; shift 2 ;;
        --genepredtogtf)  TOOL_GP2GTF="$2"; shift 2 ;;
        --genepredtobed)  TOOL_GP2BED="$2"; shift 2 ;;
        --gtftogenepred)  TOOL_GTF2GP="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        -v|--version) echo "$(basename "$0") ${VERSION}"; exit 0 ;;
        *) log_error "Unknown option: $1" ;;
    esac
done

# Validate
[[ -z "${INPUT}" ]] && log_error "Missing required option: -i/--input"
[[ -z "${MODE}" ]]  && log_error "Missing required option: -m/--mode"
[[ ! -f "${INPUT}" ]] && log_error "Input file not found: ${INPUT}"

VALID_MODES="gff2gtf gtf2bed gff2bed gtf2gp gp2bed"
echo "${VALID_MODES}" | grep -qw "${MODE}" || log_error "Invalid mode: ${MODE}. Valid modes: ${VALID_MODES}"

# Auto-generate output name
if [[ -z "${OUTPUT}" ]]; then
    BASENAME="${INPUT%%.*}"
    case "${MODE}" in
        gff2gtf) OUTPUT="${BASENAME}.gtf" ;;
        gtf2bed|gff2bed) OUTPUT="${BASENAME}.bed" ;;
        gtf2gp) OUTPUT="${BASENAME}.gp" ;;
        gp2bed) OUTPUT="${BASENAME}.bed" ;;
    esac
    log_info "Output file: ${OUTPUT}"
fi

# Resolve tools
GFF3TOGP=$(find_tool "gff3ToGenePred" "${TOOL_GFF3TOGP}")
GP2GTF=$(find_tool "genePredToGtf" "${TOOL_GP2GTF}")
GP2BED=$(find_tool "genePredToBed" "${TOOL_GP2BED}")
GTF2GP=$(find_tool "gtfToGenePred" "${TOOL_GTF2GP}")

# Execute conversion
log_info "Converting: ${MODE}"
case "${MODE}" in
    gff2gtf)
        "${GFF3TOGP}" "${INPUT}" stdout | "${GP2GTF}" file stdin "${OUTPUT}"
        ;;
    gtf2bed)
        "${GTF2GP}" "${INPUT}" stdout | "${GP2BED}" stdin "${OUTPUT}"
        ;;
    gff2bed)
        "${GFF3TOGP}" "${INPUT}" stdout | "${GP2BED}" stdin "${OUTPUT}"
        ;;
    gtf2gp)
        "${GTF2GP}" "${INPUT}" "${OUTPUT}"
        ;;
    gp2bed)
        "${GP2BED}" "${INPUT}" "${OUTPUT}"
        ;;
esac

log_info "Done: ${OUTPUT}"
