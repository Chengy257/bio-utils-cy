#!/bin/bash
#########################################################################
# File Name: extract_species_rna_rfam.sh
# Author: ChengYu
# Description: Extract species-specific RNA sequences from Rfam fasta files
# Created Time: 2026
# Changelog:
#   v1.1.0  2026-09-27
#   - FIX: a zero-match RNA class aborted the whole script (grep rc=1
#     under pipefail); empty classes now yield empty .RFids lists
#   - FIX: "grep -c ... || echo 0" produced a two-line count for empty
#     outputs ("0" then "0"); guarded with || true instead
#   - FIX: outputs were written non-atomically -- a failed seqkit left a
#     truncated file that later runs treated as complete; writes now go to
#     a temp file that is moved into place only on success
#   - FIX: species names were interpolated unescaped into seqkit's regex;
#     metacharacters ("E. coli") are now escaped so they match literally
#   - CLEAN: dropped dead zcat/xargs dependency checks and the redundant
#     second snRNA grep (grep -w snRNA never matches snoRNA anyway);
#     duplicate species in the argument list are processed once, and
#     short-name collisions between species now abort with an error
#   - FIX: seqkit resolved via buc_resolve_bin (BUC_SEQKIT_BIN); all logs
#     go to stderr; -t validated as a positive integer
#########################################################################
set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf
IFS=$'\n\t'

VERSION="1.1.0"

# ======================== Logging ============================
log()  { echo "[INFO] $*" >&2; }
warn() { echo "[WARN] $*" >&2; }
err()  { echo "[ERROR] $*" >&2; }

# ======================== Defaults ========================
INDIR=""
OUTDIR="./output_Rfam_seq"
RNATYPE_FILE=""
RNA_CLASSES=("rRNA" "tRNA" "miRNA" "snoRNA" "snRNA")
THREADS="${BUC_THREADS:-8}"

# ======================== Usage ============================
usage() {
    cat <<EOF
Usage: $(basename "$0") [options] "Species name 1" ["Species name 2" ...]

Extract species-specific RNA sequences from Rfam fasta files.
Species names are matched as substrings of the fasta headers
(regex metacharacters in names are escaped automatically).

Options:
  -i DIR     Input directory containing Rfam .fa.gz files (required)
  -a FILE    RNA type annotation file / family.txt (required)
  -o DIR     Output directory (default: ${OUTDIR})
  -t INT     Number of parallel threads (default: ${THREADS})
  -v         Show version
  -h         Show this help message

Examples:
  # Basic usage
  $(basename "$0") -i /data/rfam/15.0/fasta_files/individual \\
      -a /data/rfam/15.0/family.txt \\
      "Homo sapiens"

  # Multiple species with custom output and threads
  $(basename "$0") -i /data/rfam/15.0/fasta_files/individual \\
      -a /data/rfam/15.0/family.txt \\
      -o ./rfam_output -t 12 \\
      "Oryza sativa Japonica Group" "Zea mays"
EOF
    exit 0
}

# ======================== Functions ============================
generate_short_name() {
    local species="$1"
    echo "${species}" | awk '{
        first = toupper(substr($1, 1, 1));
        if (NF >= 2) {
            second = tolower(substr($2, 1, 3));
        } else {
            second = tolower(substr($1, 2, 3));
        }
        printf "%s%s", first, second;
    }'
}

escape_regex() {
    # Escape regex metacharacters so species names like "E. coli" are
    # matched literally. seqkit grep requires -r (regex mode) for names
    # containing spaces to match at all, so we cannot drop regex mode;
    # escaping keeps the intended plain-substring behaviour.
    printf '%s' "$1" | sed 's/[][\.|$(){}?+*^\\]/\\&/g'
}

# ======================== Parse options ============================
while getopts ":i:o:a:t:vh" opt; do
    case "${opt}" in
        i) INDIR="${OPTARG}" ;;
        o) OUTDIR="${OPTARG}" ;;
        a) RNATYPE_FILE="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        v) echo "extract_species_rna_rfam.sh version ${VERSION}"; exit 0 ;;
        h) usage ;;
        :) err "Option -${OPTARG} requires an argument."; usage ;;
        *) err "Unknown option: -${OPTARG}"; usage ;;
    esac
done
shift $((OPTIND - 1))

# ======================== Validate required args ============================
if [[ -z "${INDIR}" ]]; then
    err "Input directory (-i) is required."
    usage
fi
if [[ -z "${RNATYPE_FILE}" ]]; then
    err "RNA type annotation file (-a) is required."
    usage
fi
if [[ $# -eq 0 ]]; then
    err "No species specified."
    usage
fi
if ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -eq 0 ]]; then
    err "Thread count (-t) must be a positive integer, got '${THREADS}'."
    exit 1
fi

# ======================== Check dependencies ============================
SEQKIT_BIN="$(buc_resolve_bin BUC_SEQKIT_BIN seqkit)" || exit 1

# ======================== Validate paths ============================
if [[ ! -f "${RNATYPE_FILE}" ]]; then
    err "RNA type file not found: ${RNATYPE_FILE}"
    exit 1
fi
if [[ ! -d "${INDIR}" ]]; then
    err "Input directory not found: ${INDIR}"
    exit 1
fi
mkdir -p "${OUTDIR}"

# ======================== Normalize species list ============================
# Process each distinct species once, and fail early if two species would
# map to the same output short name (their results would overwrite each
# other silently).
mapfile -t SPECIES_LIST < <(printf '%s\n' "$@" | awk '!seen[$0]++')
declare -A _short_owner=()
for sp in "${SPECIES_LIST[@]}"; do
    _short="$(generate_short_name "${sp}")"
    if [[ -n "${_short_owner[${_short}]:-}" ]]; then
        err "Species '${sp}' and '${_short_owner[${_short}]}' share the short name '${_short}'; their outputs would be conflated. Aborting."
        exit 1
    fi
    _short_owner["${_short}"]="${sp}"
done
unset -v _short_owner

# ======================== Build RF ID lists ============================
log "Building RF ID lists from ${RNATYPE_FILE} ..."
for class in "${RNA_CLASSES[@]}"; do
    # grep returns 1 on zero matches -- tolerate it (pipefail would
    # otherwise abort the script) and let the empty list flow through.
    { grep -w "${class}" "${RNATYPE_FILE}" || true; } \
        | cut -f1 > "${OUTDIR}/${class}.RFids"
    log "  ${class}: $(wc -l < "${OUTDIR}/${class}.RFids") families"
done

# ======================== Extraction ============================
extract_species_rna() {
    local species="$1"
    local short_name
    short_name=$(generate_short_name "${species}")

    log "Processing species: ${species} (${short_name})"

    local class rfids_file outfa tmpfa fpath rf_id count
    for class in "${RNA_CLASSES[@]}"; do
        rfids_file="${OUTDIR}/${class}.RFids"
        outfa="${OUTDIR}/${short_name}_Rfam_${class}.fa"
        log "  Extracting ${class} -> ${outfa}"

        if [[ -s "${outfa}" ]]; then
            log "  SKIP: ${outfa} already exists."
            continue
        fi

        # Build list of input fasta paths that actually exist
        local input_files=()
        while IFS= read -r rf_id; do
            fpath="${INDIR}/${rf_id}.fa.gz"
            if [[ -f "${fpath}" ]]; then
                input_files+=("${fpath}")
            fi
        done < "${rfids_file}"

        if [[ ${#input_files[@]} -eq 0 ]]; then
            warn "  No fasta files found for ${class}; nothing extracted."
            continue
        fi

        # Atomic write: stage into a temp file, move into place on success
        tmpfa="${outfa}.tmp.$$"
        if "${SEQKIT_BIN}" grep -n -r -p "$(escape_regex "${species}")" -j "${THREADS}" "${input_files[@]}" \
                | "${SEQKIT_BIN}" seq -w 100 > "${tmpfa}"; then
            mv "${tmpfa}" "${outfa}"
        else
            rm -f "${tmpfa}"
            err "  seqkit failed for ${species} / ${class}; aborting."
            return 1
        fi

        count=$(grep -c "^>" "${outfa}" || true)
        log "  DONE: ${outfa} (${count} sequences)"
    done
}

# ======================== Main ============================
log "Input directory : ${INDIR}"
log "RNA type file   : ${RNATYPE_FILE}"
log "Output directory: ${OUTDIR}"
log "Threads         : ${THREADS}"
log "Species         : ${#SPECIES_LIST[@]}"

for sp in "${SPECIES_LIST[@]}"; do
    extract_species_rna "${sp}"
done

log "All done. Results saved to: ${OUTDIR}"
