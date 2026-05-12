#!/bin/bash
#########################################################################
# File Name: extract_species_rna_rfam.sh
# Author: ChengYu
# Description: Extract species-specific RNA sequences from Rfam fasta files
# Created Time: 2026
#########################################################################
set -euo pipefail
IFS=$'\n\t'

VERSION="1.0.0"

# ======================== Defaults ========================
INDIR=""
OUTDIR="./output_Rfam_seq"
RNATYPE_FILE=""
RNA_CLASSES=("rRNA" "tRNA" "miRNA" "snoRNA" "snRNA")
THREADS=8

# ======================== Usage ============================
usage() {
    cat <<EOF
Usage: $(basename "$0") [options] "Species name 1" ["Species name 2" ...]

Extract species-specific RNA sequences from Rfam fasta files.

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

# ======================== Parse options ============================
while getopts ":i:o:a:t:vh" opt; do
    case "${opt}" in
        i) INDIR="${OPTARG}" ;;
        o) OUTDIR="${OPTARG}" ;;
        a) RNATYPE_FILE="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        v) echo "extract_species_rna_rfam.sh version ${VERSION}"; exit 0 ;;
        h) usage ;;
        :) echo "[ERROR] Option -${OPTARG} requires an argument." >&2; usage ;;
        *) echo "[ERROR] Unknown option: -${OPTARG}" >&2; usage ;;
    esac
done
shift $((OPTIND - 1))

# ======================== Validate required args ============================
if [[ -z "${INDIR}" ]]; then
    echo "[ERROR] Input directory (-i) is required." >&2
    usage
fi
if [[ -z "${RNATYPE_FILE}" ]]; then
    echo "[ERROR] RNA type annotation file (-a) is required." >&2
    usage
fi
if [[ $# -eq 0 ]]; then
    echo "[ERROR] No species specified." >&2
    usage
fi
SPECIES_LIST=("$@")

# ======================== Check dependencies ============================
for cmd in seqkit zcat grep xargs; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "[ERROR] Command '${cmd}' not found. Please install it first." >&2
        exit 1
    fi
done

# ======================== Validate paths ============================
if [[ ! -f "${RNATYPE_FILE}" ]]; then
    echo "[ERROR] RNA type file not found: ${RNATYPE_FILE}" >&2
    exit 1
fi
if [[ ! -d "${INDIR}" ]]; then
    echo "[ERROR] Input directory not found: ${INDIR}" >&2
    exit 1
fi
mkdir -p "${OUTDIR}"

# ======================== Build RF ID lists ============================
echo "[INFO] Building RF ID lists from ${RNATYPE_FILE} ..."
for class in "${RNA_CLASSES[@]}"; do
    grep -w "${class}" "${RNATYPE_FILE}" | cut -f1 > "${OUTDIR}/${class}.RFids"
done
# snRNA without snoRNA
grep -w snRNA "${RNATYPE_FILE}" | grep -v snoRNA | cut -f1 > "${OUTDIR}/snRNA.RFids"

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

extract_species_rna() {
    local species="$1"
    local short_name
    short_name=$(generate_short_name "${species}")

    echo "[INFO] Processing species: ${species} (${short_name})"

    for class in "${RNA_CLASSES[@]}"; do
        local rfids_file="${OUTDIR}/${class}.RFids"
        local outfa="${OUTDIR}/${short_name}_Rfam_${class}.fa"

        echo "[INFO]   Extracting ${class} -> ${outfa}"

        if [[ -s "${outfa}" ]]; then
            echo "[SKIP]   ${outfa} already exists, skipping."
            continue
        fi

        # Build list of input fasta paths that actually exist
        local input_files=()
        while IFS= read -r rf_id; do
            local fpath="${INDIR}/${rf_id}.fa.gz"
            if [[ -f "${fpath}" ]]; then
                input_files+=("${fpath}")
            fi
        done < "${rfids_file}"

        if [[ ${#input_files[@]} -eq 0 ]]; then
            echo "[WARN]   No fasta files found for ${class}."
            continue
        fi

        seqkit grep -n -r -p "${species}" -j "${THREADS}" "${input_files[@]}" \
            | seqkit seq -w 100 > "${outfa}"

        local count
        count=$(grep -c "^>" "${outfa}" || echo 0)
        echo "[DONE]   ${outfa} (${count} sequences)"
    done
}

# ======================== Main ============================
echo "[INFO] Input directory : ${INDIR}"
echo "[INFO] RNA type file   : ${RNATYPE_FILE}"
echo "[INFO] Output directory: ${OUTDIR}"
echo "[INFO] Threads         : ${THREADS}"
echo "[INFO] Species         : ${#SPECIES_LIST[@]}"
echo ""

for sp in "${SPECIES_LIST[@]}"; do
    extract_species_rna "${sp}"
done

echo ""
echo "[ALL DONE] Results saved to: ${OUTDIR}"
