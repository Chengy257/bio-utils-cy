#!/bin/bash
#########################################################################
# File Name: getSpeciesRNA_fromRfam_251113.sh
# Author: ChengYu
# Description: Extract species-specific RNA sequences from Rfam fasta files
# Created Time: Thu 13 Nov 2025 04:07:30 PM CST
#########################################################################

set -euo pipefail
IFS=$'\n\t'

# ======================== Default ========================
INDIR="/home/chengyu/data/database/rfam/15.0/fasta_files/individual"
OUTDIR="./output_Rfam_seq" 
RNATYPE_FILE="/home/chengyu/data/database/rfam/15.0/family.txt"
RNA_CLASSES=("rRNA" "tRNA" "miRNA" "snoRNA" "snRNA")
THREADS=8

# ======================== Help ============================
usage() {
    echo "
Usage: $(basename "$0") [options] \"Species name 1\" [\"Species name 2\" ...]

Extract species-specific RNA sequences from Rfam fasta files.

Options:
  -i DIR     Input directory containing Rfam .fa.gz files  [default: ${INDIR}]
  -a FILE    RNA type annotation file (anno_RNAtype_info)  [default: ${RNATYPE_FILE}]
  -o DIR     Output directory                              [default: ${OUTDIR}]
  -t INT     Number of parallel threads                    [default: ${THREADS}]
  -h         Show this help message

Examples:
  bash getSpeciesRNA.sh \"Homo sapiens\"
  bash getSpeciesRNA.sh -t 12 -o ./out \"Oryza sativa Japonica Group\" \"Zea mays\"
"
    exit 0
}

# ======================== parameter ============================
while getopts ":i:o:a:t:h" opt; do
    case $opt in
        i) INDIR="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        a) RNATYPE_FILE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

if [[ $# -eq 0 ]]; then
    echo "[ERROR] No species specified."
    usage
fi
SPECIES_LIST=("$@")

# ======================== commands ============================
for cmd in seqkit zcat grep xargs; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "[ERROR] Command '$cmd' not found. Please install it first." >&2
        exit 1
    }
done

[[ -f "${RNATYPE_FILE}" ]] || { echo "[ERROR] RNA type file not found: ${RNATYPE_FILE}"; exit 1; }
[[ -d "${INDIR}" ]] || { echo "[ERROR] Input directory not found: ${INDIR}"; exit 1; }
mkdir -p "${OUTDIR}"

# ======================== RF id list ============================
echo "[INFO] Generating RF id lists..."
for class in "${RNA_CLASSES[@]}"; do
    grep -w "${class}" "${RNATYPE_FILE}" | cut -f1 > "${class}.RFids"
done
grep -w snRNA "${RNATYPE_FILE}" | grep -v snoRNA | cut -f1 > snRNA.RFids

# ======================== Function ============================
extract_species_rna() {
    local species="$1"
    local short_name
    # short_name=$(echo "$species" | awk '{print tolower(substr($1,1,3))}')
    short_name=$(echo "$species" | awk '{
        
        first = toupper(substr($1,1,1));
        if (NF >= 2) {
            
            second = tolower(substr($2,1,3));
        } else {
            
            second = tolower(substr($1,2,3));
        }
        
        printf "%s%s", first, second;
    }')
    
    echo "[INFO] Processing species: ${species} (${short_name})"

    for class in "${RNA_CLASSES[@]}"; do
        local outfa="${OUTDIR}/${short_name}_Rfam_${class}.fa"
        echo "[INFO]   Extracting ${class} → ${outfa}"

        if [[ -s "${outfa}" ]]; then
            echo "[SKIP]   ${outfa} already exists."
            continue
        fi

        # cat "${class}.RFids" | \
        #     xargs -P "${THREADS}" -I{} bash -c '
        #         f="'${INDIR}'/{}.fa.gz"
        #         if [[ -f "$f" ]]; then
        #             zcat "$f" | seqkit seq -w 0
        #         fi
        #     ' | \
        #     grep -A 1 -F -w "${species}" | grep -v "^--$" | seqkit seq -w 80 > "${outfa}"

        # local count
        # count=$(grep -c "^>" "${outfa}" || echo 0)
        # echo "[DONE]   ${outfa} (${count} sequences)"

        seqkit grep -n -r -p "${species}" -j "${THREADS}" \
            $(sed "s|^|${INDIR}/|; s|$|.fa.gz|" "${class}.RFids") \
            | seqkit seq -w 100 > "${outfa}"
        local count
        count=$(grep -c "^>" "${outfa}" || echo 0)
        echo "[DONE]   ${outfa} (${count} sequences)"

    done
}

# ======================== MAIN ============================
for sp in "${SPECIES_LIST[@]}"; do
    extract_species_rna "${sp}"
done

echo "[ALL DONE] Results saved to: ${OUTDIR}"


