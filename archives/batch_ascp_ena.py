#!/usr/bin/env bash

set -euo pipefail

########################################
# Default built-in config
########################################
DEFAULT_ASCP="/opt/anaconda3/envs/download/bin/ascp"
DEFAULT_SSH_KEY="/opt/anaconda3/envs/download/etc/asperaweb_id_dsa.openssh"
DEFAULT_BANDWIDTH="300m"
DEFAULT_THREADS=1

ASCP_BIN=""
SSH_KEY=""
OUTDIR="."
BANDWIDTH="$DEFAULT_BANDWIDTH"
THREADS="$DEFAULT_THREADS"

########################################
# Usage
########################################
usage() {
    echo "=================================================="
    echo "ENA Aspera FASTQ Downloader (Stable Version)"
    echo "=================================================="
    echo "Usage:"
    echo "  $0 -i SRRxxxxxxx"
    echo "  $0 -l accession_list.txt"
    echo ""
    echo "Options:"
    echo "  -i    Single accession"
    echo "  -l    Accession list file"
    echo "  -o    Output directory (default: ./)"
    echo "  -b    Bandwidth (default: 300m)"
    echo "  -t    Parallel threads (default: 1)"
    echo "  -a    ascp binary path (override default)"
    echo "  -k    aspera ssh key path (override default)"
    echo "  -h    Help"
    exit 1
}

########################################
# Parse arguments
########################################
INPUT_ID=""
LIST_FILE=""

while getopts ":i:l:o:b:t:a:k:h" opt; do
    case $opt in
        i) INPUT_ID="$OPTARG" ;;
        l) LIST_FILE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        b) BANDWIDTH="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        a) ASCP_BIN="$OPTARG" ;;
        k) SSH_KEY="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_ID" && -z "$LIST_FILE" ]]; then
    usage
fi

mkdir -p "$OUTDIR"

########################################
# Resolve ascp
########################################
if [[ -z "$ASCP_BIN" ]]; then
    ASCP_BIN="$DEFAULT_ASCP"
fi
if [[ ! -x "$ASCP_BIN" ]]; then
    echo "ERROR: ascp not executable: $ASCP_BIN"
    exit 1
fi

########################################
# Resolve ssh key
########################################
if [[ -z "$SSH_KEY" ]]; then
    SSH_KEY="$DEFAULT_SSH_KEY"
fi
if [[ ! -f "$SSH_KEY" ]]; then
    echo "ERROR: SSH key not found: $SSH_KEY"
    exit 1
fi

########################################
# Get fastq_aspera path (with retry)
########################################
get_fastq_paths() {
    local ACC="$1"
    local RETRY=5
    local COUNT=0
    local PATHS=""

    while [[ $COUNT -lt $RETRY ]]; do
        RESPONSE=$(curl -s --max-time 20 \
            "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=fastq_aspera")

        PATHS=$(echo "$RESPONSE" | awk -F'\t' '
            NR==1 {for(i=1;i<=NF;i++){if($i=="fastq_aspera"){col=i}}}
            NR==2 {print $col}
        ' | tr -d '\r')

        if [[ -n "$PATHS" && "$PATHS" != "fastq_aspera" ]]; then
            echo "$PATHS"
            return 0
        fi

        sleep 1
        COUNT=$((COUNT+1))
    done

    echo ""
}

########################################
# Download one accession
########################################
download_one() {
    local ACC="$1"
    echo ">>> Processing $ACC"

    local PATHS
    PATHS=$(get_fastq_paths "$ACC")

    if [[ -z "$PATHS" ]]; then
        echo "ERROR: Cannot retrieve fastq path for $ACC"
        return 1
    fi

    IFS=';' read -ra FILE_ARRAY <<< "$PATHS"

    for FILE in "${FILE_ARRAY[@]}"; do
        FILE=$(echo "$FILE" | xargs)  # trim spaces

        # ✅ 补全 host
        if [[ "$FILE" == fasp.sra.ebi.ac.uk:* ]]; then
            FILE="era-fasp@${FILE}"
        elif [[ "$FILE" == /vol1/* ]]; then
            FILE="era-fasp@fasp.sra.ebi.ac.uk:${FILE}"
        fi

        local BASENAME
        BASENAME=$(basename "$FILE")

        if [[ -f "$OUTDIR/$BASENAME" ]]; then
            echo "    Skip $BASENAME (exists)"
            continue
        fi

        echo "    Downloading $BASENAME"
		
        echo "$ASCP_BIN" -k1 -QT -l "$BANDWIDTH" -P33001 \
            -i "$SSH_KEY" \
            "$FILE" \
            "$OUTDIR" >> "$OUTDIR"/running_comands.sh 
        
		"$ASCP_BIN" -k1 -QT -l "$BANDWIDTH" -P33001 \
            -i "$SSH_KEY" \
            "$FILE" \
            "$OUTDIR"
    done
}


export -f get_fastq_paths
export -f download_one
export ASCP_BIN SSH_KEY OUTDIR BANDWIDTH

########################################
# Run
########################################
if [[ -n "$INPUT_ID" ]]; then
    download_one "$INPUT_ID"
fi

if [[ -n "$LIST_FILE" ]]; then
    xargs -n 1 -P "$THREADS" -I {} bash -c 'download_one "$@"' _ {} < "$LIST_FILE"
fi

echo "All done."
