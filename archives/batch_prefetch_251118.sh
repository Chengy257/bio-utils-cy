#!/bin/bash
#########################################################################
# File Name: batch_prefetch.sh
# Author: ChengYu
# Description: Download fastq files from NCBI SRA (supports single-cell mode)
# Created: 2025-11-18
#########################################################################

set -euo pipefail

jobs_number=8
SC_MODE=0    # 0: bulk; 1: single-cell
FASTQC_MODE=0
DIR=""
SRR_FILE=""

# ------------------------------
# Parse arguments
# ------------------------------
usage() {
    echo "Usage: bash batch_prefetch.sh -d <work_dir> -i <srr_id_file> [--sc] [--fastqc]"
    echo ""
    echo "  -d DIR         Working directory"
    echo "  -i FILE        SRR ID list file"
    echo "  --sc           Enable single-cell mode (use fasterq-dump)"
    echo "  --fastqc       Enable running fastqc"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--dir)
            DIR="$2"
            shift 2
            ;;
        -i|--input)
            SRR_FILE="$2"
            shift 2
            ;;
        --sc|--single-cell)
            SC_MODE=1
            shift
            ;;
        --fastqc)
            FASTQC_MODE=1
            shift
            ;;
        *)
            echo "❌ Unknown option: $1"
            usage
            ;;
    esac
done


if [[ -z "$DIR" || -z "$SRR_FILE" ]]; then
    usage
fi

SRR_FILE="$(readlink -f "$SRR_FILE")"
DIR="$(readlink -f "$DIR")"


mkdir -p "$DIR/1.rawdata"
cd "$DIR/1.rawdata"

# ------------------------------
# Prefetch + dump
# ------------------------------
runPrefetch() {
    echo "[$(date)] Building command list..."
    rm -f prefetch.cmd prefetch.cmd.Failed FailedCommands || true

    while read -r id; do
        [[ -z "$id" ]] && continue

        # ----------------------
        # Choose dump tool
        # ----------------------
        if [[ $SC_MODE -eq 1 ]]; then
            # Single-cell → use fasterq-dump
            dump_cmd="/opt/anaconda3/envs/download/bin/fasterq-dump ./${id}/${id}.sra -O ./ -e 4 --split-files --include-technical "
        else
            # Bulk RNA → fastq-dump
            dump_cmd="/opt/anaconda3/envs/download/bin/fastq-dump --gzip --split-3 --defline-qual '+' --defline-seq '@\$ac-\$si/\$ri' ./${id}/${id}.sra"
        fi

        echo "
/opt/anaconda3/envs/download/bin/prefetch -X 200G -r yes ${id} -O ./ \
&& ${dump_cmd} \
&& rm -f ./${id}/${id}.sra
" >> prefetch.cmd

    done < "$SRR_FILE"

    # Escape $ symbol
    sed -i 's/\$/\\$/g' prefetch.cmd

    echo "[$(date)] Start parallel download..."
    /home/chengyu/soft/miniconda3/bin/ParaFly -c prefetch.cmd -CPU $jobs_number

    # Retry failed
    while [[ -f "FailedCommands" ]]; do
        echo "[$(date)] Retrying failed commands..."
        mv FailedCommands prefetch.cmd.Failed
        /home/chengyu/soft/miniconda3/bin/ParaFly -c prefetch.cmd.Failed -CPU $jobs_number
        rm -f prefetch.cmd.Failed || true
    done

    echo "[$(date)] Download finished."
}

# ------------------------------
# QC Module
# ------------------------------
runFastQC() {
    echo "[$(date)] Running fastqc..."
    mkdir -p fastqc/multiqc
    /home/chengyu/soft/miniconda3/envs/common/bin/fastqc *.fastq.gz -o fastqc/ -t 24
    /home/chengyu/soft/miniconda3/envs/common/bin/multiqc fastqc/ -o fastqc/multiqc/
    echo "[$(date)] QC finished!"
}

# ------------------------------
# Main
# ------------------------------
echo "[$(date)] Program start!"
echo "Working directory: $DIR"
echo "SRR file: $SRR_FILE"
echo "Single-cell mode: $SC_MODE"

runPrefetch

# Move gz to top folder
find "$DIR/1.rawdata" -maxdepth 2 -name "*.fastq.gz" -exec mv {} "$DIR/1.rawdata/" \;

echo "[$(date)] Download ALL completed!"

# ------------------------------
# Run FASTQC if enabled
# ------------------------------
if [[ $FASTQC_MODE -eq 1 ]]; then
    runFastQC
else
    echo "[$(date)] FastQC skipped (enable via --fastqc)"
fi

echo "[$(date)] Download ALL completed!"
