#!/usr/bin/env bash
#########################################################################
# File Name: /home/chengyu/myscripts/run_featurecount_pipeline_251114.sh
# Author: ChengYu
# Description: 
# Created Time: Fri 14 Nov 2025 10:35:48 AM CST
#########################################################################

set -euo pipefail

usage() {
    echo "
Usage: $0 -b BAM_DIR -o OUT_DIR -r BED -g GTF -t THREADS

Required:
  -b   Path to BAM folder (must contain *_Aligned.sortedByCoord.out.bam)
  -o   Output directory
  -r   BED file for infer_experiment.py
  -g   GTF annotation
  -t   Threads for featureCounts
"
    exit 1
}

while getopts "b:o:r:g:t:" opt; do
    case $opt in
        b) BAM_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        r) BED="$OPTARG" ;;
        g) GTF="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${BAM_DIR:-}" ]] && usage
[[ -z "${OUT_DIR:-}" ]] && usage
[[ -z "${BED:-}" ]] && usage
[[ -z "${GTF:-}" ]] && usage
[[ -z "${THREADS:-}" ]] && usage

if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory not found: $BAM_DIR"
    exit 1
fi
if [[ ! -f "$BED" ]]; then
    echo "ERROR: BED file not found: $BED"
    exit 1
fi
if [[ ! -f "$GTF" ]]; then
    echo "ERROR: GTF file not found: $GTF"
    exit 1
fi

mkdir -p "$OUT_DIR"/expression

echo ">>> Start RNA-seq strandedness + featureCounts pipeline"

for bam in "$BAM_DIR"/*_Aligned.sortedByCoord.out.bam; do
    [[ -e "$bam" ]] || { echo "No BAM files found."; exit 1; }
    sample=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    echo -e "\n### Processing sample: $sample"

    # 0) bam index 
    bai_file="${bam}.bai"
    if [[ ! -f "$bai_file" ]]; then
        echo " - BAI index not found for $bam, generating..."
        /home/chengyu/soft/miniconda3/envs/rna-seq/bin/samtools index -@ 12 "$bam"
    else
        echo " - BAI index exists."
    fi

    # 1) Strandedness detection
    infer_out="$BAM_DIR/${sample}_infer_experiment.out"
    stranded_file="$BAM_DIR/${sample}.strandedness"

    if [[ ! -f "$stranded_file" ]]; then
        
        if [[ ! -f "$infer_out" ]]; then
            echo " - Detecting strandedness ..."
            /home/chengyu/soft/miniconda3/envs/rna-seq/bin/infer_experiment.py -r "$BED" -i "$bam" > "$infer_out"
        fi

        grep -v '^\s*$' "$infer_out" | tail -2 | awk '{
            a[NR] = $NF
        }
        END {
            total = a[1] + a[2]
            if (total < 0.5) {
                print "unstrand"
            } else {
                diff = (a[1] - a[2]) / total
                if (diff > 0.3) print "secondstrand"
                else if (diff < -0.3) print "firststrand"
                else print "unstrand"
            }
        }' > "$stranded_file"

    else
        echo " - Strandedness exists, skip."
    fi

    strandedness=$(cat "$stranded_file")
    echo "   Strand type: $strandedness"

    # 2) Convert strandedness to featureCounts strand code
    if [[ "$strandedness" == "firststrand" ]]; then
        strand=2
    elif [[ "$strandedness" == "secondstrand" ]]; then
        strand=1
    else
        strand=0
    fi
    echo "   featureCounts strand code = $strand"

    # 3) Determine PE / SE from infer_experiment.out (skip empty lines)
    first_nonempty_line=$(grep -v '^\s*$' "$infer_out" | head -1)
    if [[ "$first_nonempty_line" == *"PairEnd"* ]]; then
        isPairedEnd=True
    else
        isPairedEnd=False
    fi
    echo "   Paired-end: $isPairedEnd"

    # 4) Run featureCounts
    echo " - Running featureCounts ..."
    Rscript /share/workflows/rna-seq/scripts/run-featurecounts.R \
        -t "$THREADS" \
        -b "$bam" \
        -g "$GTF" \
        -s "$strand" \
        -i "$isPairedEnd" \
        -o "$OUT_DIR/expression/${sample}"

    echo "   Done sample: $sample"

done

echo -e "\n>>> All samples processed successfully!"
