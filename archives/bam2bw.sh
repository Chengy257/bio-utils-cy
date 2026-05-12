#!/bin/bash
#########################################################################
# File Name: /home/chengyu/myscripts/bam2bw.sh
# Author: ChengYu
# Description: 
# Created Time: Fri 15 Sep 2023 10:43:14 AM CST
#########################################################################

threads=24

function bam2bw() {
	if [ ! -f ${1}.bai ];then
		samtools index -@ $threads $id  
	fi
	/home/chengyu/soft/miniconda3/envs/deeptools/bin/bamCoverage \
		-b ${1} -o ${id%.*}.bw -of bigwig -p $threads -bs 1 \
		--normalizeUsing RPKM --effectiveGenomeSize 373128865
}

## main 
cat $1|while read id;
do
	bam2bw $id
done

