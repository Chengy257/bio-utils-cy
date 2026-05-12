#!/bin/bash
#########################################################################
# File Name: /home/chengyu/myscript/dl_from_NCBI.sh
# Author: ChengYu
# Description: Download fastq files from NCBI SRA with prefetch, parallelized using ParaFly.
# Created Time: Wed Mar  1 21:07:41 2023

#########################################################################
jobs_number=10

function runPrefech(){
    mkdir -p ${DIR}/1.rawdata
    cd ${DIR}/1.rawdata
    if [ -f prefech.command ];then
        rm -f prefech.command* FailedCommands
    fi
    cat ${1} |while read id; 
    do 
        fd="/opt/anaconda3/envs/download/bin/fastq-dump --gzip --split-3 --defline-qual '+' --defline-seq '@\$ac-\$si/\$ri'"; 
        echo "/opt/anaconda3/envs/download/bin/prefetch -X 200G -r yes $id -O ./  \
        && $fd ./${id}/${id}.sra && rm ./${id}/${id}.sra" >> ${DIR}/1.rawdata/prefech.command
    done
    sed -i 's/\$/\\$/g' ${DIR}/1.rawdata/prefech.command
    echo `date` "Start downloading..."
    /home/chengyu/soft/miniconda3/bin/ParaFly -c prefech.command -CPU $jobs_number

    while [ -f "./FailedCommands" ];
    do
        cat ./FailedCommands > prefech.command.Failed && rm ./FailedCommands
        /home/chengyu/soft/miniconda3/bin/ParaFly -c prefech.command.Failed -CPU $jobs_number
    done
    echo `date` "Download Finished!"
    echo 
    echo "Failed files:" `cat ${DIR}/1.rawdata/FailedCommands`
}
function runFastqc(){
    echo `date` "Start running fastqc..."
    cd ${DIR}/1.rawdata
    mkdir -p ./fastqc/multiqc
    /home/chengyu/soft/miniconda3/envs/common/bin/fastqc -f fastq -t 24 ./*.fastq.gz -o ./fastqc/ 
    /home/chengyu/soft/miniconda3/envs/common/bin/multiqc ./fastqc/ -o ./fastqc/multiqc/
    echo `date` "Finished QC!"
}
######################################################################### 
function dl_main(){
    echo `date` "Program start!"
    DIR=$1
    cd $DIR
    mkdir -p ${DIR}/1.rawdata
    runPrefech $2
	# runFastqc
	mv ${DIR}/1.rawdata/*/*gz ${DIR}/1.rawdata/
	echo `date` "Program finished!"
}
## call main function

dl_main $1 $2
# $1 :  Work Dir
# $2 :  SRR IDs file, one ID per line
# Usage: bash /home/chengyu/myscript/dl_from_NCBI.sh [DIR] [SRR ID file]
