#!/bin/bash
#########################################################################
# File Name: batchSangerBlast.sh
# Author: ChengYu
# Description: 
# Created Time: Sun Feb 19 18:52:23 2023
######################################################################### Usage:
#       bash batchSangerBlast.sh [DB.fa] [WorkDIR]
function getSeq(){
ls ./*zip|xargs -i unzip -O GBK {} ;
#rename.ul 程宇 chengyu ./* 
rename.ul "\报告成功" success ./*
if [ -f allseq.fa ];then rm allseq.fa; fi
ls ${1}/*success/*seq|while read id;
do
    name=`basename $id|cut -d"." -f1`
    echo ">"$name >> allseq.fa
    cat $id >> allseq.fa
    echo >> allseq.fa
done
}
function runBLAST(){
    makeblastdb -in ${1} -dbtype nucl -parse_seqids -out sanger.blastdb -logfile log
    blastn -query allseq.fa -db sanger.blastdb -out BLASTresult.format1.txt -outfmt 1 -num_threads 24 
    blastn -query allseq.fa -db sanger.blastdb -out BLASTresult.format6.txt -outfmt 6 -num_threads 24 
}
# main
function main(){
## Usage:
##       bash batchSangerBlast.sh [DB.fa] [WorkDIR]  
    cd $2
    getSeq $2
    runBLAST $1 
}

main $1 $2

