#!/bin/bash
#########################################################################
# File Name: /home/chengyu/soft/tools/download_ascp.sh
# Author: ChengYu
# Description: 
# Created Time: Fri 05 May 2023 10:17:17 AM CST
#########################################################################

ssh=/opt/anaconda3/envs/download/etc/asperaweb_id_dsa.openssh

download() {
if [ $2 == "SE" ];then

#对于单端数据
    if [ ${#1} == 9 ] ; #根据accession number的长度下载地址的规律会发生变化
    then
        ascp -k1 -QT -l 300m -P33001 -i $ssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${1:0:6}/$1/${1}.fastq.gz ./;
    else
        ascp -k1 -QT -l 300m -P33001 -i $ssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${1:0:6}/00${1:0-1}/$1/${1}.fastq.gz ./;
    fi;
elif [ $2 == "PE" ];then
#对于双端数据
    if [ ${#1} == 9 ] ;then
        ascp -k1 -QT -l 300m -P33001 -i $ssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${1:0:6}/$1/${1}_1.fastq.gz ./;
        ascp -k1 -QT -l 300m -P33001 -i $ssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${1:0:6}/$1/${1}_2.fastq.gz ./;
    else
        ascp -k1 -QT -l 300m -P33001 -i $ssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${1:0:6}/00${1:0-1}/$1/${1}_1.fastq.gz ./;
        ascp -k1 -QT -l 300m -P33001 -i $ssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${1:0:6}/00${1:0-1}/$1/${1}_2.fastq.gz ./;
    fi
fi
}

cat $1|while read id;
do

download $id "PE"

done

