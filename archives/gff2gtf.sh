#!/bin/bash
#########################################################################
# File Name: /home/chengyu/myscripts/gff2gtf.sh
# Author: ChengYu
# Description: 
# Created Time: Tue 03 Oct 2023 05:57:21 PM CST
#########################################################################

/home/chengyu/soft/ucsc-tools/gff3ToGenePred ${1} stdout| /home/chengyu/soft/ucsc-tools/genePredToGtf file stdin ${1%%.*}.gtf

