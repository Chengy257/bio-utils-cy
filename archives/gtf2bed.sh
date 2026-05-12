#!/bin/bash
#########################################################################
# File Name: /home/chengyu/myscripts/gtf2bed.sh
# Author: ChengYu
# Description: 
# Created Time: Mon 25 Sep 2023 08:23:30 PM CST
#########################################################################



cat ${1}| gtfToGenePred stdin stdout|genePredToBed stdin ${1%%.*}.bed



