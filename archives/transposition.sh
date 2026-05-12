#!/bin/bash
#########################################################################
# File Name: /home/chengyu/myscript/transposition.sh
# Author: ChengYu
# Description: 
# Created Time: Fri Apr 14 10:18:38 2023
#########################################################################


for i in `seq $(head -n 1 ${1} | awk '{print NF}')`; 
do 
    cut -f $i ${1} | tr "\n" " "| sed '$ s/$/\n/' >>  temp.transposition ; 
done

cat temp.transposition && rm temp.transposition

