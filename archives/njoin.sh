#!/bin/bash
#########################################################################
# File Name: join_n.sh
# Author: chengyu
# Description: join multi files
# Created Time: Wed Feb 23 11:12:34 2022
#########################################################################
### Usage:
###             njoin.sh [file1 file2 ... filen]
### Output:     join.out

cd `pwd`
array=($@)
a=$RANDOM
b=$RANDOM
join $1 $2 > join.out.$a
# echo $1
# echo $2
end=`expr $# - 1`
for i in `seq 2 $end `
do
	# echo ${array[$i]}
    join join.out.$a ${array[$i]} > temp.join.$b
    cat temp.join.$b > join.out.$a && rm temp.join.$b
done
sed -i 's/ /\t/g' join.out.$a
# echo [`date`]: Done, join all files in join.out.$a!
cat join.out.$a && rm join.out.$a 			    
