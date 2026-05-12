#!/usr/bin/R
#########################################################################
# File Name: /home/chengyu/myscripts/merge.R
# Author: ChengYu
# Description: 
# Created Time: Wed 18 Oct 2023 03:38:23 PM CST
#########################################################################

args <- commandArgs(T)

filelist <- read.table(args[1],header=FALSE,quote="")
outputfilename <- paste0(dirname(filelist[1,1]),"/Rmerged_output_file")

# apple(filelist,2,function(x))

merged_list <- read.table(filelist[1,1],header=T)
head(merged_list)
for(i in 2:length(filelist[,1])){
    df <- read.table(filelist[i,1],header=T)
    merged_list <- merge(merged_list,df,by=1)
}

write.table(merged_list,file=outputfilename,col.names=T,row.names=FALSE,quote=FALSE,sep="\t")

