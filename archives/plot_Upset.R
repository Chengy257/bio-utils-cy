#!/usr/bin/Rscript
#########################################################################
# File Name: /home/chengyu/myscripts/plot_Upset.R
# Author: ChengYu
# Description: 
# Created Time: Thu 14 Dec 2023 07:55:37 PM CST
#########################################################################

args <- commandArgs(T)
.libPaths("/home/chengyu/R/Rlib_4.2.3")
.libPaths()
## Load libraries
pkgs <- c('UpSetR')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})


getList <- function(file_list){
  myList <- list()
  for(i in 1:length(file_list[,1])){
    tmp <- read.table(file_list[i,1],header = F,sep=" ")
    name <- basename(file_list[i,1])
    myList[[name]] <- tmp[,1]
  }
  return(myList)
}

files <- read.table(args[1],header = F)
out <- args[2]
head(files)
llist <- getList(files)
pdf(paste0(out,"upset_plot.pdf"))
  upset(fromList(llist)) 
dev.off()
