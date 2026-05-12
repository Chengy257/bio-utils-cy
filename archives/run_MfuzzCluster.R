#!/usr/bin/env Rscript
#########################################################################
# File Name: 
# Author: ChengYu
# Description: 
# Created Time: 
#########################################################################

##
# dat <- na.omit(dat)   # filter.NA(dat) 
# std_dev <- apply(dat, 1, sd)
# dat <- dat[std_dev>0, ]
# dat <- scale(dat)
runMfuzz <- function(dat,cluster_num){

dat <- as.matrix(dat)
dat <- new('ExpressionSet',exprs = dat)
dat <- filter.NA(dat)
#dat <- fill.NA(dat) # for illustration only; rather use knn method
dat <- filter.std(dat, min.std = 0)#
dat <- standardise(dat)
#
set.seed(12345)
cl <- mfuzz(dat,c=cluster_num,m=mestimate(dat))

     return(cl)
 }

args <- commandArgs(T)
.libPaths("/home/chengyu/R/Rlib_4.2.3")
## Load libraries
pkgs <- c('Mfuzz')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})
## 
dat <- read.table(args[1],header=T,row.names=1,sep="\t")
cluster_num <- as.numeric(args[2])
outname <- as.character(args[3])

#cl <- runMfuzz(dat,cluster_num)

dat <- as.matrix(dat)
dat <- new('ExpressionSet',exprs = dat)
dat <- filter.NA(dat)
#dat <- fill.NA(dat) # for illustration only; rather use knn method
dat <- filter.std(dat, min.std = 0)#
dat <- standardise(dat)
#
set.seed(12345)
cl <- mfuzz(dat,c=cluster_num,m=mestimate(dat))


pdf(paste0(outname,"_MfuzzPlot.pdf"),width=6,height=4)
    mfuzz.plot2(dat,cl=cl,x11=FALSE,centre = TRUE)
dev.off()
print("Mfuzz cluster member numbers stat:",quote=F)
print(table(cl$cluster))

cluster <- as.data.frame(cbind(names(cl$cluster),cl$cluster))
membership <- as.data.frame(cbind(rownames(cl$membership),cl$membership))
out_merged <- merge(cluster,membership,by=1)
colnames(out_merged)[1:2] <- c("Name","Cluster")
write.table(out_merged,paste0(outname,"_Mfuzz_clusterMembership.xls"),col.names = T,row.names = F,quote = F,sep="\t")
##
