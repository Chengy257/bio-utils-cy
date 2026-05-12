#!/usr/bin/R
#########################################################################
# File Name: /home/chengyu/myscripts/calculate_EffectiveLength.R
# Author: ChengYu
# Description: calculate non-redundant exons length sum of gene
# Created Time: Wed 05 Jul 2023 10:16:47 AM CST
#########################################################################
## Usage:
##      Rscript /home/chengyu/myscripts/calculate_EffectiveLength.R [gtf]    

suppressMessages(library(GenomicFeatures))
suppressMessages(library(parallel))
args<-commandArgs(T)
cl <- makeCluster(6)  
gtf <- args[1]
output <- paste(gtf,"efflen",sep=".")
txdb <- makeTxDbFromGFF(gtf,format="gtf") 
exons_gene <- exonsBy(txdb, by = "gene") 
exons_gene_lens <- parLapply(cl,exons_gene,function(x){sum(width(reduce(x)))}) 
geneid_efflen <- data.frame(geneid=names(exons_gene_lens),efflen=as.numeric(exons_gene_lens))
write.table(geneid_efflen,output,col.names = F,row.names = F,sep="\t",quote = F)
