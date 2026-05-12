#!/usr/bin/R
#########################################################################
# File Name: /home/chengyu/myscripts/cal_expr_FPKM_TPM.R
# Author: ChengYu
# Description: calculate gene expression, using non-redundant exons length sum of gene
# Created Time: Mon Aug 21 16:25:17 CST 2023
#########################################################################
## Usage:
## Rscript /home/chengyu/myscripts/cal_expr_FPKM_TPM.R [gtf] [gene count tab file, STAR output ReadsPerGene.out.tab]   
suppressMessages(library(GenomicFeatures))
suppressMessages(library(parallel))
suppressMessages(library(getopt))
spec <- matrix(c("gtf","g",2,"character","Input gtf file path.",
                 "count","c",2,"character","Input count file path.",
                 "strand","s",2,"numeric","Library strandness type.",
                 "threads","t",1,"numeric","Threads to be used."),
                 byrow=T,ncol=5)
opt <- getopt(spec=spec)
# args<-commandArgs(T)
gtf <- opt$gtf
count_file <- opt$count
strand <- opt$strand
threads <- opt$threads
cl <- makeCluster(threads)
txdb <- makeTxDbFromGFF(gtf,format="gtf") 
exons_gene <- exonsBy(txdb, by = "gene") 
exons_gene_lens <- parLapply(cl,exons_gene,function(x){sum(width(reduce(x)))}) 
geneid_efflen <- data.frame(geneid=names(exons_gene_lens),efflen=as.numeric(exons_gene_lens))
count <- read.table(count_file,skip=4)
if(strand == 1){
    count = count[,c(1,3)]  ## 1(FeatureCounts) fr-secondstrand(Tophat) F/FR(HISAT2) yes(HTSeq) fr(StringTie) 
}else if(strand == 2){
    count = count[,c(1,4)]  ## 2(FeatureCounts) fr-firststrand(Tophat) R/RF(HISAT2) reverse(HTSeq) rf(StringTie) 
}else if(strand == 0){
    count = count[,c(1,2)]  ## unstranded
}else{
    count = count[,c(1,2)]  ## others, unstranded, defaulted value
}
## 
merged = merge(geneid_efflen,count,by=1)
colnames(merged) <- c("geneid","efflen","count")
countToExpr <- function(df){
    counts <- df[,3]
    effLen <- df[,2] 
    rate <- log(counts)-log(effLen)
    denom <- log(sum(exp(rate)))
    TPM <- exp(rate-denom+log(1e6))
    N <- sum(counts)
    FPKM <- exp(log(counts)+log(1e9)-log(effLen)-log(N))
    return(cbind(FPKM,TPM))
}
## 
out <- cbind(merge,countToExpr(merge))
write.table(out,filename=paste(basename(count_file),"_expr.xls",sep="."),header=T,row.names=FALSE,quote=FALSE,sep="\t")

# ReadsPerGene.out.tab:
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
# column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
