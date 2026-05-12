#!/usr/bin/env Rscript
#########################################################################
# File Name: 
# Author: ChengYu
# Description: 
# Created Time: 
#########################################################################
args <- commandArgs(T)
.libPaths("/home/chengyu/R/Rlib_4.2.3")
.libPaths()
## Load libraries
pkgs <- c('ComplexHeatmap')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})
filepathsList <- read.table(args[1],header=F)[,1]
groups <- read.table(args[2],header=F)[,1]  ## sampleHeader 
outname <- args[3]
# label_GeneName <- args[4]  ## add column_labels
# print(groups)
# print(outname)

set.seed("12345")
ht_list = NULL
for( i in 1:length(filepathsList)){
    filepath <- filepathsList[i]
    dat <- read.table(filepath,header = T,row.names=1) 
	# dat <- dat[,c(1,4,5,2,6,7,3,8,9)]
    dat.scaled <- t(apply(dat,1,scale))
    rownames(dat.scaled) <- rownames(dat)
    if (anyNA(dat.scaled)) {
        dat.scaled <- na.omit(dat.scaled)  
    }
    if (any(is.infinite(dat.scaled))) {
        dat.scaled[is.infinite(dat.scaled)] <- NA  # Replace infinite values with NA
        dat.scaled <- na.omit(dat.scaled)  # Remove rows with NA values
    }
    # add column name
    colnames(dat.scaled) <- colnames(dat)
    
    # groups <- apply(colnames(dat) %>% as.data.frame(),1,function(x) unlist(strsplit(x,split = "_"))[1])   
    GOcluster <- unlist(strsplit(basename(filepath),split = "\\."))[2]
    # print(GOcluster)
    col_anno = HeatmapAnnotation(groups=groups)
    row_anno = HeatmapAnnotation(GOcluster=rep(GOcluster,length(dat.scaled[,1])),which = "row")
    if( i == 1 ){
        p_single <- Heatmap(dat.scaled,cluster_rows = T,show_column_names = T,show_row_names = T,cluster_columns = F,column_names_rot = 45,top_annotation = col_anno,left_annotation = row_anno,name = "Expr. z-score",row_names_gp = gpar(fontsize=6),col=c("#3171AC","#FFFFFF","#D25536"))        
        p <- Heatmap(dat.scaled,cluster_rows = T,show_column_names = T,show_row_names = F,cluster_columns = F,column_names_rot = 45,top_annotation = col_anno,left_annotation = row_anno,name = "Expr. z-score",col=c("#3171AC","#FFFFFF","#D25536"))
        pdf(paste0(outname,"_",GOcluster,"_Single_Heatmap.pdf"),height=12 ,width=12)
            draw(p_single)
        dev.off()
    }
    else{
        p_single <- Heatmap(dat.scaled,cluster_rows = T,show_column_names = T,show_row_names = T,cluster_columns = F,column_names_rot = 45,top_annotation = col_anno,left_annotation = row_anno,name = "Expr. z-score",row_names_gp = gpar(fontsize=6),col=c("#3171AC","#FFFFFF","#D25536"))
        p <- Heatmap(dat.scaled,cluster_rows = T,show_column_names = T,show_row_names = F,cluster_columns = F,column_names_rot = 45,left_annotation = row_anno,name = "Expr. z-score",col=c("#3171AC","#FFFFFF","#D25536"))
        pdf(paste0(outname,"_",GOcluster,"_Single_Heatmap.pdf"),height=12 ,width=12)
            draw(p_single)
        dev.off()
    }
    ht_list = ht_list %v% p
}
pdf(paste0(outname,"_Multi_Heatmap.pdf"),height=12 ,width=12)
    draw(ht_list)
dev.off()
