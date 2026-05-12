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
# filepathsList <- read.table(args[1],header=F)[,1]
# groups <- read.table(args[2],header=F)[,1]  ## sampleHeader 
# outname <- args[3]
# label_GeneName <- args[4]  ## add column_labels
# print(groups)
# print(outname)
# ht_list = NULL
# for( i in 1:length(filepathsList)){
#     filepath <- filepathsList[i]
#     dat <- read.table(filepath,header = T,row.names=1) 
#     dat.scaled <- t(apply(dat,1,scale))
#     # add row name labels 
#     label_GeneName <- read.table(paste0(filepath,".Gene_Name_anno"))
#     # print(label_GeneName)
#     rownames(dat.scaled) <- label_GeneName[,2]
#     # print(rownames(dat.scaled))
#     # check data
#     if (anyNA(dat.scaled)) {
#     # Handle missing values 
#         dat.scaled <- na.omit(dat.scaled)  
#     }
#     # Check for infinite values
#     if (any(is.infinite(dat.scaled))) {
#     # Handle infinite values 
#         dat.scaled[is.infinite(dat.scaled)] <- NA  # Replace infinite values with NA
#         dat.scaled <- na.omit(dat.scaled)  # Remove rows with NA values
#     }
#     # add column name
#     colnames(dat.scaled) <- colnames(dat)
    
#     # groups <- apply(colnames(dat) %>% as.data.frame(),1,function(x) unlist(strsplit(x,split = "_"))[1])   
#     GOcluster <- unlist(strsplit(basename(filepath),split = "\\."))[2]
#     # print(GOcluster)
#     col_anno = HeatmapAnnotation(groups=groups)
#     row_anno = HeatmapAnnotation(GOcluster=rep(GOcluster,length(dat.scaled[,1])),which = "row")
#     if( i == 1 ){
#         p_single <- Heatmap(dat.scaled,cluster_rows = T,show_column_names = T,show_row_names = T,cluster_columns = T,column_names_rot = 45,top_annotation = col_anno,left_annotation = row_anno,name = "Expr. z-score",row_names_gp = gpar(fontsize=6))        
#         p <- Heatmap(dat.scaled,cluster_rows = T,show_column_names = T,show_row_names = F,cluster_columns = T,column_names_rot = 45,top_annotation = col_anno,left_annotation = row_anno,name = "Expr. z-score")
#         pdf(paste0(outname,"_",GOcluster,"_Single_Heatmap.pdf"),height=12 ,width=12)
#             draw(p_single)
#         dev.off()
#     }
#     else{
#         p_single <- Heatmap(dat.scaled,cluster_rows = T,show_column_names = T,show_row_names = T,cluster_columns = T,column_names_rot = 45,top_annotation = col_anno,left_annotation = row_anno,name = "Expr. z-score",row_names_gp = gpar(fontsize=6))
#         p <- Heatmap(dat.scaled,cluster_rows = T,show_column_names = T,show_row_names = F,cluster_columns = T,column_names_rot = 45,left_annotation = row_anno,name = "Expr. z-score")
#         pdf(paste0(outname,"_",GOcluster,"_Single_Heatmap.pdf"),height=12 ,width=12)
#             draw(p_single)
#         dev.off()
#     }
#     ht_list = ht_list %v% p
# }

# spec <- matrix(c("dat","d",2,"character","Input fasta filename.",
#                  "col_anno","c",1,"character","Type of base to be counted.",
#                  "row_anno","r",1,"numeric","Sliding window size.",
#                  "labels","l",1,"numeric","Sliding window size."),
#                  byrow=T,ncol=5)
# opt <- getopt(spec=spec)

# fasta_file <- opt$fasta
# BaseContent <- opt$base
# window_size <- opt$window

dat <- read.table(args[1],header=F)
groups <- read.table(args[2],header=F)[,1]  ## sampleHeader 
outname <- args[3]

dat <- read.table(args[1],header = T,row.names=1) 
dat.scaled <- t(apply(dat,1,scale))
# add row name labels 
# label_GeneName <- read.table(paste0(filepath,".Gene_Name_anno"))
# print(label_GeneName)
# rownames(dat.scaled) <- label_GeneName[,2]
# print(rownames(dat.scaled))
# check data
if (anyNA(dat.scaled)) {
# Handle missing values 
    dat.scaled <- na.omit(dat.scaled)  
}
# Check for infinite values
if (any(is.infinite(dat.scaled))) {
# Handle infinite values 
    dat.scaled[is.infinite(dat.scaled)] <- NA  # Replace infinite values with NA
    dat.scaled <- na.omit(dat.scaled)  # Remove rows with NA values
}
# add column name
colnames(dat.scaled) <- colnames(dat)

data_grouped <- data[, match(group_table$ID, colnames(data))] # 按照分组表的顺序重新排列列
colnames(data_grouped) <- paste0("Group_", group_table$Group) # 替换列名为分组名


pdf(paste0(outname,"_Multi_Heatmap.pdf"),height=12 ,width=12)
    draw(ht_list)
dev.off()

parseData <- function(file){
  df <- read.table(file,header = T,row.names = 1)
  df.scaled <- t(apply(df,1,scale))
  colnames(df.scaled) <- sample_info[colnames(df),1]
  df.scaled[is.infinite(df.scaled)] <- NA 
  df.scaled <- na.omit(df.scaled)
  return(df.scaled)
}

plotHeatmapSimply <- function(dat){
  col_order <- order(colnames(dat))
  col_anno = HeatmapAnnotation(groups=colnames(dat))
  p <- Heatmap(dat,cluster_rows = T,show_column_names = T,show_row_names = T,cluster_columns = F,column_names_rot = 45,top_annotation = col_anno,name = "scaled expr.",column_order = col_order,row_names_gp = gpar(fontsize=8))
  return(p)
}

onestep <- function(input){
  dat <- parseData(input);
  p <- plotHeatmapSimply(dat)
  return(p)
}

