#!/usr/bin/Rscript
#########################################################################
# File Name: /home/chengyu/myscripts/qPCR_analysis_plot.R
# Author: ChengYu
# Description: 
# Created Time: Wed 08 May 2024 09:39:36 PM CST
#########################################################################
.libPaths(c("/share/R/library/4.2.3",.libPaths()))
pkgs <- c('ggplot2','getopt','ggsci','patchwork')
lapply(pkgs, function(x){
   suppressPackageStartupMessages(library(x, character.only = T, quietly = T))})

spec <- matrix(c("dat","d",2,"character","Input data table, csv format [separated by comma].",
                 "ref","r",2,"character","Reference gene name [e.g. GAPDH/ACTIN].",
                 "control","c",2,"character","Control sample name [e.g. Control/CK/WT].",
                 "help", "h", 0, "logical", "Show this help information."),
                 byrow=T,ncol=5)
opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$dat) || is.null(opt$ref)|| is.null(opt$control)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}
# if ( is.null(opt$ref) ) {
#     opt$threads <- as.character("GAPDH") }
# if ( is.null(opt$control) ) {
#     opt$threads <- as.character("WT") }

dat <- as.character(opt$dat)
outname <- dat
# outname <- basename(dat)

dat <- read.csv(dat,skip = 34)
dat <- dat[1:(nrow(dat)-5),c(4,5,15)]
colnames(dat) <- c("sample","gene","CT")
dat$CT[which(dat$CT == "Undetermined")] <- 45  ## set as max loop number
dat$CT <- as.numeric(dat$CT)

ref_gene <- as.character(opt$ref)
sample_control <- as.character(opt$control)

## calculate delta CT
dat_genes <- data.frame()
for(sample in unique(dat$sample)){
  detect_gene <- unique(dat$gene[which(dat$sample == sample)])
  dat_ref_gene <- dat[which(dat$sample == sample & dat$gene == ref_gene),] 
  for(gene in detect_gene){
    if( gene != ref_gene ){
      dat_gene <- dat[which(dat$sample == sample & dat$gene == gene),]
      dat_gene <- cbind(dat_gene,dat_ref_gene)
      dat_genes <- rbind(dat_genes,dat_gene)
    }
  }
}
dat_genes$deltaCT <- dat_genes[,3] - dat_genes[,6]

## calculate delta delta CT
genes <- unique(dat_genes[,2])
dat_samples <- data.frame()
for( gene in genes){
  sub_dat <- dat_genes[which(dat_genes[,2]==gene),c(1,2,3,7)]
  control <- unique(sub_dat$sample[grep(sample_control,sub_dat$sample)])
  dat_control <- sub_dat[which(sub_dat$sample == control),]
  for(sample in unique(sub_dat$sample) ){
    dat_sample <- sub_dat[which(sub_dat$sample == sample),]
    dat_sample <- cbind(dat_sample,dat_control)
    dat_samples <- rbind(dat_samples,dat_sample)
  }
}
dat_samples$deltadeltaCT <- dat_samples[,4] - dat_samples[,8]
dat_samples$rel_exp <- 2^( -dat_samples$deltadeltaCT )
write.table(dat_samples,paste0(dirname(outname),"/parsed_result_table.xls"),sep="\t",quote=FALSE,col.names = TRUE,row.names =FALSE)

## plot 
# relativa exp.
p1 <- ggplot(dat_samples[,c(1,2,10)],aes(x=sample,y=(rel_exp),fill=gene))+geom_bar(position = "dodge",stat="summary",fun="mean")+geom_jitter(size=1.5)+stat_summary(fun = mean,geom = "errorbar",fun.max = function(x) mean(x) + sd(x),fun.min = function(x) mean(x) - sd(x),width=0.3)+facet_wrap(~gene,scales = "free")+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),legend.position = "null",text = element_text(family="sans",face="bold"),strip.background = element_rect(fill="white",color="white"),strip.placement = "outside",panel.border = element_rect(size = .8))+labs(x=NULL,y="Relative Expression",title = "Relative Expression")+scale_fill_npg()
# CT value
p2 <- ggplot(dat_samples[,c(1,2,3)],aes(x=sample,y=(CT),fill=gene))+geom_bar(position = "dodge",stat="summary",fun="mean")+geom_jitter(size=1.5)+stat_summary(fun = mean,geom = "errorbar",fun.max = function(x) mean(x) + sd(x),fun.min = function(x) mean(x) - sd(x),width=0.3)+facet_wrap(~gene,scales = "free")+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),legend.position = "null",text = element_text(family="sans",face="bold"),strip.background = element_rect(fill="white",color="white"),strip.placement = "outside",panel.border = element_rect(size = .8))+labs(x=NULL,y="CT value",title = "CT value")+scale_fill_npg()
# delta CT value
p3 <- ggplot(dat_samples[,c(1,2,4)],aes(x=sample,y=(deltaCT),fill=gene))+geom_bar(position = "dodge",stat="summary",fun="mean")+geom_jitter(size=1.5)+theme_bw()+stat_summary(fun = mean,geom = "errorbar",fun.max = function(x) mean(x) + sd(x),fun.min = function(x) mean(x) - sd(x),width=0.3)+facet_wrap(~gene,scales = "free")+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),legend.position = "null",text = element_text(family="sans",face="bold"),strip.background = element_rect(fill="white",color="white"),strip.placement = "outside",panel.border = element_rect(size = .8))+labs(x=NULL,y="deltaCT value",title = "deltaCT value")+scale_fill_npg()
# delta delta CT value 
# p4 <- ggplot(dat_samples[,c(1,2,9)],aes(x=sample,y=(deltadeltaCT),fill=gene))+geom_bar(position = "dodge",stat="summary",fun="mean")+geom_jitter(size=1.5)+theme_bw()+stat_summary(fun = mean,geom = "errorbar",fun.max = function(x) mean(x) + sd(x),fun.min = function(x) mean(x) - sd(x),width=0.3)+facet_wrap(~gene,scales = "free")+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),legend.position = "null",text = element_text(family="sans",face="bold"),strip.background = element_rect(fill="white",color="white"),strip.placement = "outside",panel.border = element_rect(size = .8))+labs(x=NULL,y="deltadeltaCT value",title = "deltadeltaCT value")

# Reference Gene CT value
p4 <- ggplot(dat_genes[,c(1,2,6)],aes(x=sample,y=(CT),fill=gene))+geom_bar(position = "dodge",stat="summary",fun="mean")+geom_jitter(size=1.5)+theme_bw()+stat_summary(fun = mean,geom = "errorbar",fun.max = function(x) mean(x) + sd(x),fun.min = function(x) mean(x) - sd(x),width=0.3)+facet_wrap(~gene,scales = "free")+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),legend.position = "null",text = element_text(family="sans",face="bold"),strip.background = element_rect(fill="white",color="white"),strip.placement = "outside",panel.border = element_rect(size = .8))+scale_fill_npg()+labs(x=NULL,y="CT value",title = paste0(ref_gene," CT value"))+scale_fill_npg()

(p1|p2)/(p3|p4)+plot_annotation(title = paste0("Figure. QPCR Parsed Results Barplots.      [",date(),"]"),tag_levels = 'A',caption = paste("Generated by qPCR_analysis_plot Rscript.","Author: Chengyu",sep="\n"))
ggsave(paste0(outname,".parsed_result_barplot.pdf"),width = 10,height = 10)
