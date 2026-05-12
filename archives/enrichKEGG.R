#!/usr/bin/
#########################################################################
# File Name: /home/chengyu/myscripts/enrichKEGG.R
# Author: ChengYu
# Description: 
# Created Time: Wed 24 Apr 2024 03:04:55 PM CST
#########################################################################
#!/usr/bin/Rscript
## Usage:
## /usr/bin/Rscript /home/chengyu/workflows/enrich_GO_KEGG_clusterProfiler_gProfilerGO.R [gene list] [output dir prefix]     
.libPaths("/home/chengyu/R/Rlib_4.2.3")
.libPaths()
args <- commandArgs(T)
pkgs <- c('clusterProfiler','ggplot2','aPEAR','svglite','magrittr')

genelist <- read.table(args[1])[,1]
output_dir <- args[2]
out <- paste(output_dir,"/",sep="")
prefix <- gsub(".DEGs.txt","",basename(args[1]))

lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})
R.utils::setOption("clusterProfiler.download.method",'auto')

T2G <- read.table("/share/data/reference/osa/GO/KEGG/KEGG_TranscriptID2GeneIDs",header=F,sep="\t")

tmp <- merge(x=T2G,y=genelist,by.x=2,by.y=1)
translist <- unique(tmp[,2])

ekegg <- enrichKEGG(translist, organism = "dosa",keyType = "kegg",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1)

write.table(as.data.frame(ekegg@result %>% filter( pvalue <= 0.05)),file = paste0(out,prefix,"_EnrichResult_KEGG_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)

dotplot(ekegg,showCategory =10)+labs(title = "KEGG")
ggsave(filename = paste0(out,prefix,"_KEGG_",Sys.Date(),"_dotplot.pdf"),device = "pdf",width = 6,height = 6)

aPEAR::enrichmentNetwork(ekegg@result %>%  filter(pvalue<=0.05),colorBy = 'p.adjust', colorType = 'pval', drawEllipses = FALSE,repelLabels = TRUE,verbose = F,nodeSize = "Count")
ggsave(filename = paste0(out,prefix,"_KEGG_",Sys.Date(),"_Network.pdf"),device = "pdf",width = 6,height = 6)

