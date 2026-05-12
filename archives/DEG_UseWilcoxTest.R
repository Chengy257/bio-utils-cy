# Source code from: https://rpubs.com/LiYumei/806213
# Cite: Li Y, Ge X, Peng F, Li W, Li JJ. Exaggerated false positives by popular differential expression methods when analyzing human population samples. Genome Biol 23, 79 (2022). https://doi.org/10.1186/s13059-022-02648-4

# Example:
# example.condition.tsv file; tab delimited:
# condition2	condition2	condition1	condition2	condition1	condition2	condition1
# example.countMatrix.tsv file; tab delimited:
# Gene	Sample1	Sample2	Sample3	Sample4	Sample5
# Gene1	516	349	2699	692	2630
# Gene2	0	0	3	2	1

###########################################################
suppressWarnings(library(edgeR, quietly = T))
suppressWarnings(library(getopt))

# args <- commandArgs(T)
spec <- matrix(
  c("countMatrix",  "c", 2, "character", "Raw count matrix file, tab delimited.",
    "conditions", "t", 2, "character",  "Samples conditon ordered same with countMatrix, tab delimited.",
    "help", "h", 0, "logical", "Show this help information."),byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
if( !is.null(opt$help) || is.null(opt$countMatrix) || is.null(opt$conditions)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

# Read the read count matrix file and the condition labels file.
count <- as.character(opt$countMatrix)
sample <- as.character(opt$conditions)

readCount<-read.table(count, header = T, row.names = 1, stringsAsFactors = F,check.names = F)

readCount<-readCount[rowSums(readCount)>10,]

conditions<-read.table(sample, header = F)
conditions<-factor(t(conditions))

# Count matrix preprocessing using edgeR package
y <- DGEList(counts=readCount,group=conditions)
##Remove rows consistently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
##Perform TMM normalization and transfer to CPM (Counts Per Million)
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)


# Run the Wilcoxon rank-sum test for each gene
pvalues <- sapply(1:nrow(count_norm),function(i){
     data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
     p=wilcox.test(gene~conditions, data)$p.value
     return(p)
   })
fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene

conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

# Output results based on FDR threshold
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
fdrThres=0.05
write.table(outRst[outRst$FDR<fdrThres,], file="WilcoxonTest_DEGs_output.tsv",sep="\t", quote=F,row.names = T,col.names = T)
write.table(outRst, file="WilcoxonTest_all_output.tsv",sep="\t", quote=F,row.names = T,col.names = T)
saveRDS(object = outRst, file = "WilcoxonTest.rds")
system("sed -i '1 s/^/gene_id\t/' WilcoxonTest_*.tsv")
