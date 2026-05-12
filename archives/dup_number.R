#!/mnt/c/miniconda3/envs/r4.1/bin/R
args<-commandArgs(T)
a=read.table(args[1],sep="\t", encoding = "UTF-8")
out=args[2]
a <- as.data.frame(a)
freq <- as.data.frame(table(a))
a[,2] <- 1:length(a[,1])
b <- c()
for (i in 1:length(freq[,1])) {
	b <- c(b,paste(freq[i,1],1:freq[i,2],sep = "_"))
}
a <- a[order(a[,1]),]
a[,3] <- b
b <- a[order(a[,2]),]
b <- b[,-2]
write.table(b,file=out,sep="\t",col.names=F,row.names=F,quote=F)
q(save="no")
