#!/usr/bin/Rscript
suppressMessages(library(getopt))
suppressMessages(library(Biostrings))
suppressMessages(library(ggplot2))

spec <- matrix(c("fasta","f",2,"character","Input fasta filename.",
                 "base","b",2,"character","Type of base to be counted.",
                 "window","w",2,"numeric","Sliding window size."),
                 byrow=T,ncol=5)
opt <- getopt(spec=spec)

fasta_file <- opt$fasta
BaseContent <- opt$base
window_size <- opt$window
cutoff <- 75
fasta_sequences <- readDNAStringSet(fasta_file)

for(j in 1:length(fasta_sequences)){

  name = paste(BaseContent,names(fasta_sequences)[j],"window_size",window_size,sep="_")
  # seq_len <- length(fasta_sequences[[j]])
  seq_start = round(window_size/2)  ## 从window一半处开始统计
  seq_end = length(fasta_sequences[[j]]) - round(window_size/2) + 1 ## 剩余window一半处结束统计

  base_content <- numeric()
  for (i in seq_start:seq_end) {
    window_seq <- subseq(fasta_sequences[[j]], start = i - round(window_size/2) +1, end = i + round(window_size/2) - 1)  ## 提取滑动窗口序列
    base_content[i] <- round(letterFrequency(window_seq ,BaseContent, as.prob = TRUE)[[1]] * 100,2)  ## 碱基百分比，保留两位小数
  }
  ## 
  base_content <- na.omit(base_content)
  index = seq_start:seq_end
  data <- data.frame(index = index, Content = base_content, Base = BaseContent)
  write.table(data,file = paste0(name,"_statistics_table.xls"),sep="\t",col.names = T,row.names = F)
  
  ggplot(data, aes(x = index, y = Content)) +
  geom_line()  +
  geom_hline(yintercept = cutoff,color="purple",size=1.1,lty=2) + 
  labs(title = paste(BaseContent," Content in Sliding Window of ",window_size," bp : ",name,sep=""),
        y = paste(BaseContent," Content %"))
  ggsave( filename = paste0(name,"_plot.pdf"), width = 6, height = 3, device = "pdf")
}
