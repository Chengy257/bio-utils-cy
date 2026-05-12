#!/usr/bin/Rscript
## Usage:
## 
.libPaths("/home/chengyu/R/Rlib_4.2.3")
.libPaths()

# library(knitr, quietly=TRUE)
# library(sangerseqR, quietly=TRUE)
# library(Biostrings, quietly=TRUE)
# # opts_chunk$set(tidy=TRUE, tidy.opts=list(width.cutoff=70))


runParseSeq <- function(ab1,ref,filename,SignalCutoff){
    sangerseq <- readsangerseq(ab1)
    basecalls <- makeBaseCalls(sangerseq, ratio=SignalCutoff)
    ## phase alletes based on the input reference sequence
    phaseAlleles <- setAllelePhase(basecalls, ref,trim5=20)
    ## ouput the chromatogram plot
    chromatogram(phaseAlleles,trim5=20, width=100, height=1,showcalls = "both",filename = paste(filename,"_chromatogram.pdf", sep = ""))
    ## primarySeq, Basecalls fit to reference sequence
    pairAlign_prim <- pairwiseAlignment(primarySeq(phaseAlleles),ref,type="local-global")
    ## secondarySeq, the Alleles sequence
    pairAlign_second <- pairwiseAlignment(secondarySeq(phaseAlleles),ref,type="local-global")
    ## pairwiseAlignment
    writePairwiseAlignments(pairAlign_prim,file = paste(filename,"_pairwiseAlignment_ref.txt", sep = ""))
    writePairwiseAlignments(pairAlign_second,file = paste(filename,"_pairwiseAlignment_allele.txt", sep = ""))
}

main <- function(files,ref,SignalCutoff){
    for(i in 1:length(files)){
        ab1_file <- files[i]
        filename <- strsplit(basename(files[i]),"\\.")[[1]][1]
        print(paste0("[",date(),"] Parsing file: ", files[i]))
        runParseSeq(ab1_file,ref,filename,SignalCutoff)
    }
    print(paste("[",date(),"] All done!",sep=""))
}

## call main function

pkgs <- c('knitr','sangerseqR','Biostrings')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T,quietly=TRUE))})

SignalCutoff <- 0.33
ref <- read.table("ref.txt")
ref <- as.character(ref)
if( file.exists("ab1.files.txt") == FALSE){
    system("ls *ab1 > ab1.files.txt")
}
files <- read.table("ab1.files.txt") 
files <- files[,1]
main(files,ref,SignalCutoff)


################################
# sangerseq <- readsangerseq("C6-2.(P5382)L6-F.25850929.A04.ab1")
# ## 
# # chromatogram(test, width=80, height=3, trim5=50, trim3=100,showcalls='both')
# basecalls <- makeBaseCalls(sangerseq, ratio=0.33)
# # testcalls
# # chromatogram(testcalls, width=80, height=3, trim5=50, trim3=100,showcalls='both')
# chromatogram(testalleles,width=100, height=1,showcalls = "both",filename = "ttt.pdf")
# # ref <- subseq(primarySeq(test, string=TRUE))
# phaseAlleles <- setAllelePhase(basecalls, ref)
# pa_pri <- pairwiseAlignment(primarySeq(phaseAlleles),ref,type="local-global")
# pa_sec <- pairwiseAlignment(secondarySeq(phaseAlleles),ref,type="local-global")
# writePairwiseAlignments(pa,file = "test.txt")
# writePairwiseAlignments(pa_pri)
