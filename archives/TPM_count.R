


# len <- c();
# for (i in c2$loci) {len <- c(len,genelength$length[which(i==genelength$loci)])};

# for (i in c2$) {len <- c(len,genelength$length[which(i==genelength$loci)])};
# c2$V3 <- len;



# c2 <- read.table("LR20I22AW2.raw_count.sorted")  
# c2$V3 <- c2$V1
# c2 <- c2[-1]
# names(c2) <- c("loci","counts")
# len <- c();
# for (i in c2$loci) {len <- c(len,genelength$length[which(i==genelength$loci)])};  
# c2$lenth <- len;

# c1$NL <- c1$counts/c1$lenth
# c1$TPM <- (c1$NL/(sum(c1$NL)))*(10^6)

# c6$NL <- c6$counts/c6$lenth; c6$TPM <- (c6$NL/(sum(c6$NL)))*(10^6)


countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen);
  denom <- log(sum(exp(rate)));
  exp(rate - denom + log(1e6));
};
tpm <- data.frame();
for(i in seq(3,20)){
tpm[,i] <- round(countToTpm(all[,i],all[,2]),2);
}

################
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen);
  denom <- log(sum(exp(rate)));
  exp(rate - denom + log(1e6));
};
##
countToFpkm <- function(counts, effLen)
{
  N <- sum(counts);
  exp( log(counts) + log(1e9) - log(effLen) - log(N) );
};
##
fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6));
};
##
countToEffCounts <- function(counts, len, effLen){
  counts * (len / effLen);
};
##
countToRPM <- function(counts){
  exp(log(counts) + log(1e6) - log(sum(counts)));
}
# An example
################################################################################
# cnts <- c(4250, 3300, 200, 1750, 50, 0);
# lens <- c(900, 1020, 2000, 770, 3000, 1777);
# countDf <- data.frame(count = cnts, length = lens);
 
# # assume a mean(FLD) = 203.7;
# countDf$effLength <- countDf$length - 203.7 + 1;
# countDf$tpm <- with(countDf, countToTpm(count, effLength));
# countDf$fpkm <- with(countDf, countToFpkm(count, effLength));
# with(countDf, all.equal(tpm, fpkmToTpm(fpkm)));
# countDf$effCounts <- with(countDf, countToEffCounts(count, length, effLength));

