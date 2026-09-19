#!/usr/bin/env Rscript
#########################################################################
# File Name: calculate_effective_length.R
# Author: ChengYu
# Description: Calculate non-redundant exon length per gene from GTF.
# Created Time: 2026
#########################################################################
"""Calculate non-redundant exon length per gene from GTF.

Overlapping exons within each gene are merged (reduced) before summing,
giving the effective genomic length used for TPM/FPKM calculations.
"""

suppressMessages(library(GenomicFeatures))
suppressMessages(library(parallel))
suppressMessages(library(getopt))

spec <- matrix(c(
    "gtf",     "g", 2, "character", "GTF annotation file",
    "output",  "o", 1, "character", "Output file (default: {gtf}.efflen)",
    "threads", "t", 1, "numeric",   "Threads (default: 4)",
    "help",    "h", 0, "logical",   "Show help"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (!is.null(opt$help) || is.null(opt$gtf)) {
    cat("Usage: Rscript calculate_effective_length.R -g annotation.gtf [-o output.efflen] [-t 4]\n")
    quit(status = if (!is.null(opt$help)) 0 else 1)
}

if (is.null(opt$output))  opt$output <- paste0(opt$gtf, ".efflen")
if (is.null(opt$threads)) opt$threads <- 4

if (!file.exists(opt$gtf)) stop("GTF not found: ", opt$gtf)

cat("[INFO] GTF:", opt$gtf, "\n")
cat("[INFO] Threads:", opt$threads, "\n")

# Build TxDb
cat("[INFO] Building TxDb...\n")
cl <- makeCluster(opt$threads)
txdb <- makeTxDbFromGFF(opt$gtf, format = "gtf")

# Calculate effective lengths
exons_by_gene <- exonsBy(txdb, by = "gene")
eff_lengths <- sapply(exons_by_gene, function(exons) {
    sum(width(reduce(exons)))
})
stopCluster(cl)

# Write output
result <- data.frame(gene_id = names(eff_lengths), efflen = as.numeric(eff_lengths),
                     row.names = NULL, stringsAsFactors = FALSE)
write.table(result, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

cat("[INFO]", nrow(result), "genes ->", opt$output, "\n")
