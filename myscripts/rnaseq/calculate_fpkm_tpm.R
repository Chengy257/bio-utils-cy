#!/usr/bin/env Rscript
#########################################################################
# File Name: calculate_fpkm_tpm.R
# Author: ChengYu
# Description: Calculate FPKM and TPM from STAR gene counts using
#              effective gene lengths from GTF.
# Created Time: 2026
#########################################################################
"""Calculate FPKM and TPM from STAR gene-level counts.

Reads a STAR ReadsPerGene.out.tab file and GTF annotation, computes
non-redundant exon lengths, then calculates FPKM and TPM values.
"""

suppressMessages(library(GenomicFeatures))
suppressMessages(library(parallel))
suppressMessages(library(getopt))

spec <- matrix(c(
    "gtf",     "g", 2, "character", "GTF annotation file",
    "count",   "c", 2, "character", "STAR ReadsPerGene.out.tab file",
    "strand",  "s", 2, "numeric",   "Strandness: 0=unstranded, 1=secondstrand, 2=firststrand",
    "threads", "t", 1, "numeric",   "Number of threads (default: 4)",
    "output",  "o", 1, "character", "Output file (default: {count}_expr.xls)",
    "help",    "h", 0, "logical",   "Show help"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (!is.null(opt$help) || is.null(opt$gtf) || is.null(opt$count) || is.null(opt$strand)) {
    cat("Usage: Rscript calculate_fpkm_tpm.R -g annotation.gtf -c ReadsPerGene.out.tab -s 0 [-t 4]\n")
    cat("\nStrandness codes (same as STAR):\n")
    cat("  0 = unstranded (column 2)\n")
    cat("  1 = stranded, forward first (column 3)\n")
    cat("  2 = stranded, reverse first (column 4)\n")
    quit(status = if (!is.null(opt$help)) 0 else 1)
}

if (is.null(opt$threads)) opt$threads <- 4
if (is.null(opt$output))  opt$output <- paste0(opt$count, "_expr.xls")

if (!file.exists(opt$gtf))   stop("GTF not found: ", opt$gtf)
if (!file.exists(opt$count)) stop("Count file not found: ", opt$count)

cat("[INFO] GTF:", opt$gtf, "\n")
cat("[INFO] Count:", opt$count, "\n")
cat("[INFO] Strandness:", opt$strand, "\n")

# Build TxDb and calculate effective lengths
cat("[INFO] Building TxDb from GTF...\n")
cl <- makeCluster(opt$threads)
txdb <- makeTxDbFromGFF(opt$gtf, format = "gtf")

exons_by_gene <- exonsBy(txdb, by = "gene")
eff_lengths <- sapply(exons_by_gene, function(exons) {
    sum(width(reduce(exons)))
})
stopCluster(cl)

eff_df <- data.frame(gene_id = names(eff_lengths), efflen = as.numeric(eff_lengths),
                     row.names = NULL, stringsAsFactors = FALSE)

# Read STAR counts (skip 4 header lines)
cat("[INFO] Reading STAR counts...\n")
star_data <- read.table(opt$count, skip = 4, sep = "\t", header = FALSE,
                        colClasses = c("character", "numeric", "numeric", "numeric"),
                        stringsAsFactors = FALSE)
colnames(star_data) <- c("gene_id", "unstranded", "strand_fwd", "strand_rev")

# Select correct strand column
count_col <- switch(as.character(opt$strand),
    "0" = "unstranded",
    "1" = "strand_fwd",
    "2" = "strand_rev",
    stop("Invalid strandness: ", opt$strand, ". Use 0, 1, or 2.")
)

df <- data.frame(
    gene_id = star_data$gene_id,
    count = star_data[[count_col]],
    stringsAsFactors = FALSE
)

# Merge with effective lengths
df <- merge(df, eff_df, by = "gene_id", all.x = TRUE)
df <- df[!is.na(df$efflen) & df$efflen > 0, ]

cat("[INFO]", nrow(df), "genes with valid effective lengths.\n")

# Calculate TPM and FPKM
rpk <- df$count / df$efflen
scaling_factor <- sum(rpk) / 1e6
df$TPM <- rpk / scaling_factor

total_counts <- sum(df$count)
df$FPKM <- df$count / (df$efflen / 1e3) / (total_counts / 1e6)

# Write output
df <- df[, c("gene_id", "efflen", "count", "FPKM", "TPM")]
write.table(df, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

cat("[INFO] Output:", opt$output, "\n")
cat("[INFO] Done.\n")
