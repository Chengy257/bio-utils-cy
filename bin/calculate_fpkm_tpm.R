#!/usr/bin/env Rscript
#########################################################################
# File Name: calculate_fpkm_tpm.R
# Author: ChengYu
# Description: Calculate FPKM and TPM from STAR gene counts using
#              effective gene lengths from GTF.
# Created Time: 2026
#########################################################################
# Calculate FPKM and TPM from STAR gene-level counts.
#
# Reads a STAR ReadsPerGene.out.tab file and GTF annotation, computes
# non-redundant exon lengths, then calculates FPKM and TPM values.
#
# Changelog:
#   v1.1.0  2026-09-27
#   - FIX: sum(rpk) == 0 or sum(counts) == 0 (all-zero matrix) stopped
#     with a bare division-by-zero warning and NaN output; both now stop
#     with a clear message
#   - CLEAN: the makeCluster() call was dead code — the cluster was
#     created and stopped but never used; -t is accepted for
#     compatibility and validated but does not change behaviour
#   - DOC: strand column semantics clarified (col 3 = sense / secondstrand
#     protocol, i.e. featureCounts -s 1; col 4 = antisense / firststrand,
#     i.e. featureCounts -s 2) and genes dropped for missing/zero
#     effective length are counted
#   - CHANGE: logs go to stderr (message() instead of cat(); help stays
#     on stdout)

suppressMessages(library(GenomicFeatures))
suppressMessages(library(getopt))

spec <- matrix(c(
    "gtf",     "g", 2, "character", "GTF annotation file",
    "count",   "c", 2, "character", "STAR ReadsPerGene.out.tab file",
    "strand",  "s", 2, "numeric",   "Strandness: 0=unstranded, 1=secondstrand (sense, featureCounts -s 1), 2=firststrand (dUTP, featureCounts -s 2)",
    "threads", "t", 1, "numeric",   "Accepted for compatibility; the script is single-threaded (default: 4)",
    "output",  "o", 1, "character", "Output file (default: {count}_expr.xls)",
    "help",    "h", 0, "logical",   "Show help"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (!is.null(opt$help) || is.null(opt$gtf) || is.null(opt$count) || is.null(opt$strand)) {
    cat("Usage: Rscript calculate_fpkm_tpm.R -g annotation.gtf -c ReadsPerGene.out.tab -s 0 [-t 4]\n")
    cat("\nStrandness codes (same as STAR):\n")
    cat("  0 = unstranded          -> column 2\n")
    cat("  1 = sense / secondstrand (featureCounts -s 1) -> column 3\n")
    cat("  2 = antisense / firststrand, dUTP (featureCounts -s 2) -> column 4\n")
    quit(status = if (!is.null(opt$help)) 0 else 1)
}

if (is.null(opt$threads)) opt$threads <- 4
if (is.null(opt$output))  opt$output <- paste0(opt$count, "_expr.xls")

if (opt$threads < 1) stop("--threads must be >= 1, got: ", opt$threads)
if (!opt$strand %in% c(0, 1, 2)) {
    stop("Invalid strandness: ", opt$strand, ". Use 0, 1, or 2.")
}
if (!file.exists(opt$gtf))   stop("GTF not found: ", opt$gtf)
if (!file.exists(opt$count)) stop("Count file not found: ", opt$count)

message("[INFO] GTF: ", opt$gtf)
message("[INFO] Count: ", opt$count)
message("[INFO] Strandness: ", opt$strand)

# Build TxDb and calculate effective lengths
message("[INFO] Building TxDb from GTF...")
txdb <- makeTxDbFromGFF(opt$gtf, format = "gtf")

exons_by_gene <- exonsBy(txdb, by = "gene")
eff_lengths <- sapply(exons_by_gene, function(exons) {
    sum(width(reduce(exons)))
})

eff_df <- data.frame(gene_id = names(eff_lengths), efflen = as.numeric(eff_lengths),
                     row.names = NULL, stringsAsFactors = FALSE)

# Read STAR counts (skip 4 header lines)
message("[INFO] Reading STAR counts...")
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
n_input <- nrow(df)
df <- merge(df, eff_df, by = "gene_id", all.x = TRUE)
n_dropped <- sum(is.na(df$efflen) | df$efflen <= 0)
df <- df[!is.na(df$efflen) & df$efflen > 0, ]

message("[INFO] ", nrow(df), " genes with valid effective lengths (", n_input, " in STAR table; ",
        n_dropped, " dropped: not in GTF or zero-length).")

# Calculate TPM and FPKM
rpk <- df$count / df$efflen
scaling_factor <- sum(rpk) / 1e6
if (scaling_factor <= 0) {
    stop("sum(count/efflen) is 0 — the count matrix is all zero; TPM is undefined.")
}
df$TPM <- rpk / scaling_factor

total_counts <- sum(df$count)
if (total_counts <= 0) {
    stop("sum(count) is 0 — the count matrix is all zero; FPKM is undefined.")
}
df$FPKM <- df$count / (df$efflen / 1e3) / (total_counts / 1e6)

# Write output
df <- df[, c("gene_id", "efflen", "count", "FPKM", "TPM")]
write.table(df, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

message("[INFO] Output: ", opt$output)
message("[INFO] Done.")
