#!/usr/bin/env Rscript
#########################################################################
# File Name: deg_wilcox_test.R
# Author: ChengYu
# Description: Differential expression using Wilcoxon rank-sum test
#              with TMM normalization (suitable for population data).
# Created Time: 2026
#########################################################################
# Differential expression using Wilcoxon rank-sum test.
#
# TMM-normalizes counts, then applies Wilcoxon rank-sum test per gene
# between two conditions. Suitable for population-level studies where
# negative binomial assumptions may not hold.

suppressMessages(library(edgeR))
suppressMessages(library(getopt))

spec <- matrix(c(
    "count",    "c", 2, "character", "Raw count matrix (TSV)",
    "condition","t", 2, "character", "Condition labels file (one-line TSV, same order as count columns)",
    "fdr",      "p", 1, "numeric",   "FDR threshold (default: 0.05)",
    "output",   "o", 1, "character", "Output prefix (default: WilcoxTest)",
    "help",     "h", 0, "logical",   "Show help"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (!is.null(opt$help) || is.null(opt$count) || is.null(opt$condition)) {
    cat("Usage: Rscript deg_wilcox_test.R -c counts.tsv -t conditions.tsv [-p 0.05] [-o prefix]\n")
    quit(status = if (!is.null(opt$help)) 0 else 1)
}

if (is.null(opt$fdr))    opt$fdr <- 0.05
if (is.null(opt$output)) opt$output <- "WilcoxTest"

if (!file.exists(opt$count))     stop("Count file not found: ", opt$count)
if (!file.exists(opt$condition)) stop("Condition file not found: ", opt$condition)

# Read data
cat("[INFO] Reading count matrix...\n")
readCount <- read.table(opt$count, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
conditions <- unlist(read.table(opt$condition, sep = "\t", stringsAsFactors = FALSE))

if (ncol(readCount) != length(conditions)) {
    stop("Column count (", ncol(readCount), ") != condition count (", length(conditions), ")")
}

cat("[INFO]", nrow(readCount), "genes,", ncol(readCount), "samples\n")
cat("[INFO] Conditions:", paste(unique(conditions), collapse = ", "), "\n")

if (length(unique(conditions)) != 2) {
    stop("Exactly 2 conditions required. Found: ", paste(unique(conditions), collapse = ", "))
}

# Filter lowly expressed genes
keep <- rowSums(readCount) > 10
readCount <- readCount[keep, ]
cat("[INFO] After filtering (rowSums > 10):", nrow(readCount), "genes\n")

# TMM normalization
dge <- DGEList(counts = readCount)
dge <- calcNormFactors(dge, method = "TMM")
cpm_data <- cpm(dge)

# Determine groups
cond_levels <- unique(conditions)
idx_a <- which(conditions == cond_levels[1])
idx_b <- which(conditions == cond_levels[2])

cat("[INFO] Comparing:", cond_levels[1], "(n=", length(idx_a), ") vs",
    cond_levels[2], "(n=", length(idx_b), ")\n")

# Wilcoxon test per gene
cat("[INFO] Running Wilcoxon tests...\n")
pvals <- apply(cpm_data, 1, function(row) {
    wilcox.test(row[idx_a], row[idx_b], exact = FALSE)$p.value
})

# Fold change (log2 mean)
log2fc <- log2(rowMeans(cpm_data[, idx_b, drop = FALSE]) + 1) -
          log2(rowMeans(cpm_data[, idx_a, drop = FALSE]) + 1)

# FDR adjustment
padj <- p.adjust(pvals, method = "fdr")

# Build result table
result <- data.frame(
    gene_id = rownames(readCount),
    log2FoldChange = round(log2fc, 4),
    pvalue = round(pvals, 6),
    padj = round(padj, 6),
    stringsAsFactors = FALSE
)
result <- result[order(result$padj), ]

# Write outputs
sig <- result[result$padj < opt$fdr & !is.na(result$padj), ]
write.table(result, paste0(opt$output, "_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig, paste0(opt$output, "_DEGs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("[INFO] Total:", nrow(result), "| Significant (FDR <", opt$fdr, "):", nrow(sig), "\n")
cat("[INFO] Up:", sum(sig$log2FoldChange > 0), "| Down:", sum(sig$log2FoldChange < 0), "\n")
cat("[INFO] Output:", opt$output, "_all.tsv,", opt$output, "_DEGs.tsv\n")
