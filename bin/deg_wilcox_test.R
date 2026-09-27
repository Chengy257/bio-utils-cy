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
#
# Changelog:
#   v1.1.0  2026-09-27
#   - FIX: genes whose values are fully tied across both groups get
#     p = NaN from wilcox.test; they were silently mixed into _all.tsv.
#     They are now counted and reported (kept with NA padj, excluded
#     from _DEGs.tsv), and the per-gene tie warnings are suppressed
#   - FIX: condition labels are read as character — numeric-looking
#     group names ("1", "2") were coerced to integers
#   - NEW: -r/--control names the reference condition explicitly; the
#     comparison direction used to depend on file-column order
#     (default unchanged: first condition in file order = reference)
#   - NEW: -m/--min_count parameterizes the rowSums filter (default 10,
#     unchanged). NOTE: this is a raw-count threshold, not a CPM one —
#     its stringency drifts with library size; the log2FC column is a
#     difference of log2(mean CPM + 1) (pseudo-log scale, not the DESeq2
#     shrunken log2FC)
#   - CHANGE: logs go to stderr (message() instead of cat())

suppressMessages(library(edgeR))
suppressMessages(library(getopt))

spec <- matrix(c(
    "count",    "c", 2, "character", "Raw count matrix (TSV)",
    "condition","t", 2, "character", "Condition labels file (one-line TSV, same order as count columns)",
    "control",  "r", 1, "character", "Reference condition for FC direction (default: first condition in file order)",
    "min_count","m", 1, "numeric",   "Keep genes with rowSums > this raw-count total (default: 10)",
    "fdr",      "p", 1, "numeric",   "FDR threshold (default: 0.05)",
    "output",   "o", 1, "character", "Output prefix (default: WilcoxTest)",
    "help",     "h", 0, "logical",   "Show help"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (!is.null(opt$help) || is.null(opt$count) || is.null(opt$condition)) {
    cat("Usage: Rscript deg_wilcox_test.R -c counts.tsv -t conditions.tsv [-p 0.05] [-o prefix]\n")
    cat("       [-r control] [-m 10]\n\n")
    cat("Options:\n")
    cat("  -c --count      Raw count matrix (TSV)\n")
    cat("  -t --condition  Condition labels file (one-line TSV, same order as count columns)\n")
    cat("  -r --control    Reference condition for FC direction (default: first in file order)\n")
    cat("  -m --min_count  Keep genes with rowSums > this total (default: 10)\n")
    cat("  -p --fdr        FDR threshold (default: 0.05)\n")
    cat("  -o --output     Output prefix (default: WilcoxTest)\n")
    cat("  -h --help       Show this help\n")
    quit(status = if (!is.null(opt$help)) 0 else 1)
}

if (is.null(opt$fdr))       opt$fdr <- 0.05
if (is.null(opt$output))    opt$output <- "WilcoxTest"
if (is.null(opt$min_count)) opt$min_count <- 10

if (!file.exists(opt$count))     stop("Count file not found: ", opt$count)
if (!file.exists(opt$condition)) stop("Condition file not found: ", opt$condition)

# Read data
message("[INFO] Reading count matrix...")
readCount <- read.table(opt$count, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# read labels as character: numeric-looking group names ("1", "2") must
# not silently become integers
conditions <- unlist(read.table(opt$condition, sep = "\t", stringsAsFactors = FALSE,
                                colClasses = "character"))

if (ncol(readCount) != length(conditions)) {
    stop("Column count (", ncol(readCount), ") != condition count (", length(conditions), ")")
}

message("[INFO] ", nrow(readCount), " genes, ", ncol(readCount), " samples")
message("[INFO] Conditions: ", paste(unique(conditions), collapse = ", "))

if (length(unique(conditions)) != 2) {
    stop("Exactly 2 conditions required. Found: ", paste(unique(conditions), collapse = ", "))
}
if (!is.null(opt$control) && !(opt$control %in% conditions)) {
    stop("Control condition '", opt$control, "' not found. Available: ",
         paste(unique(conditions), collapse = ", "))
}

# Filter lowly expressed genes
keep <- rowSums(readCount) > opt$min_count
readCount <- readCount[keep, ]
message("[INFO] After filtering (rowSums > ", opt$min_count, "): ", nrow(readCount), " genes")

# TMM normalization
dge <- DGEList(counts = readCount)
dge <- calcNormFactors(dge, method = "TMM")
cpm_data <- cpm(dge)

# Determine groups: control (if given) is the reference; otherwise the
# first condition in file order is (legacy behaviour)
cond_levels <- unique(conditions)
if (!is.null(opt$control)) {
    cond_levels <- c(opt$control, setdiff(cond_levels, opt$control))
}
idx_a <- which(conditions == cond_levels[1])
idx_b <- which(conditions == cond_levels[2])

message("[INFO] Comparing: ", cond_levels[2], " (n=", length(idx_b), ") vs ",
        cond_levels[1], " (n=", length(idx_a), ") [reference]")

# Wilcoxon test per gene (ties warnings are ubiquitous on CPM values and
# add no information; fully tied genes yield NaN and are reported below)
message("[INFO] Running Wilcoxon tests...")
pvals <- apply(cpm_data, 1, function(row) {
    suppressWarnings(wilcox.test(row[idx_a], row[idx_b], exact = FALSE)$p.value)
})

# Fold change (log2 mean)
log2fc <- log2(rowMeans(cpm_data[, idx_b, drop = FALSE]) + 1) -
          log2(rowMeans(cpm_data[, idx_a, drop = FALSE]) + 1)

# FDR adjustment
padj <- p.adjust(pvals, method = "fdr")

# Fully tied genes (identical values in both groups) get p = NaN from
# wilcox.test — report them instead of letting them slip through silently
n_nan <- sum(is.na(pvals))
if (n_nan > 0) {
    message("[WARN] ", n_nan, " gene(s) have undefined p-values (all values tied across groups); ",
            "they keep NA padj in _all.tsv and are excluded from _DEGs.tsv.")
}

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

message("[INFO] Total: ", nrow(result), " | Significant (FDR <", opt$fdr, "): ", nrow(sig))
message("[INFO] Up: ", sum(sig$log2FoldChange > 0, na.rm = TRUE),
        " | Down: ", sum(sig$log2FoldChange < 0, na.rm = TRUE))
message("[INFO] Output: ", opt$output, "_all.tsv, ", opt$output, "_DEGs.tsv")
