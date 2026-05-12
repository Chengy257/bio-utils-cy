#!/usr/bin/env Rscript
#########################################################################
# File Name: deseq2_multigroup.R
# Author: ChengYu
# Description: DESeq2 differential expression analysis for multi-group
#              comparisons with optional batch correction.
# Created Time: 2026
#########################################################################
"""DESeq2 multi-group differential expression analysis.

Performs DESeq2 analysis with all-vs-control pairwise comparisons.
Generates normalized counts, heatmaps, PCA plots, volcano plots,
and MA plots. Supports optional batch effect correction.
"""

pkgs <- c("DESeq2", "ggplot2", "BiocParallel", "gplots", "RColorBrewer", "amap", "getopt")
for (pkg in pkgs) {
    suppressMessages(library(pkg, character.only = TRUE))
}

spec <- matrix(c(
    "count",       "c", 2, "character", "Raw count matrix (TSV: genes x samples)",
    "sample",      "s", 2, "character", "Sample info CSV (columns: id,group[,batch])",
    "output",      "o", 2, "character", "Output prefix / directory",
    "batch",       "b", 1, "logical",   "Enable batch correction (default: FALSE)",
    "control",     "r", 1, "character", "Control group name (default: control)",
    "fdr",         "p", 1, "numeric",   "FDR threshold (default: 0.05)",
    "foldchange",  "f", 1, "numeric",   "Fold change threshold (default: 2)",
    "threads",     "t", 1, "numeric",   "Number of threads (default: 1)",
    "help",        "h", 0, "logical",   "Show help"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (!is.null(opt$help) || is.null(opt$count) || is.null(opt$sample) || is.null(opt$output)) {
    cat("Usage: Rscript deseq2_multigroup.R -c counts.tsv -s samples.csv -o result_prefix\n")
    cat("       [-b] [-r control] [-p 0.05] [-f 2] [-t 4]\n\n")
    cat("Options:\n")
    cat("  -c --count       Raw count matrix (TSV)\n")
    cat("  -s --sample      Sample info CSV (id,group[,batch])\n")
    cat("  -o --output      Output prefix\n")
    cat("  -b --batch       Enable batch correction\n")
    cat("  -r --control     Control group name (default: control)\n")
    cat("  -p --fdr         FDR threshold (default: 0.05)\n")
    cat("  -f --foldchange  Fold change threshold (default: 2)\n")
    cat("  -t --threads     Threads (default: 1)\n")
    quit(status = if (!is.null(opt$help)) 0 else 1)
}

# Defaults
if (is.null(opt$batch))      opt$batch <- FALSE
if (is.null(opt$control))    opt$control <- "control"
if (is.null(opt$fdr))        opt$fdr <- 0.05
if (is.null(opt$foldchange)) opt$foldchange <- 2
if (is.null(opt$threads))    opt$threads <- 1

logFC_threshold <- log2(opt$foldchange)

# Validate inputs
if (!file.exists(opt$count))  stop("Count file not found: ", opt$count)
if (!file.exists(opt$sample)) stop("Sample file not found: ", opt$sample)

out_dir <- paste0(opt$output, "_DESeq2_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("[INFO] Output directory:", out_dir, "\n")

# Register parallel backend
register(MulticoreParam(opt$threads))

# --- Read data ---
count_data <- read.table(opt$count, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
sample_info <- read.csv(opt$sample, stringsAsFactors = FALSE)

cat("[INFO] Count matrix:", nrow(count_data), "genes x", ncol(count_data), "samples\n")
cat("[INFO] Sample info:", nrow(sample_info), "samples\n")

# Validate sample overlap
common_samples <- intersect(colnames(count_data), sample_info$id)
if (length(common_samples) == 0) stop("No matching samples between count matrix and sample info.")
count_data <- count_data[, common_samples]
sample_info <- sample_info[sample_info$id %in% common_samples, ]
sample_info$group <- factor(sample_info$group)

if (!(opt$control %in% levels(sample_info$group))) {
    stop("Control group '", opt$control, "' not found. Available: ", paste(levels(sample_info$group), collapse = ", "))
}

cat("[INFO] Groups:", paste(levels(sample_info$group), collapse = ", "), "\n")

# --- Build DESeqDataSet ---
if (opt$batch && "batch" %in% colnames(sample_info)) {
    cat("[INFO] Using batch correction (batch column detected).\n")
    dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = sample_info,
        design = ~ batch + group
    )
} else {
    dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = sample_info,
        design = ~ group
    )
}

# --- Pre-filter low counts ---
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep, ]
cat("[INFO] Genes after filtering:", nrow(dds), "\n")

# --- Run DESeq2 ---
cat("[INFO] Running DESeq2...\n")
dds <- DESeq(dds, parallel = TRUE)

# --- Sample plots ---
cat("[INFO] Generating diagnostic plots...\n")
vsd <- vst(dds, blind = FALSE)

# Heatmap
pdf(file.path(out_dir, "vst_Pearson_heatmap.pdf"), height = 14, width = 12)
pheatmap(
    cor(assay(vsd), method = "pearson"),
    display_numbers = TRUE,
    color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
    fontsize_number = 8
)
dev.off()

# PCA
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw()
ggsave(file.path(out_dir, "vst_PCA_plot.pdf"), p, width = 8, height = 6)

# --- Pairwise comparisons (vs control) ---
treatments <- setdiff(levels(sample_info$group), opt$control)
cat("[INFO] Comparisons:", paste(treatments, "vs", opt$control, collapse = "; "), "\n")

res_dir <- file.path(out_dir, "DEG_tables")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

for (trt in treatments) {
    cat("[INFO] Analyzing:", trt, "vs", opt$control, "\n")
    contrast <- c("group", trt, opt$control)
    res <- results(dds, contrast = contrast, alpha = opt$fdr)
    res <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr")

    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df$regulation <- ifelse(
        res_df$padj < opt$fdr & res_df$log2FoldChange > logFC_threshold, "Up",
        ifelse(res_df$padj < opt$fdr & res_df$log2FoldChange < -logFC_threshold, "Down", "Unsig")
    )

    # Write results
    out_prefix <- paste0(gsub("[^a-zA-Z0-9]", "_", trt), "_vs_", gsub("[^a-zA-Z0-9]", "_", opt$control))
    write.table(res_df, file.path(res_dir, paste0(out_prefix, "_DESeq2.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Volcano plot
    res_df$neg_log10_padj <- -log10(res_df$padj)
    vp <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = regulation)) +
        geom_point(alpha = 0.6, size = 0.8) +
        scale_color_manual(values = c(Down = "blue", Up = "red", Unsig = "grey")) +
        geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "grey") +
        geom_hline(yintercept = -log10(opt$fdr), linetype = "dashed", color = "grey") +
        xlab("log2 Fold Change") + ylab("-log10 FDR") +
        ggtitle(paste(trt, "vs", opt$control)) +
        theme_bw() + theme(legend.position = "bottom")
    ggsave(file.path(res_dir, paste0(out_prefix, "_volcano.pdf")), vp, width = 6, height = 5)

    # Summary
    cat("  Up:", sum(res_df$regulation == "Up", na.rm = TRUE),
        " Down:", sum(res_df$regulation == "Down", na.rm = TRUE), "\n")
}

# --- Normalized counts ---
norm_counts <- counts(dds, normalized = TRUE)
write.table(norm_counts, file.path(out_dir, "normalized_counts.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("[INFO] Done. Results in:", out_dir, "\n")
