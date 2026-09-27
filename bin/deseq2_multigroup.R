#!/usr/bin/env Rscript
#########################################################################
# File Name: deseq2_multigroup.R
# Author: ChengYu
# Description: DESeq2 differential expression analysis for multi-group
#              comparisons with optional batch correction.
# Created Time: 2026
#########################################################################
# DESeq2 multi-group differential expression analysis.
#
# Performs DESeq2 analysis with all-vs-control pairwise comparisons.
# Generates normalized counts, heatmaps, PCA plots, volcano plots,
# and MA plots. Supports optional batch effect correction.
#
# Changelog:
#   v1.1.0  2026-09-27
#   - FIX: pheatmap is loaded at startup — the correlation heatmap call
#     crashed with "could not find function \"pheatmap\"" after the whole
#     DESeq run had completed
#   - FIX: ashr is checked at startup; lfcShrink(type="ashr") used to
#     crash mid-run on machines without it
#   - FIX: the getopt spec required an argument for the -b flag ("flag
#     \"b\" requires an argument" — -b was unusable as documented); it is
#     a plain flag now
#   - FIX: vst() errors ("less than 'nsub' rows with mean normalized
#     count > 5") on small or low-count matrices — the diagnostic plots
#     crashed after the whole DESeq run; falls back to
#     varianceStabilizingTransformation exactly when vst's own
#     applicability condition does not hold
#   - FIX: -b without a batch column in the sample table now warns
#     instead of silently falling back to ~group
#   - FIX: the sample table is validated (id/group columns required, no
#     duplicate ids) and the count matrix is checked for NA / negative /
#     non-integer values before DESeq2 sees them (DESeq2's own errors
#     are obscure)
#   - FIX: genes with NA padj (independent filtering) are counted and
#     reported; they were silently dropped from the volcano plot and
#     masked by na.rm=TRUE in the summary
#   - CLEAN: dead dependencies gplots/amap removed (never used; the
#     heatmap is drawn with pheatmap)
#   - CHANGE: logs go to stderr (message() instead of cat(); help output
#     stays on stdout)

pkgs <- c("DESeq2", "ggplot2", "BiocParallel", "pheatmap", "RColorBrewer", "ashr", "getopt")
for (pkg in pkgs) {
    suppressMessages(library(pkg, character.only = TRUE))
}

spec <- matrix(c(
    "count",       "c", 2, "character", "Raw count matrix (TSV: genes x samples)",
    "sample",      "s", 2, "character", "Sample info CSV (columns: id,group[,batch])",
    "output",      "o", 2, "character", "Output prefix / directory",
    "batch",       "b", 0, "logical",   "Enable batch correction (default: FALSE)",
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

if (opt$threads < 1) stop("--threads must be >= 1, got: ", opt$threads)

logFC_threshold <- log2(opt$foldchange)

# Validate inputs
if (!file.exists(opt$count))  stop("Count file not found: ", opt$count)
if (!file.exists(opt$sample)) stop("Sample file not found: ", opt$sample)

out_dir <- paste0(opt$output, "_DESeq2_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("[INFO] Output directory: ", out_dir)

# Register parallel backend
register(MulticoreParam(opt$threads))

# --- Read data ---
count_data <- read.table(opt$count, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
sample_info <- read.csv(opt$sample, stringsAsFactors = FALSE)

message("[INFO] Count matrix: ", nrow(count_data), " genes x ", ncol(count_data), " samples")
message("[INFO] Sample info: ", nrow(sample_info), " samples")

# Validate sample table (columns, duplicate ids) and count values up
# front — DESeq2's own errors for these are cryptic.
missing_cols <- setdiff(c("id", "group"), colnames(sample_info))
if (length(missing_cols) > 0) {
    stop("Sample table is missing required column(s): ",
         paste(missing_cols, collapse = ", "), ". Required: id,group[,batch]")
}
if (anyDuplicated(sample_info$id)) {
    stop("Duplicate sample ids in sample table: ",
         paste(unique(sample_info$id[duplicated(sample_info$id)]), collapse = ", "))
}
count_mat <- as.matrix(count_data)
if (anyNA(count_mat)) {
    stop("Count matrix contains ", sum(is.na(count_mat)),
         " NA value(s); fill or drop the affected genes first.")
}
if (any(count_mat < 0)) {
    stop("Count matrix contains ", sum(count_mat < 0), " negative value(s).")
}
if (any(count_mat != round(count_mat))) {
    stop("Count matrix contains ", sum(count_mat != round(count_mat)),
         " non-integer value(s); DESeq2 requires raw counts.")
}

# Validate sample overlap
common_samples <- intersect(colnames(count_data), sample_info$id)
if (length(common_samples) == 0) stop("No matching samples between count matrix and sample info.")
count_data <- count_data[, common_samples]
sample_info <- sample_info[sample_info$id %in% common_samples, ]
sample_info$group <- factor(sample_info$group)

if (!(opt$control %in% levels(sample_info$group))) {
    stop("Control group '", opt$control, "' not found. Available: ", paste(levels(sample_info$group), collapse = ", "))
}

message("[INFO] Groups: ", paste(levels(sample_info$group), collapse = ", "))

# --- Build DESeqDataSet ---
if (opt$batch && "batch" %in% colnames(sample_info)) {
    message("[INFO] Using batch correction (batch column detected).")
    dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = sample_info,
        design = ~ batch + group
    )
} else {
    if (opt$batch) {
        message("[WARN] Batch correction requested (-b) but the sample table has no 'batch' column; using ~group.")
    }
    dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = sample_info,
        design = ~ group
    )
}

# --- Pre-filter low counts ---
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep, ]
message("[INFO] Genes after filtering: ", nrow(dds))

# --- Run DESeq2 ---
message("[INFO] Running DESeq2...")
dds <- DESeq(dds, parallel = TRUE)

# --- Sample plots ---
message("[INFO] Generating diagnostic plots...")
# vst() requires >= nsub (1000) rows with mean normalized count > 5 and
# errors otherwise; fall back to the full varianceStabilizingTransformation
# (vst's own recommendation) exactly when that condition does not hold.
if (sum(rowMeans(counts(dds, normalized = TRUE)) > 5) >= 1000) {
    vsd <- vst(dds, blind = FALSE)
} else {
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
}

# Heatmap (invisible(): top-level dev.off() would auto-print "null device"
# to stdout and pollute the pipeline's data channel)
pdf(file.path(out_dir, "vst_Pearson_heatmap.pdf"), height = 14, width = 12)
pheatmap(
    cor(assay(vsd), method = "pearson"),
    display_numbers = TRUE,
    color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
    fontsize_number = 8
)
invisible(dev.off())

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
message("[INFO] Comparisons: ", paste(treatments, "vs", opt$control, collapse = "; "))

res_dir <- file.path(out_dir, "DEG_tables")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

for (trt in treatments) {
    message("[INFO] Analyzing: ", trt, " vs ", opt$control)
    contrast <- c("group", trt, opt$control)
    res <- results(dds, contrast = contrast, alpha = opt$fdr)
    res <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr")

    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df$regulation <- ifelse(
        res_df$padj < opt$fdr & res_df$log2FoldChange > logFC_threshold, "Up",
        ifelse(res_df$padj < opt$fdr & res_df$log2FoldChange < -logFC_threshold, "Down", "Unsig")
    )

    # Genes with NA padj (independent filtering / all-zero counts) are
    # excluded from the volcano plot — report them instead of dropping
    # them silently.
    n_na_padj <- sum(is.na(res_df$padj))
    if (n_na_padj > 0) {
        message("[WARN] ", n_na_padj, " gene(s) have NA padj and are excluded from the volcano plot / summary counts.")
    }

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
    message("[INFO]   Up: ", sum(res_df$regulation == "Up", na.rm = TRUE),
            " Down: ", sum(res_df$regulation == "Down", na.rm = TRUE),
            " Unsig: ", sum(res_df$regulation == "Unsig", na.rm = TRUE),
            " NA padj: ", n_na_padj)
}

# --- Normalized counts ---
norm_counts <- counts(dds, normalized = TRUE)
write.table(norm_counts, file.path(out_dir, "normalized_counts.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

message("[INFO] Done. Results in: ", out_dir)
