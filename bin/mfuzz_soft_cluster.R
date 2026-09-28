#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/mfuzz_soft_cluster.R
# Author: ChengYu
# Description: Soft clustering of time-series expression data using Mfuzz.
#              Reads a gene x sample expression matrix (TSV), performs
#              filtering, standardization, and fuzzy c-means clustering,
#              then outputs cluster membership table and plot.
# Created Time: 2026
# Changelog:
#   v1.1.0  2026-09-28
#   - FIX: filter.std() kept its default visu = TRUE, which draws the
#     standard-deviation histogram before any PDF device is open - a
#     stray Rplots.pdf appeared in the working directory (or an error
#     on headless devices). Plotting is now disabled.
#   - FIX: duplicated sample column names are rejected (check.names =
#     FALSE previously let identical time points through silently).
#   - FIX: -k/--clusters must be >= 2 (mfuzz requires at least two
#     clusters).
#   - FIX: mfuzz.plot2 was called with its default x11 = TRUE, which
#     opens an interactive X11 window instead of drawing to the PDF
#     device - the plot could never be produced on a headless server.
#     It now draws to the PDF (x11 = FALSE).
#   - DOCUMENTED: filter.NA drops genes with more than 25% missing
#     values (Mfuzz default threshold).
#   - CHANGE: dev.off() is silent.
#########################################################################

suppressMessages(library(Mfuzz))
suppressMessages(library(Biobase))

## ---- Argument parsing with getopt ----
library(getopt)

spec <- matrix(c(
    "input",       "i", 1, "character", "Expression matrix TSV (rows=genes, cols=samples; required)",
    "clusters",    "k", 1, "integer",   "Number of soft clusters (required)",
    "output",      "o", 1, "character", "Output prefix (required)",
    "min-std",     "s", 2, "numeric",   "Minimum standard deviation filter (default: 0)",
    "seed",        "S", 2, "integer",   "Random seed for reproducibility (default: 12345)",
    "plot-width",  "W", 2, "numeric",   "Plot width in inches (default: 6)",
    "plot-height", "H", 2, "numeric",   "Plot height in inches (default: 4)",
    "help",        "h", 0, "logical",   "Show this help message"
), ncol = 5, byrow = TRUE)

opt <- getopt(spec)

## ---- Help ----
if (!is.null(opt$help)) {
    cat("Usage: Rscript mfuzz_soft_cluster.R -i expr.tsv -k 6 -o output_prefix [options]\n\n")
    cat("Soft clustering of time-series expression data using Mfuzz.\n\n")
    cat("Options:\n")
    cat("  -i, --input       Expression matrix TSV (rows=genes, cols=samples) [required]\n")
    cat("  -k, --clusters    Number of clusters [required]\n")
    cat("  -o, --output      Output file prefix [required]\n")
    cat("  -s, --min-std     Minimum standard deviation filter (default: 0)\n")
    cat("  -S, --seed        Random seed (default: 12345)\n")
    cat("  -W, --plot-width  Plot width in inches (default: 6)\n")
    cat("  -H, --plot-height Plot height in inches (default: 4)\n")
    cat("  -h, --help        Show this help message\n\n")
    cat("Outputs:\n")
    cat("  {prefix}_MfuzzPlot.pdf           Cluster membership plot\n")
    cat("  {prefix}_Mfuzz_clusterMembership.tsv  Cluster membership table\n\n")
    cat("Example:\n")
    cat("  Rscript mfuzz_soft_cluster.R -i counts.tsv -k 8 -o results/mfuzz\n")
    quit(status = 0)
}

## ---- Validate required arguments ----
if (is.null(opt$input)) {
    stop("Error: -i/--input is required. Use -h for help.")
}
if (is.null(opt$clusters)) {
    stop("Error: -k/--clusters is required. Use -h for help.")
}
if (is.null(opt$output)) {
    stop("Error: -o/--output is required. Use -h for help.")
}

if (!file.exists(opt$input)) {
    stop(paste("Error: Input file not found:", opt$input))
}
if (length(opt$clusters) != 1 || is.na(opt$clusters) || opt$clusters < 2) {
    stop("Error: -k/--clusters must be an integer >= 2.")
}

## ---- Set defaults ----
min_std    <- if (is.null(opt[["min-std"]]))    0     else opt[["min-std"]]
seed_val   <- if (is.null(opt$seed))            12345 else opt$seed
plot_w     <- if (is.null(opt[["plot-width"]]))  6     else opt[["plot-width"]]
plot_h     <- if (is.null(opt[["plot-height"]])) 4     else opt[["plot-height"]]

## ---- Main analysis function ----
run_mfuzz <- function(input_file, cluster_num, out_prefix,
                      min.std = 0, seed = 12345,
                      plot.width = 6, plot.height = 4) {

    message("Reading expression matrix: ", input_file)
    raw_dat <- read.table(input_file, header = TRUE, row.names = 1,
                          sep = "\t", quote = "", comment.char = "",
                          check.names = FALSE)

    if (ncol(raw_dat) < 2) {
        stop("Error: Expression matrix must have at least 2 sample columns.")
    }
    if (anyDuplicated(colnames(raw_dat))) {
        stop("Error: duplicated sample column names in the expression matrix: ",
             paste(unique(colnames(raw_dat)[duplicated(colnames(raw_dat))]),
                   collapse = ", "))
    }
    if (nrow(raw_dat) < cluster_num) {
        stop("Error: Number of genes (", nrow(raw_dat),
             ") must be >= number of clusters (", cluster_num, ").")
    }

    message("Loaded ", nrow(raw_dat), " genes x ", ncol(raw_dat), " samples")

    # Convert to ExpressionSet
    expr_mat <- as.matrix(raw_dat)
    eset <- Biobase::ExpressionSet(assayData = expr_mat)

    # Filter NA values: Mfuzz default thres = 0.25 removes genes with
    # more than 25% missing values
    message("Filtering NA values (genes with > 25% missing values are dropped)...")
    eset <- Mfuzz::filter.NA(eset)

    # Filter by standard deviation (visu = FALSE: the histogram would
    # otherwise be drawn on a device that is not open yet)
    message("Filtering by minimum standard deviation: ", min.std)
    eset <- Mfuzz::filter.std(eset, min.std = min.std, visu = FALSE)

    gene_count <- nrow(Biobase::exprs(eset))
    if (gene_count < cluster_num) {
        stop("Error: Only ", gene_count,
             " genes remain after filtering; need >= ", cluster_num, ".")
    }
    message("Retained ", gene_count, " genes after filtering")

    # Standardize
    message("Standardizing expression values...")
    eset <- Mfuzz::standardise(eset)

    # Estimate fuzzifier (m parameter)
    m_est <- Mfuzz::mestimate(eset)
    message("Estimated fuzzifier (m): ", round(m_est, 4))

    # Run Mfuzz clustering
    message("Running Mfuzz clustering with k = ", cluster_num, "...")
    set.seed(seed)
    cl <- Mfuzz::mfuzz(eset, c = cluster_num, m = m_est)

    # --- Output 1: Cluster membership plot ---
    plot_file <- paste0(out_prefix, "_MfuzzPlot.pdf")
    message("Writing cluster plot: ", plot_file)
    pdf(plot_file, width = plot.width, height = plot.height)
    # x11 = FALSE: the default opens an interactive X11 window instead
    # of drawing to the open PDF device (the script could never produce
    # its plot on a headless server without this)
    Mfuzz::mfuzz.plot2(eset, cl, xlab = "Sample", ylab = "Expression",
                       x11 = FALSE)
    invisible(dev.off())

    # --- Output 2: Cluster membership table ---
    membership_file <- paste0(out_prefix, "_Mfuzz_clusterMembership.tsv")
    message("Writing cluster membership table: ", membership_file)

    membership <- cl$membership
    cluster_assign <- apply(membership, 1, which.max)
    result_df <- data.frame(
        Gene       = rownames(membership),
        Cluster    = cluster_assign,
        Membership = membership[cbind(seq_len(nrow(membership)), cluster_assign)],
        stringsAsFactors = FALSE
    )

    # Append individual cluster membership columns
    for (k_idx in seq_len(ncol(membership))) {
        colname <- paste0("Cluster_", k_idx)
        result_df[[colname]] <- membership[, k_idx]
    }

    write.table(result_df, file = membership_file, sep = "\t",
                quote = FALSE, row.names = FALSE)

    # --- Print cluster member count summary ---
    message("\n=== Cluster Member Counts (max membership assignment) ===")
    counts <- table(cluster_assign)
    for (cid in sort(as.integer(names(counts)))) {
        message(sprintf("  Cluster %d: %d genes", cid, counts[as.character(cid)]))
    }
    message(sprintf("  Total: %d genes", sum(counts)))
    message("\nDone.")

    invisible(cl)
}

## ---- Execute ----
tryCatch({
    run_mfuzz(
        input_file  = opt$input,
        cluster_num = opt$clusters,
        out_prefix  = opt$output,
        min.std     = min_std,
        seed        = seed_val,
        plot.width  = plot_w,
        plot.height = plot_h
    )
}, error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(status = 1)
})
