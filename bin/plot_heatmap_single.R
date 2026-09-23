#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/plot_heatmap_single.R
# Author: ChengYu
# Description: Generate a publication-quality heatmap from a single
#              expression matrix (genes x samples) using ComplexHeatmap.
#              Rows are z-score scaled; NA/Inf values are handled
#              gracefully. Optional sample-group annotation is displayed
#              as a top column bar.
# Created Time: 2026
#########################################################################

suppressMessages(library(getopt))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(RColorBrewer))

# ── CLI specification ──────────────────────────────────────────────────
spec <- matrix(c(
    "data",              "d", 1, "character", "Path to expression matrix (TSV, rows=genes, cols=samples, header=TRUE, row names in col 1). Required.",
    "groups",            "g", 2, "character", "Path to sample-group file (one label per line, same order as matrix columns). Optional.",
    "output",            "o", 2, "character", "Output prefix for PDF file. Default: 'heatmap'.",
    "show-row-names",    "s", 2, "logical",   "Show gene (row) names on the heatmap. Default: TRUE.",
    "row-names-fontsize","r", 2, "numeric",   "Font size for row names. Default: 6.",
    "column-cluster",    "c", 2, "logical",   "Perform hierarchical clustering on columns. Default: FALSE.",
    "height",            "H", 2, "numeric",   "Plot height in inches. Default: 10.",
    "width",             "W", 2, "numeric",   "Plot width in inches. Default: 10.",
    "help",              "h", 0, "logical",   "Print this help message and exit."
), byrow = TRUE, ncol = 5)
colnames(spec) <- c("long", "short", "argflag", "type", "help")

opt <- getopt(spec)

# ── Help ───────────────────────────────────────────────────────────────
if (!is.null(opt$help)) {
    cat("Usage: Rscript plot_heatmap_single.R -d <matrix.tsv> [options]\n\n")
    cat("Generate a z-score scaled heatmap with optional sample-group annotation.\n\n")
    cat("Options:\n")
    for (i in seq_len(nrow(spec))) {
        flag <- paste0("-", spec[i, "short"], ", --", spec[i, "long"])
        cat(sprintf("  %-30s %s\n", flag, spec[i, "help"]))
    }
    cat("\nExamples:\n")
    cat("  Rscript plot_heatmap_single.R -d expr_matrix.tsv -o result\n")
    cat("  Rscript plot_heatmap_single.R -d expr.tsv -g groups.txt --column-cluster TRUE\n")
    cat("  Rscript plot_heatmap_single.R -d expr.tsv --show-row-names FALSE --height 8\n")
    cat("\nSample-group file format (one label per line, matching column order):\n")
    cat("  Control\n")
    cat("  Control\n")
    cat("  Treatment\n")
    cat("  Treatment\n")
    quit(status = 0)
}

# ── Validate required arguments ────────────────────────────────────────
if (is.null(opt$data)) {
    stop("Error: -d/--data is required. Use -h for help.")
}
if (!file.exists(opt$data)) {
    stop(paste0("Error: data file not found: ", opt$data))
}

# ── Set defaults ───────────────────────────────────────────────────────
output_prefix    <- if (is.null(opt$output))             "heatmap"      else opt$output
show_row_names   <- if (is.null(opt$`show-row-names`))    TRUE           else opt$`show-row-names`
row_names_size   <- if (is.null(opt$`row-names-fontsize`)) 6             else opt$`row-names-fontsize`
column_cluster   <- if (is.null(opt$`column-cluster`))    FALSE          else opt$`column-cluster`
plot_height      <- if (is.null(opt$height))              10             else opt$height
plot_width       <- if (is.null(opt$width))               10             else opt$width

# ── Read data ──────────────────────────────────────────────────────────
dat <- tryCatch(
    read.table(opt$data, header = TRUE, row.names = 1, sep = "\t",
               quote = "", comment.char = "", check.names = FALSE),
    error = function(e) {
        stop(paste0("Error reading data file: ", e$message))
    }
)

if (ncol(dat) == 0) {
    stop("Error: data matrix has zero columns. Check that the file is tab-delimited with a header.")
}
if (nrow(dat) == 0) {
    stop("Error: data matrix has zero rows.")
}

original_nrow <- nrow(dat)
original_ncol <- ncol(dat)

# ── Scale rows (z-score) ──────────────────────────────────────────────
dat_scaled <- t(scale(t(as.matrix(dat))))

# Handle NA from constant rows (sd = 0 => z-score = NaN)
if (anyNA(dat_scaled)) {
    n_na <- sum(is.na(dat_scaled))
    warning(sprintf("Removed %d NA values (likely from constant rows with zero variance).", n_na))
    dat_scaled <- na.omit(dat_scaled)
}

# Handle Inf values
if (any(is.infinite(dat_scaled))) {
    n_inf <- sum(is.infinite(dat_scaled))
    warning(sprintf("Replaced %d Inf values with NA and removed affected rows.", n_inf))
    dat_scaled[is.infinite(dat_scaled)] <- NA
    dat_scaled <- na.omit(dat_scaled)
}

if (nrow(dat_scaled) == 0) {
    stop("Error: no rows remaining after NA/Inf removal.")
}

# Restore column names (scale can drop them via t())
colnames(dat_scaled) <- colnames(dat)

kept_rows <- nrow(dat_scaled)
if (kept_rows < original_nrow) {
    warning(sprintf("Retained %d of %d rows after scaling and filtering.", kept_rows, original_nrow))
}

# ── Read sample groups (optional) ─────────────────────────────────────
col_annot <- NULL
group_colors <- NULL

if (!is.null(opt$groups)) {
    if (!file.exists(opt$groups)) {
        stop(paste0("Error: group file not found: ", opt$groups))
    }
    groups <- tryCatch(
        readLines(opt$groups),
        error = function(e) {
            stop(paste0("Error reading group file: ", e$message))
        }
    )
    # Strip empty lines
    groups <- groups[groups != ""]
    if (length(groups) != ncol(dat_scaled)) {
        stop(sprintf(
            "Error: group file has %d labels but data matrix has %d columns. They must match.",
            length(groups), ncol(dat_scaled)
        ))
    }
    group_df <- data.frame(Group = factor(groups))
    rownames(group_df) <- colnames(dat_scaled)

    # Assign distinct colours per group
    n_groups <- length(unique(groups))
    if (n_groups <= 9) {
        pal <- brewer.pal(max(3, n_groups), "Set1")[seq_len(n_groups)]
    } else {
        pal <- rainbow(n_groups)
    }
    group_colors <- setNames(pal, levels(group_df$Group))

    col_annot <- HeatmapAnnotation(
        df = group_df,
        col = list(Group = group_colors),
        annotation_name_side = "left",
        annotation_legend_param = list(
            Group = list(title = "Group", title_gp = gpar(fontsize = 10))
        )
    )
}

# ── Colour palette ─────────────────────────────────────────────────────
# Low (#3171AC = blue) — mid (#FFFFFF = white) — high (#D25536 = red)
col_fun <- colorRamp2(
    c(-2, 0, 2),
    c("#3171AC", "#FFFFFF", "#D25536")
)

# ── Build heatmap ──────────────────────────────────────────────────────
hm <- Heatmap(
    dat_scaled,
    name              = "Z-score",
    col               = col_fun,
    show_row_names    = show_row_names,
    row_names_gp      = gpar(fontsize = row_names_size),
    show_column_names = TRUE,
    column_names_gp   = gpar(fontsize = 9),
    cluster_rows      = TRUE,
    cluster_columns   = column_cluster,
    show_row_dend     = TRUE,
    show_column_dend  = column_cluster,
    top_annotation    = col_annot,
    heatmap_legend_param = list(
        title = "Z-score",
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9),
        at = c(-2, -1, 0, 1, 2)
    )
)

# ── Draw to PDF ────────────────────────────────────────────────────────
pdf_path <- paste0(output_prefix, ".pdf")
pdf(file = pdf_path, height = plot_height, width = plot_width)
tryCatch(
    {
        draw(hm)
        cat(sprintf("Heatmap saved to %s\n", pdf_path))
        cat(sprintf("  Rows (genes):   %d / %d retained\n", kept_rows, original_nrow))
        cat(sprintf("  Columns (samples): %d\n", ncol(dat_scaled)))
        if (!is.null(opt$groups)) {
            cat(sprintf("  Groups: %s\n", paste(unique(groups), collapse = ", ")))
        }
    },
    error = function(e) {
        stop(paste0("Error drawing heatmap: ", e$message))
    },
    finally = {
        dev.off()
    }
)
