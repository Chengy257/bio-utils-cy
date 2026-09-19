#!/usr/bin/env Rscript
#########################################################################
# File Name: /mnt/d/myscripts/myscripts/visualization/plot_heatmap_multi.R
# Author: ChengYu
# Description: Generate multi-cluster ComplexHeatmap panels stacked
#              vertically. Each input expression matrix becomes one
#              heatmap panel (cluster), with z-score scaled rows.
#              Produces individual PDFs per panel plus a combined
#              multi-panel PDF. Cluster names are extracted from the
#              second field of each filename when split by ".".
# Created Time: 2026
#########################################################################

suppressMessages(library(getopt))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))

# -- CLI specification ---------------------------------------------------
spec <- matrix(c(
    "file-list",  "f", 1, "character", "Path to a text file listing expression matrix paths, one per line. Required.",
    "groups",     "g", 1, "character", "Path to sample-group file (one label per line, matching column order of each matrix). Required.",
    "output",     "o", 1, "character", "Output prefix for PDF files. Required.",
    "color",      "c", 2, "character", "Color scheme: 'blue_white_red' (#3171AC/#FFFFFF/#D25536) or 'viridis'. Default: blue_white_red.",
    "height",     "H", 2, "numeric",   "PDF height in inches. Default: 12.",
    "width",      "W", 2, "numeric",   "PDF width in inches. Default: 12.",
    "seed",       "S", 2, "numeric",   "Random seed for reproducible clustering. Default: 12345.",
    "help",       "h", 0, "logical",   "Print this help message and exit."
), byrow = TRUE, ncol = 5)
colnames(spec) <- c("long", "short", "argflag", "type", "help")

opt <- getopt(spec)

# -- Help ----------------------------------------------------------------
if (!is.null(opt$help)) {
    cat("Usage: Rscript plot_heatmap_multi.R -f <file_list.txt> -g <groups.txt> -o <prefix> [options]\n\n")
    cat("Generate stacked multi-cluster heatmaps with z-score scaled rows.\n")
    cat("Each input matrix produces its own panel and individual PDF, then all\n")
    cat("panels are combined into a single multi-panel PDF.\n\n")
    cat("Options:\n")
    for (i in seq_len(nrow(spec))) {
        flag <- paste0("-", spec[i, "short"], ", --", spec[i, "long"])
        cat(sprintf("  %-30s %s\n", flag, spec[i, "help"]))
    }
    cat("\nFile list format (one matrix path per line):\n")
    cat("  /path/to/expr.GO0006355.tsv\n")
    cat("  /path/to/expr.GO0006366.tsv\n")
    cat("  /path/to/expr.GO0006351.tsv\n\n")
    cat("Cluster names are extracted from the 2nd field of each filename\n")
    cat("split by '.' (e.g. 'expr.GO0006355.tsv' -> 'GO0006355').\n\n")
    cat("Sample-group file format (one label per line, matching column order):\n")
    cat("  Control\n")
    cat("  Control\n")
    cat("  Treatment\n")
    cat("  Treatment\n\n")
    cat("Examples:\n")
    cat("  Rscript plot_heatmap_multi.R -f filelist.txt -g groups.txt -o result\n")
    cat("  Rscript plot_heatmap_multi.R -f list.txt -g grp.txt -o out --color viridis --height 15\n")
    cat("  Rscript plot_heatmap_multi.R -f list.txt -g grp.txt -o out --seed 42 --width 8\n")
    quit(status = 0)
}

# -- Validate required arguments -----------------------------------------
if (is.null(opt$`file-list`)) {
    stop("Error: -f/--file-list is required. Use -h for help.")
}
if (is.null(opt$groups)) {
    stop("Error: -g/--groups is required. Use -h for help.")
}
if (is.null(opt$output)) {
    stop("Error: -o/--output is required. Use -h for help.")
}
if (!file.exists(opt$`file-list`)) {
    stop(paste0("Error: file list not found: ", opt$`file-list`))
}
if (!file.exists(opt$groups)) {
    stop(paste0("Error: groups file not found: ", opt$groups))
}

# -- Set defaults --------------------------------------------------------
color_scheme <- if (is.null(opt$color)) "blue_white_red" else opt$color
plot_height  <- if (is.null(opt$height)) 12 else opt$height
plot_width   <- if (is.null(opt$width))  12 else opt$width
seed_val     <- if (is.null(opt$seed))   12345 else opt$seed

set.seed(seed_val)

# -- Build colour function -----------------------------------------------
if (color_scheme == "viridis") {
    viridis_cols <- c("#440154", "#21908C", "#FDE725")
    col_fun <- colorRamp2(c(-2, 0, 2), viridis_cols)
} else {
    # blue_white_red (default)
    col_fun <- colorRamp2(c(-2, 0, 2), c("#3171AC", "#FFFFFF", "#D25536"))
}

# -- Read file list ------------------------------------------------------
filepaths <- tryCatch(
    read.table(opt$`file-list`, header = FALSE, stringsAsFactors = FALSE)[, 1],
    error = function(e) {
        stop(paste0("Error reading file list: ", e$message))
    }
)

if (length(filepaths) == 0) {
    stop("Error: file list is empty.")
}

# Validate all input files exist
missing <- filepaths[!file.exists(filepaths)]
if (length(missing) > 0) {
    stop(paste0("Error: the following input files were not found:\n",
                paste0("  ", missing, collapse = "\n")))
}

# -- Read groups ---------------------------------------------------------
groups <- tryCatch(
    {
        lines <- readLines(opt$groups)
        lines <- lines[lines != ""]
        lines
    },
    error = function(e) {
        stop(paste0("Error reading groups file: ", e$message))
    }
)

if (length(groups) == 0) {
    stop("Error: groups file is empty.")
}

# -- Helper: scale a matrix row-wise (z-score) --------------------------
scale_rows <- function(mat) {
    mat <- as.matrix(mat)
    scaled <- t(scale(t(mat)))

    # Replace Inf with NA
    if (any(is.infinite(scaled))) {
        n_inf <- sum(is.infinite(scaled))
        warning(sprintf("Replaced %d Inf values with NA.", n_inf))
        scaled[is.infinite(scaled)] <- NA
    }

    # Omit rows with any NA (constant rows produce NaN)
    if (anyNA(scaled)) {
        n_na <- sum(is.na(scaled))
        warning(sprintf("Removed rows containing %d NA values (likely from constant rows).", n_na))
        scaled <- na.omit(scaled)
    }

    scaled
}

# -- Helper: extract cluster name from filename -------------------------
# Splits basename by "." and returns the second field.
extract_cluster_name <- function(filepath) {
    base <- basename(filepath)
    parts <- strsplit(base, split = "\\.")[[1]]
    if (length(parts) >= 2) {
        return(parts[2])
    }
    # Fallback: use the full basename without extension
    return(tools::file_path_sans_ext(base))
}

# -- Build heatmaps panel by panel ---------------------------------------
ht_list <- NULL
panel_info <- character(0)
n_panels <- length(filepaths)

for (i in seq_len(n_panels)) {
    filepath <- filepaths[i]
    cluster_name <- extract_cluster_name(filepath)

    cat(sprintf("[%d/%d] Processing: %s (cluster: %s)\n", i, n_panels, filepath, cluster_name))

    # Read expression matrix
    dat <- tryCatch(
        read.table(filepath, header = TRUE, row.names = 1, sep = "\t",
                   quote = "", comment.char = "", check.names = FALSE),
        error = function(e) {
            stop(paste0("Error reading matrix '", filepath, "': ", e$message))
        }
    )

    if (ncol(dat) == 0) {
        warning(paste0("Skipping '", filepath, "': zero columns."))
        next
    }
    if (nrow(dat) == 0) {
        warning(paste0("Skipping '", filepath, "': zero rows."))
        next
    }

    original_nrow <- nrow(dat)

    # Validate group count matches column count
    if (length(groups) != ncol(dat)) {
        stop(sprintf(
            "Error: matrix '%s' has %d columns but groups file has %d labels. They must match.",
            filepath, ncol(dat), length(groups)
        ))
    }

    # Scale rows
    dat_scaled <- scale_rows(dat)

    if (nrow(dat_scaled) == 0) {
        warning(paste0("Skipping '", filepath, "': no rows remaining after scaling/filtering."))
        next
    }

    colnames(dat_scaled) <- colnames(dat)
    kept_rows <- nrow(dat_scaled)

    if (kept_rows < original_nrow) {
        cat(sprintf("  Retained %d of %d rows after scaling.\n", kept_rows, original_nrow))
    }

    # Annotations
    col_anno <- HeatmapAnnotation(
        groups = groups,
        annotation_legend_param = list(
            groups = list(title = "Group", title_gp = gpar(fontsize = 10))
        )
    )
    row_anno <- HeatmapAnnotation(
        GOcluster = rep(cluster_name, nrow(dat_scaled)),
        which = "row"
    )

    # Individual heatmap with row names (for single PDF)
    p_single <- Heatmap(
        dat_scaled,
        cluster_rows      = TRUE,
        show_column_names = TRUE,
        show_row_names    = TRUE,
        cluster_columns   = FALSE,
        column_names_rot  = 45,
        top_annotation    = col_anno,
        left_annotation   = row_anno,
        name              = "Expr. z-score",
        row_names_gp      = gpar(fontsize = 6),
        col               = col_fun,
        heatmap_legend_param = list(
            title = "Expr. z-score",
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 9),
            at = c(-2, -1, 0, 1, 2)
        )
    )

    # Write individual PDF
    single_pdf <- paste0(opt$output, "_", cluster_name, "_Single_Heatmap.pdf")
    tryCatch(
        {
            pdf(file = single_pdf, height = plot_height, width = plot_width)
            draw(p_single)
            dev.off()
        },
        error = function(e) {
            if (dev.cur() > 1) dev.off()
            stop(paste0("Error writing individual PDF '", single_pdf, "': ", e$message))
        }
    )
    cat(sprintf("  Individual heatmap saved: %s (%d rows x %d cols)\n",
                single_pdf, kept_rows, ncol(dat_scaled)))

    # Panel for combined plot (no row names to save space)
    # First panel keeps top_annotation with group legend; subsequent panels omit it
    # to avoid duplicate legends
    if (is.null(ht_list)) {
        p_combined <- Heatmap(
            dat_scaled,
            cluster_rows      = TRUE,
            show_column_names = TRUE,
            show_row_names    = FALSE,
            cluster_columns   = FALSE,
            column_names_rot  = 45,
            top_annotation    = col_anno,
            left_annotation   = row_anno,
            name              = "Expr. z-score",
            col               = col_fun
        )
    } else {
        p_combined <- Heatmap(
            dat_scaled,
            cluster_rows      = TRUE,
            show_column_names = TRUE,
            show_row_names    = FALSE,
            cluster_columns   = FALSE,
            column_names_rot  = 45,
            left_annotation   = row_anno,
            name              = paste0("Expr. z-score.", i),
            col               = col_fun,
            show_heatmap_legend = FALSE
        )
    }

    # Stack vertically
    if (is.null(ht_list)) {
        ht_list <- p_combined
    } else {
        ht_list <- ht_list %v% p_combined
    }

    panel_info <- c(panel_info, sprintf("  Panel %d: %s (%d rows)", i, cluster_name, kept_rows))
}

# -- Validate at least one panel was produced ----------------------------
if (is.null(ht_list)) {
    stop("Error: no valid heatmap panels were produced. Check input files and data.")
}

# -- Draw combined multi-panel PDF ---------------------------------------
multi_pdf <- paste0(opt$output, "_Multi_Heatmap.pdf")
tryCatch(
    {
        pdf(file = multi_pdf, height = plot_height, width = plot_width)
        draw(ht_list,
             column_title = "Multi-Cluster Expression Heatmap",
             column_title_gp = gpar(fontsize = 14, fontface = "bold"))
        dev.off()
    },
    error = function(e) {
        if (dev.cur() > 1) dev.off()
        stop(paste0("Error writing combined PDF '", multi_pdf, "': ", e$message))
    }
)

# -- Summary -------------------------------------------------------------
cat(sprintf("\n=== Summary ===\n"))
cat(sprintf("Panels plotted: %d\n", length(panel_info)))
cat(paste(panel_info, collapse = "\n"), "\n")
cat(sprintf("Individual PDFs: %s_*_Single_Heatmap.pdf\n", opt$output))
cat(sprintf("Combined PDF:   %s\n", multi_pdf))
cat(sprintf("Color scheme:   %s\n", color_scheme))
