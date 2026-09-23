#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/plot_upset.R
# Author: ChengYu
# Description: Generate UpSet intersection plots for multi-set
#              comparison. Reads a file listing gene-list paths (one
#              gene ID per line, no header), builds an UpSetR plot from
#              the union of all sets, and writes a PDF.
# Created Time: 2026
#########################################################################

suppressMessages(library(getopt))
suppressMessages(library(UpSetR))

# -- CLI specification ---------------------------------------------------
spec <- matrix(c(
    "input",         "i", 1, "character", "Text file listing gene-list file paths, one per line. Required.",
    "output",        "o", 1, "character", "Output prefix for PDF file. Required.",
    "n-intersects",  "n", 2, "integer",   "Number of intersections to display. Default: 15.",
    "order-by",      "O", 2, "character", "Sort intersections by 'degree' or 'freq'. Default: freq.",
    "plot-width",    "w", 2, "numeric",   "PDF width in inches. Default: 8.",
    "plot-height",   "H", 2, "numeric",   "PDF height in inches. Default: 6.",
    "help",          "h", 0, "logical",   "Print this help message and exit."
), byrow = TRUE, ncol = 5)
colnames(spec) <- c("long", "short", "argflag", "type", "help")

opt <- getopt(spec)

# -- Help ----------------------------------------------------------------
if (!is.null(opt$help)) {
    cat("Usage: Rscript plot_upset.R -i <file_list.txt> -o <prefix> [options]\n\n")
    cat("Generate UpSet intersection plots for multi-set comparison.\n")
    cat("Each gene-list file should contain one gene ID per line with no header.\n\n")
    cat("Options:\n")
    for (i in seq_len(nrow(spec))) {
        flag <- paste0("-", spec[i, "short"], ", --", spec[i, "long"])
        cat(sprintf("  %-30s %s\n", flag, spec[i, "help"]))
    }
    cat("\nInput file list format (one gene-list path per line):\n")
    cat("  /path/to/DEG_condition_A.txt\n")
    cat("  /path/to/DEG_condition_B.txt\n")
    cat("  /path/to/DEG_condition_C.txt\n\n")
    cat("Gene-list file format (one gene ID per line, no header):\n")
    cat("  GENE001\n")
    cat("  GENE002\n")
    cat("  GENE003\n\n")
    cat("Examples:\n")
    cat("  Rscript plot_upset.R -i gene_lists.txt -o result\n")
    cat("  Rscript plot_upset.R -i lists.txt -o out --n-intersects 20 --order-by degree\n")
    cat("  Rscript plot_upset.R -i lists.txt -o out --plot-width 10 --plot-height 8\n")
    quit(status = 0)
}

# -- Validate required arguments -----------------------------------------
if (is.null(opt$input)) {
    stop("Error: -i/--input is required. Use -h for help.")
}
if (is.null(opt$output)) {
    stop("Error: -o/--output is required. Use -h for help.")
}
if (!file.exists(opt$input)) {
    stop(paste0("Error: input file list not found: ", opt$input))
}

# -- Set defaults --------------------------------------------------------
n_intersects <- if (is.null(opt$`n-intersects`)) 15L else as.integer(opt$`n-intersects`)
order_by     <- if (is.null(opt$`order-by`))     "freq" else opt$`order-by`
plot_width   <- if (is.null(opt$`plot-width`))   8 else opt$`plot-width`
plot_height  <- if (is.null(opt$`plot-height`))  6 else opt$`plot-height`

# Validate order-by value
if (!order_by %in% c("degree", "freq")) {
    stop("Error: --order-by must be 'degree' or 'freq'.")
}

# -- Read file list ------------------------------------------------------
filepaths <- tryCatch(
    read.table(opt$input, header = FALSE, stringsAsFactors = FALSE)[, 1],
    error = function(e) {
        stop(paste0("Error reading input file list: ", e$message))
    }
)

if (length(filepaths) == 0) {
    stop("Error: input file list is empty.")
}

# Validate all gene-list files exist
missing <- filepaths[!file.exists(filepaths)]
if (length(missing) > 0) {
    stop(paste0("Error: the following gene-list files were not found:\n",
                paste0("  ", missing, collapse = "\n")))
}

# -- Read gene lists -----------------------------------------------------
gene_lists <- list()
empty_files <- character(0)

for (fp in filepaths) {
    set_name <- tools::file_path_sans_ext(basename(fp))
    genes <- tryCatch(
        {
            lines <- readLines(fp)
            lines <- lines[lines != ""]
            lines
        },
        error = function(e) {
            warning(paste0("Warning: could not read '", fp, "': ", e$message))
            character(0)
        }
    )

    if (length(genes) == 0) {
        empty_files <- c(empty_files, fp)
        warning(paste0("Warning: gene list '", fp, "' is empty. Skipping."))
        next
    }

    gene_lists[[set_name]] <- unique(genes)
    cat(sprintf("  Read %d unique genes from: %s\n", length(gene_lists[[set_name]]), fp))
}

if (length(gene_lists) < 2) {
    stop("Error: at least 2 non-empty gene lists are required for an UpSet plot.")
}

if (length(empty_files) > 0) {
    cat(sprintf("Warning: %d gene-list file(s) were empty and skipped.\n", length(empty_files)))
}

# -- Compute summary statistics ------------------------------------------
all_genes <- unique(unlist(gene_lists))
n_sets <- length(gene_lists)

cat(sprintf("\n=== Summary ===\n"))
cat(sprintf("Number of sets:        %d\n", n_sets))
cat(sprintf("Total unique genes:    %d\n", length(all_genes)))
cat(sprintf("Set sizes:\n"))
for (nm in names(gene_lists)) {
    cat(sprintf("  %-30s %d genes\n", nm, length(gene_lists[[nm]])))
}

# -- Generate UpSet plot -------------------------------------------------
output_pdf <- paste0(opt$output, "_upset_plot.pdf")

tryCatch(
    {
        pdf(file = output_pdf, width = plot_width, height = plot_height)
        upset(
            fromList(gene_lists),
            nsets       = n_sets,
            nintersects = n_intersects,
            order.by    = order_by,
            mb.ratio    = c(0.6, 0.4),
            text.scale  = c(1.3, 1.1, 1.0, 1.0, 1.2, 0.8),
            main.bar.color = "#4477AA",
            sets.bar.color = "#EECC66",
            point.size  = 3.5,
            line.size   = 1.0
        )
        dev.off()
    },
    error = function(e) {
        if (dev.cur() > 1) dev.off()
        stop(paste0("Error writing UpSet plot PDF '", output_pdf, "': ", e$message))
    }
)

cat(sprintf("\nUpSet plot saved: %s\n", output_pdf))
cat(sprintf("Parameters: n-intersects=%d, order-by=%s, width=%.1f, height=%.1f\n",
            n_intersects, order_by, plot_width, plot_height))
