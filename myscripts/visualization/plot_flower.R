#!/usr/bin/env Rscript
#########################################################################
# File Name: /mnt/d/myscripts/myscripts/visualization/plot_flower.R
# Author: ChengYu
# Description: Flower/petal plot for visualizing set intersection sizes
#              (e.g., Venn-like display for multi-sample gene overlaps).
#              Reads a two-column TSV (sample, count) and produces a PDF
#              with ellipses arranged radially around a central circle.
# Created Time: 2026
#########################################################################

suppressMessages(library(getopt))
suppressMessages(library(plotrix))

# -------------------------------------------------------------------------
# CLI specification
# -------------------------------------------------------------------------
spec <- matrix(c(
    "input",          "i", 1, "character", "Input TSV file (sample<TAB>count, with header), required",
    "output",         "o", 1, "character", "Output PDF file path, required",
    "ellipse-color",  "e", 2, "character", "Petal (ellipse) fill color (hex), default: #D3F2A3FF",
    "center-color",   "c", 2, "character", "Center circle fill color (hex), default: #97E196FF",
    "text-cex",       "x", 2, "numeric",   "Text size multiplier for sample labels, default: 1",
    "start-angle",    "s", 2, "numeric",   "Starting angle in degrees for first petal, default: 90",
    "help",           "h", 0, "logical",   "Show this help message and exit"
), ncol = 5, byrow = TRUE)

opt <- getopt(spec)

# -------------------------------------------------------------------------
# Help
# -------------------------------------------------------------------------
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE), "\n")
    cat("Examples:\n")
    cat("  Rscript plot_flower.R -i samples.tsv -o flower.pdf\n")
    cat("  Rscript plot_flower.R -i samples.tsv -o flower.pdf --ellipse-color '#FF9999' --center-color '#FF5555'\n")
    cat("  Rscript plot_flower.R -i samples.tsv -o flower.pdf --start-angle 0 --text-cex 0.8\n")
    quit(status = 0)
}

# -------------------------------------------------------------------------
# Validate required arguments
# -------------------------------------------------------------------------
if (is.null(opt$input)) {
    stop("Error: Required argument '-i'/'--input' is missing. Use -h for help.")
}
if (is.null(opt$output)) {
    stop("Error: Required argument '-o'/'--output' is missing. Use -h for help.")
}

# Validate input file exists
if (!file.exists(opt$input)) {
    stop(paste0("Error: Input file not found: ", opt$input))
}

# -------------------------------------------------------------------------
# Apply defaults
# -------------------------------------------------------------------------
ellipse_col <- if (is.null(opt$`ellipse-color`))  "#D3F2A3FF" else opt$`ellipse-color`
circle_col  <- if (is.null(opt$`center-color`))   "#97E196FF" else opt$`center-color`
text_cex    <- if (is.null(opt$`text-cex`))        1           else opt$`text-cex`
start_angle <- if (is.null(opt$`start-angle`))     90          else opt$`start-angle`

# -------------------------------------------------------------------------
# Read and validate input data
# -------------------------------------------------------------------------
data <- tryCatch(
    read.table(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               comment.char = "", quote = ""),
    error = function(e) {
        stop(paste0("Error reading input file: ", e$message))
    }
)

if (ncol(data) < 2) {
    stop("Error: Input file must have at least 2 columns (sample name, count).")
}

samples <- data[[1]]
values  <- as.numeric(data[[2]])

if (length(samples) < 2) {
    stop("Error: At least 2 samples are required for a flower plot.")
}

if (any(is.na(values))) {
    missing_idx <- which(is.na(values))
    stop(paste0("Error: Missing or non-numeric values found at row(s): ",
                paste(missing_idx, collapse = ", ")))
}

if (any(is.na(samples) | samples == "")) {
    stop("Error: Empty or missing sample names detected in input.")
}

n <- length(samples)

# -------------------------------------------------------------------------
# Print summary
# -------------------------------------------------------------------------
cat("=== Flower Plot Summary ===\n")
cat(sprintf("  Input file    : %s\n", opt$input))
cat(sprintf("  Output file   : %s\n", opt$output))
cat(sprintf("  Number of samples: %d\n", n))
cat(sprintf("  Petal color   : %s\n", ellipse_col))
cat(sprintf("  Center color  : %s\n", circle_col))
cat(sprintf("  Start angle   : %.1f degrees\n", start_angle))
cat(sprintf("  Text cex      : %.2f\n", text_cex))
cat("  Sample counts:\n")
for (i in seq_len(n)) {
    cat(sprintf("    %-30s : %g\n", samples[i], values[i]))
}
cat("===========================\n")

# -------------------------------------------------------------------------
# Flower plot function
# -------------------------------------------------------------------------
flower_plot <- function(sample, value, start, a, b,
                        ellipse_col = "#D3F2A3FF",
                        circle_col  = "#97E196FF",
                        circle_text_cex = 1) {

    par(bty = "n", ann = FALSE, xaxt = "n", yaxt = "n", mar = c(1, 1, 1, 1))
    plot(c(0, 10), c(0, 10), type = "n")

    n   <- length(sample)
    deg <- 360 / n

    for (t in seq_len(n)) {
        # Draw petal ellipse
        plotrix::draw.ellipse(
            x     = 5 + cos((start + deg * (t - 1)) * pi / 180),
            y     = 5 + sin((start + deg * (t - 1)) * pi / 180),
            col   = ellipse_col,
            border = ellipse_col,
            a     = a,
            b     = b,
            angle = deg * (t - 1)
        )

        # Value label inside petal
        text(
            x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
            y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
            value[t]
        )

        # Sample label outside petal (rotate for readability)
        angle_deg <- deg * (t - 1)
        if (angle_deg < 180 && angle_deg > 0) {
            text(
                x   = 5 + 3.3 * cos((start + angle_deg) * pi / 180),
                y   = 5 + 3.3 * sin((start + angle_deg) * pi / 180),
                sample[t],
                srt = angle_deg - start,
                adj = 1,
                cex = circle_text_cex
            )
        } else {
            text(
                x   = 5 + 3.3 * cos((start + angle_deg) * pi / 180),
                y   = 5 + 3.3 * sin((start + angle_deg) * pi / 180),
                sample[t],
                srt = angle_deg + start,
                adj = 0,
                cex = circle_text_cex
            )
        }
    }

    # Center circle (intersection)
    plotrix::draw.circle(x = 5, y = 5, r = 1,
                         col = circle_col, border = circle_col)
}

# -------------------------------------------------------------------------
# Determine petal dimensions based on number of samples
# -------------------------------------------------------------------------
# Scale ellipse semi-axes so petals do not overlap excessively.
# More samples -> narrower petals.
petal_a <- max(1.2, 2.5 - 0.06 * n)
petal_b <- 0.4

# -------------------------------------------------------------------------
# Produce PDF
# -------------------------------------------------------------------------
pdf_width  <- max(7, 2 + n * 0.5)
pdf_height <- pdf_width

tryCatch({
    pdf(opt$output, width = pdf_width, height = pdf_height)
    flower_plot(
        sample  = samples,
        value   = values,
        start   = start_angle,
        a       = petal_a,
        b       = petal_b,
        ellipse_col    = ellipse_col,
        circle_col     = circle_col,
        circle_text_cex = text_cex
    )
    dev.off()
    cat(sprintf("PDF written to: %s\n", opt$output))
}, error = function(e) {
    if (dev.cur() > 1) dev.off()
    stop(paste0("Error writing PDF: ", e$message))
})
