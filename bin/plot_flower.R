#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/plot_flower.R
# Author: ChengYu
# Description: Flower/petal plot for visualizing set intersection sizes
#              (e.g., Venn-like display for multi-sample gene overlaps).
#              Reads a two-column TSV (sample, count) and produces a PDF
#              with ellipses arranged radially around a central circle.
# Created Time: 2026
# Changelog:
#   v1.1.0  2026-09-28
#   - FIX: petal ellipses were rotated by the sector angle only,
#     missing the start angle, so with the default start=90 every
#     petal's long axis lay tangentially (the "flower" rendered as a
#     cross/plus). The ellipse rotation now equals the radial direction
#     (start + sector*(t-1)), so petals point outward.
#   - FIX: sample labels mixed coordinate frames (srt = sector ± start
#     with a branch on the un-rotated sector), rendering some labels
#     upside down / mirrored (e.g. S1/S3 at top and bottom). Labels now
#     sit radially beyond each petal tip, rotated along the radius and
#     flipped on the left half so text is never upside down.
#   - FIX: plot() lacked asp=1 (ellipses/circle distort on non-square
#     devices); added, and the plot range now adapts to the label
#     radius so nothing clips.
#   - FIX: value labels were drawn at a fixed radius that falls
#     outside the petals once there are more than ~20 samples. Labels
#     now sit mid-petal (radius adapts to the petal length).
#   - FIX: negative counts are rejected (intersection sizes must be
#     non-negative); text-cex must be positive.
#   - CHANGE: INFO/summary logs go to stderr via message(); silent
#     dev.off().
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

if (!is.numeric(text_cex) || text_cex <= 0) {
    stop("Error: --text-cex must be a positive number.")
}

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

if (any(values < 0)) {
    stop(paste0("Error: Negative counts are not valid intersection sizes (rows: ",
                paste(which(values < 0), collapse = ", "), ")."))
}

if (any(is.na(samples) | samples == "")) {
    stop("Error: Empty or missing sample names detected in input.")
}

n <- length(samples)

# -------------------------------------------------------------------------
# Print summary
# -------------------------------------------------------------------------
message("=== Flower Plot Summary ===")
message(sprintf("  Input file    : %s", opt$input))
message(sprintf("  Output file   : %s", opt$output))
message(sprintf("  Number of samples: %d", n))
message(sprintf("  Petal color   : %s", ellipse_col))
message(sprintf("  Center color  : %s", circle_col))
message(sprintf("  Start angle   : %.1f degrees", start_angle))
message(sprintf("  Text cex      : %.2f", text_cex))
message("  Sample counts:")
for (i in seq_len(n)) {
    message(sprintf("    %-30s : %g", samples[i], values[i]))
}
message("===========================")

# -------------------------------------------------------------------------
# Flower plot function
# -------------------------------------------------------------------------
flower_plot <- function(sample, value, start, a, b,
                        ellipse_col = "#D3F2A3FF",
                        circle_col  = "#97E196FF",
                        circle_text_cex = 1) {

    par(bty = "n", ann = FALSE, xaxt = "n", yaxt = "n", mar = c(1, 1, 1, 1))

    n   <- length(sample)
    deg <- 360 / n
    r_value <- 1 + a / 2      # mid-petal radius for the value label
    r_label <- 1 + a + 0.4    # sample label radius, beyond the petal tip
    lim <- r_label + 2.6      # keep rotated labels inside the canvas

    plot(5 + c(-lim, lim), 5 + c(-lim, lim), type = "n", asp = 1)

    for (t in seq_len(n)) {
        theta <- start + deg * (t - 1)          # radial direction (deg)
        ct <- cos(theta * pi / 180)
        st <- sin(theta * pi / 180)

        # Petal: centre on the circle edge, long axis along the radius.
        plotrix::draw.ellipse(
            x      = 5 + ct,
            y      = 5 + st,
            col    = ellipse_col,
            border = ellipse_col,
            a      = a,
            b      = b,
            angle  = theta
        )

        # Value label inside the petal
        text(x = 5 + r_value * ct, y = 5 + r_value * st, value[t])

        # Sample label beyond the petal tip, rotated along the radius;
        # flipped on the left half so glyphs stay upright. The anchor
        # sits at the outer end: when the reading direction points
        # inward (against the radius) the text must END at the anchor.
        tn <- theta %% 360
        srt <- tn
        if (tn > 90 && tn < 270) srt <- tn + 180
        rad <- pi / 180
        inward <- cos(tn * rad) * cos(srt * rad) + sin(tn * rad) * sin(srt * rad) < 0
        adj <- if (inward) 1 else 0
        text(x = 5 + r_label * ct, y = 5 + r_label * st,
             sample[t], srt = srt, adj = adj, cex = circle_text_cex)
    }

    # Center circle (intersection) — drawn last to cover petal roots
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
    invisible(dev.off())
    message(sprintf("PDF written to: %s", opt$output))
}, error = function(e) {
    if (dev.cur() > 1) invisible(dev.off())
    stop(paste0("Error writing PDF: ", e$message))
})
