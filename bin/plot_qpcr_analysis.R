#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/plot_qpcr_analysis.R
# Author: ChengYu
# Description: qPCR delta-delta CT analysis and visualization
# Created Time: 2026
#########################################################################

# qPCR Delta-Delta CT Analysis and Visualization
#
# Reads a CSV exported from a qPCR instrument, performs the delta-delta CT
# method for relative quantification, and produces:
#   - Parsed result table (TSV with .xls extension)
#   - 4-panel PDF figure (relative expression, CT, delta CT, reference gene CT)
#
# Dependencies:
#   install.packages(c("getopt", "ggplot2", "ggsci", "patchwork"))

suppressPackageStartupMessages({
  library(getopt)
})

# -------------------------------------------------------------------------
# Command-line interface
# -------------------------------------------------------------------------

spec <- matrix(c(
  "help",       "h", 0, "logical",   "Show this help message",
  "data",       "d", 1, "character", "Input CSV data exported from qPCR instrument (required)",
  "ref",        "r", 1, "character", "Reference gene name, e.g. ACTIN or GAPDH (required)",
  "control",    "c", 1, "character", "Control sample name, e.g. WT (required)",
  "output",     "o", 0, "character", "Output prefix (default: input filename without extension)",
  "skip-rows",  "s", 0, "numeric",   "Number of header rows to skip at top (default: 34)",
  "tail-rows",  "t", 0, "numeric",   "Number of rows to trim at bottom (default: 5)",
  "ct-max",     "m", 0, "numeric",   "CT value to assign Undetermined entries (default: 45)",
  "plot-width", "W", 0, "numeric",   "Plot width in inches (default: 10)",
  "plot-height","H", 0, "numeric",   "Plot height in inches (default: 10)"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

# Help
if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  cat("\n")
  cat("Description:\n")
  cat("  qPCR delta-delta CT analysis and visualization.\n")
  cat("  Reads a CSV exported from a qPCR instrument, computes delta CT,\n")
  cat("  delta-delta CT, and relative expression (2^-ddCT), then generates\n")
  cat("  a 4-panel PDF figure and a result table.\n\n")
  cat("Examples:\n")
  cat("  # Basic usage\n")
  cat("  Rscript plot_qpcr_analysis.R \\\n")
  cat("    -d qPCR_data.csv -r ACTIN -c WT\n\n")
  cat("  # Custom output and parameters\n")
  cat("  Rscript plot_qpcr_analysis.R \\\n")
  cat("    -d qPCR_data.csv -r GAPDH -c Control \\\n")
  cat("    -o my_experiment --skip-rows 30 --tail-rows 3 --ct-max 40\n\n")
  cat("  # Custom plot dimensions\n")
  cat("  Rscript plot_qpcr_analysis.R \\\n")
  cat("    -d data.csv -r ACTIN -c WT \\\n")
  cat("    --plot-width 12 --plot-height 8\n\n")
  quit(status = 0)
}

# Validate required arguments
if (is.null(opt$data) || is.null(opt$ref) || is.null(opt$control)) {
  stop("Missing required arguments (-d, -r, -c). Use -h for help.")
}

# -------------------------------------------------------------------------
# Parameters and defaults
# -------------------------------------------------------------------------

input_file  <- opt$data
ref_gene    <- opt$ref
sample_control <- opt$control
skip_rows   <- if (is.null(opt[["skip-rows"]]))  34 else opt[["skip-rows"]]
tail_rows   <- if (is.null(opt[["tail-rows"]]))   5 else opt[["tail-rows"]]
ct_max      <- if (is.null(opt[["ct-max"]]))     45 else opt[["ct-max"]]
plot_width  <- if (is.null(opt[["plot-width"]]))  10 else opt[["plot-width"]]
plot_height <- if (is.null(opt[["plot-height"]])) 10 else opt[["plot-height"]]

# Output prefix: use provided or derive from input filename
if (is.null(opt$output)) {
  out_prefix <- tools::file_path_sans_ext(basename(input_file))
} else {
  out_prefix <- opt$output
}

# -------------------------------------------------------------------------
# Load packages
# -------------------------------------------------------------------------

pkgs <- c("ggplot2", "ggsci", "patchwork")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Please install it first.", pkg))
  }
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

# -------------------------------------------------------------------------
# Read and parse input data
# -------------------------------------------------------------------------

message("[INFO] Reading input file: ", input_file)
if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

dat <- tryCatch(
  read.csv(input_file, skip = skip_rows, stringsAsFactors = FALSE),
  error = function(e) {
    stop("Failed to read CSV file: ", e$message)
  }
)

if (nrow(dat) <= tail_rows) {
  stop("Not enough data rows after skipping ", skip_rows, " header rows.")
}

# Keep columns 4 (sample), 5 (gene), 15 (CT)
if (ncol(dat) < 15) {
  stop("Expected at least 15 columns in the data (after skipping header). Found: ", ncol(dat))
}

dat <- dat[1:(nrow(dat) - tail_rows), c(4, 5, 15)]
colnames(dat) <- c("sample", "gene", "CT")

# Replace "Undetermined" with ct_max and convert to numeric
dat$CT[dat$CT == "Undetermined"] <- ct_max
dat$CT <- as.numeric(dat$CT)

# Report counts
all_samples <- unique(dat$sample)
all_genes   <- unique(dat$gene)
message("[INFO] Samples found: ", paste(all_samples, collapse = ", "))
message("[INFO] Genes found: ", paste(all_genes, collapse = ", "))
message("[INFO] Total rows: ", nrow(dat))

# Validate reference gene
if (!(ref_gene %in% all_genes)) {
  stop("Reference gene '", ref_gene, "' not found in data. Available genes: ",
       paste(all_genes, collapse = ", "))
}

# Validate control sample
control_matches <- unique(dat$sample[grep(sample_control, dat$sample)])
if (length(control_matches) == 0) {
  stop("Control sample '", sample_control, "' not found in data. Available samples: ",
       paste(all_samples, collapse = ", "))
}

# -------------------------------------------------------------------------
# Calculate delta CT (target CT - reference gene CT per sample)
# -------------------------------------------------------------------------

dat_genes <- data.frame()
for (sample_name in unique(dat$sample)) {
  detect_gene <- unique(dat$gene[dat$sample == sample_name])
  dat_ref <- dat[dat$sample == sample_name & dat$gene == ref_gene, ]
  for (gene in detect_gene) {
    if (gene != ref_gene) {
      dat_gene <- dat[dat$sample == sample_name & dat$gene == gene, ]
      dat_gene <- cbind(dat_gene, dat_ref)
      dat_genes <- rbind(dat_genes, dat_gene)
    }
  }
}
# dat_genes columns: sample, gene, CT (target), sample.1, gene.1, CT (ref)
dat_genes$deltaCT <- dat_genes[, 3] - dat_genes[, 6]

# -------------------------------------------------------------------------
# Calculate delta-delta CT and relative expression
# -------------------------------------------------------------------------

genes <- unique(dat_genes[, 2])
dat_samples <- data.frame()
for (gene in genes) {
  sub_dat <- dat_genes[dat_genes[, 2] == gene, c(1, 2, 3, 7)]
  control_name <- unique(sub_dat$sample[grep(sample_control, sub_dat$sample)])
  if (length(control_name) == 0) {
    warning("No control sample found for gene: ", gene, ". Skipping.")
    next
  }
  dat_control <- sub_dat[sub_dat$sample == control_name, ]
  for (sample_name in unique(sub_dat$sample)) {
    dat_sample <- sub_dat[sub_dat$sample == sample_name, ]
    dat_sample <- cbind(dat_sample, dat_control)
    dat_samples <- rbind(dat_samples, dat_sample)
  }
}
# dat_samples columns: sample, gene, CT, deltaCT, sample.1, gene.1, CT.1, deltaCT.1
dat_samples$deltadeltaCT <- dat_samples[, 4] - dat_samples[, 8]
dat_samples$rel_exp <- 2^(-dat_samples$deltadeltaCT)

# -------------------------------------------------------------------------
# Save result table
# -------------------------------------------------------------------------

xls_path <- paste0(out_prefix, "_parsed_result_table.xls")
write.table(
  dat_samples,
  file = xls_path,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
message("[INFO] Saved result table: ", xls_path)

# -------------------------------------------------------------------------
# Print summary
# -------------------------------------------------------------------------

message("[INFO] === Summary ===")
message("[INFO] Genes analyzed: ", paste(unique(dat_samples[, 2]), collapse = ", "))
message("[INFO] Samples: ", paste(unique(dat_samples[, 1]), collapse = ", "))
message("[INFO] Reference gene: ", ref_gene)
message("[INFO] Control sample: ", control_matches)

# -------------------------------------------------------------------------
# Shared ggplot theme for all panels
# -------------------------------------------------------------------------

shared_theme <- theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "null",
    text = element_text(family = "sans", face = "bold"),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.placement = "outside",
    panel.border = element_rect(size = 0.8)
  )

# -------------------------------------------------------------------------
# Panel 1: Relative Expression
# -------------------------------------------------------------------------

p1 <- ggplot(dat_samples[, c(1, 2, 10)], aes(x = sample, y = rel_exp, fill = gene)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean") +
  geom_jitter(size = 1.5) +
  stat_summary(
    fun = mean, geom = "errorbar",
    fun.max = function(x) mean(x) + sd(x),
    fun.min = function(x) mean(x) - sd(x),
    width = 0.3
  ) +
  facet_wrap(~gene, scales = "free") +
  shared_theme +
  labs(x = NULL, y = "Relative Expression", title = "Relative Expression") +
  scale_fill_npg()

# -------------------------------------------------------------------------
# Panel 2: CT Values
# -------------------------------------------------------------------------

p2 <- ggplot(dat_samples[, c(1, 2, 3)], aes(x = sample, y = CT, fill = gene)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean") +
  geom_jitter(size = 1.5) +
  stat_summary(
    fun = mean, geom = "errorbar",
    fun.max = function(x) mean(x) + sd(x),
    fun.min = function(x) mean(x) - sd(x),
    width = 0.3
  ) +
  facet_wrap(~gene, scales = "free") +
  shared_theme +
  labs(x = NULL, y = "CT value", title = "CT value") +
  scale_fill_npg()

# -------------------------------------------------------------------------
# Panel 3: Delta CT Values
# -------------------------------------------------------------------------

p3 <- ggplot(dat_samples[, c(1, 2, 4)], aes(x = sample, y = deltaCT, fill = gene)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean") +
  geom_jitter(size = 1.5) +
  stat_summary(
    fun = mean, geom = "errorbar",
    fun.max = function(x) mean(x) + sd(x),
    fun.min = function(x) mean(x) - sd(x),
    width = 0.3
  ) +
  facet_wrap(~gene, scales = "free") +
  shared_theme +
  labs(x = NULL, y = "deltaCT value", title = "deltaCT value") +
  scale_fill_npg()

# -------------------------------------------------------------------------
# Panel 4: Reference Gene CT Values
# -------------------------------------------------------------------------

p4 <- ggplot(dat_genes[, c(1, 2, 6)], aes(x = sample, y = CT, fill = gene)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean") +
  geom_jitter(size = 1.5) +
  stat_summary(
    fun = mean, geom = "errorbar",
    fun.max = function(x) mean(x) + sd(x),
    fun.min = function(x) mean(x) - sd(x),
    width = 0.3
  ) +
  facet_wrap(~gene, scales = "free") +
  shared_theme +
  labs(x = NULL, y = "CT value", title = paste0(ref_gene, " CT value")) +
  scale_fill_npg()

# -------------------------------------------------------------------------
# Assemble and save 4-panel figure
# -------------------------------------------------------------------------

combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = paste0("Figure. QPCR Parsed Results Barplots.      [", date(), "]"),
    tag_levels = "A",
    caption = "Generated by plot_qpcr_analysis.R\nAuthor: ChengYu"
  )

pdf_path <- paste0(out_prefix, "_parsed_result_barplot.pdf")
ggsave(pdf_path, plot = combined, width = plot_width, height = plot_height)
message("[INFO] Saved 4-panel figure: ", pdf_path)
message("[INFO] Done.")
