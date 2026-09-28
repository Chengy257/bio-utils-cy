#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/plot_qpcr_analysis.R
# Author: ChengYu
# Description: qPCR delta-delta CT analysis and visualization
# Created Time: 2026
# Changelog:
#   v1.1.0  2026-09-28
#   - FIX: --output/--skip-rows/--tail-rows/--ct-max/--plot-width/
#     --plot-height were declared as boolean flags (argflag 0) while
#     taking values, so every documented example that passed a value
#     aborted in getopt, and a bare flag was silently read as TRUE
#     (e.g. read.csv(skip = TRUE)). They now require a value; omitted
#     means the code-side default.
#   - FIX: the control sample was located with regex grep and then
#     compared with a vector equality (`sample == control_name`) that
#     recycles when several replicates match: rows paired arbitrarily
#     and the control's own replicates got relative expression != 1
#     (confirmed in a rendered figure). The control is now matched
#     exactly, or as a replicate-group prefix ("WT" -> "WT_rep1"...),
#     and delta-delta CT uses the control group's mean delta CT — the
#     standard ddCt definition — so no row pairing happens at all.
#   - FIX: Panel D ("reference gene CT") was faceted/filled by the
#     TARGET gene column (col 2) instead of the reference gene (col 5);
#     it rendered only because R partial-matched y=CT to CT.1. It now
#     uses the reference gene column explicitly.
#   - FIX: replicate count mismatches between target and reference gene
#     wells were silently recycled by cbind; merge() now warns when a
#     sample's well counts differ.
#   - FIX: legend.position "null" (invalid) -> "none"; element_rect
#     size= -> linewidth=; the run timestamp (date()) was dropped from
#     the figure title so output bytes are reproducible.
#   - CHANGE: result table extension .xls -> .tsv (always was TSV);
#     columns are now named (sample, gene, CT, deltaCT,
#     control_mean_deltaCT, deltadeltaCT, rel_exp) instead of being
#     addressed positionally.
#   - CHANGE: "Undetermined" CT handling is case-insensitive and
#     non-numeric CT leftovers are counted and reported.
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
  "output",     "o", 1, "character", "Output prefix (default: input filename without extension)",
  "skip-rows",  "s", 1, "numeric",   "Number of header rows to skip at top (default: 34)",
  "tail-rows",  "t", 1, "numeric",   "Number of rows to trim at bottom (default: 5)",
  "ct-max",     "m", 1, "numeric",   "CT value to assign Undetermined entries (default: 45)",
  "plot-width", "W", 1, "numeric",   "Plot width in inches (default: 10)",
  "plot-height","H", 1, "numeric",   "Plot height in inches (default: 10)"
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

# Replace "Undetermined" (any case) with ct_max and convert to numeric;
# any other non-numeric CT is counted and reported instead of silently
# becoming NA.
dat$CT <- as.character(dat$CT)
undetermined <- toupper(dat$CT) == "UNDETERMINED"
dat$CT[undetermined] <- as.character(ct_max)
dat$CT <- suppressWarnings(as.numeric(dat$CT))
n_bad_ct <- sum(is.na(dat$CT) & !undetermined)
if (n_bad_ct > 0) {
  warning(sprintf("%d row(s) have non-numeric CT values (set to NA): %s",
                  n_bad_ct,
                  paste(head(dat[is.na(dat$CT) & !undetermined, "sample"], n = 5),
                        collapse = ", ")))
}
if (all(is.na(dat$CT))) {
  stop("No numeric CT values found in column 15 after parsing.")
}

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

# Validate control sample: exact name first, then a replicate-group
# prefix ("WT" matches "WT_rep1", "WT_rep2", ...) as a fixed prefix.
# The old regex grep also matched unrelated names like "nWT".
if (sample_control %in% all_samples) {
  control_samples <- sample_control
} else {
  control_samples <- all_samples[startsWith(as.character(all_samples),
                                            paste0(sample_control, "_")) |
                                   startsWith(as.character(all_samples),
                                              paste0(sample_control, "-"))]
}
if (length(control_samples) == 0) {
  stop("Control sample '", sample_control, "' not found in data. Available samples: ",
       paste(all_samples, collapse = ", "))
}
message("[INFO] Control samples: ", paste(control_samples, collapse = ", "))

# -------------------------------------------------------------------------
# Calculate delta CT (target CT - reference gene CT per sample)
# -------------------------------------------------------------------------

ref_rows <- dat[dat$gene == ref_gene, c("sample", "CT")]
names(ref_rows) <- c("sample", "ref_CT")
targets <- dat[dat$gene != ref_gene, ]
if (nrow(targets) == 0) {
  stop("No target genes found besides the reference gene '", ref_gene, "'.")
}

# Warn when a sample has different numbers of target vs reference wells
# for the same gene (merge pairs all combinations; the old cbind
# silently recycled rows in that case).
n_ref <- table(ref_rows$sample)
tgt_grid <- as.data.frame(table(targets$sample, targets$gene),
                          stringsAsFactors = FALSE)
names(tgt_grid) <- c("sample", "gene", "n_wells")
tgt_grid$ref_wells <- as.integer(n_ref[as.character(tgt_grid$sample)])
bad_wells <- tgt_grid[tgt_grid$n_wells > 0 & tgt_grid$n_wells != tgt_grid$ref_wells, ]
if (nrow(bad_wells) > 0) {
  warning(paste0(
    "Sample(s) with differing target/reference well counts (all ",
    "combinations paired): ",
    paste(sprintf("%s/%s (%d target vs %d reference)",
                  bad_wells$sample, bad_wells$gene,
                  bad_wells$n_wells, bad_wells$ref_wells),
          collapse = ", ")))
}
no_ref <- setdiff(unique(targets$sample), names(n_ref))
if (length(no_ref) > 0) {
  stop("Sample(s) with target wells but no reference-gene wells: ",
       paste(no_ref, collapse = ", "))
}

dat_genes <- merge(targets, ref_rows, by = "sample")
dat_genes$deltaCT <- dat_genes$CT - dat_genes$ref_CT

# -------------------------------------------------------------------------
# Calculate delta-delta CT and relative expression
# -------------------------------------------------------------------------
# Standard ddCt: each delta CT is compared with the MEAN delta CT of the
# control group for the same gene. No row pairing is involved, so the
# control's own replicates all sit at relative expression ~ 1 (the old
# positional cbind gave them arbitrary values when several replicates
# matched).

control_mean <- aggregate(deltaCT ~ gene, data = dat_genes,
                          subset = sample %in% control_samples, FUN = mean)
names(control_mean)[2] <- "control_mean_deltaCT"

genes_no_control <- setdiff(unique(dat_genes$gene), control_mean$gene)
if (length(genes_no_control) > 0) {
  warning("No control sample found for gene(s): ",
          paste(genes_no_control, collapse = ", "), " - skipped.")
  dat_genes <- dat_genes[!(dat_genes$gene %in% genes_no_control), ]
}

dat_samples <- merge(dat_genes, control_mean, by = "gene")
dat_samples$deltadeltaCT <- dat_samples$deltaCT - dat_samples$control_mean_deltaCT
dat_samples$rel_exp <- 2^(-dat_samples$deltadeltaCT)

# -------------------------------------------------------------------------
# Save result table
# -------------------------------------------------------------------------

tsv_path <- paste0(out_prefix, "_parsed_result_table.tsv")
write.table(
  dat_samples,
  file = tsv_path,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
message("[INFO] Saved result table: ", tsv_path)

# -------------------------------------------------------------------------
# Print summary
# -------------------------------------------------------------------------

message("[INFO] === Summary ===")
message("[INFO] Genes analyzed: ", paste(unique(dat_samples[, "gene"]), collapse = ", "))
message("[INFO] Samples: ", paste(unique(dat_samples[, "sample"]), collapse = ", "))
message("[INFO] Reference gene: ", ref_gene)
message("[INFO] Control sample(s): ", paste(control_samples, collapse = ", "))

# -------------------------------------------------------------------------
# Shared ggplot theme for all panels
# -------------------------------------------------------------------------

shared_theme <- theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    text = element_text(family = "sans", face = "bold"),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.placement = "outside",
    panel.border = element_rect(linewidth = 0.8)
  )

# -------------------------------------------------------------------------
# Panel 1: Relative Expression
# -------------------------------------------------------------------------

p1 <- ggplot(dat_samples[, c("sample", "gene", "rel_exp")], aes(x = sample, y = rel_exp, fill = gene)) +
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

p2 <- ggplot(dat_samples[, c("sample", "gene", "CT")], aes(x = sample, y = CT, fill = gene)) +
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

p3 <- ggplot(dat_samples[, c("sample", "gene", "deltaCT")], aes(x = sample, y = deltaCT, fill = gene)) +
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

p4_data <- ref_rows
p4_data$gene <- ref_gene
p4 <- ggplot(p4_data, aes(x = sample, y = ref_CT, fill = gene)) +
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
    title = "Figure. QPCR Parsed Results Barplots.",
    tag_levels = "A",
    caption = "Generated by plot_qpcr_analysis.R\nAuthor: ChengYu"
  )

pdf_path <- paste0(out_prefix, "_parsed_result_barplot.pdf")
ggsave(pdf_path, plot = combined, width = plot_width, height = plot_height)
message("[INFO] Saved 4-panel figure: ", pdf_path)
message("[INFO] Done.")
