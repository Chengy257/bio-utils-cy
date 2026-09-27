#!/usr/bin/env Rscript
#########################################################################
# File Name: kegg_enrichment.R
# Author: ChengYu
# Description: KEGG pathway enrichment analysis with clusterProfiler
# Created Time: 2026
#########################################################################

# KEGG Pathway Enrichment Analysis with clusterProfiler
#
# Performs KEGG pathway enrichment on a gene list and produces:
#   - Enrichment result table (XLS)
#   - Dot plot (PDF)
#   - Enrichment network (PDF)
#
# Dependencies:
#   install.packages(c("getopt", "magrittr", "dplyr", "R.utils"))
#   BiocManager::install(c("clusterProfiler", "aPEAR"))
#
# Notes:
#   - The T2G mapping file must be a 2-column TSV: column 1 = KEGG
#     transcript id, column 2 = your gene id.
#   - enrichKEGG is run WITHOUT a custom universe: the background is all
#     KEGG genes of the organism, not your detected-gene set. Provide a
#     universe explicitly (enrichKEGG(universe = ...)) if you need that
#     statistics.
#
# Changelog:
#   v1.1.0  2026-09-27
#   - FIX: the getopt spec declared pvalue/qvalue/padjust/showcat/libpath
#     as flags (3rd column 0) — passing e.g. "-p 0.01" crashed at parse
#     time ("\"0.01\" is not a valid option"), so every documented
#     threshold invocation was unusable; they take values now
#   - FIX: the network plot was never saved — enrichmentNetwork()'s
#     return value was discarded and ggsave() re-used the dot plot (or
#     failed on last_plot()); the ggplot object is now captured and saved
#   - FIX: the package check omitted dplyr and R.utils (both used
#     unconditionally) and required svglite which is never used
#   - FIX: the output-prefix regex ".DEGs.txt$" did not escape its dots
#     (matched e.g. "XDEGs.txt")
#   - FIX: the T2G file must have >= 2 columns — merge() silently produced
#     wrong mappings otherwise; the expected column order is documented
#   - CHANGE: @result slot access replaced with as.data.frame(); the
#     no-universe background choice is documented in the header/help

suppressPackageStartupMessages({
  library(getopt)
})

# -------------------------------------------------------------------------
# Command-line interface
# -------------------------------------------------------------------------

spec <- matrix(c(
  "help",      "h", 0, "logical",   "Show this help message",
  "input",     "i", 1, "character", "Input gene list file (one gene per line)",
  "output",    "o", 1, "character", "Output directory",
  "t2g",       "t", 1, "character", "Transcript-to-Gene mapping file (2-column TSV: col1=KEGG id, col2=gene id)",
  "organism",  "g", 1, "character", "KEGG organism code (e.g., dosa, hsa, mmu)",
  "pvalue",    "p", 1, "numeric",   "P-value cutoff (default: 0.05)",
  "qvalue",    "q", 1, "numeric",   "Q-value cutoff (default: 0.05)",
  "padjust",   "a", 1, "character", "P-adjust method (default: BH)",
  "showcat",   "n", 1, "numeric",   "Number of categories in dot plot (default: 10)",
  "libpath",   "l", 1, "character", "Custom R library path"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

# Help
if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  cat("\n")
  cat("Description:\n")
  cat("  KEGG pathway enrichment analysis using clusterProfiler.\n")
  cat("  Reads a gene list, maps transcripts to KEGG gene IDs,\n")
  cat("  runs enrichKEGG, and outputs an XLS table, dot plot,\n")
  cat("  and enrichment network.\n\n")
  cat("Examples:\n")
  cat("  # Basic usage with rice (dosa)\n")
  cat("  Rscript kegg_enrichment.R \\\n")
  cat("    -i DEGs.txt -o ./kegg_out \\\n")
  cat("    -t KEGG_TranscriptID2GeneIDs \\\n")
  cat("    -g dosa\n\n")
  cat("  # Human KEGG with custom thresholds\n")
  cat("  Rscript kegg_enrichment.R \\\n")
  cat("    -i gene_list.txt -o ./kegg_human \\\n")
  cat("    -t tx2gene.tsv -g hsa \\\n")
  cat("    -p 0.01 -n 15\n\n")
  cat("  # With custom library path\n")
  cat("  Rscript kegg_enrichment.R \\\n")
  cat("    -i genes.txt -o ./out \\\n")
  cat("    -t map.tsv -g mmu \\\n")
  cat("    -l /home/user/R/Rlib_4.2.3\n\n")
  quit(status = 0)
}

# Validate required arguments
if (is.null(opt$input) || is.null(opt$output) || is.null(opt$t2g) || is.null(opt$organism)) {
  stop("Missing required arguments. Use -h for help.")
}

# Defaults
pvalue_cutoff <- if (is.null(opt$pvalue)) 0.05 else opt$pvalue
qvalue_cutoff <- if (is.null(opt$qvalue)) 0.05 else opt$qvalue
padjust_method <- if (is.null(opt$padjust)) "BH" else opt$padjust
show_category <- if (is.null(opt$showcat)) 10 else opt$showcat

input_file   <- opt$input
output_dir   <- opt$output
t2g_file     <- opt$t2g
organism     <- opt$organism

# -------------------------------------------------------------------------
# Library setup
# -------------------------------------------------------------------------

if (!is.null(opt$libpath)) {
  .libPaths(c(opt$libpath, .libPaths()))
}

pkgs <- c("clusterProfiler", "ggplot2", "aPEAR", "dplyr", "R.utils", "magrittr")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Please install it first.", pkg))
  }
}
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
  library(aPEAR)
  library(magrittr)
})

# Use auto download method for clusterProfiler
R.utils::setOption("clusterProfiler.download.method", "auto")

# -------------------------------------------------------------------------
# Prepare output directory
# -------------------------------------------------------------------------

out_prefix <- gsub("\\.DEGs\\.txt$", "", basename(input_file))
out_dir_slash <- paste0(output_dir, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Read data
# -------------------------------------------------------------------------

message("[INFO] Reading gene list: ", input_file)
gene_list <- read.table(input_file, stringsAsFactors = FALSE)[, 1]
message("[INFO] Gene count: ", length(gene_list))

message("[INFO] Reading transcript-to-gene mapping: ", t2g_file)
T2G <- read.table(t2g_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
if (ncol(T2G) < 2) {
  stop("T2G mapping file must have >= 2 columns (col1 = KEGG transcript id, ",
       "col2 = your gene id); got ", ncol(T2G), " column(s) in ", t2g_file)
}

# -------------------------------------------------------------------------
# Map genes to KEGG transcript IDs
# -------------------------------------------------------------------------

tmp <- merge(x = T2G, y = data.frame(V1 = gene_list), by.x = 2, by.y = 1)
trans_list <- unique(tmp[, 2])
message("[INFO] Mapped transcript count: ", length(trans_list))

if (length(trans_list) == 0) {
  stop("No transcripts mapped. Check your gene list and T2G mapping file.")
}

# -------------------------------------------------------------------------
# KEGG enrichment
# -------------------------------------------------------------------------

message("[INFO] Running enrichKEGG for organism: ", organism)
ekegg <- enrichKEGG(
  gene         = trans_list,
  organism     = organism,
  keyType      = "kegg",
  pvalueCutoff = pvalue_cutoff,
  pAdjustMethod = padjust_method,
  qvalueCutoff = qvalue_cutoff
)

if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
  warning("No significant KEGG pathways found.")
  quit(status = 0)
}

# -------------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------------

date_str <- Sys.Date()

# Enrichment table (filtered by p-value)
sig_results <- subset(as.data.frame(ekegg), pvalue <= pvalue_cutoff)
xls_path <- paste0(out_dir_slash, out_prefix, "_EnrichResult_KEGG_", date_str, ".xls")
write.table(
  sig_results,
  file = xls_path,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)
message("[INFO] Saved enrichment table: ", xls_path)

# Dot plot
dot_path <- paste0(out_dir_slash, out_prefix, "_KEGG_", date_str, "_dotplot.pdf")
p <- dotplot(ekegg, showCategory = show_category) + labs(title = "KEGG")
ggsave(filename = dot_path, plot = p, device = "pdf", width = 6, height = 6)
message("[INFO] Saved dot plot: ", dot_path)

# Enrichment network
sig_for_network <- subset(as.data.frame(ekegg), pvalue <= pvalue_cutoff)
if (nrow(sig_for_network) > 0) {
  net_path <- paste0(out_dir_slash, out_prefix, "_KEGG_", date_str, "_Network.pdf")
  net <- aPEAR::enrichmentNetwork(
    sig_for_network,
    colorBy = "p.adjust",
    colorType = "pval",
    drawEllipses = FALSE,
    repelLabels = TRUE,
    verbose = FALSE,
    nodeSize = "Count"
  )
  ggsave(filename = net_path, plot = net, device = "pdf", width = 6, height = 6)
  message("[INFO] Saved network plot: ", net_path)
}

message("[INFO] KEGG enrichment analysis complete.")
