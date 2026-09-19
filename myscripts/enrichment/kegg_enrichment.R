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
#   install.packages(c("getopt", "magrittr"))
#   BiocManager::install(c("clusterProfiler", "aPEAR"))
#   install.packages(c("ggplot2", "svglite"))

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
  "t2g",       "t", 1, "character", "Transcript-to-Gene mapping file (TSV)",
  "organism",  "g", 1, "character", "KEGG organism code (e.g., dosa, hsa, mmu)",
  "pvalue",    "p", 0, "numeric",   "P-value cutoff (default: 0.05)",
  "qvalue",    "q", 0, "numeric",   "Q-value cutoff (default: 0.05)",
  "padjust",   "a", 0, "character", "P-adjust method (default: BH)",
  "showcat",   "n", 0, "numeric",   "Number of categories in dot plot (default: 10)",
  "libpath",   "l", 0, "character", "Custom R library path"
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

pkgs <- c("clusterProfiler", "ggplot2", "aPEAR", "svglite", "magrittr")
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

out_prefix <- gsub(".DEGs.txt$", "", basename(input_file))
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

if (is.null(ekegg) || nrow(ekegg@result) == 0) {
  warning("No significant KEGG pathways found.")
  quit(status = 0)
}

# -------------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------------

date_str <- Sys.Date()

# Enrichment table (filtered by p-value)
sig_results <- ekegg@result %>% dplyr::filter(pvalue <= pvalue_cutoff)
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
sig_for_network <- ekegg@result %>% dplyr::filter(pvalue <= pvalue_cutoff)
if (nrow(sig_for_network) > 0) {
  net_path <- paste0(out_dir_slash, out_prefix, "_KEGG_", date_str, "_Network.pdf")
  aPEAR::enrichmentNetwork(
    sig_for_network,
    colorBy = "p.adjust",
    colorType = "pval",
    drawEllipses = FALSE,
    repelLabels = TRUE,
    verbose = FALSE,
    nodeSize = "Count"
  )
  ggsave(filename = net_path, device = "pdf", width = 6, height = 6)
  message("[INFO] Saved network plot: ", net_path)
}

message("[INFO] KEGG enrichment analysis complete.")
