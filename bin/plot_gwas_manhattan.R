#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/plot_gwas_manhattan.R
# Author: ChengYu
# Description: GWAS Manhattan plot, QQ plot, and regional Manhattan plot
#              with gene structure overlay from EMMAX association results.
#              Reads EMMAX output + SNP map file, produces genome-wide
#              Manhattan and QQ plots in PDF, and optionally a regional
#              Manhattan plot with Gviz gene tracks when a BED12 file is
#              supplied.
# Created Time: 2026
#########################################################################

suppressMessages(library(getopt))
suppressMessages(library(ggplot2))
suppressMessages(library(qqman))
suppressMessages(library(data.table))

# -------------------------------------------------------------------------
# CLI argument specification
# -------------------------------------------------------------------------
spec <- matrix(
  c(
    "emmax",       "e", 1, "character", "EMMAX association result file (required)",
    "map",         "m", 1, "character", "Map file with CHR/SNP/BP columns (required)",
    "output",      "o", 1, "character", "Output file prefix (required)",
    "bed",         "b", 2, "character", "BED12 gene file for regional plot (optional)",
    "extension",   "x", 2, "numeric",   "Extension in bp around BED region (default: 50000)",
    "signif-line", "s", 2, "numeric",   "Genome-wide significance threshold (default: 5e-8)",
    "suggest-line","u", 2, "numeric",   "Suggestive significance threshold (default: 1e-5)",
    "width",       "w", 2, "numeric",   "Manhattan/QQ plot width in inches (default: 12)",
    "height",      "H", 2, "numeric",   "Manhattan/QQ plot height in inches (default: 6)",
    "help",        "h", 0, "logical",   "Print this help message and exit"
  ),
  ncol = 5, byrow = TRUE
)

colnames(spec) <- c("long", "short", "argflag", "type", "help")

# Parse arguments
opts <- tryCatch(
  getopt(spec, opt = commandArgs(trailingOnly = TRUE)),
  error = function(e) {
    message("Error parsing arguments: ", e$message)
    quit(status = 1)
  }
)

# Help handler
if (!is.null(opts$help)) {
  cat("Usage: Rscript plot_gwas_manhattan.R -e <emmax> -m <map> -o <prefix> [options]\n\n")
  cat("GWAS Manhattan plot, QQ plot, and optional regional Manhattan plot\n")
  cat("with gene structure from EMMAX results.\n\n")
  cat("Required arguments:\n")
  cat("  -e, --emmax         EMMAX result file (SNP Beta SE P, no header)\n")
  cat("  -m, --map           Map file (CHR SNP C BP, no header)\n")
  cat("  -o, --output        Output file prefix (e.g. mygwas)\n\n")
  cat("Optional arguments:\n")
  cat("  -b, --bed           BED12 gene file for regional Manhattan plot\n")
  cat("  -x, --extension     Extension in bp around BED region (default: 50000)\n")
  cat("  -s, --signif-line   Genome-wide significance threshold (default: 5e-8)\n")
  cat("  -u, --suggest-line  Suggestive significance threshold (default: 1e-5)\n")
  cat("  -w, --width         Plot width in inches (default: 12)\n")
  cat("  -H, --height        Plot height in inches (default: 6)\n")
  cat("  -h, --help          Show this help message\n\n")
  cat("Output files:\n")
  cat("  {prefix}_manhattan.pdf   Genome-wide Manhattan plot\n")
  cat("  {prefix}_qq.pdf          QQ plot of p-values\n")
  cat("  {prefix}_regional.pdf    Regional Manhattan + gene track (if --bed given)\n\n")
  cat("Examples:\n")
  cat("  # Basic Manhattan + QQ\n")
  cat("  Rscript plot_gwas_manhattan.R -e emmax.out -m snp.map -o result\n\n")
  cat("  # With regional plot around a gene region\n")
  cat("  Rscript plot_gwas_manhattan.R -e emmax.out -m snp.map -o result \\\n")
  cat("      -b gene_region.bed -x 100000\n\n")
  cat("  # Custom significance thresholds\n")
  cat("  Rscript plot_gwas_manhattan.R -e emmax.out -m snp.map -o result \\\n")
  cat("      --signif-line 1e-6 --suggest-line 1e-4\n")
  quit(status = 0)
}

# -------------------------------------------------------------------------
# Validate required arguments
# -------------------------------------------------------------------------
if (is.null(opts$emmax)) {
  stop("Missing required argument: -e/--emmax (EMMAX result file). Use -h for help.")
}
if (is.null(opts$map)) {
  stop("Missing required argument: -m/--map (map file). Use -h for help.")
}
if (is.null(opts$output)) {
  stop("Missing required argument: -o/--output (output prefix). Use -h for help.")
}
if (!file.exists(opts$emmax)) {
  stop("EMMAX file not found: ", opts$emmax)
}
if (!file.exists(opts$map)) {
  stop("Map file not found: ", opts$map)
}

# Apply defaults for optional arguments
extension    <- if (is.null(opts$extension))     50000  else opts$extension
signif_line  <- if (is.null(opts$`signif-line`)) 5e-8   else opts$`signif-line`
suggest_line <- if (is.null(opts$`suggest-line`))1e-5   else opts$`suggest-line`
plot_width   <- if (is.null(opts$width))         12     else opts$width
plot_height  <- if (is.null(opts$height))        6      else opts$height
bed_file     <- opts$bed

if (!is.null(bed_file) && !file.exists(bed_file)) {
  stop("BED file not found: ", bed_file)
}

# -------------------------------------------------------------------------
# Helper: read and merge EMMAX + map data
# -------------------------------------------------------------------------
read_data <- function(emmax_file, map_file) {
  message("[INFO] Reading EMMAX file: ", emmax_file)
  emmax_data <- tryCatch(
    fread(emmax_file, header = FALSE),
    error = function(e) stop("Failed to read EMMAX file: ", e$message)
  )
  if (ncol(emmax_data) < 4) {
    stop("EMMAX file must have at least 4 columns (SNP Beta SE P). Found: ", ncol(emmax_data))
  }
  emmax_data <- emmax_data[, 1:4]
  setnames(emmax_data, c("SNP", "Beta", "SE", "P"))

  message("[INFO] Reading map file: ", map_file)
  map_data <- tryCatch(
    fread(map_file, header = FALSE),
    error = function(e) stop("Failed to read map file: ", e$message)
  )
  if (ncol(map_data) < 4) {
    stop("Map file must have at least 4 columns (CHR SNP C BP). Found: ", ncol(map_data))
  }
  map_data <- map_data[, 1:4]
  setnames(map_data, c("CHR", "SNP", "C", "BP"))

  # Normalise chromosome: strip common prefixes (Chr, chr, Chr0, chr0) to bare integer
  map_data$CHR <- as.numeric(gsub("^[Cc]hr0*", "", map_data$CHR))
  map_data$BP  <- as.numeric(map_data$BP)

  merged_data <- merge(emmax_data, map_data, by = "SNP")
  merged_data$P <- as.numeric(merged_data$P)

  # Remove rows with missing or non-finite p-values
  merged_data <- merged_data[!is.na(merged_data$P) & is.finite(merged_data$P), ]

  if (nrow(merged_data) == 0) {
    stop("No valid SNPs remaining after merging EMMAX and map data.")
  }

  return(merged_data)
}

# -------------------------------------------------------------------------
# Helper: normalise chromosome from BED chrom field
# -------------------------------------------------------------------------
normalise_chr <- function(chrom_str) {
  as.numeric(gsub("^[Cc]hr0*", "", chrom_str))
}

# -------------------------------------------------------------------------
# Helper: plot genome-wide Manhattan plot (PDF)
# -------------------------------------------------------------------------
plot_manhattan <- function(data, output_file, signif_thresh, suggest_thresh) {
  message("[INFO] Generating Manhattan plot: ", output_file)
  pdf(output_file, width = plot_width, height = plot_height)
  manhattan(
    data,
    chr  = "CHR",
    bp   = "BP",
    p    = "P",
    snp  = "SNP",
    genomewideline   = -log10(signif_thresh),
    suggestiveline   = -log10(suggest_thresh),
    main = "Manhattan Plot"
  )
  dev.off()
  message("[INFO] Manhattan plot saved.")
}

# -------------------------------------------------------------------------
# Helper: plot QQ plot (PDF)
# -------------------------------------------------------------------------
plot_qq <- function(data, output_file) {
  message("[INFO] Generating QQ plot: ", output_file)
  pdf(output_file, width = plot_height, height = plot_height)
  qq(data$P, main = "QQ Plot")
  dev.off()
  message("[INFO] QQ plot saved.")
}

# -------------------------------------------------------------------------
# Helper: plot regional Manhattan + gene structure (PDF)
# -------------------------------------------------------------------------
plot_regional_manhattan <- function(data, bed_file, extension_val, output_file) {
  message("[INFO] Generating regional Manhattan plot: ", output_file)

  # Load Gviz only when needed
  suppressMessages(library(Gviz))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(cowplot))

  # Read BED12 file (first row)
  bed_data <- tryCatch(
    fread(bed_file, header = FALSE),
    error = function(e) stop("Failed to read BED file: ", e$message)
  )
  if (ncol(bed_data) < 12) {
    stop("BED file must be in BED12 format (12 columns). Found: ", ncol(bed_data))
  }
  bed_row <- bed_data[1, ]
  setnames(bed_row, c("chrom", "start", "end", "name", "score", "strand",
                       "thickStart", "thickEnd", "itemRgb", "blockCount",
                       "blockSizes", "blockStarts"))

  chrom_raw   <- as.character(bed_row$chrom)
  region_start <- as.numeric(bed_row$start)
  region_end   <- as.numeric(bed_row$end)
  region_name  <- as.character(bed_row$name)

  chr_num <- normalise_chr(chrom_raw)
  ext_start <- max(0, region_start - extension_val)
  ext_end   <- region_end + extension_val

  message("[INFO] Region: ", chrom_raw, ":", ext_start, "-", ext_end,
          " (gene: ", region_name, ")")

  # Subset data for the region
  region_data <- data[CHR == chr_num & BP >= ext_start & BP <= ext_end, ]
  if (nrow(region_data) == 0) {
    stop("No data points found in the specified region: ",
         chrom_raw, ":", ext_start, "-", ext_end)
  }
  message("[INFO] SNPs in region: ", nrow(region_data))

  # Regional Manhattan scatter plot
  p1 <- ggplot(region_data, aes(x = BP, y = -log10(P))) +
    geom_point(aes(color = factor(CHR)), alpha = 0.7, size = 1.2) +
    geom_hline(yintercept = -log10(signif_line), color = "red",
               linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = -log10(suggest_line), color = "blue",
               linetype = "dotted", linewidth = 0.5) +
    scale_color_manual(values = rep(c("steelblue", "firebrick"), 22)) +
    theme_bw() +
    labs(
      title = paste0("Regional Manhattan: ", chrom_raw, ":",
                     format(ext_start, big.mark = ","), "-",
                     format(ext_end, big.mark = ",")),
      x = "Base Pair Position",
      y = expression(-log[10](italic(p)))
    ) +
    theme(legend.position = "none")

  # Build Gviz GeneRegionTrack
  # Parse BED12 block structure for gene model
  block_sizes   <- as.numeric(strsplit(as.character(bed_row$blockSizes), ",")[[1]])
  block_starts  <- as.numeric(strsplit(as.character(bed_row$blockStarts), ",")[[1]])
  block_sizes   <- block_sizes[!is.na(block_sizes) & block_sizes > 0]
  block_starts  <- block_starts[!is.na(block_starts)]

  # Construct GRanges with exon structure
  if (length(block_sizes) > 0 && length(block_starts) > 0) {
    exon_starts <- region_start + block_starts
    exon_ends   <- exon_starts + block_sizes - 1
    gr <- GRanges(
      seqnames = chrom_raw,
      ranges   = IRanges(start = exon_starts, end = exon_ends),
      gene     = region_name,
      exon     = paste0("exon_", seq_along(exon_starts))
    )
  } else {
    # Fallback: use the whole region as a single feature
    gr <- GRanges(
      seqnames = chrom_raw,
      ranges   = IRanges(start = region_start, end = region_end),
      gene     = region_name,
      exon     = "exon_1"
    )
  }

  gene_track <- GeneRegionTrack(
    gr,
    name     = "Genes",
    fill     = "darkblue",
    col      = "black",
    transcriptAnnotation = "gene"
  )

  # Render gene track to temporary PDF, then capture as image
  tmp_gene_file <- tempfile(fileext = ".pdf")
  pdf(tmp_gene_file, width = plot_width, height = 3)
  plotTracks(
    list(gene_track),
    from       = ext_start,
    to         = ext_end,
    chromosome = chrom_raw,
    main       = paste("Gene:", region_name)
  )
  dev.off()

  # Read gene track PDF as image via cowplot
  p2 <- ggdraw() + draw_image(tmp_gene_file)

  # Combine plots
  combined_plot <- plot_grid(p1, p2, ncol = 1, rel_heights = c(2, 1))
  ggsave(output_file, plot = combined_plot, width = plot_width, height = 10)

  # Clean up temp file
  unlink(tmp_gene_file)

  message("[INFO] Regional Manhattan plot saved.")
}

# -------------------------------------------------------------------------
# Main logic
# -------------------------------------------------------------------------
main <- function() {
  # Read and merge data
  data <- read_data(opts$emmax, opts$map)
  total_snps <- nrow(data)

  # Count significant hits
  signif_count <- sum(data$P < signif_line, na.rm = TRUE)
  suggest_count <- sum(data$P < suggest_line, na.rm = TRUE)

  message("[INFO] Total SNPs: ", total_snps)
  message("[INFO] Significant hits (P < ", signif_line, "): ", signif_count)
  message("[INFO] Suggestive hits (P < ", suggest_line, "): ", suggest_count)

  # Generate Manhattan plot
  manhattan_file <- paste0(opts$output, "_manhattan.pdf")
  plot_manhattan(data, manhattan_file, signif_line, suggest_line)

  # Generate QQ plot
  qq_file <- paste0(opts$output, "_qq.pdf")
  plot_qq(data, qq_file)

  # Generate regional Manhattan plot if BED file provided
  if (!is.null(bed_file)) {
    regional_file <- paste0(opts$output, "_regional.pdf")
    plot_regional_manhattan(data, bed_file, extension, regional_file)
  }

  # Print summary
  message("========================================")
  message("Summary")
  message("========================================")
  message("Total SNPs analysed:          ", total_snps)
  message("Genome-wide significant (P < ", signif_line, "): ", signif_count)
  message("Suggestive (P < ", suggest_line, "):       ", suggest_count)
  message("Manhattan plot:  ", manhattan_file)
  message("QQ plot:         ", qq_file)
  if (!is.null(bed_file)) {
    message("Regional plot:   ", paste0(opts$output, "_regional.pdf"))
  }
  message("========================================")
  message("Done.")
}

# Run main
tryCatch(
  main(),
  error = function(e) {
    message("ERROR: ", e$message)
    quit(status = 1)
  }
)
