#!/usr/bin/
#########################################################################
# File Name: plot.R
# Author: ChengYu
# Description: 
# Created Time: Fri 07 Jun 2024 04:58:54 PM CST
#########################################################################
# Load necessary libraries
library(ggplot2)
library(qqman)
library(data.table)
library(Gviz)
library(optparse)
library(cowplot)

# Function to read EMMAX result and map files
read_data <- function(emmax_file, map_file) {
  emmax_data <- fread(emmax_file, header=FALSE)
  setnames(emmax_data, c("SNP", "Beta", "SE", "P"))
  
  map_data <- fread(map_file, header=FALSE)
  setnames(map_data, c("CHR", "SNP", "C", "BP"))
  
  merged_data <- merge(emmax_data, map_data, by="SNP")
  return(merged_data)
}

# Function to plot QQ plot
plot_qq <- function(data, output_file) {
  png(output_file)
  qq(data$P, main = "QQ Plot")
  dev.off()
}

# Function to plot Manhattan plot
plot_manhattan <- function(data, output_file) {
  png(output_file)
  manhattan(data, chr="CHR", bp="BP", p="P", snp="SNP", main = "Manhattan Plot")
  dev.off()
}

plot_regional_manhattan <- function(data, bed_file, extension, output_file) {
  # Read bed file for gene structure
  bed_data <- fread(bed_file, header=FALSE)
  setnames(bed_data, c("chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"))
  
  # Extend the region by specified upstream and downstream distances
  bed_data$extended_start <- bed_data$start - extension
  bed_data$extended_end <- bed_data$end + extension

  # Debugging: print bed_data
  print(bed_data)

  # Subset data for the specified region
  chr_num <- as.numeric(gsub("Chr", "", bed_data$chrom))
  region_data <- subset(data, CHR == chr_num & BP >= bed_data$extended_start & BP <= bed_data$extended_end)

  # Debugging: print region_data
  print(region_data)

  if (nrow(region_data) == 0) {
    stop("No data points found in the specified region.")
  }

  # Create regional Manhattan plot using ggplot2
  p1 <- ggplot(region_data, aes(x = BP, y = -log10(P))) +
    geom_point(aes(color = factor(CHR))) +
    scale_color_manual(values = rep(c("blue", "red"), 22)) +
    theme_bw() +
    labs(title = paste("Regional Manhattan Plot: Chr", bed_data$chrom, ":", bed_data$extended_start, "-", bed_data$extended_end),
         x = "Base Pair Position",
         y = "-log10(p-value)")

  # Convert bed_data to GRanges object
  gr <- GRanges(seqnames = bed_data$chrom,
                ranges = IRanges(start = bed_data$start, end = bed_data$end),
                gene = bed_data$name)

  # Create a GeneRegionTrack
  gene_track <- GeneRegionTrack(gr, name = "Genes")

  # Plot the gene structure using Gviz and save as temporary file
  tmp_gene_file <- tempfile(fileext = ".png")
  png(tmp_gene_file, width = 800, height = 400)
  plotTracks(list(gene_track), from=bed_data$extended_start, to=bed_data$extended_end, chromosome=bed_data$chrom, main="Gene Structure")
  dev.off()

  # Read the gene structure plot
  p2 <- ggdraw() + draw_image(tmp_gene_file)

  # Combine the regional Manhattan plot and gene structure plot
  combined_plot <- plot_grid(p1, p2, ncol = 1, rel_heights = c(2, 1))

  # Save the combined plot as PDF
  ggsave(output_file, plot = combined_plot, width = 10, height = 8)
}

# Main function to execute the plotting
main <- function(emmax_file, map_file, bed_file, extension, qq_output, manhattan_output, regional_output) {
  data <- read_data(emmax_file, map_file)
  
  # Plot QQ plot
  plot_qq(data, qq_output)
  
  # Plot Manhattan plot
  plot_manhattan(data, manhattan_output)
  
  # Plot regional Manhattan plot with gene structure
  plot_regional_manhattan(data, bed_file, extension, regional_output)
}

# Command line argument parsing
option_list <- list(
  make_option(c("-e", "--emmax"), type="character", default=NULL, help="EMMAX result file in .ps format", metavar="character"),
  make_option(c("-m", "--map"), type="character", default=NULL, help="Map file", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default=NULL, help="BED12 file for gene structure", metavar="character"),
  make_option(c("-x", "--extension"), type="integer", default=0, help="Extension distance for upstream and downstream", metavar="integer"),
  make_option(c("-q", "--qq_output"), type="character", default="qq_plot.png", help="Output file for QQ plot", metavar="character"),
  make_option(c("-o", "--manhattan_output"), type="character", default="manhattan_plot.png", help="Output file for Manhattan plot", metavar="character"),
  make_option(c("-r", "--regional_output"), type="character", default="regional_manhattan_plot.png", help="Output file for regional Manhattan plot with gene structure", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$emmax) || is.null(opt$map) || is.null(opt$bed)) {
  print_help(opt_parser)
  stop("EMMAX file, map file, and BED12 file are required.\n", call.=FALSE)
}

# Execute main function with parsed arguments
main(opt$emmax, opt$map, opt$bed, opt$extension, opt$qq_output, opt$manhattan_output, opt$regional_output)
