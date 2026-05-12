

# Load necessary libraries
library(ggplot2)
library(rtracklayer)
library(dplyr)

# Function to parse GTF file and plot gene structure for all genes
plot_all_genes_structure <- function(gtf_file) {
  # Read GTF file
  gtf_data <- import(gtf_file)
  
  # Extract relevant features: exons and CDS
  exons <- gtf_data[gtf_data$type == "exon",]
  cds <- gtf_data[gtf_data$type == "CDS",]
  
  # Create a data frame for ggplot
  exons_df <- data.frame(
    gene_id = exons$gene_id,
    start = start(exons),
    end = end(exons),
    type = "exon",
    strand = as.character(strand(exons))
  )
  
  cds_df <- data.frame(
    gene_id = cds$gene_id,
    start = start(cds),
    end = end(cds),
    type = "CDS",
    strand = as.character(strand(cds))
  )
  
  # Calculate introns by finding gaps between exons
  introns_df <- exons_df %>%
    group_by(gene_id, strand) %>%
    arrange(start) %>%
    mutate(
      next_start = lead(start),
      intron_start = end,
      intron_end = next_start
    ) %>%
    filter(!is.na(intron_end)) %>%
    select(gene_id, intron_start, intron_end, strand) %>%
    rename(start = intron_start, end = intron_end) %>%
    mutate(type = "intron")
  
  plot_data <- bind_rows(exons_df, cds_df, introns_df)
  
  # Create a summary data frame for gene labels
  gene_labels <- plot_data %>%
    group_by(gene_id) %>%
    summarise(gene_start = min(start), gene_end = max(end)) %>%
    mutate(midpoint = (gene_start + gene_end) / 2)
  
  # Define arrow types based on strand
  plot_data$arrow_type <- ifelse(plot_data$strand == "+", "last", "first")
  plot_data_intron <- subset(plot_data, type == "intron")
  # Plot gene structures
  p <- ggplot() +
    geom_segment(data = plot_data_intron,
                 aes(x = start, xend = end, y = gene_id, yend = gene_id),
                 arrow = arrow(length = unit(0.1, "inches"), ends = ifelse(plot_data_intron$strand == "+","first","last"), type = "open"), color = "black") +
    
    geom_rect(data = subset(plot_data, type == "exon"),
              aes(xmin = start, xmax = end, ymin = as.numeric(factor(gene_id)) - 0.1, ymax = as.numeric(factor(gene_id)) + 0.1),
              fill = "white", color = "black") +
    geom_rect(data = subset(plot_data, type == "CDS"),
              aes(xmin = start, xmax = end, ymin = as.numeric(factor(gene_id)) - 0.1, ymax = as.numeric(factor(gene_id)) + 0.1),
              fill = "blue") +
    geom_text(data = gene_labels,
              aes(x = midpoint, y = gene_id, label = gene_id),
              vjust = 1.5, hjust = 0.5, size = 3,) +
    theme_void() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14)) 
  print(p)
}

# Example usage
gtf_file <- "~/test.gtf"
plot_all_genes_structure(gtf_file)
