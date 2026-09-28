#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/plot_gene_structure.R
# Author: ChengYu
# Description: Gene structure visualization from GTF files. Plots exons,
#              CDS regions, and introns (with strand-direction arrows)
#              for specified genes in a clean publication-ready PDF.
# Created Time: 2026
# Changelog:
#   v1.1.0  2026-09-28
#   - FIX: a GTF without a "gene_id" attribute crashed with the cryptic
#     "arguments imply differing number of rows" while building the
#     exon data frame. The attribute is now checked up front with a
#     clear message.
#   - FIX: GTFs with no CDS features crashed the same way (a constant
#     "CDS" column cannot be mixed with zero-length columns in the
#     empty data.frame); the empty case now builds all zero-row
#     columns.
#   - FIX: --genes entries are trimmed and empty tokens dropped, so
#     "G1, G2" no longer fails to match " G2".
#   - FIX: the multi-chromosome check fired whenever two different
#     genes sat on different chromosomes (harmless - rows are
#     independent) and skipped the real hazard. It now warns only when
#     a single gene's own exons span more than one chromosome, which
#     genuinely mixes coordinates on one row.
#   - CHANGE: intron arrows are drawn in two vectorised layers
#     (plus/minus strand) instead of one ggplot layer per intron.
#   - CHANGE: dev.off() is silent.
#########################################################################

suppressMessages(library(ggplot2))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
suppressMessages(library(getopt))

# -------------------------------------------------------------------------
# CLI argument specification
# -------------------------------------------------------------------------
spec <- matrix(
  c(
    "gtf",        "g", 1, "character", "Path to input GTF file (required)",
    "output",     "o", 1, "character", "Path to output PDF file (required)",
    "genes",      "s", 2, "character", "Comma-separated gene IDs to plot (default: all genes)",
    "exon-color", "e", 2, "character", "Fill color for exon rectangles (default: white)",
    "cds-color",  "c", 2, "character", "Fill color for CDS rectangles (default: blue)",
    "intron-color","n", 2, "character", "Color for intron lines (default: black)",
    "plot-width", "w", 2, "numeric",   "Plot width in inches (default: 10)",
    "plot-height","H", 2, "numeric",   "Plot height in inches (default: auto from gene count)",
    "help",       "h", 0, "logical",   "Print this help message and exit"
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
  cat("Usage: Rscript plot_gene_structure.R -g <gtf> -o <output.pdf> [options]\n\n")
  cat("Gene structure visualization from GTF files.\n\n")
  cat("Plots exons, CDS regions, and introns (with strand-direction arrows)\n")
  cat("for specified genes in a clean publication-ready PDF.\n\n")
  cat("Required arguments:\n")
  cat("  -g, --gtf         GTF file path\n")
  cat("  -o, --output      Output PDF file path\n\n")
  cat("Optional arguments:\n")
  cat("  -s, --genes       Comma-separated gene IDs to plot (default: all)\n")
  cat("  -e, --exon-color  Exon fill color (default: white)\n")
  cat("  -c, --cds-color   CDS fill color (default: blue)\n")
  cat("  -n, --intron-color Intron line color (default: black)\n")
  cat("  -w, --plot-width  Plot width in inches (default: 10)\n")
  cat("  -H, --plot-height Plot height in inches (default: auto)\n")
  cat("  -h, --help        Show this help message\n\n")
  cat("Examples:\n")
  cat("  Rscript plot_gene_structure.R -g annotation.gtf -o genes.pdf\n")
  cat("  Rscript plot_gene_structure.R -g annotation.gtf -o genes.pdf \\\n")
  cat("      --genes AT1G01010,AT1G01020 --cds-color steelblue\n")
  quit(status = 0)
}

# Validate required arguments
if (is.null(opts$gtf)) {
  stop("Missing required argument: -g/--gtf (GTF file path). Use -h for help.")
}
if (is.null(opts$output)) {
  stop("Missing required argument: -o/--output (output PDF path). Use -h for help.")
}

# Validate GTF file exists
if (!file.exists(opts$gtf)) {
  stop("GTF file not found: ", opts$gtf)
}

# Apply defaults for optional arguments
gene_filter   <- if (is.null(opts$genes)) NULL else {
                  gs <- trimws(strsplit(opts$genes, ",")[[1]])
                  gs[gs != ""]
                }
exon_color    <- if (is.null(opts$`exon-color`)) "white" else opts$`exon-color`
cds_color     <- if (is.null(opts$`cds-color`))  "blue"  else opts$`cds-color`
intron_color  <- if (is.null(opts$`intron-color`)) "black" else opts$`intron-color`
plot_width    <- if (is.null(opts$`plot-width`))  10      else opts$`plot-width`
plot_height   <- opts$`plot-height`

# -------------------------------------------------------------------------
# Helper: compute intron regions from exon coordinates per gene
# -------------------------------------------------------------------------
compute_introns <- function(exons_df) {
  introns <- exons_df %>%
    group_by(gene_id) %>%
    arrange(start) %>%
    mutate(
      next_start = lead(start),
      intron_start = end,
      intron_end   = next_start
    ) %>%
    filter(!is.na(intron_end)) %>%
    select(gene_id, intron_start, intron_end, strand) %>%
    rename(start = intron_start, end = intron_end) %>%
    mutate(type = "intron") %>%
    ungroup()

  return(introns)
}

# -------------------------------------------------------------------------
# Helper: build the gene structure plot
# -------------------------------------------------------------------------
build_gene_structure_plot <- function(plot_data, gene_labels,
                                      exon_color, cds_color, intron_color) {
  # Separate feature types
  exons_subset   <- subset(plot_data, type == "exon")
  cds_subset     <- subset(plot_data, type == "CDS")
  intron_subset  <- subset(plot_data, type == "intron")

  # Build gene factor for consistent y-axis ordering
  all_gene_ids <- unique(plot_data$gene_id)
  gene_factor  <- factor(all_gene_ids, levels = all_gene_ids)

  # Map gene_id to numeric y position
  gene_y_map <- setNames(seq_along(all_gene_ids), all_gene_ids)

  # Prepare intron data with numeric y
  intron_plot <- intron_subset %>%
    mutate(y = gene_y_map[gene_id])

  # Determine arrow direction based on strand
  # "+" strand: arrow on last segment end (right); "-" strand: first end (left)
  intron_plot$arrow_ends <- ifelse(intron_plot$strand == "+", "last", "first")

  # Prepare exon/CDS data with numeric y
  exons_plot <- exons_subset %>%
    mutate(
      y     = gene_y_map[gene_id],
      ymin  = y - 0.1,
      ymax  = y + 0.1
    )

  cds_plot <- cds_subset %>%
    mutate(
      y     = gene_y_map[gene_id],
      ymin  = y - 0.1,
      ymax  = y + 0.1
    )

  # Prepare label data
  labels_plot <- gene_labels %>%
    mutate(y = gene_y_map[gene_id])

  # Build ggplot
  p <- ggplot()

  # Introns as arrowed segments: two vectorised layers (plus strand
  # arrows point downstream, minus strand upstream) instead of one
  # layer per intron.
  if (nrow(intron_plot) > 0) {
    for (dir in c("last", "first")) {
      seg <- intron_plot[intron_plot$arrow_ends == dir, , drop = FALSE]
      if (nrow(seg) == 0) next
      p <- p + geom_segment(
        data = seg,
        aes(x = .data$start, xend = .data$end, y = .data$y, yend = .data$y),
        arrow = arrow(length = unit(0.1, "inches"), ends = dir, type = "open"),
        color = intron_color,
        linewidth = 0.4
      )
    }
  }

  # Exons as outlined rectangles
  if (nrow(exons_plot) > 0) {
    p <- p + geom_rect(
      data    = exons_plot,
      aes(xmin = .data$start, xmax = .data$end,
          ymin = .data$ymin,  ymax = .data$ymax),
      fill  = exon_color,
      color = "black",
      linewidth = 0.3
    )
  }

  # CDS as filled rectangles (overlaid on exons)
  if (nrow(cds_plot) > 0) {
    p <- p + geom_rect(
      data    = cds_plot,
      aes(xmin = .data$start, xmax = .data$end,
          ymin = .data$ymin,  ymax = .data$ymax),
      fill  = cds_color,
      color = "black",
      linewidth = 0.3
    )
  }

  # Gene labels
  if (nrow(labels_plot) > 0) {
    p <- p + geom_text(
      data  = labels_plot,
      aes(x = .data$midpoint, y = .data$y, label = .data$gene_id),
      vjust = 1.5, hjust = 0.5, size = 3
    )
  }

  # Theme and scales
  p <- p +
    scale_y_continuous(
      breaks = seq_along(all_gene_ids),
      labels = all_gene_ids
    ) +
    theme_void() +
    theme(
      panel.grid       = element_blank(),
      plot.title        = element_text(hjust = 0.5, size = 14),
      axis.text.y       = element_text(size = 10),
      axis.ticks.y      = element_line(linewidth = 0.3),
      axis.line.y.left  = element_line(linewidth = 0.3)
    )

  return(p)
}

# -------------------------------------------------------------------------
# Main logic
# -------------------------------------------------------------------------
main <- function() {
  message("Reading GTF file: ", opts$gtf)

  gtf_data <- tryCatch(
    import(opts$gtf),
    error = function(e) {
      stop("Failed to import GTF file: ", e$message)
    }
  )

  if (!"gene_id" %in% colnames(mcols(gtf_data))) {
    stop("GTF file has no 'gene_id' attribute; features cannot be grouped by gene.")
  }

  # Extract exons and CDS
  exons <- gtf_data[gtf_data$type == "exon", ]
  cds   <- gtf_data[gtf_data$type == "CDS", ]

  if (length(exons) == 0) {
    stop("No exon features found in GTF file.")
  }

  # Build data frames
  exons_df <- data.frame(
    gene_id = exons$gene_id,
    start   = start(exons),
    end     = end(exons),
    type    = "exon",
    strand  = as.character(strand(exons)),
    stringsAsFactors = FALSE
  )

  cds_df <- if (length(cds) > 0) {
    data.frame(
      gene_id = cds$gene_id,
      start   = start(cds),
      end     = end(cds),
      type    = "CDS",
      strand  = as.character(strand(cds)),
      stringsAsFactors = FALSE
    )
  } else {
    # all zero-row columns: mixing a constant with them would fail
    data.frame(gene_id = character(0), start = numeric(0), end = numeric(0),
               type = character(0), strand = character(0),
               stringsAsFactors = FALSE)
  }

  # Filter to requested genes
  if (!is.null(gene_filter)) {
    exons_df <- exons_df[exons_df$gene_id %in% gene_filter, ]
    cds_df   <- cds_df[cds_df$gene_id %in% gene_filter, ]

    if (nrow(exons_df) == 0) {
      stop("None of the specified gene IDs were found in the GTF exon features.")
    }

    missing <- setdiff(gene_filter, unique(exons_df$gene_id))
    if (length(missing) > 0) {
      warning("The following gene IDs were not found: ", paste(missing, collapse = ", "))
    }
  }

  # Warn only about genes whose OWN exons span several chromosomes:
  # different chromosomes across different genes are fine (each gene
  # gets its own row), but one gene mixing chromosomes genuinely
  # garbles its coordinates on the shared x axis.
  genes_in_plot <- unique(exons_df$gene_id)
  n_chr_per_gene <- tapply(as.character(seqnames(exons)), exons$gene_id,
                           function(x) length(unique(x)))
  multi_chr <- intersect(names(n_chr_per_gene)[n_chr_per_gene > 1],
                         genes_in_plot)
  if (length(multi_chr) > 0) {
    warning("Gene(s) span multiple chromosomes and will be mixed on one axis: ",
            paste(multi_chr, collapse = ", "))
  }

  # Compute introns
  introns_df <- compute_introns(exons_df)

  # Combine all features
  plot_data <- bind_rows(exons_df, cds_df, introns_df)

  if (nrow(plot_data) == 0) {
    stop("No plottable features found after processing.")
  }

  # Compute gene label positions (midpoint of each gene)
  gene_labels <- plot_data %>%
    group_by(gene_id) %>%
    summarise(
      gene_start = min(start),
      gene_end   = max(end),
      .groups    = "drop"
    ) %>%
    mutate(midpoint = (gene_start + gene_end) / 2)

  # Auto-calculate plot height if not specified
  n_genes <- length(unique(plot_data$gene_id))
  if (is.null(plot_height)) {
    plot_height <- max(3, 1.2 * n_genes + 1)
  }

  message("Plotting ", n_genes, " gene(s) ...")

  # Build plot
  p <- build_gene_structure_plot(
    plot_data, gene_labels,
    exon_color, cds_color, intron_color
  )

  # Write PDF
  message("Writing output to: ", opts$output)
  pdf(opts$output, width = plot_width, height = plot_height)
  print(p)
  invisible(dev.off())

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
