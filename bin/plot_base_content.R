#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/plot_base_content.R
# Author: ChengYu
# Description: Sliding window base composition analysis and plotting
#              for FASTA sequences.
# Created Time: 2026
# Changelog:
#   v1.1.0  2026-09-28
#   - FIX: --cutoff/--output/--plot-width/--plot-height were declared
#     as boolean flags (getopt argflag 0) while taking values, so
#     "-c 50" aborted with a getopt error and a bare "--cutoff" was
#     silently read as TRUE (reference line drawn at 1). They now
#     require a value; when omitted, the code-side defaults apply and
#     a bare flag is a loud getopt error.
#   - FIX: odd window sizes analysed window_size - 1 bp (half_win =
#     floor(w/2) built a 2*floor(w/2)-bp window). Odd windows now span
#     the full w, centred on the position; even windows keep their
#     previous anchoring, so existing even-window output is unchanged.
#   - FIX: with multiple sequences and an explicit -o/--output prefix,
#     every sequence wrote to the same files and only the last one
#     survived. The sequence name is appended to the prefix whenever
#     the FASTA holds more than one sequence.
#   - CHANGE: the base-type split inside the sliding-window loop is
#     hoisted out (computed once per sequence).
#   - CHANGE: the statistics file extension is now .tsv (the file was
#     always tab-delimited, never Excel).
#   - DOCUMENTED: the denominator includes N and other ambiguous
#     codes, i.e. the value is the fraction of the whole window.
#   - CHANGE: INFO/WARN logs go to stderr via message().
#########################################################################

# Sliding Window Base Composition Analysis
#
# Reads a FASTA file, computes the percentage of a specified base (or
# base combination such as GC or AT) in a sliding window across each
# sequence, writes a per-sequence statistics TSV, and generates a
# per-sequence line plot in PDF format.
#
# Dependencies:
#   install.packages("getopt")
#   BiocManager::install(c("Biostrings", "ggplot2"))

suppressPackageStartupMessages({
  library(getopt)
  library(Biostrings)
  library(ggplot2)
})

# -------------------------------------------------------------------------
# Command-line interface
# -------------------------------------------------------------------------

spec <- matrix(c(
  "help",        "h", 0, "logical",   "Show this help message",
  "fasta",       "f", 2, "character", "Input FASTA file (required)",
  "base",        "b", 2, "character", "Base type to count: A, T, G, C, GC, or AT (required)",
  "window",      "w", 2, "numeric",   "Sliding window size in bp (required)",
  "cutoff",      "c", 1, "numeric",   "Horizontal reference line threshold (default: 75)",
  "output",      "o", 1, "character", "Output prefix (default: auto from base and sequence name)",
  "plot-width",  "p", 1, "numeric",   "Plot width in inches (default: 6)",
  "plot-height", "q", 1, "numeric",   "Plot height in inches (default: 3)"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

# Help
if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  cat("\n")
  cat("Description:\n")
  cat("  Sliding window base composition analysis for FASTA sequences.\n")
  cat("  For each sequence, computes the percentage of the requested base(s)\n")
  cat("  in a sliding window, writes a TSV statistics table, and produces\n")
  cat("  a line plot in PDF format.\n\n")
  cat("Examples:\n")
  cat("  # GC content with 100 bp window\n")
  cat("  Rscript plot_base_content.R -f genome.fa -b GC -w 100\n\n")
  cat("  # AT content with custom cutoff and output prefix\n")
  cat("  Rscript plot_base_content.R -f seqs.fa -b AT -w 200 \\\n")
  cat("    --cutoff 50 -o at_content\n\n")
  cat("  # Single base with larger plot\n")
  cat("  Rscript plot_base_content.R -f contigs.fa -b G -w 500 \\\n")
  cat("    --plot-width 10 --plot-height 4\n\n")
  quit(status = 0)
}

# Validate required arguments
if (is.null(opt$fasta) || is.null(opt$base) || is.null(opt$window)) {
  stop("Missing required arguments (-f, -b, -w). Use -h for help.")
}

# -------------------------------------------------------------------------
# Parameters and validation
# -------------------------------------------------------------------------

fasta_file  <- opt$fasta
base_type   <- toupper(opt$base)
window_size <- opt$window
cutoff_val  <- if (is.null(opt$cutoff)) 75 else opt$cutoff
plot_width  <- if (is.null(opt$`plot-width`)) 6 else opt$`plot-width`
plot_height <- if (is.null(opt$`plot-height`)) 3 else opt$`plot-height`

# Validate FASTA file
if (!file.exists(fasta_file)) {
  stop("FASTA file not found: ", fasta_file)
}

# Validate base type
valid_bases <- c("A", "T", "G", "C", "GC", "AT")
if (!(base_type %in% valid_bases)) {
  stop("Invalid base type '", opt$base, "'. Must be one of: ",
       paste(valid_bases, collapse = ", "))
}

# Validate window size
if (window_size < 1 || window_size != as.integer(window_size)) {
  stop("Window size must be a positive integer.")
}

# Guard against a bare "--cutoff" (getopt then yields TRUE, not a number)
if (!is.numeric(cutoff_val)) {
  stop("--cutoff needs a numeric value, e.g. --cutoff 50")
}
if (!is.numeric(plot_width) || plot_width <= 0) {
  stop("--plot-width needs a positive number")
}
if (!is.numeric(plot_height) || plot_height <= 0) {
  stop("--plot-height needs a positive number")
}

# -------------------------------------------------------------------------
# Read FASTA
# -------------------------------------------------------------------------

fasta_sequences <- tryCatch(
  readDNAStringSet(fasta_file),
  error = function(e) {
    stop("Failed to read FASTA file: ", e$message)
  }
)

n_seq <- length(fasta_sequences)
if (n_seq == 0) {
  message("[INFO] FASTA file contains no sequences. Nothing to do.")
  quit(status = 0)
}

message("[INFO] Loaded ", n_seq, " sequence(s) from ", fasta_file)
message("[INFO] Base type: ", base_type, " | Window size: ", window_size,
        " bp | Cutoff: ", cutoff_val, "%")

# -------------------------------------------------------------------------
# Process each sequence
# -------------------------------------------------------------------------

half_win   <- as.integer(floor(window_size / 2))
# Window is [i - left, i + right]: for even w this keeps the historical
# anchoring ([i - w/2, i + w/2 - 1]); for odd w it is truly centred and
# spans the full w bp (the old floor(w/2) scheme analysed w - 1 bp).
if (window_size %% 2 == 1) {
  win_left  <- half_win
  win_right <- half_win
} else {
  win_left  <- half_win
  win_right <- half_win - 1L
}
processed  <- 0
skipped    <- 0

for (j in seq_len(n_seq)) {
  seq_name <- names(fasta_sequences)[j]
  seq_len  <- length(fasta_sequences[[j]])
  # With an explicit -o prefix and multiple sequences, every sequence
  # used to write to the same files (only the last survived); append
  # the sequence name to keep each sequence's outputs distinct.
  out_prefix <- if (!is.null(opt$output)) {
    if (n_seq > 1) {
      paste0(opt$output, "_",
             gsub("[^A-Za-z0-9._-]+", "_", seq_name))
    } else {
      opt$output
    }
  } else {
    paste0(base_type, "_", gsub("[[:space:]]+", "_", seq_name),
           "_window", window_size)
  }

  # Edge case: window larger than sequence
  if (window_size > seq_len) {
    message("[WARN] Skipping '", seq_name, "' (length ", seq_len,
          " bp): window size (", window_size, " bp) exceeds sequence length.")
    skipped <- skipped + 1
    next
  }

  message("[INFO] Processing sequence ", j, "/", n_seq, ": ", seq_name,
        " (", seq_len, " bp)")

  # Compute sliding window base content
  seq_start <- win_left + 1
  seq_end   <- seq_len - win_right
  positions <- seq_start:seq_end
  target_bases <- if (base_type %in% c("GC", "AT")) strsplit(base_type, "")[[1]] else base_type

  base_content <- numeric(seq_len)  # pre-allocate full length, subset later
  for (i in positions) {
    win_start <- i - win_left
    win_end   <- i + win_right
    window_seq <- subseq(fasta_sequences[[j]], start = win_start, end = win_end)

    if (base_type %in% c("GC", "AT")) {
      # For compound base types, sum individual base frequencies
      freq_sum <- 0
      for (b in target_bases) {
        freq_sum <- freq_sum + letterFrequency(window_seq, b, as.prob = TRUE)[[1]]
      }
      base_content[i] <- round(freq_sum * 100, 2)
    } else {
      base_content[i] <- round(
        letterFrequency(window_seq, base_type, as.prob = TRUE)[[1]] * 100, 2
      )
    }
  }

  # Extract only the computed positions
  base_content <- base_content[positions]
  base_content <- na.omit(base_content)

  # Build data frame
  data <- data.frame(
    index   = positions,
    Content = base_content,
    Base    = base_type,
    stringsAsFactors = FALSE
  )

  # Write statistics table
  tsv_file <- paste0(out_prefix, "_statistics_table.tsv")
  write.table(data, file = tsv_file, sep = "\t",
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  message("[INFO]   Wrote statistics: ", tsv_file)

  # Generate line plot
  p <- ggplot(data, aes(x = index, y = Content)) +
    geom_line(linewidth = 0.4) +
    geom_hline(yintercept = cutoff_val, color = "purple",
               linewidth = 0.8, linetype = "dashed") +
    labs(
      title = paste0(base_type, " Content in Sliding Window of ",
                     window_size, " bp: ", seq_name),
      x = "Position (bp)",
      y = paste0(base_type, " Content (%)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 9),
      axis.text  = element_text(size = 8)
    )

  pdf_file <- paste0(out_prefix, "_plot.pdf")
  ggsave(filename = pdf_file, plot = p, width = plot_width,
         height = plot_height, device = "pdf")
  message("[INFO]   Wrote plot: ", pdf_file)

  processed <- processed + 1
}

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

message("")
message("[INFO] ===== Summary =====")
message("[INFO] Total sequences in FASTA: ", n_seq)
message("[INFO] Successfully processed:    ", processed)
if (skipped > 0) {
  message("[WARN] Skipped (window > seq):    ", skipped)
}
message("[INFO] Done.")
