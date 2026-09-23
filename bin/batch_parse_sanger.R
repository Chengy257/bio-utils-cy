#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/batch_parse_sanger.R
# Author: ChengYu
# Description: Batch parse Sanger sequencing .ab1 trace files using
#              sangerseqR. For each file: read trace, call bases, phase
#              alleles against reference, generate chromatogram PDF and
#              pairwise alignment text files.
# Created Time: 2026
#########################################################################

suppressMessages(library(sangerseqR))
suppressMessages(library(Biostrings))

## ---- Argument parsing with getopt ----
library(getopt)

spec <- matrix(c(
    "input-dir",     "i", 2, "character", "Directory containing .ab1 files (default: current dir)",
    "reference",     "r", 1, "character", "Reference sequence file (plain text; required)",
    "signal-cutoff", "s", 2, "numeric",   "Signal ratio cutoff for base calling (default: 0.33)",
    "trim5",         "t", 2, "integer",   "Bases to trim from 5-prime end (default: 20)",
    "output-dir",    "o", 2, "character", "Output directory (default: current dir)",
    "help",          "h", 0, "logical",   "Show this help message"
), ncol = 5, byrow = TRUE)

opt <- getopt(spec)

## ---- Help ----
if (!is.null(opt$help)) {
    cat("Usage: Rscript batch_parse_sanger.R -r ref.txt [-i ab1_dir] [options]\n\n")
    cat("Batch parse Sanger .ab1 sequencing trace files.\n\n")
    cat("Options:\n")
    cat("  -i, --input-dir     Directory with .ab1 files (default: current dir)\n")
    cat("  -r, --reference     Reference sequence file, plain text [required]\n")
    cat("  -s, --signal-cutoff Base call signal ratio cutoff (default: 0.33)\n")
    cat("  -t, --trim5         Bases to trim from 5-prime end (default: 20)\n")
    cat("  -o, --output-dir    Output directory (default: current dir)\n")
    cat("  -h, --help          Show this help message\n\n")
    cat("Outputs per .ab1 file:\n")
    cat("  {basename}_chromatogram.pdf\n")
    cat("  {basename}_pairwiseAlignment_ref.txt\n")
    cat("  {basename}_pairwiseAlignment_allele.txt\n\n")
    cat("Example:\n")
    cat("  Rscript batch_parse_sanger.R -r reference.txt -i traces/ -o results/\n")
    quit(status = 0)
}

## ---- Validate required arguments ----
if (is.null(opt$reference)) {
    stop("Error: -r/--reference is required. Use -h for help.")
}
if (!file.exists(opt$reference)) {
    stop(paste("Error: Reference file not found:", opt$reference))
}

## ---- Set defaults ----
input_dir   <- if (is.null(opt[["input-dir"]]))    "."      else opt[["input-dir"]]
signal_cut  <- if (is.null(opt[["signal-cutoff"]])) 0.33    else opt[["signal-cutoff"]]
trim5_val   <- if (is.null(opt$trim5))              20      else opt$trim5
output_dir  <- if (is.null(opt[["output-dir"]]))    "."      else opt[["output-dir"]]

if (!dir.exists(input_dir)) {
    stop(paste("Error: Input directory not found:", input_dir))
}

## ---- Create output directory if needed ----
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
}

## ---- Read reference sequence ----
read_ref_sequence <- function(ref_file) {
    lines <- readLines(ref_file, warn = FALSE)
    # Remove header lines starting with '>' and blank lines
    seq_lines <- lines[!grepl("^>", lines) & nzchar(trimws(lines))]
    seq <- paste(seq_lines, collapse = "")
    seq <- gsub("\\s+", "", seq)
    if (nchar(seq) == 0) {
        stop("Error: Reference file contains no sequence data.")
    }
    return(toupper(seq))
}

ref_seq <- tryCatch({
    read_ref_sequence(opt$reference)
}, error = function(e) {
    stop(paste("Error reading reference file:", conditionMessage(e)))
})
message("Reference sequence loaded: ", nchar(ref_seq), " bases")

## ---- Find .ab1 files ----
ab1_files <- list.files(input_dir, pattern = "\\.ab1$", full.names = TRUE,
                        ignore.case = TRUE)

if (length(ab1_files) == 0) {
    stop("Error: No .ab1 files found in: ", input_dir)
}
message("Found ", length(ab1_files), " .ab1 file(s) to process\n")

## ---- Parse single .ab1 file ----
parse_single_sanger <- function(ab1_path, ref_seq, signal_cutoff,
                                trim5, out_dir) {
    basename_str <- tools::file_path_sans_ext(basename(ab1_path))
    message("Processing: ", basename_str)

    # Read Sanger sequencing data
    sangerseq_obj <- sangerseqR::readsangerseq(ab1_path)

    # Make base calls with specified signal ratio
    basecalls <- sangerseqR::makeBaseCalls(sangerseq_obj, ratio = signal_cutoff)

    # Phase alleles using reference sequence
    phase_alleles <- tryCatch({
        sangerseqR::setAllelePhase(basecalls, ref_seq, trim5 = trim5)
    }, error = function(e) {
        message("  Warning: Could not phase alleles for ", basename_str,
                " - ", conditionMessage(e))
        return(NULL)
    })

    # Generate chromatogram PDF
    chromo_file <- file.path(out_dir, paste0(basename_str, "_chromatogram.pdf"))
    message("  Writing chromatogram: ", chromo_file)
    tryCatch({
        pdf(chromo_file, width = 100, height = 1)
        sangerseqR::chromatogram(basecalls, trim5 = trim5,
                                 width = 100, height = 1,
                                 showcalls = "both")
        dev.off()
    }, error = function(e) {
        message("  Warning: Chromatogram failed for ", basename_str,
                " - ", conditionMessage(e))
        if (dev.cur() > 1) dev.off()
    })

    if (!is.null(phase_alleles)) {
        # Pairwise alignment: primary sequence vs reference
        prim_seq  <- sangerseqR::primarySeq(phase_alleles)
        sec_seq   <- sangerseqR::secondarySeq(phase_alleles)

        align_ref_file <- file.path(out_dir,
                            paste0(basename_str, "_pairwiseAlignment_ref.txt"))
        align_allele_file <- file.path(out_dir,
                            paste0(basename_str, "_pairwiseAlignment_allele.txt"))

        tryCatch({
            message("  Writing pairwise alignment (primary vs ref): ",
                    align_ref_file)
            pair_prim <- Biostrings::pairwiseAlignment(
                prim_seq, ref_seq, type = "local-global")
            Biostrings::writePairwiseAlignments(pair_prim,
                file = align_ref_file)
        }, error = function(e) {
            message("  Warning: Primary alignment failed - ",
                    conditionMessage(e))
            writeLines(paste("Primary alignment failed:",
                             conditionMessage(e)),
                       con = align_ref_file)
        })

        tryCatch({
            message("  Writing pairwise alignment (secondary vs ref): ",
                    align_allele_file)
            pair_sec <- Biostrings::pairwiseAlignment(
                sec_seq, ref_seq, type = "local-global")
            Biostrings::writePairwiseAlignments(pair_sec,
                file = align_allele_file)
        }, error = function(e) {
            message("  Warning: Secondary alignment failed - ",
                    conditionMessage(e))
            writeLines(paste("Secondary alignment failed:",
                             conditionMessage(e)),
                       con = align_allele_file)
        })
    } else {
        # No phasing possible - still try basic alignments
        message("  Skipping allele-specific alignments (phasing failed)")
    }

    message("  Done: ", basename_str, "\n")
    invisible(TRUE)
}

## ---- Process all files ----
success_count <- 0
for (i in seq_along(ab1_files)) {
    message(sprintf("[%d/%d]", i, length(ab1_files)))
    tryCatch({
        parse_single_sanger(
            ab1_path      = ab1_files[i],
            ref_seq       = ref_seq,
            signal_cutoff = signal_cut,
            trim5         = trim5_val,
            out_dir       = output_dir
        )
        success_count <- success_count + 1
    }, error = function(e) {
        message("  Error processing ", basename(ab1_files[i]), ": ",
                conditionMessage(e))
    })
}

## ---- Summary ----
message("=== Summary ===")
message(sprintf("Processed: %d / %d files successfully",
                success_count, length(ab1_files)))
message("Output directory: ", normalizePath(output_dir))
