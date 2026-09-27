#!/usr/bin/env Rscript
#########################################################################
# File Name: calculate_effective_length.R
# Author: ChengYu
# Description: Calculate non-redundant exon length per gene from GTF.
# Created Time: 2026
#########################################################################
# Calculate non-redundant exon length per gene from GTF.
#
# Overlapping exons within each gene are merged (reduced) before summing,
# giving the effective genomic length used for TPM/FPKM calculations.
#
# Changelog:
#   v1.1.0  2026-09-27
#   - NEW: -F/--format exposes the makeTxDbFromGFF format (default gtf;
#     "gff3"/"gtf" accepted — previously a GFF3 file was force-parsed as
#     GTF)
#   - FIX: a GTF that yields no genes now stops with a clear message
#   - CLEAN: the makeCluster() call was dead code — the cluster was
#     created and stopped but never used; -t is accepted for
#     compatibility and validated but does not change behaviour
#   - DOC: exons on opposite strands of the same gene are NOT merged
#     (GenomicRanges reduce(ignore.strand=FALSE) default) — stated here
#   - CHANGE: logs go to stderr (message() instead of cat(); help stays
#     on stdout)

suppressMessages(library(GenomicFeatures))
suppressMessages(library(getopt))

spec <- matrix(c(
    "gtf",     "g", 2, "character", "GTF annotation file",
    "format",  "F", 1, "character", "Annotation format: gtf or gff3 (default: gtf)",
    "output",  "o", 1, "character", "Output file (default: {gtf}.efflen)",
    "threads", "t", 1, "numeric",   "Accepted for compatibility; the script is single-threaded (default: 4)",
    "help",    "h", 0, "logical",   "Show help"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (!is.null(opt$help) || is.null(opt$gtf)) {
    cat("Usage: Rscript calculate_effective_length.R -g annotation.gtf [-o output.efflen] [-F gtf|gff3] [-t 4]\n")
    quit(status = if (!is.null(opt$help)) 0 else 1)
}

if (is.null(opt$output))  opt$output <- paste0(opt$gtf, ".efflen")
if (is.null(opt$threads)) opt$threads <- 4
if (is.null(opt$format))  opt$format  <- "gtf"

if (!opt$format %in% c("gtf", "gff3")) {
    stop("--format must be gtf or gff3, got: ", opt$format)
}
if (opt$threads < 1) stop("--threads must be >= 1, got: ", opt$threads)
if (!file.exists(opt$gtf)) stop("GTF not found: ", opt$gtf)

message("[INFO] GTF: ", opt$gtf)

# Build TxDb
message("[INFO] Building TxDb...")
txdb <- makeTxDbFromGFF(opt$gtf, format = opt$format)

# Calculate effective lengths (exons on opposite strands of the same gene
# are NOT merged — GenomicRanges reduce() with ignore.strand=FALSE)
exons_by_gene <- exonsBy(txdb, by = "gene")
eff_lengths <- sapply(exons_by_gene, function(exons) {
    sum(width(reduce(exons)))
})

if (length(eff_lengths) == 0) {
    stop("No genes found in ", opt$gtf,
         " — check the file and the --format value (gtf vs gff3).")
}

# Write output
result <- data.frame(gene_id = names(eff_lengths), efflen = as.numeric(eff_lengths),
                     row.names = NULL, stringsAsFactors = FALSE)
write.table(result, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

message("[INFO] ", nrow(result), " genes -> ", opt$output)
