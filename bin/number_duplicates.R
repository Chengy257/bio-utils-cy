#!/usr/bin/env Rscript
#########################################################################
# File Name: number_duplicates.R
# Author: ChengYu
# Description: Append sequential numbers to duplicate values in a
#              single-column data file.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-25
#   - FIX: values are read verbatim as strings (colClasses="character",
#     comment.char="", quote=""); previously numeric-looking values were
#     type-converted (leading zeros lost) and '#' or quotes truncated them.
#   - FIX: calling the script with no arguments printed usage and then
#     crashed with an R error (NA in a conditional); it now exits cleanly.
#   - CHANGE: CLI moved to getopt (-i/-o) per project convention; the old
#     positional form `number_duplicates.R in.tsv out.tsv` is gone.
#   - DOC: help now states that EVERY value receives a _N suffix (its
#     occurrence index), guaranteeing uniqueness (unchanged behavior).
#   - CLEANUP: removed unused freq table; summary line goes to stderr;
#     occurrence counter vectorized.
#########################################################################

spec <- matrix(c(
    "input",  "i", 1, "character",
    "output", "o", 1, "character",
    "help",   "h", 0, "logical"
), byrow = TRUE, ncol = 4)

usage <- function() {
    cat("Usage: Rscript number_duplicates.R -i <input.tsv> -o <output.tsv>\n")
    cat("\nAppend sequential occurrence numbers to the values in column 1.\n")
    cat("Every value receives a _N suffix (1st occurrence _1, 2nd _2, ...),\n")
    cat("so the numbered column is guaranteed unique. Row order is preserved.\n")
    cat("\nInput:  single-column TSV file\n")
    cat("Output: two-column TSV: original_value <tab> numbered_value\n")
    cat("\nOptions:\n")
    cat("  -i, --input    input TSV file (required)\n")
    cat("  -o, --output   output TSV file (required)\n")
    cat("  -h, --help     show this help\n")
}

opt <- getopt::getopt(spec)
if (!is.null(opt$help)) {
    usage()
    quit(status = 0)
}
if (is.null(opt$input) || is.null(opt$output)) {
    usage()
    quit(status = 1)
}

if (!file.exists(opt$input)) {
    message("ERROR: input file not found: ", opt$input)
    quit(status = 1)
}

data <- tryCatch(
    read.table(opt$input, sep = "\t", header = FALSE,
               colClasses = "character", comment.char = "", quote = "",
               encoding = "UTF-8"),
    error = function(e) NULL
)
if (is.null(data) || nrow(data) < 1) {
    message("ERROR: input file has no data rows: ", opt$input)
    quit(status = 1)
}

values <- data[, 1]

# occurrence index per value, in input order (replaces the old for-loop)
occurrence <- ave(rep(1L, length(values)), values, FUN = cumsum)
numbered <- paste0(values, "_", occurrence)

result <- data.frame(original = values, numbered = numbered,
                     stringsAsFactors = FALSE)
write.table(result, file = opt$output, sep = "\t", col.names = FALSE,
            row.names = FALSE, quote = FALSE)

message("Processed ", length(values), " rows -> ", opt$output)
