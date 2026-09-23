#!/usr/bin/env Rscript
#########################################################################
# File Name: number_duplicates.R
# Author: ChengYu
# Description: Append sequential numbers to duplicate values in a
#              single-column data file.
# Created Time: 2026
#########################################################################
# Append sequential numbers to duplicate values.
#
# For each unique value in the input column, appends _1, _2, _3, etc.
# to duplicates, preserving the original row order.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || args[1] == "--help" || args[1] == "-h") {
    cat("Usage: Rscript number_duplicates.R <input.tsv> <output.tsv>\n")
    cat("\nAppend sequential numbers to duplicate values in column 1.\n")
    cat("Each unique value gets _1, _2, _3, ... appended to duplicates.\n")
    cat("\nInput:  Single-column TSV file\n")
    cat("Output: Two-column TSV: original_value, numbered_value\n")
    quit(status = if (args[1] == "--help" || args[1] == "-h") 0 else 1)
}

input_file <- args[1]
output_file <- args[2]

# Validate input
if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
}

# Read data
data <- read.table(input_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")
if (ncol(data) < 1) {
    stop("Input file has no columns.")
}

values <- data[, 1]
freq <- table(values)

# Build numbered vector
numbered <- character(length(values))
freq_counter <- rep(0, length(unique(values)))
names(freq_counter) <- unique(values)

for (i in seq_along(values)) {
    v <- values[i]
    freq_counter[v] <- freq_counter[v] + 1
    numbered[i] <- paste0(v, "_", freq_counter[v])
}

# Write output (original order)
result <- data.frame(original = values, numbered = numbered, stringsAsFactors = FALSE)
write.table(result, file = output_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

cat("Processed", length(values), "rows ->", output_file, "\n")
