#!/usr/bin/env Rscript
#########################################################################
# File Name: calculate_tpm.R
# Author: ChengYu
# Description: Expression unit conversion functions: counts to TPM,
#              FPKM, RPM, and inter-conversions.
# Created Time: 2026
#########################################################################
"""Expression unit conversion utility.

Provides functions for converting between expression units:
  - Counts -> TPM, FPKM, RPM, effective counts
  - FPKM -> TPM

Usage:
  source("calculate_tpm.R")

  # Convert count vector to TPM
  tpm_values <- countToTpm(counts, effLengths)

  # Convert count vector to FPKM
  fpkm_values <- countToFpkm(counts, effLengths)
"""

#' Convert counts to TPM
#' @param counts Numeric vector of raw counts
#' @param effLen Numeric vector of effective lengths (same length as counts)
#' @return Numeric vector of TPM values
countToTpm <- function(counts, effLen) {
    rate <- counts / effLen
    rate / sum(rate) * 1e6
}

#' Convert counts to FPKM
#' @param counts Numeric vector of raw counts
#' @param effLen Numeric vector of effective lengths
#' @return Numeric vector of FPKM values
countToFpkm <- function(counts, effLen) {
    N <- sum(counts)
    counts / (effLen / 1e3) / (N / 1e6)
}

#' Convert FPKM to TPM
#' @param fpkm Numeric vector of FPKM values
#' @return Numeric vector of TPM values
fpkmToTpm <- function(fpkm) {
    fpkm / sum(fpkm) * 1e6
}

#' Convert counts to RPM (reads per million)
#' @param counts Numeric vector of raw counts
#' @return Numeric vector of RPM values
countToRpm <- function(counts) {
    counts / sum(counts) * 1e6
}

#' Convert counts to effective counts
#' @param counts Numeric vector of raw counts
#' @param len Numeric vector of actual transcript lengths
#' @param effLen Numeric vector of effective lengths
#' @return Numeric vector of effective counts
countToEffCounts <- function(counts, len, effLen) {
    counts * (len / effLen)
}

# Print info when sourced
cat("[INFO] calculate_tpm.R loaded. Available functions:\n")
cat("  countToTpm(counts, effLen)\n")
cat("  countToFpkm(counts, effLen)\n")
cat("  fpkmToTpm(fpkm)\n")
cat("  countToRpm(counts)\n")
cat("  countToEffCounts(counts, len, effLen)\n")
