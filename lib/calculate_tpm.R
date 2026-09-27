#########################################################################
# File Name: calculate_tpm.R
# Author: ChengYu
# Description: Expression unit conversion functions (counts to TPM,
#              FPKM, RPM, effective counts) — shared library.
# Created Time: 2026
#########################################################################
# calculate_tpm.R — expression unit conversion functions (shared library).
#
# Pure-function library, sourced by other scripts — NOT a CLI tool:
#   source("<repo>/lib/calculate_tpm.R")
#   tpm_values  <- countToTpm(counts, effLengths)
#   fpkm_values <- countToFpkm(counts, effLengths)
#
# Conventions:
#   TPM  = (count / effLen) / sum(count / effLen) * 1e6
#   FPKM = count / (effLen / 1e3) / (sum(count) / 1e6)
#   RPM  = count / sum(count) * 1e6
#
# Changelog:
#   v1.1.0  2026-09-27
#   - MOVE: bin/ -> lib/ (non-CLI shared code belongs in lib/ per README)
#   - FIX: zero/negative effective lengths and all-zero count vectors now
#     stop with a clear message instead of returning NaN/Inf silently
#   - FIX: input vectors must be non-empty and of equal length
#   - CLEAN: Python-style docstring replaced with R comments; the banner
#     cat()s on source() are gone (they polluted every caller's stdout);
#     the shebang is gone (this file is sourced, not executed)

#' Convert counts to TPM
#' @param counts Numeric vector of raw counts
#' @param effLen Numeric vector of effective lengths (same length as counts)
#' @return Numeric vector of TPM values
countToTpm <- function(counts, effLen) {
    counts <- as.numeric(counts); effLen <- as.numeric(effLen)
    if (length(counts) == 0 || length(counts) != length(effLen)) {
        stop("countToTpm: counts and effLen must be non-empty and of equal length (",
             length(counts), " vs ", length(effLen), ")")
    }
    if (anyNA(effLen) || any(effLen <= 0)) {
        stop("countToTpm: effective lengths must be positive and non-NA")
    }
    rate <- counts / effLen
    total <- sum(rate)
    if (is.na(total) || total <= 0) {
        stop("countToTpm: sum(count/effLen) is ", format(total),
             " — counts are all zero or non-numeric; TPM is undefined")
    }
    rate / total * 1e6
}

#' Convert counts to FPKM
#' @param counts Numeric vector of raw counts
#' @param effLen Numeric vector of effective lengths
#' @return Numeric vector of FPKM values
countToFpkm <- function(counts, effLen) {
    counts <- as.numeric(counts); effLen <- as.numeric(effLen)
    if (length(counts) == 0 || length(counts) != length(effLen)) {
        stop("countToFpkm: counts and effLen must be non-empty and of equal length (",
             length(counts), " vs ", length(effLen), ")")
    }
    if (anyNA(effLen) || any(effLen <= 0)) {
        stop("countToFpkm: effective lengths must be positive and non-NA")
    }
    n_total <- sum(counts)
    if (is.na(n_total) || n_total <= 0) {
        stop("countToFpkm: sum(counts) is ", format(n_total),
             " — FPKM is undefined for an all-zero count vector")
    }
    counts / (effLen / 1e3) / (n_total / 1e6)
}

#' Convert FPKM to TPM
#' @param fpkm Numeric vector of FPKM values
#' @return Numeric vector of TPM values
fpkmToTpm <- function(fpkm) {
    fpkm <- as.numeric(fpkm)
    if (length(fpkm) == 0) {
        stop("fpkmToTpm: fpkm must be non-empty")
    }
    total <- sum(fpkm)
    if (is.na(total) || total <= 0) {
        stop("fpkmToTpm: sum(fpkm) is ", format(total), " — TPM is undefined")
    }
    fpkm / total * 1e6
}

#' Convert counts to RPM (reads per million)
#' @param counts Numeric vector of raw counts
#' @return Numeric vector of RPM values
countToRpm <- function(counts) {
    counts <- as.numeric(counts)
    if (length(counts) == 0) {
        stop("countToRpm: counts must be non-empty")
    }
    total <- sum(counts)
    if (is.na(total) || total <= 0) {
        stop("countToRpm: sum(counts) is ", format(total), " — RPM is undefined")
    }
    counts / total * 1e6
}

#' Convert counts to effective counts
#' @param counts Numeric vector of raw counts
#' @param len Numeric vector of actual transcript lengths
#' @param effLen Numeric vector of effective lengths
#' @return Numeric vector of effective counts
countToEffCounts <- function(counts, len, effLen) {
    counts <- as.numeric(counts); len <- as.numeric(len); effLen <- as.numeric(effLen)
    if (length(counts) == 0 || length(counts) != length(len) ||
        length(counts) != length(effLen)) {
        stop("countToEffCounts: counts/len/effLen must be non-empty and of equal length (",
             length(counts), "/", length(len), "/", length(effLen), ")")
    }
    if (any(effLen <= 0)) {
        stop("countToEffCounts: effective lengths must be positive")
    }
    counts * (len / effLen)
}
