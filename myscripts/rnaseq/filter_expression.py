#!/usr/bin/env python3
#########################################################################
# File Name: filter_expression.py
# Author: ChengYu
# Description: Filter low-expression and low-variability genes from
#              an expression matrix.
# Created Time: 2026
#########################################################################
"""Filter genes from an expression matrix by expression level and MAD.

Applies two filters:
1. Expression threshold: keep genes with expression >= min_expr in at
   least min_samples samples.
2. MAD percentile: keep genes in the top N% by median absolute deviation
   (i.e., most variable genes).
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

__version__ = "1.0.0"


def filter_by_expression(
    df: pd.DataFrame,
    min_expr: float = 1.0,
    min_samples: int = 2,
) -> pd.DataFrame:
    """Filter genes by minimum expression in at least N samples.

    Args:
        df: Expression DataFrame (genes x samples).
        min_expr: Minimum expression value.
        min_samples: Minimum number of samples meeting threshold.

    Returns:
        Filtered DataFrame.
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    passing = (df[numeric_cols] >= min_expr).sum(axis=1) >= min_samples
    result = df[passing]
    logging.info(
        "Expression filter (>= %.1f in >= %d samples): %d -> %d genes (%d removed).",
        min_expr, min_samples, len(df), len(result), len(df) - len(result),
    )
    return result


def filter_by_mad(
    df: pd.DataFrame,
    mad_percentile: float = 75.0,
) -> pd.DataFrame:
    """Filter genes by MAD percentile (keep top N%).

    Args:
        df: Expression DataFrame.
        mad_percentile: Keep genes above this percentile of MAD.

    Returns:
        Filtered DataFrame.
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    mad_values = df[numeric_cols].apply(lambda row: np.median(np.abs(row - np.median(row))), axis=1)
    threshold = np.percentile(mad_values, mad_percentile)
    result = df[mad_values >= threshold]
    logging.info(
        "MAD filter (top %.0f%%, threshold=%.4f): %d -> %d genes.",
        100 - mad_percentile, threshold, len(df), len(result),
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Filter low-expression and low-variability genes from an expression matrix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
The input TSV should have genes as rows and samples as columns.
The first column is used as the gene ID index.

Filters are applied sequentially:
1. Expression filter: keep genes with expression >= threshold in >= N samples
2. MAD filter: keep genes in the top X% by variability

examples:
  python filter_expression.py -i expr.tsv -o filtered.tsv
  python filter_expression.py -i expr.tsv -o filtered.tsv --min_expr 5 --min_samples 3 --mad_percentile 80
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input expression matrix TSV.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output filtered TSV.")
    parser.add_argument(
        "--min_expr", type=float, default=1.0,
        help="Minimum expression threshold (default: 1.0).",
    )
    parser.add_argument(
        "--min_samples", type=int, default=2,
        help="Minimum samples with expression >= threshold (default: 2).",
    )
    parser.add_argument(
        "--mad_percentile", type=float, default=75.0,
        help="Keep top N%% by MAD (default: 75, i.e. top 25%%).",
    )
    parser.add_argument(
        "--skip-mad", action="store_true",
        help="Skip MAD variability filter.",
    )
    parser.add_argument(
        "--log-level", type=str, default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO).",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main() -> None:
    """Entry point."""
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=getattr(logging, args.log_level),
    )

    if not Path(args.input).is_file():
        logging.error("Input file not found: %s", args.input)
        sys.exit(1)

    df = pd.read_csv(args.input, sep="\t", index_col=0)
    logging.info("Input: %d genes x %d samples.", len(df), len(df.columns))

    # Expression filter
    df = filter_by_expression(df, args.min_expr, args.min_samples)

    # MAD filter
    if not args.skip_mad:
        df = filter_by_mad(df, args.mad_percentile)

    if df.empty:
        logging.warning("All genes filtered out! Consider relaxing thresholds.")

    df.to_csv(args.output, sep="\t")
    logging.info("Output: %d genes -> %s", len(df), args.output)


if __name__ == "__main__":
    main()
