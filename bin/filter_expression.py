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
2. MAD percentile: keep genes whose median absolute deviation (raw MAD,
   NOT scaled by 1.4826) is >= the mad_percentile-th percentile of all
   MADs. With the default 75 this keeps the most variable 25% of genes.

Changelog:
  v1.1.0  2026-09-27
  - FIX: the --mad_percentile help was self-contradictory ("Keep top N%
    (default: 75, i.e. top 25%)"); the semantics are: keep genes with
    MAD >= this percentile, so 75 keeps the top 25%
  - FIX: a single gene containing NaN made the MAD threshold NaN
    (np.percentile does not skip NaN), which failed the comparison for
    EVERY gene and silently emptied the output; the percentile is now
    NaN-aware (nanpercentile) and NaN-MAD genes are counted, warned
    about, and dropped
  - FIX: matrices without any numeric column now fail with a clear
    message instead of filtering everything away
  - FIX: --min_samples < 1 and --mad_percentile outside [0, 100] are
    rejected instead of producing empty/nonsense output
  - CHANGE: the row-wise MAD is computed vectorized (apply(axis=1) was
    O(n) Python calls) with skipna=False, preserving the original
    per-row semantics
  - DOC: raw MAD is used deliberately (not scaled by 1.4826); NaN in the
    expression filter counts as a missing measurement (the row can still
    pass on its valid samples)
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

__version__ = "1.1.0"


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

    Rows containing NaN cannot pass the comparison and are counted in a
    warning (they are dropped like any non-passing row).
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) == 0:
        raise ValueError("No numeric columns in the input matrix — nothing to filter.")
    values = df[numeric_cols]
    passing = (values >= min_expr).sum(axis=1) >= min_samples
    n_nan = int(values.isna().any(axis=1).sum())
    if n_nan:
        logging.warning(
            "%d gene(s) contain NaN and cannot pass the expression filter (dropped).",
            n_nan,
        )
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
    """Filter genes by MAD percentile (keep genes with MAD >= the percentile).

    Args:
        df: Expression DataFrame.
        mad_percentile: Percentile (0-100) of the MAD distribution; genes
            with MAD >= this value are kept. 75 keeps the most variable 25%.

    Returns:
        Filtered DataFrame.

    Raw MAD is used (median absolute deviation, not scaled by 1.4826) —
    percentile ranking is unaffected by the constant anyway.
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) == 0:
        raise ValueError("No numeric columns in the input matrix — nothing to filter.")
    values = df[numeric_cols]
    # vectorized row MAD: median(|x - median(x)|), skipna=False so that a
    # row containing NaN gets an undefined MAD (the pre-vectorization
    # behaviour) instead of a silently recomputed one. The percentile
    # threshold ignores NaN rows (np.percentile would return NaN and
    # filter the whole matrix away).
    centers = values.median(axis=1, skipna=False)
    mad_values = values.sub(centers, axis=0).abs().median(axis=1, skipna=False)
    n_nan = int(mad_values.isna().sum())
    if n_nan:
        logging.warning(
            "%d gene(s) have undefined MAD (NaN in some samples) and are dropped.",
            n_nan,
        )
    threshold = np.nanpercentile(mad_values, mad_percentile)
    result = df[mad_values >= threshold]
    logging.info(
        "MAD filter (MAD >= %.0fth percentile, threshold=%.4g, i.e. the most variable %.0f%%): "
        "%d -> %d genes.",
        mad_percentile, threshold, 100 - mad_percentile, len(df), len(result),
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
        help="Keep genes with MAD >= this percentile, 0-100 (default: 75, "
             "i.e. the most variable 25%%). Raw MAD, not scaled by 1.4826.",
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

    if args.min_samples < 1:
        logging.error("--min_samples must be >= 1, got: %d", args.min_samples)
        sys.exit(1)
    if not 0 <= args.mad_percentile <= 100:
        logging.error("--mad_percentile must be within [0, 100], got: %s", args.mad_percentile)
        sys.exit(1)

    df = pd.read_csv(args.input, sep="\t", index_col=0)
    logging.info("Input: %d genes x %d samples.", len(df), len(df.columns))

    try:
        # Expression filter
        df = filter_by_expression(df, args.min_expr, args.min_samples)

        # MAD filter
        if not args.skip_mad:
            df = filter_by_mad(df, args.mad_percentile)
    except ValueError as exc:
        logging.error("%s", exc)
        sys.exit(1)

    if df.empty:
        logging.warning("All genes filtered out! Consider relaxing thresholds.")

    df.to_csv(args.output, sep="\t")
    logging.info("Output: %d genes -> %s", len(df), args.output)


if __name__ == "__main__":
    main()
