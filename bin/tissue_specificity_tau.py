#!/usr/bin/env python3
#########################################################################
# File Name: tissue_specificity_tau.py
# Author: ChengYu
# Description: Calculate tissue specificity index (tau) for genes
#              from an expression matrix.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-25
#   - FIX: headerless input no longer produces a bogus "0  tau" header
#     line in the output.
#   - FIX: gene IDs are written as-is (strings); purely numeric IDs no
#     longer depend on pandas type inference.
#   - FIX: negative expression values fail with a clean error (tau is
#     undefined); rows containing NaN are counted and warned about
#     instead of being silently lumped into the all-zero statistic.
#   - PERF: tau is computed vectorized over the whole matrix instead of
#     one apply() call per gene.
#########################################################################
"""Calculate tissue specificity index (tau) for genes.

Tau measures tissue specificity: 0 = ubiquitously expressed,
1 = expressed in a single tissue. Input is a TSV with gene IDs in
column 1 and expression values in subsequent columns.
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

__version__ = "1.1.0"


def compute_tau(input_path: str, output_path: str, has_header: bool = False) -> None:
    """Compute tau for all genes in an expression matrix.

    Args:
        input_path: Input TSV file (gene ID in col 1, expressions in cols 2+).
        output_path: Output TSV file.
        has_header: Whether input has a header row.
    """
    header = 0 if has_header else None
    df = pd.read_csv(input_path, sep="\t", header=header)

    # First column is gene ID; keep it verbatim (as strings)
    df.iloc[:, 0] = df.iloc[:, 0].astype(str)
    gene_col = df.columns[0]
    expr_cols = df.columns[1:]

    if len(expr_cols) == 0:
        logging.error("No expression columns found (need at least 2 columns).")
        sys.exit(1)

    values = df[expr_cols].to_numpy(dtype=float)
    n_tissues = values.shape[1]

    n_negative = int((values < 0).sum())
    if n_negative:
        logging.error(
            "Found %d negative expression value(s); tau is undefined for negative values.",
            n_negative,
        )
        sys.exit(1)

    n_nan_rows = int(np.isnan(values).any(axis=1).sum())
    if n_nan_rows:
        logging.warning("%d gene(s) contain missing values; their tau is NA.", n_nan_rows)

    logging.info("Computing tau for %d genes across %d tissues.", len(df), n_tissues)

    # tau = sum(1 - x/x_max) / (n - 1); all-zero rows -> NaN
    max_expr = np.nanmax(values, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = np.where(max_expr > 0, max_expr, np.nan)
        tau = (n_tissues - values.sum(axis=1) / denom) / (n_tissues - 1)
    tau[values.shape[1] <= 1] = np.nan

    result = pd.DataFrame({gene_col: df.iloc[:, 0], "tau": tau})
    result.to_csv(output_path, sep="\t", index=False, header=has_header)

    valid = int(np.isfinite(tau).sum())
    logging.info(
        "Results: %d genes (%d valid, %d all-zero/NA) -> %s",
        len(result), valid, len(result) - valid, output_path,
    )


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate tissue specificity index (tau) from an expression matrix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Tau ranges from 0 (ubiquitous expression) to 1 (single tissue).
Input format: gene_id <tab> tissue1 <tab> tissue2 <tab> ...

examples:
  python tissue_specificity_tau.py -i expression.tsv -o tau.tsv
  python tissue_specificity_tau.py -i expression.tsv -o tau.tsv --header
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input expression matrix TSV.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV with gene ID and tau.")
    parser.add_argument(
        "--header", action="store_true",
        help="Input file has a header row (output then carries the same header).",
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

    compute_tau(args.input, args.output, args.header)


if __name__ == "__main__":
    main()
