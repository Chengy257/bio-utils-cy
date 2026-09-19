#!/usr/bin/env python3
#########################################################################
# File Name: tissue_specificity_tau.py
# Author: ChengYu
# Description: Calculate tissue specificity index (tau) for genes
#              from an expression matrix.
# Created Time: 2026
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

__version__ = "1.0.0"


def calc_tau(values: np.ndarray) -> float:
    """Calculate tissue specificity tau for a single gene.

    Args:
        values: Array of expression values across tissues.

    Returns:
        Tau value (0-1), or NaN if all values are zero.
    """
    values = np.array(values, dtype=float)
    max_expr = np.max(values)
    if max_expr == 0:
        return np.nan
    n = len(values)
    if n <= 1:
        return np.nan
    return np.sum(1 - (values / max_expr)) / (n - 1)


def compute_tau(input_path: str, output_path: str, has_header: bool = False) -> None:
    """Compute tau for all genes in an expression matrix.

    Args:
        input_path: Input TSV file (gene ID in col 1, expressions in cols 2+).
        output_path: Output TSV file.
        has_header: Whether input has a header row.
    """
    header = 0 if has_header else None
    df = pd.read_csv(input_path, sep="\t", header=header)

    # First column is gene ID
    gene_col = df.columns[0]
    expr_cols = df.columns[1:]

    if len(expr_cols) == 0:
        logging.error("No expression columns found (need at least 2 columns).")
        sys.exit(1)

    logging.info("Computing tau for %d genes across %d tissues.", len(df), len(expr_cols))

    df["tau"] = df[expr_cols].apply(calc_tau, axis=1)
    result = df[[gene_col, "tau"]]
    result.to_csv(output_path, sep="\t", index=False)

    valid = result["tau"].notna().sum()
    logging.info("Results: %d genes (%d valid, %d all-zero) -> %s", len(result), valid, len(result) - valid, output_path)


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
        help="Input file has a header row.",
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
