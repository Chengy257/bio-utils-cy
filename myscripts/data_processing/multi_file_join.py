#!/usr/bin/env python3
#########################################################################
# File Name: multi_file_join.py
# Author: ChengYu
# Description: Join multiple TSV/CSV files by a common key column
#              using outer merge.
# Created Time: 2026
#########################################################################
"""Join multiple tabular files by a common key column.

Reads a list of file paths, extracts specified columns from each,
and merges them by an ID column using outer join (union of all keys).
"""

import argparse
import csv
import logging
import os
import sys
from pathlib import Path

import pandas as pd

__version__ = "1.0.0"


def join_files(
    file_list_path: str,
    columns: str,
    output_path: str,
    separator: str = "\t",
    join_how: str = "outer",
    na_rep: str = "NA",
) -> None:
    """Join multiple files on a common key column.

    Args:
        file_list_path: File containing one path per line.
        columns: Comma-separated 1-based column indices (e.g., '1,3').
        output_path: Output merged file.
        separator: Column separator (default: tab).
        join_how: Join method (outer, inner, left, right).
        na_rep: String representation for missing values.
    """
    col_indices = [int(x) - 1 for x in columns.split(",")]
    if len(col_indices) != 2:
        logging.error("Must specify exactly 2 columns (ID, value). Got: %s", columns)
        sys.exit(1)

    with open(file_list_path, "r") as fh:
        file_paths = [line.strip() for line in fh if line.strip()]

    if not file_paths:
        logging.error("No file paths found in %s", file_list_path)
        sys.exit(1)

    dfs = []
    for fpath in file_paths:
        if not os.path.isfile(fpath):
            logging.warning("File not found, skipping: %s", fpath)
            continue

        try:
            df = pd.read_csv(
                fpath, sep=separator, quoting=csv.QUOTE_NONE, quotechar=None,
            )
            selected = df.iloc[:, col_indices]
            name = os.path.basename(fpath)
            selected.columns = ["ID", name]
            dfs.append(selected)
            logging.debug("Loaded: %s (%d rows)", name, len(selected))
        except Exception as e:
            logging.warning("Failed to read %s: %s", fpath, e)

    if not dfs:
        logging.error("No valid files processed.")
        sys.exit(1)

    merged = dfs[0]
    for df in dfs[1:]:
        merged = pd.merge(merged, df, on="ID", how=join_how)

    merged.to_csv(output_path, index=False, sep=separator, na_rep=na_rep)
    logging.info("Merged %d files (%d rows, %d columns) -> %s", len(dfs), len(merged), len(merged.columns), output_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Join multiple TSV files by a common key column.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python multi_file_join.py -i file_list.txt -c 1,3 -o merged.tsv
  python multi_file_join.py -i file_list.txt -c 1,2 -o merged.tsv --how inner
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="File with one input path per line.")
    parser.add_argument("-c", "--columns", type=str, required=True, help="1-based column indices for ID,value (e.g., 1,3).")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output merged file.")
    parser.add_argument(
        "--how", type=str, default="outer",
        choices=["outer", "inner", "left", "right"],
        help="Join method (default: outer).",
    )
    parser.add_argument("--na-rep", type=str, default="NA", help="Missing value representation (default: NA).")
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
        logging.error("File list not found: %s", args.input)
        sys.exit(1)

    join_files(args.input, args.columns, args.output, join_how=args.how, na_rep=args.na_rep)


if __name__ == "__main__":
    main()
