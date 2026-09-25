#!/usr/bin/env python3
#########################################################################
# File Name: multi_file_join.py
# Author: ChengYu
# Description: Join multiple TSV/CSV files by a common key column
#              using outer merge.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-25
#   - FIX: all columns are read as strings, so numeric-looking IDs join
#     correctly across files whose ID columns pandas would otherwise type
#     differently (int64 vs object made matching rows silently disappear).
#   - FIX: files whose basenames collide no longer produce duplicate column
#     names (later files get a numbered suffix instead of a merge crash).
#   - FIX: invalid --columns values fail with a clean error instead of a
#     traceback.
#   - NEW: --no-header flag for headerless input files (previously the first
#     data row was always consumed as a header); --sep selects the column
#     separator (default tab; previously accepted only by the internal API).
#   - Duplicate IDs within a file are warned about (they multiply rows in
#     the join); files that fail to load are counted and reported, and the
#     exit code is 1 when any listed file was lost.
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

__version__ = "1.1.0"


def columns_type(spec: str) -> list:
    """Parse/validate the 1-based 'ID,value' column spec for argparse."""
    parts = spec.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"expected exactly 2 comma-separated column indices, got {spec!r}"
        )
    try:
        indices = [int(p) for p in parts]
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"column indices must be integers, got {spec!r}"
        ) from None
    if any(i < 1 for i in indices):
        raise argparse.ArgumentTypeError(
            f"column indices are 1-based, got {spec!r}"
        )
    return [i - 1 for i in indices]


def join_files(
    file_list_path: str,
    col_indices: list,
    output_path: str,
    separator: str = "\t",
    join_how: str = "outer",
    na_rep: str = "NA",
    no_header: bool = False,
) -> None:
    """Join multiple files on a common key column.

    Args:
        file_list_path: File containing one path per line.
        col_indices: 0-based column indices [ID, value].
        output_path: Output merged file.
        separator: Column separator (default: tab).
        join_how: Join method (outer, inner, left, right).
        na_rep: String representation for missing values.
        no_header: Input files have no header row.
    """
    with open(file_list_path, "r") as fh:
        file_paths = [line.strip() for line in fh if line.strip()]

    if not file_paths:
        logging.error("No file paths found in %s", file_list_path)
        sys.exit(1)

    dfs = []
    failed = []
    used_names = set()
    for fpath in file_paths:
        if not os.path.isfile(fpath):
            logging.warning("File not found, skipping: %s", fpath)
            failed.append(fpath)
            continue

        try:
            df = pd.read_csv(
                fpath,
                sep=separator,
                header=None if no_header else 0,
                dtype=str,
                quoting=csv.QUOTE_NONE,
                quotechar=None,
            )
            need = max(col_indices) + 1
            if df.shape[1] < need:
                raise ValueError(f"file has {df.shape[1]} column(s), need {need}")
            selected = df.iloc[:, col_indices]
            name = os.path.basename(fpath)
            if name in used_names:
                stem, ext = os.path.splitext(name)
                n = 2
                while f"{stem}_{n}{ext}" in used_names:
                    n += 1
                name = f"{stem}_{n}{ext}"
                logging.warning(
                    "Duplicate basename in file list; using column name %r for %s",
                    name, fpath,
                )
            used_names.add(name)
            selected.columns = ["ID", name]
            dup = selected["ID"].duplicated()
            if dup.any():
                sample = ", ".join(selected.loc[dup, "ID"].head(3))
                logging.warning(
                    "%s: %d duplicate ID value(s) (%s%s) will multiply rows in the join",
                    fpath, int(dup.sum()), sample, "..." if dup.sum() > 3 else "",
                )
            dfs.append(selected)
            logging.debug("Loaded: %s (%d rows)", name, len(selected))
        except Exception as e:
            logging.warning("Failed to read %s: %s", fpath, e)
            failed.append(fpath)

    if not dfs:
        logging.error("No valid files processed.")
        sys.exit(1)

    merged = dfs[0]
    for df in dfs[1:]:
        merged = pd.merge(merged, df, on="ID", how=join_how)

    merged.to_csv(output_path, index=False, sep=separator, na_rep=na_rep)
    logging.info(
        "Merged %d files (%d rows, %d columns) -> %s",
        len(dfs), len(merged), len(merged.columns), output_path,
    )

    if failed:
        logging.error(
            "%d of %d listed file(s) failed to load and are missing from the output",
            len(failed), len(file_paths),
        )
        sys.exit(1)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Join multiple TSV files by a common key column.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python multi_file_join.py -i file_list.txt -c 1,3 -o merged.tsv
  python multi_file_join.py -i file_list.txt -c 1,2 -o merged.tsv --how inner
  python multi_file_join.py -i file_list.txt -c 1,2 -o merged.tsv --no-header
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="File with one input path per line.")
    parser.add_argument("-c", "--columns", type=columns_type, required=True, metavar="ID,VALUE",
                        help="1-based column indices for ID,value (e.g., 1,3).")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output merged file.")
    parser.add_argument(
        "--how", type=str, default="outer",
        choices=["outer", "inner", "left", "right"],
        help="Join method (default: outer).",
    )
    parser.add_argument("--sep", type=str, default="\t",
                        help="Column separator of input/output files (default: tab).")
    parser.add_argument("--no-header", action="store_true",
                        help="Input files have NO header row; treat the first row as data "
                             "(default: the first row is treated as a header and dropped).")
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

    join_files(
        args.input, args.columns, args.output,
        separator=args.sep, join_how=args.how, na_rep=args.na_rep,
        no_header=args.no_header,
    )


if __name__ == "__main__":
    main()
