#!/usr/bin/env python3
#########################################################################
# File Name: merge_featurecounts.py
# Author: ChengYu
# Description: Merge featureCounts results from multiple samples into
#              unified count, FPKM, and TPM matrices.
# Created Time: 2026
#########################################################################
"""Merge featureCounts results into unified expression matrices.

Reads individual sample .count files (with id, effLength, counts, fpkm,
tpm columns) and merges them into count.matrix.tsv, FPKM, and TPM tables.
Also merges assignment log files.
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

__version__ = "1.0.0"

REQUIRED_COUNT_COLS = ["id", "effLength", "counts", "fpkm", "tpm"]


def read_count_file(path: str) -> pd.DataFrame:
    """Read a single featureCounts output file.

    Args:
        path: Path to .count TSV file.

    Returns:
        DataFrame with columns: id, effLength, counts, fpkm, tpm.
    """
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COUNT_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")
    return df


def merge_metric(count_files: list, metric: str, output_path: str) -> None:
    """Merge a specific metric across samples.

    Args:
        count_files: List of paths to .count files.
        metric: Column name to merge (counts, fpkm, tpm).
        output_path: Output TSV path.
    """
    merged = None
    for fpath in count_files:
        sample_name = Path(fpath).stem
        df = read_count_file(fpath)
        col = df[["id", metric]].rename(columns={metric: sample_name})
        if merged is None:
            merged = col
        else:
            merged = pd.merge(merged, col, on="id", how="outer")

    merged.to_csv(output_path, sep="\t", index=False)
    logging.info("Merged %s (%d genes, %d samples) -> %s", metric, len(merged), len(count_files), output_path)


def merge_log_files(log_files: list, output_path: str) -> None:
    """Merge featureCounts log files.

    Args:
        log_files: List of paths to .log files.
        output_path: Output TSV path.
    """
    rows = {}
    for fpath in log_files:
        sample = Path(fpath).stem
        with open(fpath) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    key = parts[0]
                    val = parts[1]
                    if key not in rows:
                        rows[key] = {}
                    rows[key][sample] = val

    df = pd.DataFrame.from_dict(rows, orient="index")
    df.index.name = "metric"
    df.to_csv(output_path, sep="\t")
    logging.info("Merged logs (%d metrics, %d samples) -> %s", len(df), len(log_files), output_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Merge featureCounts results from multiple samples.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Input directory should contain *.count and *.log files from featureCounts.
Each .count file must have columns: id, effLength, counts, fpkm, tpm.

examples:
  python merge_featurecounts.py -i featurecounts_output/ -o merged/
""",
    )
    parser.add_argument("-i", "--input-dir", type=str, required=True, help="Directory with *.count and *.log files.")
    parser.add_argument("-o", "--output-dir", type=str, default=None, help="Output directory (default: same as input).")
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

    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        logging.error("Input directory not found: %s", args.input_dir)
        sys.exit(1)

    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    count_files = sorted(input_dir.glob("*.count"))
    log_files = sorted(input_dir.glob("*.log"))

    if not count_files:
        logging.error("No *.count files found in %s", args.input_dir)
        sys.exit(1)

    logging.info("Found %d count files, %d log files.", len(count_files), len(log_files))

    # Extract effective lengths from first file
    first_df = read_count_file(str(count_files[0]))
    first_df[["id", "effLength"]].to_csv(
        output_dir / "effLength.txt", sep="\t", index=False
    )

    # Merge metrics
    merge_metric(count_files, "counts", output_dir / "count.matrix.tsv")
    merge_metric(count_files, "fpkm", output_dir / "GeneExpression_FPKM.xls")
    merge_metric(count_files, "tpm", output_dir / "GeneExpression_TPM.xls")

    # Merge logs
    if log_files:
        merge_log_files(log_files, output_dir / "GeneCount_Assigned_logs.xls")


if __name__ == "__main__":
    main()
