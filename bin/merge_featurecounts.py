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

Changelog:
  v1.1.0  2026-09-27
  - FIX: gene sets that differ across samples produced NaN cells in the
    counts matrix (downstream DESeq2 refuses non-integer input); the
    counts matrix now fills missing gene-sample pairs with 0, while
    FPKM/TPM keep NaN and warn with the affected gene count
  - FIX: effective lengths are now audited across all .count files —
    differing values for the same gene warn (first file still wins)
  - FIX: duplicate gene ids in any input file now fail with a clear
    message instead of silently cartesian-expanding the outer merge;
    unreadable files report the offending path instead of a bare
    pandas traceback
  - DOC: the .count contract is stated in the help text — this tool only
    merges the run-featurecounts.R output (featurecounts_pipeline.sh
    --fc-script); native featureCounts output (*.fc.tsv) is detected and
    rejected with guidance instead of "No *.count files found"
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

__version__ = "1.1.0"

REQUIRED_COUNT_COLS = ["id", "effLength", "counts", "fpkm", "tpm"]


def read_count_file(path: str) -> pd.DataFrame:
    """Read a single featureCounts output file.

    Args:
        path: Path to .count TSV file.

    Returns:
        DataFrame with columns: id, effLength, counts, fpkm, tpm.
    """
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:
        raise ValueError(f"Failed to read {path}: {exc}") from exc
    missing = [c for c in REQUIRED_COUNT_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")
    dup_ids = df.loc[df["id"].duplicated(), "id"].tolist()
    if dup_ids:
        shown = ", ".join(map(str, dup_ids[:5])) + ("..." if len(dup_ids) > 5 else "")
        raise ValueError(f"Duplicate gene ids in {path}: {shown}")
    return df


def merge_metric(count_files: list, metric: str, output_path: str) -> None:
    """Merge a specific metric across samples.

    Args:
        count_files: List of paths to .count files.
        metric: Column name to merge (counts, fpkm, tpm).
        output_path: Output TSV path.

    Gene sets may differ across samples (outer merge). For counts, missing
    gene-sample pairs are filled with 0 (a gene absent from a sample's
    .count file has no mapped reads there). For FPKM/TPM the NaN cells are
    kept and warned about — filling them would fabricate expression values.
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

    sample_cols = [c for c in merged.columns if c != "id"]
    n_missing = int(merged[sample_cols].isna().sum().sum())
    if n_missing:
        if metric == "counts":
            filled = merged[sample_cols].fillna(0)
            if bool((filled % 1 == 0).all().all()):
                filled = filled.astype(int)  # keep the matrix integer for DESeq2
            merged[sample_cols] = filled
            logging.warning(
                "Gene sets differ across samples: %d missing gene-sample pairs in the counts matrix were filled with 0.",
                n_missing,
            )
        else:
            genes = merged.loc[merged[sample_cols].isna().any(axis=1), "id"].tolist()
            shown = ", ".join(map(str, genes[:5])) + ("..." if len(genes) > 5 else "")
            logging.warning(
                "%s matrix keeps %d missing values for genes absent from some samples (e.g. %s); "
                "downstream tools must handle NaN.",
                metric, n_missing, shown,
            )

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
Input directory should contain *.count files produced by run-featurecounts.R
(i.e. featurecounts_pipeline.sh --fc-script), each with the columns:
id, effLength, counts, fpkm, tpm.

Note: native featureCounts output (*.fc.tsv, the pipeline's default branch)
is NOT supported here — convert it to the .count format first.

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
    fc_native = sorted(input_dir.glob("*.fc.tsv"))

    if fc_native and not count_files:
        logging.error(
            "Found %d native featureCounts output file(s) (*.fc.tsv) but no *.count files in %s. "
            "This tool merges the run-featurecounts.R .count format — rerun "
            "featurecounts_pipeline.sh with --fc-script, or convert the .fc.tsv files first.",
            len(fc_native), args.input_dir,
        )
        sys.exit(1)
    if fc_native:
        logging.warning("Ignoring %d native *.fc.tsv file(s) (only *.count files are merged).", len(fc_native))

    if not count_files:
        logging.error("No *.count files found in %s", args.input_dir)
        sys.exit(1)

    logging.info("Found %d count files, %d log files.", len(count_files), len(log_files))

    # Effective lengths: take them from the first file, but audit every
    # other file so annotation-version mismatches are not silent.
    first_df = read_count_file(str(count_files[0]))
    eff_ref = first_df.set_index("id")["effLength"]
    for fpath in count_files[1:]:
        df = read_count_file(str(fpath))
        eff_other = df.set_index("id")["effLength"]
        common = eff_ref.index.intersection(eff_other.index)
        n_diff = int((eff_ref.loc[common] != eff_other.loc[common]).sum())
        if n_diff:
            logging.warning(
                "effLength differs for %d shared gene(s) in %s; keeping the values from %s.",
                n_diff, fpath, count_files[0],
            )
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
