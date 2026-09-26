#!/usr/bin/env python3
#########################################################################
# File Name: filter_long_introns.py
# Author: ChengYu
# Description: Filter transcripts with abnormally long introns
#              from GTF annotation files.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-26
#   - REWRITE: pure Python GTF line filtering, no BCBio.GFF/tqdm.
#     BCBio round-trip rewriting of records (attribute reordering/
#     reformatting) is gone: retained lines are written verbatim.
#   - FIX: flat GTFs (gene and exon lines at the same level) made the
#     BCBio sub_feature tree empty, so nothing was ever filtered; the
#     tool now groups features by transcript_id and fails loudly when
#     the input carries no transcript_id at all instead of silently
#     passing everything through.
#   - gene lines whose transcripts were all removed are dropped; other
#     lines without transcript_id (comments etc.) are kept.
#   - exit 1 when the input carries no transcript_id at all (not a
#     transcript_id-based GTF) instead of silently passing everything
#     through.
#########################################################################
"""Filter transcripts with abnormally long introns from GTF files.

Groups feature lines by transcript_id, computes intron lengths from
consecutive exons (0-based half-open: next.start - prev.end), and
removes every transcript containing an intron longer than the
threshold. Retained lines are written verbatim (no reformatting).

gene lines (no transcript_id) are kept only if at least one of their
transcripts survives; comments and blank lines are always kept.

Requires: nothing beyond the Python standard library.
"""

import argparse
import logging
import sys
from pathlib import Path

__version__ = "1.1.0"


def _gtf_attr(attributes: str, key: str):
    """Return the value of a GTF attribute, or None if absent."""
    for attr in attributes.strip().split(";"):
        attr = attr.strip()
        if not attr:
            continue
        if " " in attr:
            k, _, v = attr.partition(" ")
        else:
            k, _, v = attr.partition("=")
        if k == key:
            return v.strip().strip('"')
    return None


def _classify_line(line: str):
    """Classify a GTF line.

    Returns:
        ("feature", transcript_id, gene_id, start, end) for lines with a
        transcript_id, ("gene", gene_id) for gene-level lines, or
        ("other", None) for comments/blank/unparseable lines.
    """
    if line.startswith("#") or not line.strip():
        return ("other", None)
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 9:
        return ("other", None)
    gene_id = _gtf_attr(fields[8], "gene_id")
    transcript_id = _gtf_attr(fields[8], "transcript_id")
    if transcript_id is not None:
        try:
            start, end = int(fields[3]), int(fields[4])
        except ValueError:
            return ("other", None)
        return ("feature", (transcript_id, gene_id, fields[2], start, end))
    return ("gene", gene_id)


def compute_max_intron(exons) -> int:
    """Max intron length from sorted exon (start, end) tuples; 0 if <2 exons."""
    exons = sorted(exons)
    return max(
        (exons[i + 1][0] - exons[i][1] for i in range(len(exons) - 1)),
        default=0,
    )


def filter_gtf(input_file: str, output_file: str, max_intron: int) -> None:
    """Filter transcripts with introns longer than threshold.

    Args:
        input_file: Input GTF file.
        output_file: Output filtered file.
        max_intron: Maximum intron length in bp.

    Raises:
        ValueError: if the input contains no transcript_id-bearing lines.
    """
    lines = []                # (text, classification)
    exons = {}                # transcript_id -> [(start, end)]
    gene_transcripts = {}     # gene_id -> set of transcript_ids
    seen_tx = set()
    n_unparseable = 0

    with open(input_file, "r") as fh:
        for line in fh:
            kind, payload = _classify_line(line)
            lines.append((line, kind, payload))
            if kind == "feature":
                tx_id, gene_id, ftype, start, end = payload
                seen_tx.add(tx_id)
                if gene_id:
                    gene_transcripts.setdefault(gene_id, set()).add(tx_id)
                if ftype == "exon":
                    exons.setdefault(tx_id, []).append((start, end))
            elif kind == "other" and line.strip() and not line.startswith("#"):
                n_unparseable += 1

    if not seen_tx:
        raise ValueError(
            "no transcript_id found in input; is this a transcript_id-based GTF?"
        )
    if n_unparseable:
        logging.warning("%d lines could not be parsed and are kept verbatim.", n_unparseable)

    kept_tx = {
        tx_id for tx_id in seen_tx
        if compute_max_intron(exons.get(tx_id, [])) <= max_intron
    }
    filtered_tx = len(seen_tx) - len(kept_tx)
    kept_genes = {
        gene_id for gene_id, txs in gene_transcripts.items() if txs & kept_tx
    }

    with open(output_file, "w") as out_fh:
        for line, kind, payload in lines:
            if kind == "feature":
                if payload[0] in kept_tx:
                    out_fh.write(line)
            elif kind == "gene":
                if payload in kept_genes:
                    out_fh.write(line)
            else:
                out_fh.write(line)

    logging.info(
        "Total transcripts: %d, Filtered: %d, Retained: %d",
        len(seen_tx), filtered_tx, len(kept_tx),
    )
    logging.info("Output: %s", output_file)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Filter transcripts with abnormally long introns from GTF files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python filter_long_introns.py -i annotation.gtf -o filtered.gtf -m 20000

notes:
  Lines are written verbatim; only whole transcripts (and gene lines
  left without any retained transcript) are removed.
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input GTF file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output filtered file.")
    parser.add_argument(
        "-m", "--max-intron", type=int, default=20000,
        help="Maximum intron length in bp (default: 20000).",
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

    if args.max_intron < 0:
        parser.error("--max-intron must be >= 0")

    try:
        filter_gtf(args.input, args.output, args.max_intron)
    except ValueError as e:
        logging.error("%s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
