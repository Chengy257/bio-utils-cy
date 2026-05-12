#!/usr/bin/env python3
#########################################################################
# File Name: filter_long_introns.py
# Author: ChengYu
# Description: Filter transcripts with abnormally long introns
#              from GTF/GFF annotation files.
# Created Time: 2026
#########################################################################
"""Filter transcripts with abnormally long introns from GTF/GFF files.

Removes any transcript whose intron length exceeds the specified
threshold. Retains all other features unchanged.

Requires: BCBio.GFF (biopython), tqdm
"""

import argparse
import logging
import sys
from pathlib import Path

__version__ = "1.0.0"


def filter_gtf(input_file: str, output_file: str, max_intron: int) -> None:
    """Filter transcripts with introns longer than threshold.

    Args:
        input_file: Input GTF/GFF file.
        output_file: Output filtered file.
        max_intron: Maximum intron length in bp.
    """
    from BCBio import GFF
    from tqdm import tqdm

    total_tx = 0
    filtered_tx = 0

    with open(input_file, "r") as in_fh, open(output_file, "w") as out_fh:
        for rec in tqdm(GFF.parse(in_fh), desc="Processing", unit="seq"):
            new_feats = []
            for feat in rec.features:
                total_tx += 1
                exons = sorted(
                    [f.location for f in feat.sub_features if f.type == "exon"],
                    key=lambda x: x.start,
                )
                if len(exons) > 1:
                    intron_lengths = [
                        exons[i + 1].start - exons[i].end
                        for i in range(len(exons) - 1)
                    ]
                    if any(length > max_intron for length in intron_lengths):
                        filtered_tx += 1
                        continue
                new_feats.append(feat)
            rec.features = new_feats
            GFF.write([rec], out_fh)

    retained = total_tx - filtered_tx
    logging.info("Total: %d, Filtered: %d, Retained: %d", total_tx, filtered_tx, retained)
    logging.info("Output: %s", output_file)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Filter transcripts with abnormally long introns from GTF/GFF files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python filter_long_introns.py -i annotation.gtf -o filtered.gtf -m 20000
  python filter_long_introns.py -i annotation.gff -o filtered.gff -m 50000
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input GTF/GFF file.")
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

    try:
        from BCBio import GFF  # noqa: F401
    except ImportError:
        logging.error("Missing dependency: BCBio.GFF. Install with: pip install bcbio-gff")
        sys.exit(1)

    filter_gtf(args.input, args.output, args.max_intron)


if __name__ == "__main__":
    main()
