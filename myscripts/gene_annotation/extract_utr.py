#!/usr/bin/env python3
#########################################################################
# File Name: extract_utr.py
# Author: ChengYu
# Description: Extract 5' and 3' UTR sequences from GTF annotations
#              and a genome FASTA file.
# Created Time: 2026
#########################################################################
"""Extract 5' and 3' UTR sequences from GTF annotations.

Parses GTF to identify UTR regions (exon portions outside CDS),
extracts sequences from the genome FASTA, and outputs both a FASTA
file with UTR sequences and a CSV with UTR lengths.

Requires: biopython
"""

import argparse
import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path

__version__ = "1.0.0"


def parse_gtf(gtf_file: str) -> dict:
    """Parse GTF file into transcript structures.

    Args:
        gtf_file: Path to GTF annotation file.

    Returns:
        Dict mapping transcript_id -> {'exons': [...], 'cds': [...]}.
    """
    transcripts = defaultdict(lambda: {"exons": [], "cds": []})

    with open(gtf_file, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature_type = fields[0], fields[1], fields[2]
            start, end = int(fields[3]), int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            # Extract transcript_id
            transcript_id = None
            for attr in attributes.split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    parts = attr.split()
                    if len(parts) >= 2:
                        transcript_id = parts[1].strip('"')
                        break

            if transcript_id is None:
                continue

            if feature_type == "exon":
                transcripts[transcript_id]["exons"].append((chrom, start, end, strand))
            elif feature_type == "CDS":
                transcripts[transcript_id]["cds"].append((chrom, start, end, strand))

    logging.info("Parsed %d transcripts from GTF.", len(transcripts))
    return transcripts


def find_utr_regions(transcripts: dict) -> dict:
    """Identify 5' UTR and 3' UTR regions for each transcript.

    Args:
        transcripts: Output from parse_gtf().

    Returns:
        Dict mapping transcript_id -> {'5UTR': [...], '3UTR': [...]}.
    """
    utr_info = defaultdict(lambda: {"5UTR": [], "3UTR": []})

    for tx_id, data in transcripts.items():
        exons = sorted(data["exons"], key=lambda x: x[1])
        cds = sorted(data["cds"], key=lambda x: x[1])

        if not cds or not exons:
            continue

        cds_start = cds[0][1]
        cds_end = cds[-1][2]

        for chrom, exon_start, exon_end, strand in exons:
            if strand == "+":
                if exon_start < cds_start:
                    utr_info[tx_id]["5UTR"].append(
                        (chrom, exon_start, min(exon_end, cds_start - 1))
                    )
                if exon_end > cds_end:
                    utr_info[tx_id]["3UTR"].append(
                        (chrom, max(exon_start, cds_end + 1), exon_end)
                    )
            elif strand == "-":
                if exon_end > cds_end:
                    utr_info[tx_id]["5UTR"].append(
                        (chrom, max(exon_start, cds_end + 1), exon_end)
                    )
                if exon_start < cds_start:
                    utr_info[tx_id]["3UTR"].append(
                        (chrom, exon_start, min(exon_end, cds_start - 1))
                    )

    return utr_info


def extract_sequences(
    utr_info: dict,
    genome_fasta: str,
    output_fasta: str,
    output_csv: str,
) -> None:
    """Extract UTR sequences from genome and write outputs.

    Args:
        utr_info: Output from find_utr_regions().
        genome_fasta: Path to genome FASTA file.
        output_fasta: Output FASTA file path.
        output_csv: Output CSV with UTR lengths.
    """
    from Bio import SeqIO

    logging.info("Loading genome: %s", genome_fasta)
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    logging.info("Loaded %d chromosomes.", len(genome))

    n_5utr = 0
    n_3utr = 0

    with open(output_fasta, "w") as fasta_out, open(output_csv, "w", newline="") as csv_out:
        writer = csv.writer(csv_out)
        writer.writerow(["transcript_id", "5UTR_length", "3UTR_length"])

        for tx_id, utrs in utr_info.items():
            lengths = {"5UTR": 0, "3UTR": 0}

            for utr_type in ("5UTR", "3UTR"):
                for chrom, start, end in utrs[utr_type]:
                    if chrom not in genome:
                        logging.warning("Chromosome %s not in genome FASTA.", chrom)
                        continue
                    seq = genome[chrom].seq[start - 1:end]  # 1-based to 0-based
                    lengths[utr_type] += len(seq)
                    fasta_out.write(f">{tx_id}_{utr_type}\n{seq}\n")

            if lengths["5UTR"] > 0:
                n_5utr += 1
            if lengths["3UTR"] > 0:
                n_3utr += 1

            writer.writerow([tx_id, lengths["5UTR"], lengths["3UTR"]])

    logging.info("Extracted: %d transcripts with 5'UTR, %d with 3'UTR -> %s", n_5utr, n_3utr, output_fasta)
    logging.info("Length summary -> %s", output_csv)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Extract 5' and 3' UTR sequences from GTF and genome FASTA.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python extract_utr.py --gtf annotation.gtf --genome genome.fa -f utr.fa -c utr_lengths.csv
""",
    )
    parser.add_argument("--gtf", type=str, required=True, help="Input GTF annotation file.")
    parser.add_argument("--genome", type=str, required=True, help="Genome FASTA file.")
    parser.add_argument("-f", "--output-fasta", type=str, required=True, help="Output UTR FASTA file.")
    parser.add_argument("-c", "--output-csv", type=str, required=True, help="Output UTR length CSV file.")
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

    for path_label, path_val in [("GTF", args.gtf), ("Genome", args.genome)]:
        if not Path(path_val).is_file():
            logging.error("%s file not found: %s", path_label, path_val)
            sys.exit(1)

    transcripts = parse_gtf(args.gtf)
    utr_info = find_utr_regions(transcripts)
    extract_sequences(utr_info, args.genome, args.output_fasta, args.output_csv)


if __name__ == "__main__":
    main()
