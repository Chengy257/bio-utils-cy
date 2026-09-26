#!/usr/bin/env python3
#########################################################################
# File Name: extract_utr.py
# Author: ChengYu
# Description: Extract 5' and 3' UTR sequences from GTF annotations
#              and a genome FASTA file.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-26
#   - FIX: negative-strand UTRs are now reverse-complemented; output
#     sequences always match the transcript strand (previously the
#     forward genomic slice was written, i.e. the wrong strand).
#   - FIX: multi-exon UTR fragments are joined in transcript order
#     (5'->3') into a single FASTA record per transcript/UTR type,
#     eliminating duplicate FASTA IDs for spliced UTRs.
#   - FIX: genome FASTA is now lazy-indexed (SeqIO.index) instead of
#     fully loaded into memory.
#   - FIX: malformed GTF lines (non-integer coordinates, unknown strand)
#     are skipped with a warning instead of raising a traceback.
#   - FIX: transcript_id attribute parsing handles unquoted values and
#     no longer mismatches keys like transcript_id_version.
#   - Transcripts with CDS/exon records on different chromosomes or
#     mixed strands are skipped with a warning.
#########################################################################
"""Extract 5' and 3' UTR sequences from GTF annotations.

Parses GTF to identify UTR regions (exon portions outside CDS),
extracts sequences from the genome FASTA, and outputs both a FASTA
file with UTR sequences and a CSV with UTR lengths.

UTR fragments spanning multiple exons are joined in transcript order
(5'->3') into a single FASTA record named {transcript_id}_{5UTR|3UTR};
negative-strand sequences are reverse-complemented so output always
matches the transcript strand. Transcripts without CDS records
(non-coding) are skipped.

Requires: biopython
"""

import argparse
import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

__version__ = "1.1.0"


def _gtf_attr(attributes: str, key: str) -> Optional[str]:
    """Return the value of a GTF attribute, or None if absent.

    Handles both 'key "value";' and 'key=value;' spellings; the key must
    match exactly (transcript_id does not match transcript_id_version).
    """
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


def parse_gtf(gtf_file: str) -> dict:
    """Parse GTF file into transcript structures.

    Args:
        gtf_file: Path to GTF annotation file.

    Returns:
        Dict mapping transcript_id -> {'exons': [...], 'cds': [...]},
        each record a (chrom, start, end, strand) tuple.
    """
    transcripts = defaultdict(lambda: {"exons": [], "cds": []})
    n_malformed = 0

    with open(gtf_file, "r") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                n_malformed += 1
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                n_malformed += 1
                continue
            if end < start or fields[6] not in ("+", "-"):
                n_malformed += 1
                continue

            chrom, feature_type, strand = fields[0], fields[2], fields[6]
            transcript_id = _gtf_attr(fields[8], "transcript_id")
            if transcript_id is None:
                continue

            if feature_type == "exon":
                transcripts[transcript_id]["exons"].append((chrom, start, end, strand))
            elif feature_type == "CDS":
                transcripts[transcript_id]["cds"].append((chrom, start, end, strand))

    if n_malformed:
        logging.warning("Skipped %d malformed GTF lines.", n_malformed)
    logging.info("Parsed %d transcripts from GTF.", len(transcripts))
    return transcripts


def find_utr_regions(transcripts: dict) -> dict:
    """Identify 5' UTR and 3' UTR regions for each transcript.

    UTR regions are returned in ascending genomic order; the caller is
    responsible for transcript-order joining and strand orientation.

    Args:
        transcripts: Output from parse_gtf().

    Returns:
        Dict mapping transcript_id -> {'5UTR': [...], '3UTR': [...],
        'strand': '+'|'-'}.
    """
    utr_info = {}
    n_no_cds = 0
    n_inconsistent = 0

    for tx_id, data in transcripts.items():
        exons = sorted(data["exons"], key=lambda x: x[1])
        cds = sorted(data["cds"], key=lambda x: x[1])

        if not cds or not exons:
            n_no_cds += 1
            continue

        if len({r[3] for r in exons} | {r[3] for r in cds}) != 1:
            logging.warning("Transcript %s has mixed strands; skipped.", tx_id)
            n_inconsistent += 1
            continue
        if {r[0] for r in exons} != {r[0] for r in cds}:
            logging.warning(
                "Transcript %s has exons and CDS on different chromosomes; skipped.",
                tx_id,
            )
            n_inconsistent += 1
            continue

        strand = exons[0][3]
        cds_start = cds[0][1]
        cds_end = cds[-1][2]

        left, right = [], []
        for chrom, exon_start, exon_end, _ in exons:
            if exon_start < cds_start:
                left.append((chrom, exon_start, min(exon_end, cds_start - 1), strand))
            if exon_end > cds_end:
                right.append((chrom, max(exon_start, cds_end + 1), exon_end, strand))

        if strand == "+":
            utr_info[tx_id] = {"5UTR": left, "3UTR": right, "strand": strand}
        else:
            utr_info[tx_id] = {"5UTR": right, "3UTR": left, "strand": strand}

    if n_no_cds:
        logging.info(
            "Skipped %d transcripts without both exon and CDS records (non-coding "
            "or incomplete).",
            n_no_cds,
        )
    if n_inconsistent:
        logging.warning("Skipped %d inconsistent transcripts.", n_inconsistent)
    return utr_info


def extract_sequences(
    utr_info: dict,
    genome_fasta: str,
    output_fasta: str,
    output_csv: str,
) -> None:
    """Extract UTR sequences from genome and write outputs.

    Fragments of one UTR are joined in transcript order (5'->3');
    negative-strand sequences are reverse-complemented. A UTR whose
    chromosome is missing from the genome FASTA is omitted entirely
    (no partial sequences) with a warning.

    Args:
        utr_info: Output from find_utr_regions().
        genome_fasta: Path to genome FASTA file.
        output_fasta: Output FASTA file path.
        output_csv: Output CSV with UTR lengths.
    """
    from Bio import SeqIO

    logging.info("Indexing genome: %s", genome_fasta)
    genome = SeqIO.index(genome_fasta, "fasta")

    n_5utr = 0
    n_3utr = 0

    with open(output_fasta, "w") as fasta_out, open(output_csv, "w", newline="") as csv_out:
        writer = csv.writer(csv_out)
        writer.writerow(["transcript_id", "5UTR_length", "3UTR_length"])

        for tx_id, utrs in utr_info.items():
            strand = utrs["strand"]
            lengths = {"5UTR": 0, "3UTR": 0}

            for utr_type in ("5UTR", "3UTR"):
                fragments = utrs[utr_type]
                if not fragments:
                    continue
                # Transcript order: ascending genomic on +, descending on -.
                fragments = sorted(
                    fragments, key=lambda f: f[1], reverse=(strand == "-")
                )
                parts = []
                for chrom, start, end, _ in fragments:
                    if chrom not in genome:
                        logging.warning(
                            "Chromosome %s not in genome FASTA; %s %s omitted.",
                            chrom, tx_id, utr_type,
                        )
                        parts = None
                        break
                    seq = genome[chrom].seq[start - 1:end]  # 1-based inclusive
                    if strand == "-":
                        seq = seq.reverse_complement()
                    parts.append(str(seq))
                if parts:
                    seq_str = "".join(parts)
                    lengths[utr_type] = len(seq_str)
                    fasta_out.write(f">{tx_id}_{utr_type}\n{seq_str}\n")

            if lengths["5UTR"] > 0:
                n_5utr += 1
            if lengths["3UTR"] > 0:
                n_3utr += 1

            writer.writerow([tx_id, lengths["5UTR"], lengths["3UTR"]])

    logging.info(
        "Extracted: %d transcripts with 5'UTR, %d with 3'UTR -> %s",
        n_5utr, n_3utr, output_fasta,
    )
    logging.info("Length summary -> %s", output_csv)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Extract 5' and 3' UTR sequences from GTF and genome FASTA.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python extract_utr.py --gtf annotation.gtf --genome genome.fa -f utr.fa -c utr_lengths.csv

notes:
  UTR fragments spanning multiple exons are joined in transcript order
  (5'->3'); negative-strand sequences are reverse-complemented. One
  FASTA record per transcript and UTR type ({tx_id}_5UTR/{tx_id}_3UTR).
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
