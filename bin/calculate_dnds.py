#!/usr/bin/env python3
#########################################################################
# File Name: calculate_dnds.py
# Author: ChengYu
# Description: Calculate dN/dS ratios for all species pairs from
#              multi-species aligned DNA FASTA files.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-26
#   - FIX: identical sequences wrote ratio "Inf" (dn=ds=0); they now
#     write "NA".
#   - FIX: skipped pairs (protein length mismatch after stop-codon
#     truncation, translation errors) were logged at DEBUG and silently
#     missing from the table; they are now reported at INFO level and
#     counted, with a per-run summary.
#   - FIX: duplicate sequence IDs in one file silently overwrote each
#     other; a warning is emitted.
#   - DOC: the metric is documented as a crude per-codon difference
#     ratio (counts of differing codons classified by amino-acid
#     effect). It is NOT a Nei-Gojobori-style dN/dS estimate: there are
#     no synonymous/non-synonymous site denominators, no multi-base
#     codon-path decomposition, and no substitution-model correction.
#     Values close to 1.0 for divergent pairs are expected with this
#     method and do not mean neutral evolution.
#########################################################################
"""Calculate dN/dS-style ratios for all species pairs from aligned DNA FASTA.

Processes a directory of ``*.dna.fa`` files, each containing aligned CDS
from multiple species. For every species pair, codons are compared
column-wise and classified by amino-acid effect:

- dN = number of differing codons that change the amino acid
- dS = number of differing codons that keep the amino acid
- ratio = dN / dS

IMPORTANT — semantics of this metric: this is a *per-codon difference
ratio*, not a substitution-based dN/dS estimate. It has no synonymous/
non-synonymous site denominators, no multi-base pathway decomposition,
and no correction for multiple substitutions at the same site. It is a
quick screen, not a selection test; use PAML/HyPhy for inference.

Edge cases: identical pairs report ratio NA; pairs where every codon
difference is non-synonymous report Inf; pairs whose translated
proteins differ in length (e.g. internal stop codons) are skipped and
reported; a trailing partial codon is ignored.

Requires: biopython
"""

import argparse
import logging
import os
import sys
from itertools import combinations
from pathlib import Path

from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq
from Bio.SeqIO import parse as seqio_parse

__version__ = "1.1.0"


def calculate_dn_ds(dna_seq1: str, dna_seq2: str) -> tuple:
    """Count synonymous/non-synonymous codon differences between two
    aligned DNA sequences.

    Args:
        dna_seq1: First DNA sequence string (aligned, may contain gaps).
        dna_seq2: Second DNA sequence string.

    Returns:
        Tuple of (dN, dS, dN/dS ratio). ratio is float("inf") when
        dS == 0 and dN > 0, and None when both are 0 (identical
        sequences).

    Raises:
        ValueError: if the sequences have different lengths.
    """
    if len(dna_seq1) != len(dna_seq2):
        raise ValueError("Sequences have different lengths.")

    dn = 0
    ds = 0

    for i in range(0, len(dna_seq1) - len(dna_seq1) % 3, 3):
        codon1 = dna_seq1[i:i + 3]
        codon2 = dna_seq2[i:i + 3]
        if codon1 == codon2:
            continue
        if "-" in codon1 or "-" in codon2:
            continue

        aa1 = standard_dna_table.forward_table.get(codon1)
        aa2 = standard_dna_table.forward_table.get(codon2)
        if aa1 is None or aa2 is None:
            continue

        if aa1 == aa2:
            ds += 1
        else:
            dn += 1

    if ds > 0:
        ratio = dn / ds
    elif dn > 0:
        ratio = float("inf")
    else:
        ratio = None
    return dn, ds, ratio


def process_file(dna_file: str) -> tuple:
    """Process a single multi-species DNA FASTA file.

    Args:
        dna_file: Path to DNA FASTA file with aligned sequences.

    Returns:
        Tuple of (results, n_skipped): results a list of
        (sp1, sp2, dN, dS, dN/dS or None), n_skipped the number of
        pairs skipped for protein-length mismatch or translation
        errors.
    """
    records = {}
    for rec in seqio_parse(dna_file, "fasta"):
        if rec.id in records:
            logging.warning("Duplicate sequence ID %r in %s; keeping the last.", rec.id, dna_file)
        records[rec.id] = str(rec.seq).upper()

    results = []
    n_skipped = 0

    for sp1, sp2 in combinations(records.keys(), 2):
        seq1 = records[sp1]
        seq2 = records[sp2]

        # Translate (gap-free) and compare protein lengths; internal stop
        # codons truncate the translation, so pseudogeneized copies fail
        # this check and are skipped.
        try:
            prot1 = Seq(seq1.replace("-", "")).translate(to_stop=True)
            prot2 = Seq(seq2.replace("-", "")).translate(to_stop=True)
        except Exception as e:
            logging.warning("Translation error for %s/%s: %s", sp1, sp2, e)
            n_skipped += 1
            continue
        if len(prot1) != len(prot2):
            logging.info(
                "Protein length mismatch for %s vs %s (%d vs %d aa), skipping pair.",
                sp1, sp2, len(prot1), len(prot2),
            )
            n_skipped += 1
            continue

        try:
            dn, ds, ratio = calculate_dn_ds(seq1, seq2)
            results.append((sp1, sp2, dn, ds, ratio))
        except ValueError as e:
            logging.warning("dN/dS failed for %s/%s: %s", sp1, sp2, e)
            n_skipped += 1

    return results, n_skipped


def batch_process(input_dir: str, output_file: str, suffix: str = ".dna.fa") -> None:
    """Batch process all DNA FASTA files in a directory.

    Args:
        input_dir: Directory with .dna.fa files.
        output_file: Output TSV path.
        suffix: File suffix to match.
    """
    files = sorted(f for f in os.listdir(input_dir) if f.endswith(suffix))
    if not files:
        logging.error("No files matching '*%s' in %s", suffix, input_dir)
        sys.exit(1)

    logging.info("Processing %d files.", len(files))

    total_pairs = 0
    total_skipped = 0
    with open(output_file, "w") as fh:
        fh.write("Gene\tSpecies1\tSpecies2\tdN\tdS\tdN_dS\n")

        for fname in files:
            fpath = os.path.join(input_dir, fname)
            # Strip the full suffix ("g1.dna.fa" -> "g1", not "g1.dna").
            gene_name = fname[: -len(suffix)] if suffix and fname.endswith(suffix) else os.path.splitext(fname)[0]
            logging.info("Processing: %s", fname)

            results, n_skipped = process_file(fpath)
            total_skipped += n_skipped
            for sp1, sp2, dn, ds, ratio in results:
                if ratio is None:
                    ratio_str = "NA"
                elif ratio == float("inf"):
                    ratio_str = "Inf"
                else:
                    ratio_str = f"{ratio:.4f}"
                fh.write(f"{gene_name}\t{sp1}\t{sp2}\t{dn}\t{ds}\t{ratio_str}\n")
            total_pairs += len(results)

    if total_skipped:
        logging.info("Skipped %d pairs (protein length mismatch or translation error).", total_skipped)
    logging.info("Completed: %d pairs -> %s", total_pairs, output_file)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate dN/dS-style ratios for all species pairs from aligned DNA FASTA files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python calculate_dnds.py -i gene_alignments/ -o dnds_results.tsv
  python calculate_dnds.py -i gene_alignments/ -o dnds_results.tsv --suffix .fa

notes:
  The ratio is a per-codon difference ratio (differing-codon counts by
  amino-acid effect), NOT a substitution-model dN/dS estimate: values
  near 1.0 for divergent pairs are expected with this method and do not
  indicate neutral evolution. Use PAML/HyPhy for selection inference.
  Identical pairs report NA; pairs with protein-length mismatches (e.g.
  internal stops) are skipped and reported.
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input directory with multi-species DNA FASTA files.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV file.")
    parser.add_argument("--suffix", type=str, default=".dna.fa", help="File suffix to match (default: .dna.fa).")
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

    if not Path(args.input).is_dir():
        logging.error("Input directory not found: %s", args.input)
        sys.exit(1)

    batch_process(args.input, args.output, args.suffix)


if __name__ == "__main__":
    main()
