#!/usr/bin/env python3
#########################################################################
# File Name: calculate_dnds.py
# Author: ChengYu
# Description: Calculate dN/dS ratios for all species pairs from
#              multi-species aligned DNA FASTA files.
# Created Time: 2026
#
# Changelog:
#   v1.2.0  2026-09-29
#   - NEW: --method nei-gojobori (now the default): a real Nei-Gojobori
#     (1986) dN/dS — weighted synonymous/non-synonymous site counts,
#     unweighted pathway averaging for multi-step codon differences, and
#     a Jukes-Cantor correction. The dN/dS columns now carry these
#     per-site rates; the former differing-codon counts moved to new
#     dn_count/ds_count columns. Saturated components (p >= 0.75)
#     report NA.
#   - NEW: --method diff-ratio reproduces the previous per-codon
#     difference-count ratio byte-for-byte (same columns, same values).
#   - NOTE: NG86 assumes equal substitution rates (Jukes-Cantor) and
#     unweighted pathway averaging; for transition/transversion-aware or
#     maximum-likelihood estimates use PAML/HyPhy.
#   v1.1.0  2026-09-26
#   - FIX: identical sequences wrote ratio "Inf" (dn=ds=0); they now
#     write "NA".
#   - FIX: skipped pairs (protein length mismatch after stop-codon
#     truncation, translation errors) were logged at DEBUG and silently
#     missing from the table; they are now reported at INFO level and
#     counted, with a per-run summary.
#   - FIX: duplicate sequence IDs in one file silently overwrote each
#     other; a warning is emitted.
#   - DOC: (historical) the pre-1.2 metric was a per-codon difference
#     ratio with no site denominators, no pathway decomposition, and no
#     saturation correction; preserved verbatim as --method diff-ratio.
#########################################################################
"""Calculate dN/dS between all species pairs from aligned DNA FASTA files.

Processes a directory of ``*.dna.fa`` files, each containing aligned CDS
from multiple species, and reports dN, dS and dN/dS for every pair.

Methods (--method):

- ``nei-gojobori`` (default): Nei & Gojobori (1986). Synonymous and
  non-synonymous sites are counted per codon by weighting each
  position's single-step mutants; multi-step codon differences are
  decomposed over all unweighted mutational pathways; pS and pN get the
  Jukes-Cantor correction. dN/dS are per-site rates (typically 0.0x-1);
  saturated components (p >= 0.75) report NA. The former differing-codon
  counts are kept in the extra dn_count/ds_count columns.
- ``diff-ratio`` (legacy): counts of differing codons classified by
  amino-acid effect (dN = amino-acid-changing codons, dS = synonymous
  codons). No site denominators, no pathway decomposition, no
  saturation correction — a quick screen only.

Edge cases (both methods): identical pairs report ratio NA; pairs where
every difference is non-synonymous (with zero synonymous divergence)
report Inf; pairs whose translated proteins differ in length (e.g.
internal stop codons) are skipped and reported; a trailing partial
codon is ignored; codons containing gaps or stop/ambiguous bases are
skipped.

NG86 assumes equal substitution rates (Jukes-Cantor). For selection
inference use PAML/HyPhy.

Requires: biopython
"""

import argparse
import logging
import math
import os
import sys
from itertools import combinations, permutations
from pathlib import Path

from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq
from Bio.SeqIO import parse as seqio_parse

__version__ = "1.2.0"

METHODS = ("nei-gojobori", "diff-ratio")


def _aa(codon: str):
    """Amino acid of a codon via the standard table (None for stops/ambiguous)."""
    return standard_dna_table.forward_table.get(codon)


def count_codon_diffs(dna_seq1: str, dna_seq2: str) -> tuple:
    """Count differing codons classified by amino-acid effect (the legacy
    diff-ratio metric).

    Args:
        dna_seq1: First DNA sequence string (aligned, may contain gaps).
        dna_seq2: Second DNA sequence string.

    Returns:
        Tuple of (dn_count, ds_count): numbers of differing codons that
        change / keep the amino acid. Gap-containing and stop/ambiguous
        codons are skipped; a trailing partial codon is ignored.

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

        aa1 = _aa(codon1)
        aa2 = _aa(codon2)
        if aa1 is None or aa2 is None:
            continue

        if aa1 == aa2:
            ds += 1
        else:
            dn += 1

    return dn, ds


def calculate_dn_ds(dna_seq1: str, dna_seq2: str) -> tuple:
    """Legacy per-codon difference ratio (selected with --method diff-ratio).

    Returns:
        Tuple of (dN, dS, dN/dS ratio). ratio is float("inf") when
        dS == 0 and dN > 0, and None when both are 0 (identical
        sequences).
    """
    dn, ds = count_codon_diffs(dna_seq1, dna_seq2)
    if ds > 0:
        ratio = dn / ds
    elif dn > 0:
        ratio = float("inf")
    else:
        ratio = None
    return dn, ds, ratio


# --- Nei-Gojobori (1986) -----------------------------------------------------

_SYN_SITES_CACHE = {}
_PAIR_DIFF_CACHE = {}


def codon_syn_sites(codon: str) -> float:
    """Number of synonymous sites in a codon (Nei-Gojobori 1986): each
    position contributes the fraction of its single-step mutants that are
    synonymous; the three positions sum to a value in [0, 3]."""
    cached = _SYN_SITES_CACHE.get(codon)
    if cached is not None:
        return cached
    aa = _aa(codon)
    if aa is None:
        s = 0.0
    else:
        syn_mutants = 0.0
        for pos in range(3):
            for base in "ACGT":
                if base != codon[pos] and _aa(codon[:pos] + base + codon[pos + 1:]) == aa:
                    syn_mutants += 1.0
        s = syn_mutants / 3.0
    _SYN_SITES_CACHE[codon] = s
    return s


def codon_pair_diffs(codon1: str, codon2: str) -> tuple:
    """Synonymous/non-synonymous substitution counts between two codons,
    averaged over all unweighted mutational pathways (Nei-Gojobori 1986).

    Single-step differences are unambiguous. An intermediate stop codon
    is treated as an amino-acid state (a stop->stop step counts as
    synonymous). The result is symmetric in its arguments.

    Returns:
        Tuple of (S_d, N_d) summing to the number of nucleotide
        differences.
    """
    key = (codon1, codon2) if codon1 <= codon2 else (codon2, codon1)
    cached = _PAIR_DIFF_CACHE.get(key)
    if cached is not None:
        return cached

    diff_positions = [i for i in range(3) if codon1[i] != codon2[i]]
    if not diff_positions:
        result = (0.0, 0.0)
    else:
        n_paths = math.factorial(len(diff_positions))
        syn_total = 0.0
        cur = list(codon1)
        for perm in permutations(diff_positions):
            cur[:] = codon1
            aa_prev = _aa(codon1)
            for pos in perm:
                cur[pos] = codon2[pos]
                aa_new = _aa("".join(cur))
                if aa_new == aa_prev:
                    syn_total += 1.0
                aa_prev = aa_new
        s_d = syn_total / n_paths
        result = (s_d, float(len(diff_positions)) - s_d)

    _PAIR_DIFF_CACHE[key] = result
    return result


def jukes_cantor(p):
    """Jukes-Cantor correction d = -3/4 ln(1 - 4p/3).

    Returns None (NA) when p is undefined or saturated (p >= 0.75), and
    0.0 for p == 0.
    """
    if p is None:
        return None
    if p == 0:
        return 0.0
    if p >= 0.75:
        return None
    return -0.75 * math.log(1.0 - 4.0 * p / 3.0)


def calculate_dn_ds_ng(dna_seq1: str, dna_seq2: str) -> tuple:
    """Nei-Gojobori (1986) dN/dS between two aligned DNA sequences.

    Synonymous/nonsynonymous sites are accumulated per codon site
    (pair-averaged), multi-step codon differences are pathway-averaged,
    and pS/pN get the Jukes-Cantor correction.

    Returns:
        Tuple of (dN, dS, dN/dS ratio): per-site rates. ratio is None
        (NA) when either component is undefined (saturated p >= 0.75, no
        comparable sites) or when there are no differences at all, and
        float("inf") when dS == 0 while dN > 0.

    Raises:
        ValueError: if the sequences have different lengths.
    """
    if len(dna_seq1) != len(dna_seq2):
        raise ValueError("Sequences have different lengths.")

    s_sites = 0.0
    n_compared = 0
    s_diffs = 0.0
    n_diffs = 0.0

    for i in range(0, len(dna_seq1) - len(dna_seq1) % 3, 3):
        codon1 = dna_seq1[i:i + 3]
        codon2 = dna_seq2[i:i + 3]
        if "-" in codon1 or "-" in codon2:
            continue
        if _aa(codon1) is None or _aa(codon2) is None:
            continue
        s_pair = (codon_syn_sites(codon1) + codon_syn_sites(codon2)) / 2.0
        s_sites += s_pair
        n_compared += 1
        s_d, n_d = codon_pair_diffs(codon1, codon2)
        s_diffs += s_d
        n_diffs += n_d

    n_sites = 3.0 * n_compared - s_sites

    p_s = s_diffs / s_sites if s_sites > 0 else None
    p_n = n_diffs / n_sites if n_sites > 0 else None
    d_s = jukes_cantor(p_s)
    d_n = jukes_cantor(p_n)

    if d_s is None or d_n is None:
        ratio = None
    elif d_s == 0 and d_n == 0:
        ratio = None
    elif d_s == 0:
        ratio = float("inf")
    else:
        ratio = d_n / d_s
    return d_n, d_s, ratio


def process_file(dna_file: str, method: str = "nei-gojobori") -> tuple:
    """Process a single multi-species DNA FASTA file.

    Args:
        dna_file: Path to DNA FASTA file with aligned sequences.
        method: "nei-gojobori" (Nei & Gojobori 1986 rates) or
            "diff-ratio" (legacy differing-codon counts).

    Returns:
        Tuple of (results, n_skipped): results a list of
        (sp1, sp2, dN, dS, dN/dS or None, dn_count, ds_count) where
        dN/dS are NG per-site rates (nei-gojobori method) or the raw
        differing-codon counts (diff-ratio), and dn_count/ds_count
        always carry the raw counts; n_skipped the number of pairs
        skipped for protein-length mismatch or translation errors.
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
            dn_count, ds_count = count_codon_diffs(seq1, seq2)
            if method == "nei-gojobori":
                d_n, d_s, ratio = calculate_dn_ds_ng(seq1, seq2)
            else:
                d_n, d_s, ratio = calculate_dn_ds(seq1, seq2)
            results.append((sp1, sp2, d_n, d_s, ratio, dn_count, ds_count))
        except ValueError as e:
            logging.warning("dN/dS failed for %s/%s: %s", sp1, sp2, e)
            n_skipped += 1

    return results, n_skipped


def batch_process(input_dir: str, output_file: str, suffix: str = ".dna.fa",
                  method: str = "nei-gojobori") -> None:
    """Batch process all DNA FASTA files in a directory.

    Args:
        input_dir: Directory with .dna.fa files.
        output_file: Output TSV path.
        suffix: File suffix to match.
        method: "nei-gojobori" (default; dN/dS columns are NG per-site
            rates, raw counts in dn_count/ds_count) or "diff-ratio"
            (legacy layout: dN/dS columns ARE the raw counts).
    """
    files = sorted(f for f in os.listdir(input_dir) if f.endswith(suffix))
    if not files:
        logging.error("No files matching '*%s' in %s", suffix, input_dir)
        sys.exit(1)

    logging.info("Processing %d files (method: %s).", len(files), method)

    def fmt_num(v):
        if v is None:
            return "NA"
        if isinstance(v, float) and v == float("inf"):
            return "Inf"
        return f"{v:.4f}" if isinstance(v, float) else str(v)

    total_pairs = 0
    total_skipped = 0
    with open(output_file, "w") as fh:
        if method == "nei-gojobori":
            fh.write("Gene\tSpecies1\tSpecies2\tdN\tdS\tdN_dS\tdn_count\tds_count\n")
        else:
            fh.write("Gene\tSpecies1\tSpecies2\tdN\tdS\tdN_dS\n")

        for fname in files:
            fpath = os.path.join(input_dir, fname)
            # Strip the full suffix ("g1.dna.fa" -> "g1", not "g1.dna").
            gene_name = fname[: -len(suffix)] if suffix and fname.endswith(suffix) else os.path.splitext(fname)[0]
            logging.info("Processing: %s", fname)

            results, n_skipped = process_file(fpath, method)
            total_skipped += n_skipped
            for sp1, sp2, d_n, d_s, ratio, dn_count, ds_count in results:
                if method == "nei-gojobori":
                    fh.write(
                        f"{gene_name}\t{sp1}\t{sp2}\t{fmt_num(d_n)}\t{fmt_num(d_s)}\t"
                        f"{fmt_num(ratio)}\t{dn_count}\t{ds_count}\n"
                    )
                else:
                    fh.write(f"{gene_name}\t{sp1}\t{sp2}\t{d_n}\t{d_s}\t{fmt_num(ratio)}\n")
            total_pairs += len(results)

    if total_skipped:
        logging.info("Skipped %d pairs (protein length mismatch or translation error).", total_skipped)
    logging.info("Completed: %d pairs -> %s", total_pairs, output_file)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate dN/dS between all species pairs from aligned DNA FASTA files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python calculate_dnds.py -i gene_alignments/ -o dnds_results.tsv
  python calculate_dnds.py -i gene_alignments/ -o dnds_results.tsv --suffix .fa
  python calculate_dnds.py -i gene_alignments/ -o legacy.tsv --method diff-ratio

notes:
  Default method is nei-gojobori (Nei & Gojobori 1986): per-site rates
  with pathway decomposition and a Jukes-Cantor correction; saturated
  components (p >= 0.75) report NA. The legacy per-codon difference
  ratio (differing-codon counts by amino-acid effect — no site
  denominators, no saturation correction) is available as
  --method diff-ratio and reproduces the pre-1.2 output exactly.
  NG86 assumes equal substitution rates; use PAML/HyPhy for selection
  inference. Identical pairs report NA; pairs with protein-length
  mismatches (e.g. internal stops) are skipped and reported.
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input directory with multi-species DNA FASTA files.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV file.")
    parser.add_argument("--suffix", type=str, default=".dna.fa", help="File suffix to match (default: .dna.fa).")
    parser.add_argument(
        "--method", type=str, default="nei-gojobori", choices=list(METHODS),
        help="dN/dS method (default: nei-gojobori). 'nei-gojobori' reports "
             "Jukes-Cantor-corrected per-site rates (raw differing-codon "
             "counts in the extra dn_count/ds_count columns); 'diff-ratio' "
             "reproduces the legacy count-based ratio with the original "
             "columns.",
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

    if not Path(args.input).is_dir():
        logging.error("Input directory not found: %s", args.input)
        sys.exit(1)

    batch_process(args.input, args.output, args.suffix, method=args.method)


if __name__ == "__main__":
    main()
