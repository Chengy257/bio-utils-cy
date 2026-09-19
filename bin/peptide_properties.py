#!/usr/bin/env python3
"""
Calculate peptide physicochemical properties (MW, charge, pI, hydropathy,
instability, etc.) from a FASTA file.

Author: ChengYu
Created Time: 2026
"""

__version__ = "1.0.0"

import argparse
import csv
import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

logger = logging.getLogger(__name__)

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def analyze_sequence(sequence: str, record_id: str) -> Optional[Dict[str, Any]]:
    """Compute physicochemical properties for a single peptide sequence.

    Parameters
    ----------
    sequence : str
        Amino-acid sequence (one-letter codes).
    record_id : str
        FASTA header identifier.

    Returns
    -------
    dict or None
        Property dictionary, or *None* if the sequence is invalid.
    """
    if not sequence or not all(c in VALID_AA for c in sequence):
        logger.warning("Invalid sequence for '%s'; skipping.", record_id)
        return None

    analysis = ProteinAnalysis(sequence)

    ss_frac = analysis.secondary_structure_fraction()
    molar_ext = analysis.molar_extinction_coefficient()

    return {
        "id": record_id,
        "length": len(sequence),
        "molecular_weight": round(analysis.molecular_weight(), 4),
        "net_charge_pH7": round(analysis.charge_at_pH(7.0), 4),
        "isoelectric_point": round(analysis.isoelectric_point(), 4),
        "average_hydropathy": round(analysis.gravy(), 4),
        "aromaticity": round(analysis.aromaticity(), 4),
        "instability_index": round(analysis.instability_index(), 4),
        "helix_fraction": round(ss_frac[0], 4),
        "turn_fraction": round(ss_frac[1], 4),
        "sheet_fraction": round(ss_frac[2], 4),
        "molar_ext_reduced": molar_ext[0],
        "molar_ext_cystine": molar_ext[1],
        "amino_acid_composition": analysis.get_amino_acids_percent(),
    }


def flatten_for_csv(results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Expand nested amino-acid-composition dict into flat columns."""
    flat: List[Dict[str, Any]] = []
    for r in results:
        row = {k: v for k, v in r.items() if k != "amino_acid_composition"}
        comp = r.get("amino_acid_composition", {})
        for aa in sorted(VALID_AA):
            row[f"AA_{aa}"] = round(comp.get(aa, 0.0), 6)
        flat.append(row)
    return flat


def visualize_properties(results: List[Dict[str, Any]], output_prefix: str) -> None:
    """Plot distribution histograms of key peptide properties."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mw = [r["molecular_weight"] / 1000 for r in results]
    charges = [r["net_charge_pH7"] for r in results]
    pi = [r["isoelectric_point"] for r in results]
    gravy = [r["average_hydropathy"] for r in results]
    aroma = [r["aromaticity"] for r in results]
    instab = [r["instability_index"] for r in results]

    helix = [r["helix_fraction"] for r in results]
    turn = [r["turn_fraction"] for r in results]
    sheet = [r["sheet_fraction"] for r in results]

    # Aggregate amino-acid composition
    aa_totals: Dict[str, float] = {}
    for r in results:
        comp = r.get("amino_acid_composition", {})
        for aa, pct in comp.items():
            aa_totals[aa] = aa_totals.get(aa, 0.0) + pct

    fig, axes = plt.subplots(3, 3, figsize=(15, 12))

    axes[0, 0].hist(mw, bins=30, color="skyblue")
    axes[0, 0].set_title("Molecular Weight (kDa)")
    axes[0, 0].grid(axis="y", alpha=0.75)

    axes[0, 1].hist(charges, bins=30, color="lightgreen")
    axes[0, 1].set_title("Net Charge (pH 7)")
    axes[0, 1].grid(axis="y", alpha=0.75)

    axes[0, 2].hist(pi, bins=30, color="salmon")
    axes[0, 2].set_title("Isoelectric Point")
    axes[0, 2].grid(axis="y", alpha=0.75)

    axes[1, 0].hist(instab, bins=30, color="grey")
    axes[1, 0].set_title("Instability Index")
    axes[1, 0].grid(axis="y", alpha=0.75)

    axes[1, 1].bar(sorted(aa_totals.keys()), [aa_totals[a] for a in sorted(aa_totals.keys())], color="lightblue")
    axes[1, 1].set_title("Amino Acid Composition")
    axes[1, 1].grid(axis="y", alpha=0.75)

    axes[1, 2].bar(
        ["Helix", "Turn", "Sheet"],
        [sum(helix), sum(turn), sum(sheet)],
        color=["gold", "lightgreen", "lightcoral"],
    )
    axes[1, 2].set_title("Secondary Structure Fractions")
    axes[1, 2].grid(axis="y", alpha=0.75)

    axes[2, 0].hist(gravy, bins=30, color="gold")
    axes[2, 0].set_title("Average Hydropathy (GRAVY)")
    axes[2, 0].grid(axis="y", alpha=0.75)

    axes[2, 1].hist(aroma, bins=30, color="lightcoral")
    axes[2, 1].set_title("Aromaticity")
    axes[2, 1].grid(axis="y", alpha=0.75)

    fig.delaxes(axes[2, 2])

    plt.tight_layout()
    out_path = f"{output_prefix}_properties.pdf"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Property plots saved to: %s", out_path)


def process_fasta(
    fasta_file: str,
    output_csv: str,
    visualize: bool,
    threads: int,
) -> None:
    """Read FASTA, compute properties in parallel, write CSV."""
    fasta_path = Path(fasta_file)
    if not fasta_path.is_file():
        logger.error("FASTA file not found: %s", fasta_file)
        sys.exit(1)

    sequences: List[Tuple[str, str]] = [
        (rec.id, str(rec.seq)) for rec in SeqIO.parse(str(fasta_path), "fasta")
    ]
    logger.info("Loaded %d sequences from %s", len(sequences), fasta_file)

    results: List[Dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(analyze_sequence, seq, rid): rid
            for rid, seq in sequences
        }
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                results.append(result)

    if not results:
        logger.error("No valid sequences to process.")
        sys.exit(1)

    flat = flatten_for_csv(results)
    fieldnames = list(flat[0].keys())

    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(flat)
    logger.info("Results written to: %s (%d entries)", output_csv, len(flat))

    if visualize:
        visualize_properties(results, output_csv.rsplit(".", 1)[0])


# -----------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="peptide_properties.py",
        description="Calculate peptide physicochemical properties from a FASTA file.",
        epilog=(
            "Examples:\n"
            "  %(prog)s -i peptides.fasta -o properties.tsv\n"
            "  %(prog)s -i peptides.fasta -o props.tsv --visualize\n"
            "  %(prog)s -i input.fa -o output.tsv --threads 8 --log-level DEBUG\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file.")
    parser.add_argument("--visualize", action="store_true", help="Generate property distribution plots.")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads (default: 4).")
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    process_fasta(args.input, args.output, args.visualize, args.threads)


if __name__ == "__main__":
    main()
