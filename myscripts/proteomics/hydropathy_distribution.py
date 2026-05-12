#!/usr/bin/env python3
"""
Calculate hydropathy distribution using sliding window analysis (Kyte-Doolittle).

Author: ChengYu
Created Time: 2026
"""

__version__ = "1.0.0"

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from Bio import SeqIO

logger = logging.getLogger(__name__)

# Kyte-Doolittle hydropathy scale
KYTE_DOOLITTLE: Dict[str, float] = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5,
    "F": 2.8, "G": -0.4, "H": -3.2, "I": 4.5,
    "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5,
    "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8,
    "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}


def calculate_hydropathy_distribution(
    sequence: str,
    window_size: int = 10,
) -> List[float]:
    """Calculate per-window mean hydropathy using the Kyte-Doolittle scale.

    Parameters
    ----------
    sequence : str
        Amino-acid sequence (one-letter codes).
    window_size : int
        Sliding-window size (must be >= 1).

    Returns
    -------
    list[float]
        Mean hydropathy score for each window position.
    """
    scores: List[float] = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i : i + window_size]
        score = float(np.mean([KYTE_DOOLITTLE.get(aa, 0.0) for aa in window]))
        scores.append(score)
    return scores


def plot_hydropathy(
    record_id: str,
    scores: List[float],
    window_size: int,
    output_dir: str,
    fmt: str,
) -> str:
    """Plot and save hydropathy distribution for one sequence.

    Returns the output file path.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    positions = range(len(scores))
    ax.plot(positions, scores, marker="o", markersize=2, linestyle="-", label=record_id)
    ax.axhline(0, color="red", linestyle="--", linewidth=0.8, label="Score = 0")
    ax.set_title(f"Hydropathy Distribution: {record_id} (window={window_size})")
    ax.set_xlabel("Window Position")
    ax.set_ylabel("Kyte-Doolittle Score")
    ax.legend()
    ax.grid(alpha=0.3)

    ext = fmt if fmt in ("pdf", "png", "svg") else "pdf"
    out_path = str(Path(output_dir) / f"hydropathy_{record_id}.{ext}")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def process_fasta(
    fasta_file: str,
    output_csv: str,
    window_size: int,
    visualize: bool,
    plot_dir: str,
    plot_fmt: str,
) -> None:
    """Read FASTA, compute hydropathy distributions, write CSV, optionally plot."""
    fasta_path = Path(fasta_file)
    if not fasta_path.is_file():
        logger.error("FASTA file not found: %s", fasta_file)
        sys.exit(1)

    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        logger.error("No sequences found in %s", fasta_file)
        sys.exit(1)

    logger.info("Loaded %d sequences from %s", len(records), fasta_file)

    rows: List[Dict[str, str]] = []
    for rec in records:
        seq = str(rec.seq).upper()
        scores = calculate_hydropathy_distribution(seq, window_size)
        rows.append({
            "id": rec.id,
            "length": str(len(seq)),
            "window_size": str(window_size),
            "num_windows": str(len(scores)),
            "mean_score": f"{np.mean(scores):.4f}" if scores else "NA",
            "hydropathy_scores": ",".join(f"{s:.4f}" for s in scores),
        })

        if visualize:
            out_plot = plot_hydropathy(rec.id, scores, window_size, plot_dir, plot_fmt)
            logger.info("Plot saved: %s", out_plot)

    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Results written to: %s (%d entries)", output_csv, len(rows))


# -----------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hydropathy_distribution.py",
        description="Calculate hydropathy distribution using sliding-window "
                    "analysis (Kyte-Doolittle) from a FASTA file.",
        epilog=(
            "Examples:\n"
            "  %(prog)s -i proteins.fasta -o hydropathy.csv\n"
            "  %(prog)s -i input.fa -o out.csv --window-size 15 --visualize\n"
            "  %(prog)s -i input.fa -o out.csv --plot-dir ./plots --plot-fmt png\n"
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
    parser.add_argument("-o", "--output", default="hydropathy.csv", help="Output CSV file (default: hydropathy.csv).")
    parser.add_argument(
        "-w", "--window-size", type=int, default=10,
        help="Sliding-window size (default: 10).",
    )
    parser.add_argument("--visualize", action="store_true", help="Generate per-sequence hydropathy plots.")
    parser.add_argument("--plot-dir", default="hydropathy_plots", help="Directory for plots (default: hydropathy_plots).")
    parser.add_argument(
        "--plot-fmt", default="pdf", choices=["pdf", "png", "svg"],
        help="Plot file format (default: pdf).",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    process_fasta(
        args.input, args.output, args.window_size,
        args.visualize, args.plot_dir, args.plot_fmt,
    )


if __name__ == "__main__":
    main()
