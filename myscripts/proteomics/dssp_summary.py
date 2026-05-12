#!/usr/bin/env python3
"""
Batch analyze DSSP output files and summarize secondary structure statistics.

Author: ChengYu
Created Time: 2026
"""

__version__ = "1.0.0"

import argparse
import logging
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)

# Maximum accessible surface area reference (Tien et al. 2013)
MAX_ASA: Dict[str, float] = {
    "A": 106.0, "R": 248.0, "N": 157.0, "D": 163.0, "C": 135.0,
    "Q": 198.0, "E": 194.0, "G": 84.0,  "H": 184.0, "I": 169.0,
    "L": 164.0, "K": 205.0, "M": 188.0, "F": 197.0, "P": 136.0,
    "S": 130.0, "T": 142.0, "W": 227.0, "Y": 222.0, "V": 142.0,
}

# Map raw DSSP codes to three-state secondary structure
STRUCTURE_MAP: Dict[str, str] = {
    "H": "Helix", "G": "Helix", "I": "Helix",
    "E": "Sheet", "B": "Sheet",
    "T": "Turn",  "S": "Turn",
}

ALL_RAW_TYPES = ["H", "G", "I", "E", "B", "T", "S", " ", "-"]


def parse_dssp_file(filepath: str) -> Optional[Dict[str, Any]]:
    """Parse a single DSSP output file and extract secondary-structure statistics.

    Parameters
    ----------
    filepath : str
        Path to the DSSP output file.

    Returns
    -------
    dict or None
        Summary statistics, or *None* if parsing fails.
    """
    with open(filepath, "r") as fh:
        lines = fh.readlines()

    # Locate the data section
    data_start: Optional[int] = None
    for i, line in enumerate(lines):
        if line.startswith("  #  RESIDUE AA STRUCTURE"):
            data_start = i + 1
            break

    if data_start is None:
        logger.warning("Could not find data header in: %s", filepath)
        return None

    structure_counts: Dict[str, int] = {"Helix": 0, "Sheet": 0, "Turn": 0, "Coil": 0}
    raw_counter: Counter = Counter()
    relative_asa_list: List[float] = []
    total_residues = 0

    for line in lines[data_start:]:
        if len(line) < 38:
            continue

        aa = line[13].strip()
        struct = line[16].strip()
        acc_str = line[34:38].strip()

        if not aa or aa == "!":
            continue

        max_asa = MAX_ASA.get(aa.upper())
        if max_asa is None:
            continue

        try:
            acc = float(acc_str)
            relative_asa = acc / max_asa
            relative_asa_list.append(relative_asa)
        except ValueError:
            continue

        raw_counter[struct] += 1
        struct_class = STRUCTURE_MAP.get(struct, "Coil")
        structure_counts[struct_class] += 1
        total_residues += 1

    if total_residues == 0:
        logger.warning("No residues parsed in: %s", filepath)
        return None

    structure_ratios = {k: v / total_residues for k, v in structure_counts.items()}
    avg_relative_asa = (
        sum(relative_asa_list) / len(relative_asa_list) if relative_asa_list else 0.0
    )

    result: Dict[str, Any] = {
        "Filename": os.path.basename(filepath),
        "TotalResidues": total_residues,
        "AvgRelativeASA": round(avg_relative_asa, 4),
    }

    for k in ("Helix", "Sheet", "Turn", "Coil"):
        result[f"{k}_Count"] = structure_counts[k]
        result[f"{k}_Ratio"] = round(structure_ratios[k], 4)

    for raw_type in ALL_RAW_TYPES:
        label = f"Raw_{raw_type}_Count" if raw_type.strip() else "Raw_space_Count"
        result[label] = raw_counter.get(raw_type, 0)

    return result


def batch_process(folder_path: str, output_file: str) -> None:
    """Process all .dssp files in *folder_path* and write summary CSV.

    Parameters
    ----------
    folder_path : str
        Directory containing .dssp files.
    output_file : str
        Path for the output CSV summary.
    """
    logger.info("Analyzing directory: %s", folder_path)

    results: List[Dict[str, Any]] = []
    dssp_dir = Path(folder_path)

    for dssp_file in sorted(dssp_dir.glob("*.dssp")):
        logger.info("Processing: %s", dssp_file.name)
        result = parse_dssp_file(str(dssp_file))
        if result is not None:
            results.append(result)

    if not results:
        logger.error("No DSSP files were successfully parsed.")
        sys.exit(1)

    df = pd.DataFrame(results).fillna(0)
    df.to_csv(output_file, index=False)
    logger.info("Summary written to: %s (%d entries)", output_file, len(df))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dssp_summary.py",
        description="Batch analyze DSSP output files and summarize "
                    "secondary-structure statistics.",
        epilog=(
            "Examples:\n"
            "  %(prog)s -i ./dssp_output -o summary.csv\n"
            "  %(prog)s -i ./dssp_output -o results.csv --log-level DEBUG\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input directory containing .dssp files.",
    )
    parser.add_argument(
        "-o", "--output", default="dssp_summary.csv",
        help="Output CSV file path (default: dssp_summary.csv).",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    batch_process(args.input, args.output)


if __name__ == "__main__":
    main()
