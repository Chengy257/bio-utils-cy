#!/usr/bin/env python3
"""
Generate 3D structure images with hydropathy/secondary structure coloring
using DSSP and matplotlib.

Author: ChengYu
Created Time: 2026
"""

__version__ = "1.0.0"

import argparse
import glob
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)

# Kyte-Doolittle hydropathy scale (3-letter code -> score)
HYDROPATHY: Dict[str, float] = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5,
    "CYS": 2.5, "GLN": -3.5, "GLU": -3.5, "GLY": -0.4,
    "HIS": -3.2, "ILE": 4.5, "LEU": 3.8, "LYS": -3.9,
    "MET": 1.9, "PHE": 2.8, "PRO": -1.6, "SER": -0.8,
    "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}

# Map raw DSSP codes to colour-friendly labels
SS_MAP: Dict[str, str] = {
    "H": "Helix", "G": "Helix", "I": "Helix",
    "E": "Sheet", "B": "Sheet",
    "T": "Turn",  "S": "Turn",
}

SS_COLORS: Dict[str, str] = {
    "Helix": "#3366CC",
    "Sheet": "#E6A817",
    "Turn":  "#BFBFBF",
    "Coil":  "#BFBFBF",
}

# Three-letter -> one-letter
AA3TO1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


# -----------------------------------------------------------------------
# PDB parsing
# -----------------------------------------------------------------------

def parse_pdb(filepath: str) -> List[Dict[str, Any]]:
    """Extract per-residue coordinate and metadata from a PDB file.

    Returns a list of dicts with keys: resnum, resn, x, y, z.
    """
    residues: Dict[int, Dict[str, Any]] = {}
    with open(filepath) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resn = line[17:20].strip()
            chain = line[21].strip()
            resnum = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            if resnum not in residues:
                residues[resnum] = {"resnum": resnum, "resn": resn, "chain": chain, "x": x, "y": y, "z": z}
    return sorted(residues.values(), key=lambda r: r["resnum"])


# -----------------------------------------------------------------------
# DSSP integration
# -----------------------------------------------------------------------

def detect_dssp() -> str:
    """Try to locate the DSSP executable."""
    for name in ("mkdssp", "dssp"):
        path = shutil.which(name)
        if path:
            return path
    return ""


def run_dssp(pdb_path: str, dssp_bin: str) -> Dict[int, str]:
    """Run DSSP on a PDB file and return {resnum: ss_code}."""
    fd, tmp_out = tempfile.mkstemp(suffix=".dssp")
    try:
        os.close(fd)
        subprocess.run(
            [dssp_bin, "--output-format", "dssp", pdb_path, tmp_out],
            check=True, capture_output=True,
        )
        ss_map: Dict[int, str] = {}
        with open(tmp_out) as fh:
            started = False
            for line in fh:
                if line.startswith("  #  RESIDUE AA STRUCTURE"):
                    started = True
                    continue
                if not started or len(line) < 38:
                    continue
                aa = line[13].strip()
                if aa == "!" or not aa:
                    continue
                resnum = int(line[5:10].strip())
                ss_map[resnum] = line[16].strip()
        return ss_map
    finally:
        os.unlink(tmp_out)


# -----------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------

def plot_structure(
    residues: List[Dict[str, Any]],
    ss_map: Dict[int, str],
    output_path: str,
    color_mode: str,
    width: int,
    height: int,
    dpi: int,
) -> None:
    """Generate a 3D scatter/line plot of the protein structure.

    Parameters
    ----------
    residues : list[dict]
        Per-residue coordinate data.
    ss_map : dict[int, str]
        DSSP secondary-structure codes keyed by residue number.
    output_path : str
        Destination file for the figure.
    color_mode : str
        ``hydro`` for hydropathy coloring, ``ss`` for secondary-structure.
    """
    fig = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")

    xs = np.array([r["x"] for r in residues])
    ys = np.array([r["y"] for r in residues])
    zs = np.array([r["z"] for r in residues])
    resnums = [r["resnum"] for r in residues]

    # Backbone trace
    ax.plot(xs, ys, zs, color="gray", linewidth=0.8, alpha=0.5)

    if color_mode == "hydro":
        values = [HYDROPATHY.get(r["resn"], 0.0) for r in residues]
        vmin, vmax = -4.5, 4.5
        cmap = plt.cm.coolwarm
        sc = ax.scatter(xs, ys, zs, c=values, cmap=cmap, vmin=vmin, vmax=vmax, s=20)
        fig.colorbar(sc, ax=ax, shrink=0.6, label="Kyte-Doolittle hydropathy")
    else:
        colors = []
        for r in residues:
            ss_code = ss_map.get(r["resnum"], " ")
            label = SS_MAP.get(ss_code, "Coil")
            colors.append(SS_COLORS.get(label, "#BFBFBF"))
        ax.scatter(xs, ys, zs, c=colors, s=20)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    title = "Hydropathy" if color_mode == "hydro" else "Secondary Structure"
    ax.set_title(f"3D Structure ({title})")

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved: %s", output_path)


# -----------------------------------------------------------------------
# Batch processing
# -----------------------------------------------------------------------

def process_files(
    input_folder: str,
    output_folder: str,
    dssp_bin: str,
    color_mode: str,
    width: int,
    height: int,
    dpi: int,
    pattern: str,
) -> None:
    os.makedirs(output_folder, exist_ok=True)
    pdb_files = sorted(glob.glob(os.path.join(input_folder, pattern)))
    if not pdb_files:
        logger.error("No structure files found in %s", input_folder)
        sys.exit(1)

    for pdb_file in pdb_files:
        name = os.path.splitext(os.path.basename(pdb_file))[0]
        try:
            residues = parse_pdb(pdb_file)
            if not residues:
                logger.warning("No residues parsed from %s", pdb_file)
                continue

            ss_map: Dict[int, str] = {}
            if dssp_bin and color_mode == "ss":
                try:
                    ss_map = run_dssp(pdb_file, dssp_bin)
                except Exception as exc:
                    logger.warning("DSSP failed for %s: %s", pdb_file, exc)

            out_path = os.path.join(output_folder, f"{name}.png")
            plot_structure(residues, ss_map, out_path, color_mode, width, height, dpi)

        except Exception as exc:
            logger.error("Error processing %s: %s", pdb_file, exc)

    logger.info("Batch rendering complete. Output in: %s", output_folder)


# -----------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plot_structure_3d.py",
        description="Generate 3D structure images with hydropathy or "
                    "secondary-structure coloring using DSSP + matplotlib.",
        epilog=(
            "Examples:\n"
            "  %(prog)s -i ./pdb_files -o ./images --color hydro\n"
            "  %(prog)s -i ./pdb_files -o ./images --color ss\n"
            "  %(prog)s -i ./pdbs -o ./out --dssp /usr/bin/mkdssp --color ss --dpi 300\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument("-i", "--input", required=True, help="Input folder with PDB/CIF files.")
    parser.add_argument("-o", "--output", default="structure_images", help="Output folder for PNG images.")
    parser.add_argument(
        "--dssp", default=None,
        help="Path to DSSP executable (default: auto-detect).",
    )
    parser.add_argument(
        "--color", default="ss", choices=["hydro", "ss"],
        help="Color mode: 'hydro' (hydropathy) or 'ss' (secondary structure, default: ss).",
    )
    parser.add_argument("--width", type=int, default=2400, help="Image width in pixels (default: 2400).")
    parser.add_argument("--height", type=int, default=2400, help="Image height in pixels (default: 2400).")
    parser.add_argument("--dpi", type=int, default=600, help="Image DPI (default: 600).")
    parser.add_argument("--pattern", default="*.pdb", help="Glob pattern for structure files (default: '*.pdb').")
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    dssp_bin = args.dssp or detect_dssp()
    if args.color == "ss" and not dssp_bin:
        logger.warning(
            "DSSP not found. Secondary-structure coloring may not work. "
            "Install DSSP or supply --dssp <path>."
        )

    process_files(
        args.input, args.output, dssp_bin, args.color,
        args.width, args.height, args.dpi, args.pattern,
    )


if __name__ == "__main__":
    main()
