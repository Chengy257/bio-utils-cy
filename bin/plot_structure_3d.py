#!/usr/bin/env python3
"""
Generate 3D structure images with hydropathy/secondary structure coloring
using DSSP and matplotlib.

Author: ChengYu
Created Time: 2026

Changelog:
  v1.1.0  2026-09-27
  - FIX: residues were keyed by residue number only, so multi-chain
    structures lost every same-numbered residue across chains and the
    backbone trace zig-zagged between chains; keys now include the chain.
  - FIX: residue coordinates came from the first ATOM record (usually N);
    CA coordinates are now preferred when present.
  - FIX: HETATM records (waters, ligands, ions) were plotted as residues;
    only ATOM records are parsed now.
  - FIX: the DSSP assignment map was keyed by residue number only and
    collided across chains; keyed by (chain, residue number) now.
  - FIX: 'ss' mode with DSSP missing produced an all-gray plot with only
    a warning; DSSP is now required for ss mode (hard error), and any
    per-file failure makes the batch exit non-zero.
  - FIX: DSSP detection honours $BUC_MKDSSP_BIN; unknown residues warn
    instead of silently scoring 0.0 hydropathy.
  - DOC: help no longer claims CIF support (fixed-column PDB parsing only).
  - CLEAN: dead Path/Tuple imports removed.
"""

__version__ = "1.1.0"

import argparse
import glob
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
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

    Only ATOM records are parsed (HETATM waters/ligands/ions would pollute
    the plot). Residues are keyed by (chain, residue number) so multi-chain
    structures keep their same-numbered residues; CA coordinates are
    preferred over the first backbone atom of the residue.

    Returns a list of dicts with keys: chain, resnum, resn, x, y, z,
    sorted by (chain, residue number).
    """
    residues: Dict[Any, Dict[str, Any]] = {}
    with open(filepath) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            resn = line[17:20].strip()
            chain = line[21].strip()
            resnum = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            atom_name = line[12:16].strip()
            key = (chain, resnum)
            if key not in residues:
                residues[key] = {
                    "chain": chain, "resnum": resnum, "resn": resn,
                    "x": x, "y": y, "z": z, "is_ca": atom_name == "CA",
                }
            elif atom_name == "CA" and not residues[key]["is_ca"]:
                residues[key].update({"x": x, "y": y, "z": z, "is_ca": True})
    parsed = list(residues.values())
    for r in parsed:
        r.pop("is_ca", None)
    return sorted(parsed, key=lambda r: (r["chain"], r["resnum"]))


# -----------------------------------------------------------------------
# DSSP integration
# -----------------------------------------------------------------------

def detect_dssp() -> str:
    """Locate the DSSP executable ($BUC_MKDSSP_BIN first, then PATH)."""
    env_path = os.environ.get("BUC_MKDSSP_BIN", "")
    if env_path and os.path.isfile(env_path) and os.access(env_path, os.X_OK):
        return env_path
    for name in ("mkdssp", "dssp"):
        path = shutil.which(name)
        if path:
            return path
    return ""


def run_dssp(pdb_path: str, dssp_bin: str) -> Dict[Any, str]:
    """Run DSSP on a PDB file and return {(chain, resnum): ss_code}."""
    fd, tmp_out = tempfile.mkstemp(suffix=".dssp")
    try:
        os.close(fd)
        subprocess.run(
            [dssp_bin, "--output-format", "dssp", pdb_path, tmp_out],
            check=True, capture_output=True,
        )
        ss_map: Dict[Any, str] = {}
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
                chain = line[11].strip()
                resnum = int(line[5:10].strip())
                ss_map[(chain, resnum)] = line[16].strip()
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
        values = []
        warned = set()
        for r in residues:
            if r["resn"] not in HYDROPATHY and r["resn"] not in warned:
                warned.add(r["resn"])
                logger.warning("Unknown residue '%s' scored as hydropathy 0.0", r["resn"])
            values.append(HYDROPATHY.get(r["resn"], 0.0))
        vmin, vmax = -4.5, 4.5
        cmap = plt.cm.coolwarm
        sc = ax.scatter(xs, ys, zs, c=values, cmap=cmap, vmin=vmin, vmax=vmax, s=20)
        fig.colorbar(sc, ax=ax, shrink=0.6, label="Kyte-Doolittle hydropathy")
    else:
        colors = []
        for r in residues:
            ss_code = ss_map.get((r["chain"], r["resnum"]), " ")
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
) -> int:
    """Render every matching structure; returns the number of failures."""
    os.makedirs(output_folder, exist_ok=True)
    pdb_files = sorted(glob.glob(os.path.join(input_folder, pattern)))
    if not pdb_files:
        logger.error("No structure files found in %s", input_folder)
        sys.exit(1)

    n_failed = 0
    for pdb_file in pdb_files:
        name = os.path.splitext(os.path.basename(pdb_file))[0]
        try:
            residues = parse_pdb(pdb_file)
            if not residues:
                logger.error("No residues parsed from %s", pdb_file)
                n_failed += 1
                continue

            ss_map: Dict[Any, str] = {}
            if color_mode == "ss":
                try:
                    ss_map = run_dssp(pdb_file, dssp_bin)
                except Exception as exc:
                    logger.error("DSSP failed for %s: %s", pdb_file, exc)
                    n_failed += 1
                    continue

            out_path = os.path.join(output_folder, f"{name}.png")
            plot_structure(residues, ss_map, out_path, color_mode, width, height, dpi)

        except Exception as exc:
            logger.error("Error processing %s: %s", pdb_file, exc)
            n_failed += 1

    logger.info("Batch rendering complete. Output in: %s", output_folder)
    return n_failed


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
    parser.add_argument("-i", "--input", required=True,
                        help="Input folder with PDB files (fixed-column parsing; CIF is not supported).")
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
        logger.error(
            "DSSP is required for --color ss but was not found. "
            "Install mkdssp, set BUC_MKDSSP_BIN, or supply --dssp <path>."
        )
        sys.exit(1)

    n_failed = process_files(
        args.input, args.output, dssp_bin, args.color,
        args.width, args.height, args.dpi, args.pattern,
    )
    if n_failed:
        logger.error("%d structure file(s) failed to render.", n_failed)
        sys.exit(1)


if __name__ == "__main__":
    main()
