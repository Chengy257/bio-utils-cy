#!/usr/bin/env python3
"""
Batch render PDB files with PyMOL using AlphaFold pLDDT coloring.

Author: ChengYu
Created Time: 2026
"""

__version__ = "1.0.0"

import argparse
import glob
import logging
import os
import sys
import time
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


def setup_alphafold_colors() -> None:
    """Register AlphaFold pLDDT color palette in PyMOL."""
    from pymol import cmd

    cmd.set_color("af_very_high", [0.051, 0.341, 0.827])   # >90  deep blue
    cmd.set_color("af_confident", [0.416, 0.796, 0.945])   # 70-90 light blue
    cmd.set_color("af_low",       [0.996, 0.851, 0.212])   # 50-70 yellow
    cmd.set_color("af_very_low",  [0.992, 0.490, 0.302])   # <50   orange


def color_by_plddt(obj_name: str) -> None:
    """Color a PyMOL object by B-factor (pLDDT) ranges."""
    from pymol import cmd

    try:
        cmd.color("af_very_high", f"{obj_name} and b > 90")
        cmd.color("af_confident", f"{obj_name} and b > 70 and b <= 90")
        cmd.color("af_low",       f"{obj_name} and b > 50 and b <= 70")
        cmd.color("af_very_low",  f"{obj_name} and b <= 50")
    except Exception:
        logger.warning("pLDDT coloring failed for '%s'; falling back to spectrum.", obj_name)
        cmd.spectrum("b", "blue_white_red", obj_name, minimum=50, maximum=90)


def render_single(
    pdb_path: str,
    output_dir: str,
    width: int,
    height: int,
    dpi: int,
) -> bool:
    """Load one PDB, color by pLDDT, render, and return success flag."""
    from pymol import cmd

    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]

    try:
        cmd.delete("all")
        cmd.load(pdb_path, pdb_name)
        time.sleep(0.3)

        cmd.hide("everything")
        cmd.show("cartoon", pdb_name)
        color_by_plddt(pdb_name)

        cmd.bg_color("white")
        cmd.set("antialias", 1)
        cmd.center(pdb_name)
        cmd.orient(pdb_name)
        cmd.zoom(pdb_name, buffer=2)
        cmd.refresh()
        time.sleep(0.2)

        output_path = os.path.join(output_dir, f"{pdb_name}_plddt.png")
        cmd.png(output_path, width=width, height=height, dpi=dpi, ray=1)
        logger.info("Rendered: %s", output_path)
        return True

    except Exception as exc:
        logger.error("Failed to render %s: %s", pdb_path, exc)
        return False


def batch_render(
    pdb_folder: str,
    output_folder: str,
    width: int,
    height: int,
    dpi: int,
    pattern: str,
) -> None:
    """Discover PDB files and render each with AlphaFold pLDDT coloring."""
    import pymol
    from pymol import cmd

    os.makedirs(output_folder, exist_ok=True)

    pdb_files = sorted(glob.glob(os.path.join(pdb_folder, pattern)))
    if not pdb_files:
        logger.error("No PDB files matching '%s' found in %s", pattern, pdb_folder)
        sys.exit(1)

    logger.info("Found %d PDB files in %s", len(pdb_files), pdb_folder)

    # Launch PyMOL in headless mode
    pymol.finish_launching()
    time.sleep(1)
    setup_alphafold_colors()

    success = 0
    for pdb_file in pdb_files:
        if render_single(pdb_file, output_folder, width, height, dpi):
            success += 1

    logger.info("Rendered %d / %d structures.", success, len(pdb_files))
    cmd.quit()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="batch_pymol_render.py",
        description="Batch render PDB files with PyMOL using AlphaFold pLDDT coloring.",
        epilog=(
            "Examples:\n"
            "  %(prog)s -i ./pdb_files -o ./images\n"
            "  %(prog)s -i ./pdbs -o ./out --width 1200 --height 1200 --dpi 300\n"
            "  %(prog)s -i ./pdbs -o ./out --pattern '*.pdb'\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument("-i", "--input", required=True, help="Input directory with PDB files.")
    parser.add_argument("-o", "--output", default="alphafold_images", help="Output directory for images.")
    parser.add_argument("--width", type=int, default=800, help="Image width in pixels (default: 800).")
    parser.add_argument("--height", type=int, default=600, help="Image height in pixels (default: 600).")
    parser.add_argument("--dpi", type=int, default=150, help="Image DPI (default: 150).")
    parser.add_argument("--pattern", default="*.pdb", help="Glob pattern for PDB files (default: '*.pdb').")
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    batch_render(args.input, args.output, args.width, args.height, args.dpi, args.pattern)


if __name__ == "__main__":
    # The `pymol` python module ships with the PyMOL installation. When the
    # current interpreter lacks it, re-exec under the configured PyMOL binary
    # (BUC_PYMOL_BIN from config/env.sh); arguments after "--" are forwarded.
    try:
        import pymol  # noqa: F401
    except ImportError:
        _pymol_bin = os.environ.get("BUC_PYMOL_BIN", "")
        if _pymol_bin and os.path.isfile(_pymol_bin) and os.access(_pymol_bin, os.X_OK):
            os.execv(_pymol_bin, [_pymol_bin, "-cq", os.path.abspath(__file__), "--", *sys.argv[1:]])
        # PyMOL not configured: fall through and let tool functions report the gap.
    main()
