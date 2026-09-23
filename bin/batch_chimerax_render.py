#!/usr/bin/env python3
"""
Batch render PDB files using ChimeraX command scripting.

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
from typing import List, Optional

logger = logging.getLogger(__name__)


def detect_chimerax() -> str:
    """Locate the ChimeraX executable: MYS_CHIMERAX_BIN (config/env.sh) first, then PATH."""
    env_path = os.environ.get("MYS_CHIMERAX_BIN", "")
    if env_path and os.path.isfile(env_path) and os.access(env_path, os.X_OK):
        return env_path
    for name in ("chimerax", "ChimeraX"):
        path = shutil.which(name)
        if path:
            return path
    return ""


def build_chimerax_script(
    pdb_files: List[str],
    output_folder: str,
    width: int,
    height: int,
    color_scheme: str,
) -> str:
    """Generate a ChimeraX .cxc script for batch rendering.

    Parameters
    ----------
    pdb_files : list[str]
        Absolute paths to PDB files.
    output_folder : str
        Directory to write rendered PNGs.
    width, height : int
        Image dimensions in pixels.
    color_scheme : str
        One of ``bychain``, ``bfactor``, ``ss``, ``rainbow``.

    Returns
    -------
    str
        Contents of the .cxc script.
    """
    lines = [
        "set bgColor white",
        "lighting soft",
        "graphics silhouettes true",
        "",
    ]

    for pdb_file in pdb_files:
        pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
        out_path = os.path.join(output_folder, f"{pdb_name}.png")

        lines.append(f"# --- {pdb_name} ---")
        lines.append(f"open {pdb_file}")
        lines.append("hide atoms")
        lines.append("show cartoons")

        if color_scheme == "bychain":
            lines.append("color bychain")
        elif color_scheme == "bfactor":
            lines.append("color byattribute bfactor")
        elif color_scheme == "ss":
            lines.append("color ss")
        elif color_scheme == "rainbow":
            lines.append("rainbow chain")
        else:
            lines.append("color bychain")

        lines.append("view")
        lines.append("lighting soft")
        lines.append(f"save {out_path} width {width} height {height} supersample 2")
        lines.append("close all")
        lines.append("")

    lines.append("# done")
    return "\n".join(lines)


def batch_render(
    pdb_folder: str,
    output_folder: str,
    chimerax_bin: str,
    width: int,
    height: int,
    color_scheme: str,
    pattern: str,
) -> None:
    """Discover PDB files and render them through ChimeraX."""
    os.makedirs(output_folder, exist_ok=True)

    pdb_files = sorted(glob.glob(os.path.join(pdb_folder, pattern)))
    if not pdb_files:
        logger.error("No PDB files matching '%s' found in %s", pattern, pdb_folder)
        sys.exit(1)

    logger.info("Found %d PDB files in %s", len(pdb_files), pdb_folder)

    # Use absolute paths so the script works regardless of cwd
    pdb_abs = [os.path.abspath(p) for p in pdb_files]
    out_abs = os.path.abspath(output_folder)

    script_text = build_chimerax_script(pdb_abs, out_abs, width, height, color_scheme)

    # Write to a temporary script file
    fd, script_path = tempfile.mkstemp(suffix=".cxc")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(script_text)

        cmd = [chimerax_bin, "--nogui", "--exit", script_path]
        logger.info("Launching ChimeraX: %s", " ".join(cmd))

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error("ChimeraX stderr:\n%s", result.stderr)
            sys.exit(1)

        logger.info("ChimeraX rendering complete. Output in: %s", output_folder)

    finally:
        if os.path.exists(script_path):
            os.unlink(script_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="batch_chimerax_render.py",
        description="Batch render PDB files using ChimeraX command scripting.",
        epilog=(
            "Examples:\n"
            "  %(prog)s -i ./pdb_files -o ./images\n"
            "  %(prog)s -i ./pdbs -o ./out --chimerax /usr/bin/chimerax --color bfactor\n"
            "  %(prog)s -i ./pdbs -o ./out --width 1200 --height 1200 --color ss\n"
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
    parser.add_argument("-o", "--output", default="structure_images", help="Output directory for images.")
    parser.add_argument(
        "--chimerax", default=None,
        help="Path to ChimeraX executable (default: auto-detect on PATH).",
    )
    parser.add_argument("--width", type=int, default=600, help="Image width in pixels (default: 600).")
    parser.add_argument("--height", type=int, default=600, help="Image height in pixels (default: 600).")
    parser.add_argument(
        "--color", default="bychain",
        choices=["bychain", "bfactor", "ss", "rainbow"],
        help="Coloring scheme (default: bychain).",
    )
    parser.add_argument("--pattern", default="*.pdb", help="Glob pattern for PDB files (default: '*.pdb').")
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    # Resolve ChimeraX binary
    chimerax_bin = args.chimerax or detect_chimerax()
    if not chimerax_bin:
        logger.error(
            "ChimeraX not found on PATH. Install ChimeraX or supply --chimerax <path>."
        )
        sys.exit(1)

    batch_render(
        args.input, args.output, chimerax_bin,
        args.width, args.height, args.color, args.pattern,
    )


if __name__ == "__main__":
    main()
