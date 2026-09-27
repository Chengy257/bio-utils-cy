#!/usr/bin/env python3
"""
Batch render PDB files with PyMOL using AlphaFold pLDDT coloring.

The B-factor column is interpreted as pLDDT (AlphaFold-style PDBs);
crystal structures will be colored by whatever their B-factors encode.

Author: ChengYu
Created Time: 2026

Changelog:
  v1.1.0  2026-09-27
  - FIX: when PyMOL was not importable and BUC_PYMOL_BIN was not
    configured, the script crashed with a bare ImportError traceback; it
    now exits with a clear error (and still re-execs under
    BUC_PYMOL_BIN when configured). --help/--version work without PyMOL.
  - FIX: the exit code was always 0 -- PyMOL swallows SystemExit and
    cmd.quit() terminates with status 0, so per-file render failures were
    never reflected. The script now terminates via os._exit with the real
    exit code (non-zero when any file failed).
  - FIX: pLDDT color selections used '<=', which is not valid PyMOL
    selection syntax -- every coloring raised and silently fell back to a
    generic spectrum; ranges are now expressed with 'not b > n'.
  - FIX: files PyMOL could not parse loaded as empty objects and were
    "rendered" as blank images; the atom count is now checked after load.
  - FIX: the BUC_PYMOL_BIN re-exec never actually ran main(): PyMOL
    executes scripts with __name__ == "pymol", so the __main__ guard
    failed and the child exited silently with status 0. main() is now
    invoked under both names.
  - DOC: B-factor = pLDDT assumption (AlphaFold PDBs) stated in the help.
  - CLEAN: dead Path import removed.
"""

__version__ = "1.1.0"

import argparse
import glob
import logging
import os
import sys
import time
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
    """Color a PyMOL object by B-factor (pLDDT) ranges.

    PyMOL selection syntax has no '<=' operator, so upper bounds are
    expressed as 'not b > n'.
    """
    from pymol import cmd

    cmd.color("af_very_high", f"{obj_name} and b > 90")
    cmd.color("af_confident", f"{obj_name} and b > 70 and not b > 90")
    cmd.color("af_low",       f"{obj_name} and b > 50 and not b > 70")
    cmd.color("af_very_low",  f"{obj_name} and not b > 50")


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
        if cmd.count_atoms(pdb_name) == 0:
            raise ValueError("no atoms loaded -- invalid or empty PDB file")
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
) -> int:
    """Discover PDB files and render each with AlphaFold pLDDT coloring.

    Returns the process exit code (0 = all rendered, 1 = any failure).
    """
    import pymol
    from pymol import cmd

    os.makedirs(output_folder, exist_ok=True)

    pdb_files = sorted(glob.glob(os.path.join(pdb_folder, pattern)))
    if not pdb_files:
        logger.error("No PDB files matching '%s' found in %s", pattern, pdb_folder)
        return 1

    logger.info("Found %d PDB files in %s", len(pdb_files), pdb_folder)

    # Launch PyMOL in headless mode
    pymol.finish_launching()
    time.sleep(1)
    setup_alphafold_colors()

    success = 0
    for pdb_file in pdb_files:
        if render_single(pdb_file, output_folder, width, height, dpi):
            success += 1

    failed = len(pdb_files) - success
    logger.info("Rendered %d / %d structures.", success, len(pdb_files))
    if failed:
        logger.error("%d file(s) failed to render.", failed)
    return 1 if failed else 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="batch_pymol_render.py",
        description="Batch render PDB files with PyMOL using AlphaFold pLDDT "
                    "coloring (B-factor = pLDDT, i.e. AlphaFold-style PDBs).",
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

    code = batch_render(args.input, args.output, args.width, args.height, args.dpi, args.pattern)

    # PyMOL swallows SystemExit and cmd.quit() always terminates with
    # status 0 -- os._exit is the only way to propagate a real exit code
    # (cmd.png writes synchronously, so nothing is pending).
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(code)


def _maybe_reexec_under_pymol() -> None:
    """Re-exec under the configured PyMOL binary when `import pymol` fails.

    PyMOL forwards everything after '--' to sys.argv. --help/--version are
    answered by argparse directly so they work on machines without PyMOL.
    """
    try:
        import pymol  # noqa: F401
        return
    except ImportError:
        pass

    if any(a in ("-h", "--help", "--version") for a in sys.argv[1:]):
        return

    pymol_bin = os.environ.get("BUC_PYMOL_BIN", "")
    if pymol_bin and os.path.isfile(pymol_bin) and os.access(pymol_bin, os.X_OK):
        os.execv(pymol_bin, [pymol_bin, "-cq", os.path.abspath(__file__), "--", *sys.argv[1:]])

    logging.basicConfig(level=logging.ERROR, format="[%(levelname)s] %(message)s")
    logger.error(
        "The PyMOL python module is not importable with this interpreter and "
        "BUC_PYMOL_BIN is not set or not executable. Source config/env.sh or "
        "install PyMOL (open-source PyMOL provides the `pymol` module)."
    )
    sys.exit(1)


if __name__ in ("__main__", "pymol"):
    # When re-exec'd under `pymol -cq script.py`, PyMOL runs the file with
    # __name__ == "pymol" (verified with PyMOL 2.5); the plain __main__
    # guard would never fire and main() would silently never run.
    _maybe_reexec_under_pymol()
    main()
