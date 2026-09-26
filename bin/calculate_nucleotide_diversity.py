#!/usr/bin/env python3
#########################################################################
# File Name: calculate_nucleotide_diversity.py
# Author: ChengYu
# Description: Calculate nucleotide diversity (Pi) for genomic
#              intervals using PLINK + VCFTools.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-26
#   - FIX: coordinate off-by-one at interval starts. BED starts are
#     0-based half-open while PLINK --extract range wants 1-based
#     inclusive positions; the start was passed through verbatim, so
#     every interval was shifted by one base. start+1 is now passed.
#   - FIX: a real BED6 file crashed pandas with a ParserError (only 4
#     column names were provided). Files with 3-6+ columns are now
#     accepted; the 4th column is the name when present, otherwise a
#     chrom:start-end name is generated.
#   - FIX: a vcftools failure was silently recorded as pi=0.0 —
#     indistinguishable from a genuinely invariant interval. Tool
#     failures now yield an empty (NA) row, a warning, and exit 1 at
#     the end; snp_count=0/pi=0.0 is reserved for real zero diversity.
#   - FIX: an unusable --plink-bin/--vcftools-bin or missing tool
#     raised a bare FileNotFoundError traceback; it is now a clean
#     error with exit 1.
#   - FIX: duplicate interval names shared one temp prefix (races under
#     the thread pool); unique per-row prefixes are used instead.
#   - FIX: plink/vcftools stdout+stderr now go to the debug log instead
#     of leaking into the terminal.
#   - --threads defaults to BUC_THREADS (config/env.sh) when set and is
#     validated; the vcftools stderr is included in failure messages.
#########################################################################
"""Calculate nucleotide diversity (Pi) for genomic intervals.

For each interval in a BED file, extracts SNPs using PLINK and computes
average site-level Pi using VCFTools. Supports multi-threaded processing.

Input is a BED file (3-6 columns; chrom, start [0-based, half-open],
end, name...). Positions are converted to PLINK's 1-based inclusive
coordinates internally.

Results TSV columns: name, snp_count, pi_value. A row with empty
snp_count/pi_value means the tools failed for that interval (the script
then exits 1); snp_count=0 with pi_value=0.0 means the interval simply
contains no variant sites.

Requires: plink, vcftools (PATH, BUC_PLINK_BIN/BUC_VCFTOOLS_BIN, or
--plink-bin/--vcftools-bin)
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

__version__ = "1.1.0"


def sanitize_filename(name: str) -> str:
    """Sanitize a string for safe use as a filename."""
    return re.sub(r"[^\w.-]", "_", name)


def find_tool(tool_name: str, explicit_path: str = None) -> str:
    """Find a tool executable.

    Args:
        tool_name: Tool name (e.g., 'plink', 'vcftools').
        explicit_path: Explicit path if provided by user.

    Returns:
        Path to the tool executable.

    Raises:
        FileNotFoundError: If tool not found.
    """
    if explicit_path:
        if os.path.isfile(explicit_path) and os.access(explicit_path, os.X_OK):
            return explicit_path
        raise FileNotFoundError(f"Specified {tool_name} path not executable: {explicit_path}")

    # Config-provided path (BUC_<TOOL>_BIN from config/env.sh) beats PATH
    env_path = os.environ.get(f"BUC_{tool_name.upper()}_BIN", "")
    if env_path:
        if os.path.isfile(env_path) and os.access(env_path, os.X_OK):
            return env_path
        raise FileNotFoundError(
            f"BUC_{tool_name.upper()}_BIN is set but not executable: {env_path}"
        )

    found = shutil.which(tool_name)
    if found:
        return found
    raise FileNotFoundError(
        f"{tool_name} not found in PATH. Install it or specify path via --{tool_name}-bin."
    )


def read_bed(bed_file: str):
    """Read a BED file with 3-6+ columns into (chrom, start, end, name) rows.

    Args:
        bed_file: Path to BED file.

    Returns:
        List of (chrom, start, end, name) tuples; names generated when
        the 4th column is absent. Malformed lines are skipped with a
        warning.

    Raises:
        ValueError: if no valid intervals remain.
    """
    rows = []
    n_malformed = 0
    with open(bed_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                n_malformed += 1
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                n_malformed += 1
                continue
            if end < start:
                n_malformed += 1
                continue
            name = fields[3] if len(fields) >= 4 and fields[3] else f"{fields[0]}:{start}-{end}"
            rows.append((fields[0], start, end, name))
    if n_malformed:
        logging.warning("Skipped %d malformed BED lines.", n_malformed)
    if not rows:
        raise ValueError(f"no valid BED intervals found in {bed_file}")
    return rows


def calculate_interval_pi(
    plink_prefix: str,
    plink_bin: str,
    vcftools_bin: str,
    chrom: str,
    start: int,
    end: int,
    name: str,
    temp_dir: str,
    unique_tag: str,
    threads: int = 1,
    memory: int = 4000,
) -> tuple:
    """Calculate average Pi for a single genomic interval.

    Args:
        plink_prefix: PLINK fileset prefix.
        plink_bin: Path to plink executable.
        vcftools_bin: Path to vcftools executable.
        chrom: Chromosome.
        start: BED start (0-based half-open).
        end: BED end (0-based exclusive).
        name: Interval name.
        temp_dir: Directory for temporary files.
        unique_tag: Unique string making temp prefixes collision-free.
        threads: PLINK threads.
        memory: PLINK memory in MB.

    Returns:
        Tuple of (name, snp_count, pi_value); snp_count/pi_value are
        None when a tool failed (as opposed to a genuinely empty
        interval, which yields (0, 0.0)).
    """
    # BED [start, end) -> PLINK 1-based inclusive [start+1, end].
    plink_start = start + 1

    prefix = os.path.join(temp_dir, unique_tag)
    range_file = f"{prefix}_range.txt"
    with open(range_file, "w") as fh:
        fh.write(f"{chrom}\t{plink_start}\t{end}\t{sanitize_filename(name)}\n")

    try:
        plink_cmd = [
            plink_bin,
            "--bfile", plink_prefix,
            "--threads", str(threads),
            "--memory", str(memory),
            "--extract", "range", range_file,
            "--recode", "vcf-iid",
            "--out", prefix,
        ]
        proc = subprocess.run(plink_cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            # PLINK exits nonzero when the range matches no variant at
            # all ("All variants excluded") — a legitimate empty interval.
            if "All variants excluded" in (proc.stderr or ""):
                logging.info("No variant sites in interval %s.", name)
                return name, 0, 0.0
            logging.error("PLINK failed for %s (rc=%d): %s", name, proc.returncode,
                          (proc.stderr or "").strip()[:500])
            return name, None, None
        logging.debug("plink stdout: %s", (proc.stdout or "")[:500])

        pi_prefix = f"{prefix}_pi"
        proc = subprocess.run(
            [vcftools_bin, "--vcf", f"{prefix}.vcf", "--site-pi", "--out", pi_prefix],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            logging.error("vcftools failed for %s (rc=%d): %s", name, proc.returncode,
                          (proc.stderr or "").strip()[:500])
            return name, None, None
        logging.debug("vcftools stderr: %s", (proc.stderr or "")[:500])

        pi_file = f"{pi_prefix}.sites.pi"
        with open(pi_file, "r") as fh:
            lines = fh.readlines()[1:]  # Skip header
            values = [float(line.strip().split()[-1]) for line in lines if line.strip()]
            if values:
                avg_pi = sum(values) / len(values)
                return name, len(values), avg_pi
            logging.info("No variant sites in interval %s.", name)
            return name, 0, 0.0
    except (OSError, ValueError) as e:
        logging.error("Interval %s failed: %s", name, e)
        return name, None, None


def process_bed(
    bed_file: str,
    plink_prefix: str,
    output_file: str,
    threads: int = 4,
    plink_threads: int = 1,
    memory: int = 4000,
    plink_bin: str = None,
    vcftools_bin: str = None,
) -> int:
    """Process all intervals in a BED file.

    Args:
        bed_file: Input BED file.
        plink_prefix: PLINK fileset prefix.
        output_file: Output TSV file.
        threads: Number of parallel workers.
        plink_threads: Threads per PLINK run.
        memory: PLINK memory in MB.
        plink_bin: Path to plink.
        vcftools_bin: Path to vcftools.

    Returns:
        Number of intervals whose tools failed (NA rows).
    """
    plink_path = find_tool("plink", plink_bin)
    vcftools_path = find_tool("vcftools", vcftools_bin)

    bed_rows = read_bed(bed_file)
    logging.info("Processing %d intervals with %d workers.", len(bed_rows), threads)

    results = []
    n_failed = 0
    with tempfile.TemporaryDirectory(prefix="pi_calc_") as temp_dir:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = [
                executor.submit(
                    calculate_interval_pi,
                    plink_prefix, plink_path, vcftools_path,
                    chrom, start, end, name,
                    temp_dir, f"interval_{i:06d}", plink_threads, memory,
                )
                for i, (chrom, start, end, name) in enumerate(bed_rows)
            ]
            for future in futures:
                try:
                    results.append(future.result())
                except Exception as e:
                    logging.error("Interval processing failed: %s", e)
                    results.append(("", None, None))

    results_df = pd.DataFrame(results, columns=["name", "snp_count", "pi_value"])
    results_df.to_csv(output_file, sep="\t", index=False)
    n_failed = int(results_df["pi_value"].isna().sum())
    if n_failed:
        logging.warning("%d/%d intervals failed (NA rows).", n_failed, len(results_df))
    logging.info("Results (%d intervals) -> %s", len(results_df), output_file)
    return n_failed


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate nucleotide diversity (Pi) for genomic intervals using PLINK + VCFTools.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python calculate_nucleotide_diversity.py -b regions.bed -p sample_data -o pi_results.tsv
  python calculate_nucleotide_diversity.py -b regions.bed -p sample_data -o pi.tsv -t 8 --memory 8000

notes:
  BED coordinates are 0-based half-open and are converted to PLINK's
  1-based inclusive range internally. Intervals where a tool fails get
  empty (NA) result cells and the script exits 1; snp_count=0 with
  pi_value=0.0 means genuinely no variant sites.
""",
    )
    parser.add_argument("-b", "--bed", type=str, required=True, help="Input BED file (3-6 columns).")
    parser.add_argument("-p", "--plink", type=str, required=True, help="PLINK fileset prefix (.bed/.bim/.fam).")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV file.")
    parser.add_argument(
        "-t", "--threads", type=int,
        default=int(os.environ.get("BUC_THREADS", "4") or 4),
        help="Number of parallel workers (default: BUC_THREADS or 4).",
    )
    parser.add_argument("--plink-threads", type=int, default=1, help="Threads per PLINK invocation (default: 1).")
    parser.add_argument("--memory", type=int, default=4000, help="PLINK memory limit in MB (default: 4000).")
    parser.add_argument("--plink-bin", type=str, default=None, help="Explicit path to plink executable.")
    parser.add_argument("--vcftools-bin", type=str, default=None, help="Explicit path to vcftools executable.")
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

    if args.threads < 1 or args.plink_threads < 1:
        parser.error("--threads and --plink-threads must be >= 1")

    if not Path(args.bed).is_file():
        logging.error("BED file not found: %s", args.bed)
        sys.exit(1)

    for ext in (".bed", ".bim", ".fam"):
        if not Path(args.plink + ext).is_file():
            logging.error("PLINK file not found: %s", args.plink + ext)
            sys.exit(1)

    try:
        n_failed = process_bed(
            bed_file=args.bed,
            plink_prefix=args.plink,
            output_file=args.output,
            threads=args.threads,
            plink_threads=args.plink_threads,
            memory=args.memory,
            plink_bin=args.plink_bin,
            vcftools_bin=args.vcftools_bin,
        )
    except (FileNotFoundError, ValueError) as e:
        logging.error("%s", e)
        sys.exit(1)

    if n_failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
