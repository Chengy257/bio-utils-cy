#!/usr/bin/env python3
#########################################################################
# File Name: calculate_nucleotide_diversity.py
# Author: ChengYu
# Description: Calculate nucleotide diversity (Pi) for genomic
#              intervals using PLINK + VCFTools.
# Created Time: 2026
#########################################################################
"""Calculate nucleotide diversity (Pi) for genomic intervals.

For each interval in a BED6 file, extracts SNPs using PLINK and computes
average site-level Pi using VCFTools. Supports multi-threaded processing.

Requires: plink, vcftools (must be in PATH or specified via --plink-bin / --vcftools-bin)
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

__version__ = "1.0.0"


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


def calculate_interval_pi(
    plink_prefix: str,
    plink_bin: str,
    vcftools_bin: str,
    chrom: str,
    start: int,
    end: int,
    name: str,
    temp_dir: str,
    threads: int = 1,
    memory: int = 4000,
) -> tuple:
    """Calculate average Pi for a single genomic interval.

    Args:
        plink_prefix: PLINK fileset prefix.
        plink_bin: Path to plink executable.
        vcftools_bin: Path to vcftools executable.
        chrom: Chromosome.
        start: Start position.
        end: End position.
        name: Interval name.
        temp_dir: Directory for temporary files.
        threads: PLINK threads.
        memory: PLINK memory in MB.

    Returns:
        Tuple of (name, snp_count, pi_value).
    """
    safe_name = sanitize_filename(name)
    prefix = os.path.join(temp_dir, safe_name)

    # Write range file
    range_file = f"{prefix}_range.txt"
    with open(range_file, "w") as fh:
        fh.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    # Extract SNPs with plink
    vcf_prefix = f"{prefix}"
    plink_cmd = [
        plink_bin,
        "--bfile", plink_prefix,
        "--threads", str(threads),
        "--memory", str(memory),
        "--extract", "range", range_file,
        "--recode", "vcf-iid",
        "--out", vcf_prefix,
    ]
    try:
        subprocess.run(plink_cmd, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.warning("PLINK failed for %s: %s", name, e.stderr.decode() if e.stderr else str(e))
        return name, 0, 0.0

    # Calculate Pi with vcftools
    vcf_file = f"{prefix}.vcf"
    pi_prefix = f"{prefix}_pi"
    vcftools_cmd = [
        vcftools_bin,
        "--vcf", vcf_file,
        "--site-pi",
        "--out", pi_prefix,
    ]
    subprocess.run(vcftools_cmd, capture_output=True)

    # Parse site Pi
    pi_file = f"{pi_prefix}.sites.pi"
    try:
        with open(pi_file, "r") as fh:
            lines = fh.readlines()[1:]  # Skip header
            if lines:
                values = [float(line.strip().split()[-1]) for line in lines if line.strip()]
                avg_pi = sum(values) / len(values) if values else 0.0
                return name, len(values), avg_pi
            else:
                logging.debug("No SNPs for interval %s.", name)
                return name, 0, 0.0
    except FileNotFoundError:
        logging.debug("Pi file missing for interval %s.", name)
        return name, 0, 0.0


def process_bed(
    bed_file: str,
    plink_prefix: str,
    output_file: str,
    threads: int = 4,
    plink_threads: int = 1,
    memory: int = 4000,
    plink_bin: str = None,
    vcftools_bin: str = None,
) -> None:
    """Process all intervals in a BED6 file.

    Args:
        bed_file: Input BED6 file.
        plink_prefix: PLINK fileset prefix.
        output_file: Output TSV file.
        threads: Number of parallel workers.
        plink_threads: Threads per PLINK run.
        memory: PLINK memory in MB.
        plink_bin: Path to plink.
        vcftools_bin: Path to vcftools.
    """
    plink_path = find_tool("plink", plink_bin)
    vcftools_path = find_tool("vcftools", vcftools_bin)

    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "name"])
    if bed_df.empty:
        logging.error("BED file is empty: %s", bed_file)
        sys.exit(1)

    logging.info("Processing %d intervals with %d workers.", len(bed_df), threads)

    results = []
    with tempfile.TemporaryDirectory(prefix="pi_calc_") as temp_dir:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = [
                executor.submit(
                    calculate_interval_pi,
                    plink_prefix, plink_path, vcftools_path,
                    row.chrom, row.start, row.end, row.name,
                    temp_dir, plink_threads, memory,
                )
                for row in bed_df.itertuples(index=False)
            ]
            for future in futures:
                try:
                    results.append(future.result())
                except Exception as e:
                    logging.error("Interval processing failed: %s", e)

    results_df = pd.DataFrame(results, columns=["name", "snp_count", "pi_value"])
    results_df.to_csv(output_file, sep="\t", index=False)
    logging.info("Results (%d intervals) -> %s", len(results_df), output_file)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate nucleotide diversity (Pi) for genomic intervals using PLINK + VCFTools.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python calculate_nucleotide_diversity.py -b regions.bed -p sample_data -o pi_results.tsv
  python calculate_nucleotide_diversity.py -b regions.bed -p sample_data -o pi.tsv -t 8 --memory 8000
""",
    )
    parser.add_argument("-b", "--bed", type=str, required=True, help="Input BED6 file with genomic intervals.")
    parser.add_argument("-p", "--plink", type=str, required=True, help="PLINK fileset prefix (.bed/.bim/.fam).")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV file.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of parallel workers (default: 4).")
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

    if not Path(args.bed).is_file():
        logging.error("BED file not found: %s", args.bed)
        sys.exit(1)

    for ext in (".bed", ".bim", ".fam"):
        if not Path(args.plink + ext).is_file():
            logging.error("PLINK file not found: %s", args.plink + ext)
            sys.exit(1)

    process_bed(
        bed_file=args.bed,
        plink_prefix=args.plink,
        output_file=args.output,
        threads=args.threads,
        plink_threads=args.plink_threads,
        memory=args.memory,
        plink_bin=args.plink_bin,
        vcftools_bin=args.vcftools_bin,
    )


if __name__ == "__main__":
    main()
