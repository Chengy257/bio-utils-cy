#!/usr/bin/env python3
# ==============================================================================
# File Name:    ena_ascp_download.py
# Author:       ChengYu
# Description:  Download FASTQ files from ENA (European Nucleotide Archive) using
#               Aspera ascp for high-speed transfer. Resolves download URLs via
#               the ENA API, supports parallel downloads, bandwidth control,
#               and automatic retry on failure.
# Created Time: 2026
# ==============================================================================

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import time
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Tuple

__version__ = "2.0.0"

logger = logging.getLogger("ena_ascp_download")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_ENA_HOST = "era-fasp@fasp.sra.ebi.ac.uk:"
DEFAULT_ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
DEFAULT_ASP_KEY = None  # will auto-detect
DEFAULT_BANDWIDTH = "500M"
DEFAULT_RETRIES = 3
DEFAULT_THREADS = 4
DEFAULT_TIMEOUT = 600  # seconds per download attempt


# ---------------------------------------------------------------------------
# Helper: find ascp binary and key
# ---------------------------------------------------------------------------
def find_ascp(path: Optional[str] = None) -> str:
    """Locate the ascp executable.

    Parameters
    ----------
    path : str or None
        Explicit path supplied by the user. If *None*, search the system PATH
        and common Aspera install locations.

    Returns
    -------
    str
        Absolute path to the ascp binary.

    Raises
    ------
    FileNotFoundError
        If ascp cannot be located.
    """
    if path:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return os.path.abspath(path)
        raise FileNotFoundError(f"Provided ascp path is not executable: {path}")

    # Try PATH first
    found = shutil.which("ascp")
    if found:
        return found

    # Common locations
    candidates = [
        os.path.expanduser("~/.aspera/connect/bin/ascp"),
        "/opt/aspera/connect/bin/ascp",
        "/usr/local/bin/ascp",
        "/opt/aspera/cli/bin/ascp",
    ]
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c

    raise FileNotFoundError(
        "Cannot find ascp. Install Aspera Connect/CLI or provide --ascp-bin."
    )


def find_ascp_key(ascp_bin: str) -> str:
    """Locate the Aspera private key that ships with the ascp installation.

    Parameters
    ----------
    ascp_bin : str
        Path to the ascp binary; the key is expected in the same ``bin/``
        directory or one level up under ``etc/``.

    Returns
    -------
    str
        Absolute path to the private key file.

    Raises
    ------
    FileNotFoundError
        If no suitable key file is found.
    """
    bindir = os.path.dirname(os.path.abspath(ascp_bin))
    rootdir = os.path.dirname(bindir)

    candidates = [
        os.path.join(bindir, "asperaworkshop.pem"),
        os.path.join(bindir, "aspera_id_rsa"),
        os.path.join(bindir, "aspera.openssh"),
        os.path.join(rootdir, "etc", "asperaworkshop.pem"),
        os.path.join(rootdir, "etc", "aspera_id_rsa"),
        os.path.join(rootdir, "etc", "aspera.openssh"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c

    raise FileNotFoundError(
        "Cannot find Aspera private key near ascp binary. "
        "Provide --ascp-key explicitly."
    )


# ---------------------------------------------------------------------------
# ENA API helpers
# ---------------------------------------------------------------------------
def resolve_ena_urls(
    accessions: List[str],
    api_url: str = DEFAULT_ENA_API,
    timeout: int = 30,
) -> List[Tuple[str, str, str]]:
    """Query the ENA file report API to retrieve FTP/ascp URLs for accessions.

    Parameters
    ----------
    accessions : list[str]
        List of SRR/ERR/DRR accession identifiers.
    api_url : str
        Base URL of the ENA file report API endpoint.
    timeout : int
        HTTP request timeout in seconds.

    Returns
    -------
    list[tuple[str, str, str]]
        Each tuple is ``(accession, fastq_ftp_url, fastq_aspera_url)``.
        URLs are semicolon-separated when paired-end.

    Raises
    ------
    RuntimeError
        If the API response cannot be parsed.
    """
    params = "&".join(
        [
            "accession=" + ",".join(accessions),
            "result=read_run",
            "fields=run_accession,fastq_ftp,fastq_aspera",
            "format=json",
            "limit=0",
        ]
    )
    url = f"{api_url}?{params}"
    logger.info("Querying ENA API: %s", url)

    req = urllib.request.Request(url)
    req.add_header("User-Agent", "ena_ascp_download/2.0")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode())
    except (urllib.error.URLError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"ENA API request failed: {exc}") from exc

    results: List[Tuple[str, str, str]] = []
    for entry in data:
        acc = entry.get("run_accession", "")
        ftp = entry.get("fastq_ftp", "")
        aspera = entry.get("fastq_aspera", "")
        if not acc or not aspera:
            logger.warning("Skipping %s: no aspera URL available", acc)
            continue
        results.append((acc, ftp, aspera))
    return results


def parse_aspera_urls(
    aspera_field: str,
) -> List[str]:
    """Split a semicolon-delimited aspera URL field into individual URLs.

    Parameters
    ----------
    aspera_field : str
        Raw ``fastq_aspera`` value from ENA, e.g.
        ``fasp.sra.ebi.ac.uk:/vol1/.../file_1.fastq.gz;fasp.../file_2.fastq.gz``.

    Returns
    -------
    list[str]
    """
    return [u.strip() for u in aspera_field.split(";") if u.strip()]


# ---------------------------------------------------------------------------
# Download via ascp
# ---------------------------------------------------------------------------
def run_ascp(
    url: str,
    outdir: str,
    ascp_bin: str,
    ascp_key: str,
    host: str,
    bandwidth: str,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
) -> str:
    """Execute a single ascp download with retry logic.

    Parameters
    ----------
    url : str
        The aspera URL (e.g. ``fasp.sra.ebi.ac.uk:/vol1/fastq/SRR...``).
    outdir : str
        Local directory to receive the file.
    ascp_bin : str
        Path to the ascp binary.
    ascp_key : str
        Path to the Aspera private key.
    host : str
        Remote host prefix, e.g. ``era-fasp@fasp.sra.ebi.ac.uk:``.
    bandwidth : str
        Bandwidth limit, e.g. ``500M``.
    timeout : int
        Per-attempt timeout in seconds.
    retries : int
        Maximum number of retry attempts.

    Returns
    -------
    str
        Path to the downloaded file.

    Raises
    ------
    RuntimeError
        If all retry attempts fail.
    """
    # Derive local filename from URL
    filename = url.rsplit("/", maxsplit=1)[-1]
    dest = os.path.join(outdir, filename)

    if os.path.isfile(dest):
        logger.info("File already exists, skipping: %s", dest)
        return dest

    # Build the remote source using configured host
    remote_path = url.split(":", 1)[-1] if ":" in url else url
    remote_src = f"{host}{remote_path}"

    cmd = [
        ascp_bin,
        "-T",               # disable encryption for speed
        "-l", bandwidth,
        "-P", "33001",      # UDP port
        "-i", ascp_key,
        "-Q",               # adaptive rate
        remote_src,
        outdir,
    ]

    logger.debug("ascp command: %s", " ".join(cmd))

    last_exc: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        logger.info(
            "Downloading %s (attempt %d/%d)", filename, attempt, retries
        )
        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=timeout,
                check=True,
            )
            if result.stderr:
                logger.debug("ascp stderr: %s", result.stderr.decode(errors="replace"))
            if os.path.isfile(dest):
                logger.info("Downloaded: %s", dest)
                return dest
            raise RuntimeError(
                f"ascp exited 0 but file not found: {dest}"
            )
        except subprocess.TimeoutExpired:
            logger.warning("Timeout on attempt %d for %s", attempt, filename)
            last_exc = RuntimeError(f"Timeout after {timeout}s for {filename}")
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
            logger.warning(
                "ascp failed (rc=%d) attempt %d: %s", exc.returncode, attempt, stderr
            )
            last_exc = exc

        if attempt < retries:
            wait = min(2**attempt, 30)
            logger.info("Retrying in %ds ...", wait)
            time.sleep(wait)

    raise RuntimeError(
        f"All {retries} attempts failed for {filename}: {last_exc}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="ena_ascp_download.py",
        description=(
            "Download FASTQ files from ENA using Aspera ascp. "
            "Accession URLs are resolved via the ENA API."
        ),
        epilog=(
            "Examples:\n"
            "  # Download a single run\n"
            "  %(prog)s -i SRR1234567 -o ./fastq\n"
            "\n"
            "  # Batch download from a file, 8 parallel, 1 Gbps limit\n"
            "  %(prog)s -i accession_list.txt -o ./fastq -t 8 -b 1G\n"
            "\n"
            "  # Specify ascp and key paths\n"
            "  %(prog)s -i SRR1234567 -o ./fastq --ascp-bin /opt/aspera/bin/ascp \\\n"
            "       --ascp-key /opt/aspera/etc/aspera_id_rsa\n"
            "\n"
            "  # Custom ENA host\n"
            "  %(prog)s -i SRR1234567 -o ./fastq --host era-fasp@fasp.sra.ebi.ac.uk:\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Single accession ID or path to a file with one accession per line")
    parser.add_argument("-o", "--outdir", default="./ena_download",
                        help="Output directory (default: ./ena_download)")
    parser.add_argument("-b", "--bandwidth", default=DEFAULT_BANDWIDTH,
                        help=f"Bandwidth limit for ascp (default: {DEFAULT_BANDWIDTH})")
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS,
                        help=f"Number of parallel downloads (default: {DEFAULT_THREADS})")
    parser.add_argument("--retries", type=int, default=DEFAULT_RETRIES,
                        help=f"Retry attempts per file (default: {DEFAULT_RETRIES})")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                        help=f"Per-download timeout in seconds (default: {DEFAULT_TIMEOUT})")
    parser.add_argument("--ascp-bin", default=None,
                        help="Path to ascp binary (default: auto-detect)")
    parser.add_argument("--ascp-key", default=None,
                        help="Path to Aspera private key (default: auto-detect)")
    parser.add_argument("--host", default=DEFAULT_ENA_HOST,
                        help=f"ENA Aspera host (default: {DEFAULT_ENA_HOST})")
    parser.add_argument("--api-url", default=DEFAULT_ENA_API,
                        help="ENA file report API URL")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Set logging verbosity (default: INFO)")
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {__version__}")
    return parser


def setup_logging(level: str) -> None:
    """Configure the root logger.

    Parameters
    ----------
    level : str
        Logging level string, e.g. ``'INFO'``.
    """
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s [%(levelname)-8s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_accessions(input_arg: str) -> List[str]:
    """Read accessions from a file or treat input as a single accession.

    Parameters
    ----------
    input_arg : str
        File path or a bare accession identifier.

    Returns
    -------
    list[str]
        Cleaned list of accession identifiers.
    """
    path = Path(input_arg)
    if path.is_file():
        accs = [
            line.strip()
            for line in path.read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
        logger.info("Read %d accession(s) from %s", len(accs), input_arg)
        return accs
    # Assume it is a single accession or comma-separated list
    accs = [a.strip() for a in input_arg.split(",") if a.strip()]
    return accs


def main() -> None:
    """Entry point for the ENA Aspera downloader."""
    parser = build_parser()
    args = parser.parse_args()
    setup_logging(args.log_level)

    logger.info("=== ena_ascp_download.py %s ===", __version__)

    # Resolve tool paths
    try:
        ascp_bin = find_ascp(args.ascp_bin)
        ascp_key = args.ascp_key or find_ascp_key(ascp_bin)
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        sys.exit(1)

    logger.info("ascp binary : %s", ascp_bin)
    logger.info("ascp key    : %s", ascp_key)
    logger.info("Host        : %s", args.host)
    logger.info("Bandwidth   : %s", args.bandwidth)

    # Read accessions
    accessions = read_accessions(args.input)
    if not accessions:
        logger.error("No accessions to process.")
        sys.exit(1)

    # Resolve URLs from ENA API
    try:
        records = resolve_ena_urls(accessions, api_url=args.api_url)
    except RuntimeError as exc:
        logger.error("Failed to resolve ENA URLs: %s", exc)
        sys.exit(1)

    if not records:
        logger.error("No downloadable URLs found for the given accessions.")
        sys.exit(1)

    # Prepare output directory
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    # Build download tasks
    tasks: List[Tuple[str, str]] = []  # (url, outdir)
    for acc, ftp, aspera in records:
        urls = parse_aspera_urls(aspera)
        for u in urls:
            tasks.append((u, outdir))
        logger.info("  %s: %d file(s)", acc, len(urls))

    total_tasks = len(tasks)
    logger.info("Total files to download: %d  |  Threads: %d", total_tasks, args.threads)

    # Execute downloads in parallel
    success = 0
    failed: List[str] = []
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(
                run_ascp,
                url=u,
                outdir=d,
                ascp_bin=ascp_bin,
                ascp_key=ascp_key,
                host=args.host,
                bandwidth=args.bandwidth,
                timeout=args.timeout,
                retries=args.retries,
            ): u
            for u, d in tasks
        }
        for future in as_completed(futures):
            url = futures[future]
            try:
                dest = future.result()
                success += 1
            except Exception as exc:
                logger.error("FAILED %s: %s", url, exc)
                failed.append(url)

    # Summary
    logger.info("=" * 50)
    logger.info("Download complete.")
    logger.info("  Successful : %d / %d", success, total_tasks)
    if failed:
        logger.warning("  Failed     : %d", len(failed))
        fail_log = os.path.join(outdir, "failed_downloads.txt")
        with open(fail_log, "w") as fh:
            for f in failed:
                fh.write(f + "\n")
        logger.warning("  Failed URLs written to %s", fail_log)
    logger.info("  Output dir : %s", outdir)


if __name__ == "__main__":
    main()
