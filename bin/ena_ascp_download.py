#!/usr/bin/env python3
# ==============================================================================
# File Name:    ena_ascp_download.py
# Author:       ChengYu
# Description:  Download FASTQ files from ENA (European Nucleotide Archive) using
#               Aspera ascp for high-speed transfer. Resolves download URLs via
#               the ENA API, supports parallel downloads, bandwidth control,
#               and automatic retry on failure.
# Created Time: 2026
# Changelog:
#   v2.1.0  2026-09-27
#   - FIX: batch with failures still exited 0 -- now exits non-zero and
#     keeps writing failed_downloads.txt
#   - FIX: 600 s per-attempt timeout was too short for multi-GB FASTQ;
#     default raised to 4 h, --timeout 0 disables the cap entirely
#   - FIX: truncated left-over files were treated as "already downloaded";
#     file sizes from the API (fastq_bytes) are now compared before skipping
#   - FIX: ENA API resolved all accessions in a single unchunked GET and
#     had no retry; now batched (500/batch) with 3 retry attempts each
#   - FIX: aspera URLs that already carry a user@host are used as-is
#     instead of being redirected to --host
#   - FIX: -P 33001 comment (TCP/SSH control port, not UDP)
#   - CLEAN: dropped dead ftp field / ET import / DEFAULT_ASP_KEY;
#     read_accessions dedupes and warns on non-run accessions;
#     --threads honours $BUC_THREADS
# ==============================================================================

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Tuple

__version__ = "2.1.0"

logger = logging.getLogger("ena_ascp_download")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_ENA_HOST = "era-fasp@fasp.sra.ebi.ac.uk:"
DEFAULT_ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
DEFAULT_BANDWIDTH = "500M"
DEFAULT_RETRIES = 3
DEFAULT_THREADS = 4
DEFAULT_TIMEOUT = 14400  # seconds per download attempt (4 h; 0 = unlimited)
ENA_API_CHUNK = 500      # accessions per ENA portal API request

_RUN_ACC_RE = re.compile(r"^[SED]RR\d+$")


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

    # Config-provided path (BUC_ASCP_BIN from config/env.sh) beats PATH
    env_path = os.environ.get("BUC_ASCP_BIN", "")
    if env_path:
        if os.path.isfile(env_path) and os.access(env_path, os.X_OK):
            return os.path.abspath(env_path)
        raise FileNotFoundError(f"BUC_ASCP_BIN is set but not executable: {env_path}")

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
def _fetch_json_with_retry(url: str, timeout: int, retries: int):
    """GET *url* and parse the JSON body, retrying transient failures."""
    last_exc: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url)
            req.add_header("User-Agent", f"ena_ascp_download/{__version__}")
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except (urllib.error.URLError, json.JSONDecodeError, OSError) as exc:
            last_exc = exc
            logger.warning("ENA API attempt %d/%d failed: %s", attempt, retries, exc)
            if attempt < retries:
                wait = min(2 ** attempt, 30)
                logger.info("Retrying ENA API in %ds ...", wait)
                time.sleep(wait)
    raise RuntimeError(
        f"ENA API request failed after {retries} attempts: {last_exc}"
    ) from last_exc


def resolve_ena_urls(
    accessions: List[str],
    api_url: str = DEFAULT_ENA_API,
    timeout: int = 30,
    chunk_size: int = ENA_API_CHUNK,
    retries: int = DEFAULT_RETRIES,
) -> List[Tuple[str, str, str]]:
    """Query the ENA file report API to retrieve aspera URLs for accessions.

    Accessions are sent in batches of *chunk_size* (a single GET with
    thousands of accessions exceeds URL limits); each batch is retried.

    Parameters
    ----------
    accessions : list[str]
        List of SRR/ERR/DRR accession identifiers.
    api_url : str
        Base URL of the ENA file report API endpoint.
    timeout : int
        HTTP request timeout in seconds.
    chunk_size : int
        Accessions per API request.
    retries : int
        Attempts per API request.

    Returns
    -------
    list[tuple[str, str, str]]
        Each tuple is ``(accession, fastq_aspera, fastq_bytes)``; the aspera
        URL and byte-size fields are semicolon-separated when paired-end.
        Entries without an aspera URL are skipped with a warning.

    Raises
    ------
    RuntimeError
        If an API request fails after all retries.
    """
    results: List[Tuple[str, str, str]] = []
    total = len(accessions)
    for start in range(0, total, chunk_size):
        chunk = accessions[start:start + chunk_size]
        params = "&".join(
            [
                "accession=" + ",".join(chunk),
                "result=read_run",
                "fields=run_accession,fastq_aspera,fastq_bytes",
                "format=json",
                "limit=0",
            ]
        )
        url = f"{api_url}?{params}"
        logger.info(
            "Querying ENA API (accessions %d-%d of %d)",
            start + 1, min(start + chunk_size, total), total,
        )
        logger.debug("ENA API URL: %s", url)
        data = _fetch_json_with_retry(url, timeout=timeout, retries=retries)

        for entry in data:
            acc = entry.get("run_accession", "")
            aspera = entry.get("fastq_aspera", "")
            if not acc or not aspera:
                logger.warning("Skipping %s: no aspera URL available", acc)
                continue
            results.append((acc, aspera, entry.get("fastq_bytes", "")))
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
def run_ascp(
    url: str,
    outdir: str,
    ascp_bin: str,
    ascp_key: str,
    host: str,
    bandwidth: str,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    expected_size: Optional[int] = None,
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
        Remote host prefix, e.g. ``era-fasp@fasp.sra.ebi.ac.uk:``. Only
        applied to URLs that do not already carry a ``user@host`` prefix.
    bandwidth : str
        Bandwidth limit, e.g. ``500M``.
    timeout : int
        Per-attempt timeout in seconds; 0 or negative disables the cap.
    retries : int
        Maximum number of retry attempts.
    expected_size : int, optional
        Expected file size in bytes (from the ENA API). When given, an
        existing local file is only treated as complete if its size matches;
        truncated left-overs are re-downloaded.

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
        if expected_size is None or os.path.getsize(dest) == expected_size:
            logger.info("File already exists, skipping: %s", dest)
            return dest
        logger.warning(
            "Existing file %s is %d bytes (expected %d), re-downloading",
            dest, os.path.getsize(dest), expected_size,
        )

    # Build the remote source. ENA returns URLs without a user prefix
    # ("fasp.sra.ebi.ac.uk:/vol1/..."), which need the anonymous user from
    # --host; URLs that already contain user@host are used verbatim.
    if "@" in url.split(":", 1)[0]:
        remote_src = url
    else:
        remote_path = url.split(":", 1)[-1] if ":" in url else url
        remote_src = f"{host}{remote_path}"

    cmd = [
        ascp_bin,
        "-T",               # disable encryption for speed
        "-l", bandwidth,
        "-P", "33001",      # TCP/SSH control port (UDP port is -O)
        "-i", ascp_key,
        "-Q",               # adaptive rate
        remote_src,
        outdir,
    ]

    logger.debug("ascp command: %s", " ".join(cmd))
    timeout_arg = timeout if timeout and timeout > 0 else None

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
                timeout=timeout_arg,
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
            last_exc = RuntimeError(f"Timeout after {timeout_arg}s for {filename}")
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
    parser.add_argument("-t", "--threads", type=int,
                        default=int(os.environ.get("BUC_THREADS", "") or DEFAULT_THREADS),
                        help=f"Number of parallel downloads (default: {DEFAULT_THREADS}, "
                             "or $BUC_THREADS)")
    parser.add_argument("--retries", type=int, default=DEFAULT_RETRIES,
                        help=f"Retry attempts per file (default: {DEFAULT_RETRIES})")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                        help=f"Per-download timeout in seconds (default: {DEFAULT_TIMEOUT}); "
                             "0 disables the cap")
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

    Duplicates are removed (order preserved). Entries that do not look like
    run accessions (SRR/ERR/DRR + digits) are kept but warned about.

    Parameters
    ----------
    input_arg : str
        File path or a bare accession identifier (comma-separated lists are
        also accepted).

    Returns
    -------
    list[str]
        Cleaned, de-duplicated list of accession identifiers.
    """
    path = Path(input_arg)
    if path.is_file():
        accs = [
            line.strip()
            for line in path.read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
        logger.info("Read %d accession(s) from %s", len(accs), input_arg)
    else:
        # Assume it is a single accession or comma-separated list
        accs = [a.strip() for a in input_arg.split(",") if a.strip()]

    deduped: List[str] = []
    seen = set()
    for acc in accs:
        if acc not in seen:
            seen.add(acc)
            deduped.append(acc)
    if len(deduped) < len(accs):
        logger.info("Removed %d duplicate accession(s)", len(accs) - len(deduped))
    for acc in deduped:
        if not _RUN_ACC_RE.match(acc):
            logger.warning(
                "'%s' does not look like a run accession (SRR/ERR/DRR + digits); "
                "the ENA API may return no data for it", acc,
            )
    return deduped


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

    # Build download tasks: (url, outdir, expected_size or None)
    tasks: List[Tuple[str, str, Optional[int]]] = []
    for acc, aspera, sizes in records:
        urls = parse_aspera_urls(aspera)
        byte_fields = parse_aspera_urls(sizes) if sizes else []
        for idx, u in enumerate(urls):
            size: Optional[int] = None
            if idx < len(byte_fields) and byte_fields[idx].isdigit():
                size = int(byte_fields[idx])
            tasks.append((u, outdir, size))
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
                expected_size=size,
            ): u
            for u, d, size in tasks
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
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
