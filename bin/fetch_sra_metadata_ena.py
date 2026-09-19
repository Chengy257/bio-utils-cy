#!/usr/bin/env python3
"""
File Name: fetch_sra_metadata_ena.py
Author: ChengYu
Description: Fetch SRA metadata from ENA filereport API.
             Accepts a list of SRA accessions (SRR/ERR/DRR) and retrieves
             specified fields from the ENA API, outputting a merged TSV.
Created Time: 2026
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from pathlib import Path
from typing import List, Optional, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

__version__ = "1.0.0"

ENA_FILEREPORT_BASE_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

DEFAULT_FIELDS = (
    "run_accession,experiment_accession,sample_accession,"
    "study_accession,experiment_title,instrument_model,"
    "library_layout,library_strategy,library_source,"
    "read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes"
)

DEFAULT_TIMEOUT = 60
DEFAULT_DELAY = 0.5
MAX_ACCESSIONS_PER_REQUEST = 500

logger = logging.getLogger(__name__)


def build_url(accessions: List[str], fields: str) -> str:
    """Build the ENA filereport API URL for the given accessions and fields.

    Args:
        accessions: List of SRA run accessions (e.g. SRR123456).
        fields: Comma-separated field names to request.

    Returns:
        Fully constructed ENA filereport URL string.
    """
    params = {
        "accession": ",".join(accessions),
        "result": "read_run",
        "fields": fields,
        "format": "tsv",
        "download": "txt",
    }
    return f"{ENA_FILEREPORT_BASE_URL}?{urlencode(params)}"


def fetch_tsv(
    url: str,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = 3,
    backoff: float = 2.0,
) -> str:
    """Fetch TSV content from the given URL with retry logic.

    Args:
        url: The URL to fetch.
        timeout: Request timeout in seconds.
        retries: Maximum number of retry attempts.
        backoff: Exponential backoff multiplier between retries.

    Returns:
        Decoded text content of the response.

    Raises:
        RuntimeError: If all retry attempts are exhausted.
    """
    last_exception: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            logger.debug("Fetching URL (attempt %d/%d): %s", attempt, retries, url)
            req = Request(url)
            with urlopen(req, timeout=timeout) as resp:
                charset = resp.headers.get_content_charset() or "utf-8"
                return resp.read().decode(charset)
        except HTTPError as exc:
            last_exception = exc
            if exc.code == 404:
                logger.warning("HTTP 404 for URL: %s", url)
                return ""
            logger.warning(
                "HTTPError %s on attempt %d/%d: %s", exc.code, attempt, retries, exc
            )
        except URLError as exc:
            last_exception = exc
            logger.warning(
                "URLError on attempt %d/%d: %s", attempt, retries, exc
            )
        except Exception as exc:
            last_exception = exc
            logger.warning(
                "Unexpected error on attempt %d/%d: %s", attempt, retries, exc
            )

        if attempt < retries:
            sleep_time = backoff ** (attempt - 1)
            logger.info("Retrying in %.1f seconds ...", sleep_time)
            time.sleep(sleep_time)

    raise RuntimeError(
        f"Failed to fetch data after {retries} attempts. Last error: {last_exception}"
    )


def read_accessions(path: str) -> List[str]:
    """Read accessions from a file, one per line.

    Blank lines and lines starting with '#' are skipped.

    Args:
        path: Path to the input file.

    Returns:
        List of stripped accession strings.
    """
    accessions: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            acc = line.strip()
            if acc and not acc.startswith("#"):
                accessions.append(acc)
    logger.info("Read %d accessions from %s", len(accessions), path)
    return accessions


def chunk_list(items: List[str], size: int) -> Sequence[List[str]]:
    """Split a list into chunks of the given size.

    Args:
        items: The list to chunk.
        size: Maximum items per chunk.

    Returns:
        List of sub-lists.
    """
    return [items[i : i + size] for i in range(0, len(items), size)]


def parse_tsv_text(text: str) -> List[dict]:
    """Parse TSV text into a list of dictionaries.

    Args:
        text: Raw TSV text with a header row.

    Returns:
        List of dicts keyed by the header column names.
    """
    rows: List[dict] = []
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    for row in reader:
        rows.append(dict(row))
    return rows


def main() -> None:
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Fetch SRA metadata from the ENA filereport API.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  %(prog)s -i accessions.txt -o metadata.tsv
  %(prog)s -i accessions.txt -o out.tsv -f "run_accession,sample_accession,fastq_ftp"
  %(prog)s -i ids.txt --timeout 120 --retries 5
  %(prog)s -i ids.txt --log-level DEBUG
""",
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Path to file with one accession per line."
    )
    parser.add_argument(
        "-o", "--output", default="-", help="Output TSV file path (default: stdout)."
    )
    parser.add_argument(
        "-f",
        "--fields",
        default=DEFAULT_FIELDS,
        help=(
            "Comma-separated ENA fields to retrieve. "
            f"(default: {DEFAULT_FIELDS})"
        ),
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help=f"HTTP request timeout in seconds (default: {DEFAULT_TIMEOUT}).",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of retry attempts per request (default: 3).",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=DEFAULT_DELAY,
        help=f"Delay in seconds between batch requests (default: {DEFAULT_DELAY}).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=MAX_ACCESSIONS_PER_REQUEST,
        help=f"Maximum accessions per API request (default: {MAX_ACCESSIONS_PER_REQUEST}).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    accessions = read_accessions(args.input)
    if not accessions:
        logger.error("No accessions found in %s", args.input)
        sys.exit(1)

    all_rows: List[dict] = []
    field_names: List[str] = []
    chunks = chunk_list(accessions, args.batch_size)
    logger.info(
        "Fetching metadata in %d batch(es) of up to %d accessions each.",
        len(chunks),
        args.batch_size,
    )

    for idx, chunk in enumerate(chunks, start=1):
        logger.info("Processing batch %d/%d (%d accessions)", idx, len(chunks), len(chunk))
        url = build_url(chunk, args.fields)
        tsv_text = fetch_tsv(url, timeout=args.timeout, retries=args.retries)
        if not tsv_text.strip():
            logger.warning("Batch %d returned empty response; skipping.", idx)
            continue
        rows = parse_tsv_text(tsv_text)
        if rows and not field_names:
            field_names = list(rows[0].keys())
        all_rows.extend(rows)
        if idx < len(chunks):
            time.sleep(args.delay)

    if not all_rows:
        logger.error("No metadata retrieved for any accession.")
        sys.exit(1)

    if not field_names:
        field_names = args.fields.split(",")

    if args.output == "-":
        writer = csv.DictWriter(sys.stdout, fieldnames=field_names, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)
    else:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=field_names, delimiter="\t")
            writer.writeheader()
            writer.writerows(all_rows)
        logger.info("Wrote %d rows to %s", len(all_rows), args.output)


if __name__ == "__main__":
    main()
