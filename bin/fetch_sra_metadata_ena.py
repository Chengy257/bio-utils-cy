#!/usr/bin/env python3
"""
File Name: fetch_sra_metadata_ena.py
Author: ChengYu
Description: Fetch SRA metadata from ENA filereport API.
             Accepts a list of SRA accessions (SRR/ERR/DRR) and retrieves
             specified fields from the ENA API, outputting a merged TSV.
Created Time: 2026

Changelog:
  v1.1.0  2026-09-27
  - FIX: HTTP 404 silently dropped the whole batch (up to 500 accessions)
    with a single warning; 404 batches are now bisected to isolate the
    unknown accession(s), and every failed accession is recorded to
    <output>.failed.txt
  - FIX: ENA error bodies (HTTP 200, non-TSV text) were parsed into
    garbage rows; they are now detected and the batch recorded as failed
  - FIX: the TSV header was taken from the first batch only -- later
    batches with a different column set crashed DictWriter; fieldnames
    are now the union across all batches
  - FIX: input accessions are deduplicated; non-run accessions warned
  - FIX: any failed accession makes the script exit non-zero (the merged
    TSV would otherwise be silently incomplete)
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

__version__ = "1.1.0"

ENA_FILEREPORT_BASE_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

DEFAULT_FIELDS = (
    "run_accession,experiment_accession,sample_accession,"
    "study_accession,experiment_title,instrument_model,"
    "library_layout,library_strategy,library_source,"
    "read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes"
)

DEFAULT_TIMEOUT = 60
DEFAULT_DELAY = 0.5
DEFAULT_BACKOFF = 2.0
MAX_ACCESSIONS_PER_REQUEST = 500

_RUN_ACC_RE = re.compile(r"^[SED]RR\d+$")

logger = logging.getLogger(__name__)


class EnaNotFoundError(RuntimeError):
    """HTTP 404 from the ENA API: at least one accession in the batch is unknown."""


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
            if exc.code == 404:
                # Deterministic: at least one accession in the batch does
                # not exist. The caller bisects to isolate it -- retrying
                # or returning "" would silently drop the whole batch.
                raise EnaNotFoundError(
                    f"HTTP 404 from ENA (unknown accession in batch): {url}"
                ) from exc
            last_exception = exc
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

    Blank lines and lines starting with '#' are skipped; duplicates are
    removed (order preserved); tokens that do not look like run accessions
    are kept but warned about.

    Args:
        path: Path to the input file.

    Returns:
        List of stripped, de-duplicated accession strings.
    """
    accessions: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            acc = line.strip()
            if acc and not acc.startswith("#") and acc not in accessions:
                accessions.append(acc)
    for acc in accessions:
        if not _RUN_ACC_RE.match(acc):
            logger.warning(
                "'%s' does not look like a run accession (SRR/ERR/DRR + digits)",
                acc,
            )
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
        List of dicts keyed by the header column names (a possible
        DictReader restkey from ragged rows is dropped).
    """
    rows: List[dict] = []
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    for row in reader:
        rows.append({k: v for k, v in row.items() if k is not None})
    return rows


def fetch_chunk_rows(
    chunk: List[str],
    fields: str,
    timeout: int,
    retries: int,
    failed: List[str],
    backoff: float = DEFAULT_BACKOFF,
) -> List[dict]:
    """Fetch one batch of accessions; on 404, bisect to isolate unknown ones.

    Args:
        chunk: Accessions in this batch.
        fields: Comma-separated ENA field names.
        timeout: HTTP timeout per request.
        retries: Retry attempts per request.
        failed: List to append accessions that could not be retrieved.
        backoff: Retry backoff multiplier.

    Returns:
        Rows for the accessions that were found.
    """
    url = build_url(chunk, fields)
    try:
        text = fetch_tsv(url, timeout=timeout, retries=retries, backoff=backoff)
    except EnaNotFoundError:
        if len(chunk) == 1:
            logger.warning("Accession %s not found in ENA; recorded as failed.", chunk[0])
            failed.append(chunk[0])
            return []
        mid = len(chunk) // 2
        logger.info("404 for batch of %d; splitting to isolate unknown accession(s).", len(chunk))
        return (
            fetch_chunk_rows(chunk[:mid], fields, timeout, retries, failed, backoff)
            + fetch_chunk_rows(chunk[mid:], fields, timeout, retries, failed, backoff)
        )

    if not text.strip():
        logger.warning("Batch of %d returned an empty response; recorded as failed.", len(chunk))
        failed.extend(chunk)
        return []

    # ENA reports problems (e.g. invalid field names) as HTTP 200 with a
    # plain-text error body; a real TSV has a header + tab separators.
    if "\t" not in text and len(fields.split(",")) > 1:
        logger.error(
            "ENA returned a non-TSV error body for a batch of %d: %s",
            len(chunk), text[:200],
        )
        failed.extend(chunk)
        return []

    return parse_tsv_text(text)


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
    failed: List[str] = []
    chunks = chunk_list(accessions, args.batch_size)
    logger.info(
        "Fetching metadata in %d batch(es) of up to %d accessions each.",
        len(chunks),
        args.batch_size,
    )

    for idx, chunk in enumerate(chunks, start=1):
        logger.info("Processing batch %d/%d (%d accessions)", idx, len(chunks), len(chunk))
        all_rows.extend(
            fetch_chunk_rows(chunk, args.fields, args.timeout, args.retries, failed)
        )
        if idx < len(chunks):
            time.sleep(args.delay)

    if not all_rows:
        logger.error("No metadata retrieved for any accession.")
        if failed:
            logger.error("Accessions not found or failed: %s", ", ".join(failed))
        sys.exit(1)

    # Fieldnames: union across all batches (later batches may differ from
    # the first one); first-seen order preserved.
    seen_fields: Dict[str, None] = {}
    for row in all_rows:
        for key in row:
            seen_fields.setdefault(key, None)
    field_names = list(seen_fields) or args.fields.split(",")

    if failed:
        if args.output == "-":
            logger.warning(
                "Accessions not found or failed (%d): %s",
                len(failed), ", ".join(failed),
            )
        else:
            failed_path = Path(args.output + ".failed.txt")
            failed_path.write_text("\n".join(failed) + "\n", encoding="utf-8")
            logger.warning(
                "Accessions not found or failed (%d); list written to %s",
                len(failed), failed_path,
            )

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

    if failed:
        # The merged TSV is incomplete; make that visible to pipelines.
        sys.exit(1)


if __name__ == "__main__":
    main()
