#!/usr/bin/env python3
#########################################################################
# File Name: blast_align_analysis.py
# Author: ChengYu
# Description: Analyze DNA FASTA files using BLAST for pairwise
#              similarity metrics between species (nucleotide + protein).
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-26
#   - FIX: protein queries passed to blastp kept their terminal stop
#     codon as '*'; translation now stops at the first stop codon
#     (to_stop=True).
#   - NEW: --task option (blastn-short default, preserving previous
#     behaviour) so long sequences can use the blastn/megablast tasks;
#     blastn-short is tuned for very short queries and changes e-value
#     semantics for longer ones (documented in --help).
#   - DOC: only the first HSP of the top alignment is reported, and
#     nucleotide identity uses the BLAST alignment length (including
#     gap columns) as denominator — the same convention as outfmt6
#     pident. Both are now stated in the docstring and --help.
#   - --threads defaults to BUC_THREADS (config/env.sh) when set and is
#     validated.
#########################################################################
"""Analyze DNA FASTA files using BLAST for inter-species similarity.

For each multi-species FASTA file, uses BLAST to compute nucleotide and
protein similarity between a reference species and all others. Outputs
a CSV with e-values, bit scores, and identity percentages.

Per pair, only the first HSP of the top BLAST alignment is reported.
Nucleotide identity uses the alignment length (gap columns included) as
denominator, matching BLAST outfmt6 pident. The default nucleotide task
is blastn-short (tuned for very short queries); for CDS-sized or longer
sequences pass --task blastn (or megablast) — e-value scales differ
between tasks.

Requires: blastn, blastp (NCBI BLAST+)
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq

__version__ = "1.1.0"

VALID_BLASTN_TASKS = ("blastn", "blastn-short", "megablast", "dc-megablast")

# blastn-short is kept as the default for backwards compatibility; it is
# only appropriate for very short queries (<~30 bp).

BLASTP_PARAMS = ["-matrix", "PAM30", "-seg", "no", "-gapopen", "9", "-gapextend", "1"]


def resolve_tool(name: str) -> str:
    """Resolve an external tool: BUC_<NAME>_BIN (config/env.sh) first, then PATH."""
    env_path = os.environ.get(f"BUC_{name.upper()}_BIN", "")
    if env_path:
        if os.path.isfile(env_path) and os.access(env_path, os.X_OK):
            return env_path
        raise FileNotFoundError(f"BUC_{name.upper()}_BIN is set but not executable: {env_path}")
    found = shutil.which(name)
    if not found:
        raise FileNotFoundError(
            f"{name} not found in PATH; set BUC_{name.upper()}_BIN in config/env.local.sh"
        )
    return found


def validate_fasta(file_path: str) -> bool:
    """Validate that a file is in FASTA format with non-empty sequences.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        True if valid, False otherwise.
    """
    try:
        records = list(SeqIO.parse(file_path, "fasta"))
        if not records:
            logging.error("No sequences found in %s", file_path)
            return False
        for record in records:
            if not record.seq or not str(record.seq).replace("-", "").isalpha():
                logging.error("Invalid sequence in %s: %s", file_path, record.id)
                return False
        return True
    except Exception as e:
        logging.error("Error parsing FASTA file %s: %s", file_path, e)
        return False


def clean_sequence(seq: Seq) -> Seq:
    """Remove gaps and trim to codon boundary.

    Args:
        seq: Input sequence (may contain gaps).

    Returns:
        Cleaned sequence trimmed to a multiple of 3, or empty Seq.
    """
    cleaned = str(seq).replace("-", "").upper()
    if len(cleaned) < 3:
        return Seq("")
    remainder = len(cleaned) % 3
    if remainder != 0:
        cleaned = cleaned[:-remainder]
    return Seq(cleaned)


def run_blast(
    query_seq: Seq,
    target_seq: Seq,
    is_protein: bool = False,
    blast_task: str = "blastn-short",
) -> dict:
    """Run BLAST between two sequences and return similarity metrics.

    Args:
        query_seq: Query sequence.
        target_seq: Target/subject sequence.
        is_protein: Whether to use blastp instead of blastn.
        blast_task: blastn task (ignored for blastp).

    Returns:
        Dict with 'e-value', 'bit_score', 'identity' keys (None on failure).
    """
    if len(query_seq) < 3 or len(target_seq) < 3:
        logging.debug("Sequence too short for BLAST (query=%d, target=%d).", len(query_seq), len(target_seq))
        return {"e-value": None, "bit_score": None, "identity": None}

    query_file = NamedTemporaryFile(delete=False, suffix=".fa")
    target_file = NamedTemporaryFile(delete=False, suffix=".fa")
    result_file = NamedTemporaryFile(delete=False, suffix=".xml")

    try:
        query_file.write(f">query\n{query_seq}\n".encode())
        target_file.write(f">target\n{target_seq}\n".encode())
        query_file.close()
        target_file.close()

        program = resolve_tool("blastp" if is_protein else "blastn")
        if is_protein:
            extra_params = BLASTP_PARAMS
        else:
            # Dust filtering stays off for every task: this tool compares
            # orthologous (often GC-rich) sequences, and DUST can mask
            # those almost completely.
            extra_params = ["-task", blast_task, "-dust", "no"]
            if blast_task == "blastn-short":
                # short-query tuning (legacy default): lenient gaps
                extra_params += ["-gapopen", "4", "-gapextend", "2"]

        cmd = [
            program,
            "-query", query_file.name,
            "-subject", target_file.name,
            "-outfmt", "5",
            "-out", result_file.name,
        ] + extra_params

        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with open(result_file.name) as fh:
            blast_record = NCBIXML.read(fh)

        if blast_record.alignments:
            hsp = blast_record.alignments[0].hsps[0]
            return {
                "e-value": hsp.expect,
                "bit_score": hsp.bits,
                "identity": hsp.identities / hsp.align_length if hsp.align_length > 0 else None,
            }
        return {"e-value": None, "bit_score": None, "identity": None}

    except subprocess.CalledProcessError as e:
        logging.error("BLAST failed: %s", e.stderr.decode() if e.stderr else str(e))
        return {"e-value": None, "bit_score": None, "identity": None}
    finally:
        for f in (query_file.name, target_file.name, result_file.name):
            if os.path.exists(f):
                os.unlink(f)


def extract_species(record_id: str) -> str:
    """Extract species name from record ID (text before first dot)."""
    return record_id.split(".")[0]


def process_file(file_path: str, reference_species: str, blast_task: str = "blastn-short") -> list:
    """Process a single multi-species FASTA file.

    Args:
        file_path: Path to the DNA FASTA file.
        reference_species: Species name to use as reference.

    Returns:
        List of result dicts with similarity metrics.
    """
    if not validate_fasta(file_path):
        return []

    records = list(SeqIO.parse(file_path, "fasta"))

    ref_record = None
    for record in records:
        if extract_species(record.id) == reference_species:
            ref_record = record
            break

    if ref_record is None:
        logging.warning("Reference species '%s' not found in %s.", reference_species, file_path)
        return []

    ref_seq = clean_sequence(ref_record.seq)
    if len(ref_seq) < 3:
        logging.warning("Reference sequence too short after cleaning in %s.", file_path)
        return []

    results = []
    for record in records:
        if record.id == ref_record.id:
            continue

        target_seq = clean_sequence(record.seq)
        if len(target_seq) < 3:
            logging.debug("Target sequence too short: %s in %s", record.id, file_path)
            continue

        nuc_blast = run_blast(ref_seq, target_seq, is_protein=False, blast_task=blast_task)
        prot_blast = run_blast(ref_seq.translate(to_stop=True), target_seq.translate(to_stop=True), is_protein=True)

        results.append({
            "file": os.path.basename(file_path),
            "reference": ref_record.id,
            "target": record.id,
            "nucleotide_evalue": nuc_blast["e-value"],
            "nucleotide_bitscore": nuc_blast["bit_score"],
            "nucleotide_identity": nuc_blast["identity"],
            "protein_identity": prot_blast["identity"],
        })

    return results


def analyze_directory(
    input_dir: str,
    reference_species: str,
    output_file: str,
    threads: int = 4,
    suffix: str = ".dna.fa",
    blast_task: str = "blastn-short",
) -> None:
    """Analyze all matching FASTA files in a directory.

    Args:
        input_dir: Directory containing FASTA files.
        reference_species: Reference species name.
        output_file: Output CSV path.
        threads: Number of parallel workers.
        suffix: File suffix to match (default: .dna.fa).
        blast_task: blastn task for nucleotide comparisons.
    """
    files = sorted(
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith(suffix)
    )

    if not files:
        logging.error("No files matching '*%s' found in %s", suffix, input_dir)
        sys.exit(1)

    logging.info("Processing %d files with %d threads.", len(files), threads)

    results = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_map = {
            executor.submit(process_file, fp, reference_species, blast_task): fp
            for fp in files
        }
        for future in as_completed(future_map):
            fp = future_map[future]
            try:
                file_results = future.result()
                results.extend(file_results)
                logging.info("Completed: %s (%d pairs)", os.path.basename(fp), len(file_results))
            except Exception as e:
                logging.error("Failed to process %s: %s", fp, e)

    if not results:
        logging.warning("No similarity results produced.")
        return

    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    logging.info("Results (%d rows) written to %s", len(df), output_file)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Analyze DNA FASTA files using BLAST for inter-species similarity metrics.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python blast_align_analysis.py -i alignments/ -r human -o results.csv
  python blast_align_analysis.py -i alignments/ -r human -o results.csv -t 8 --suffix .fa
""",
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Input directory containing multi-species DNA FASTA files.",
    )
    parser.add_argument(
        "-r", "--reference", type=str, required=True,
        help="Reference species name (matched against record ID before first dot).",
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True,
        help="Output CSV file path.",
    )
    parser.add_argument(
        "-t", "--threads", type=int,
        default=int(os.environ.get("BUC_THREADS", "4") or 4),
        help="Number of parallel threads (default: BUC_THREADS or 4).",
    )
    parser.add_argument(
        "--task", type=str, default="blastn-short", choices=VALID_BLASTN_TASKS,
        help="blastn task for nucleotide comparisons (default: blastn-short, "
             "the previous behaviour; tuned for very short queries). Use "
             "'blastn' or 'megablast' for CDS-sized or longer sequences — "
             "e-value scales differ between tasks.",
    )
    parser.add_argument(
        "--suffix", type=str, default=".dna.fa",
        help="File suffix to match in input directory (default: .dna.fa).",
    )
    parser.add_argument(
        "--log-level", type=str, default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO).",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    return parser


def main() -> None:
    """Entry point."""
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=getattr(logging, args.log_level),
    )

    input_dir = Path(args.input)
    if not input_dir.is_dir():
        logging.error("Input directory not found: %s", args.input)
        sys.exit(1)

    if args.threads < 1:
        parser.error("--threads must be >= 1")

    # Check BLAST+ availability (and resolve configured binaries)
    try:
        resolve_tool("blastn")
        resolve_tool("blastp")
        subprocess.run([resolve_tool("blastn"), "-version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        logging.error("BLAST+ not available: %s", exc)
        sys.exit(1)

    analyze_directory(
        input_dir=str(input_dir),
        reference_species=args.reference,
        output_file=args.output,
        threads=args.threads,
        suffix=args.suffix,
        blast_task=args.task,
    )


if __name__ == "__main__":
    main()
