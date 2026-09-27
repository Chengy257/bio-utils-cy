#!/usr/bin/env python3
"""
Visualize MS2 spectra with theoretical b/y ion annotations from mzID+mzML/MGF.

Author: ChengYu
Created Time: 2026

Changelog:
  v1.1.0  2026-09-27
  - FIX: top-level lxml/pyteomics imports crashed even `--help` (both
    packages are absent from the runtime environment); the mzIdentML
    parser now uses the stdlib ElementTree (namespace-agnostic local-name
    matching, works for mzid 1.1/1.2), fragment masses use a built-in
    monoisotopic residue table, and MGF is parsed by a built-in reader.
    Runtime deps are now numpy + matplotlib only; mzML stays optional
    behind pymzml with a clean error message.
  - FIX: the b/y ion formula was systematically wrong by ~19.02 Da
    (pyteomics' std_ion_comp offsets are relative to the neutral peptide
    including water, but were applied to the bare residue sum -- and the
    float + Composition arithmetic would have raised TypeError anyway).
    Standard singly-charged values are now used: b_i = sum(res) + proton,
    y_i = sum(res) + H2O + proton.
  - FIX: modification placement off-by-ones: N-terminal (location=0) mods
    were never included in b ions, and a mod on residue i was excluded
    from b_i / wrongly included from y_i (position >= len-i). Now:
    b includes mods with location <= i, y with location >= len-i+1
    (C-terminal location = len+1 lands in every y ion).
  - FIX: protein matching looked for ProteinAccession/DBSequence elements
    inside each SpectrumIdentificationItem, but those live in the
    SequenceCollection section -- every PSM was filtered out and the tool
    always reported "No PSMs found". Proteins are now resolved through
    PeptideEvidenceRef -> PeptideEvidence -> DBSequence accessions.
  - FIX: modifications were read from under each SpectrumIdentificationItem,
    but in the mzIdentML schema they live inside the Peptide element --
    every PSM silently had no modifications. They are now parsed from the
    referenced Peptide (monoisotopicMassDelta with avgMassDelta fallback).
  - FIX: mzML scan lookup used MS:1000016, which is retention time, not
    the scan number; spectra are now matched by the scan= field of the
    spectrum ID.
  - FIX: MGF "SCANS=lo-hi" ranges never matched, and the whole MGF was
    re-read for every PSM (O(n*m)); one pass now resolves all scans.
  - FIX: output directories are created on demand, and exiting 0 with
    zero plotted spectra became an error (exit 1).
  - CLEAN: dead MS-GF score field removed.
"""

__version__ = "1.1.0"

import argparse
import logging
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import pymzml
except ImportError:
    pymzml = None

logger = logging.getLogger(__name__)

ION_OFFSET = 0.1
COLORS = {"b": "#FF3030", "y": "#4876FF"}

# Monoisotopic residue masses (amino acid minus H2O), Da. Values match
# pyteomics.mass.std_aa_mass / Unimod.
STD_AA_MASS = {
    "G": 57.021464, "A": 71.037114, "S": 87.032029, "P": 97.052764,
    "V": 99.068414, "T": 101.047679, "C": 103.009185, "L": 113.084064,
    "I": 113.084064, "N": 114.042927, "D": 115.026943, "Q": 128.058578,
    "K": 128.094963, "E": 129.042593, "M": 131.040485, "H": 137.058912,
    "F": 147.068414, "R": 156.101111, "Y": 163.063329, "W": 186.079313,
}
PROTON_MASS = 1.00727646688     # monoisotopic proton, Da
H2O_MASS = 18.010564684         # monoisotopic water, Da


# ---------------------------------------------------------------------------
# mzIdentML parsing (stdlib ElementTree, namespace-agnostic)
# ---------------------------------------------------------------------------

def _local(tag: str) -> str:
    """Strip the '{namespace}' prefix from an ElementTree tag."""
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def _iter_elems(root: ET.Element, name: str) -> Iterator[ET.Element]:
    """Yield all descendant elements with the given local tag name."""
    for elem in root.iter():
        if _local(elem.tag) == name:
            yield elem


def _first_child(elem: ET.Element, name: str) -> Optional[ET.Element]:
    """Return the first direct child with the given local tag name."""
    for child in elem:
        if _local(child.tag) == name:
            return child
    return None


def _extract_scan(spectrum_id: str) -> Optional[int]:
    """Extract the scan number from a spectrumID like 'scan=123' or '123'."""
    match = re.search(r"scan=(\d+)", spectrum_id)
    if match:
        return int(match.group(1))
    try:
        return int(spectrum_id)
    except ValueError:
        return None


def _primary_accession(target: str) -> str:
    """Extract the accession to search for from a user-supplied protein ID.

    UniProt-style IDs ("sp|P12345|PROT_HUMAN") reduce to the second
    pipe-separated field; anything else reduces to its first
    '. '/'_'/'|'-separated token, upper-cased.
    """
    if target.count("|") >= 2:
        parts = target.split("|")
        return parts[1].strip().upper()
    parts = [p for p in re.split(r"[.\s_|]", target) if p]
    return parts[0].upper() if parts else ""


def _protein_matches(accession: str, target: str) -> bool:
    """Fuzzy protein match: any token of *accession* equals *target*.

    Handles full UniProt-style accessions ("sp|P12345|PROT_HUMAN") as well
    as bare accessions; comparing only the first token never matches
    "sp|..." strings.
    """
    return any(
        token == target
        for token in re.split(r"[.\s_|]", accession)
        if token
    )


def parse_mzid(
    mzid_path: str,
    target_protein_id: str,
) -> List[Dict[str, Any]]:
    """Parse an mzIdentML file and return PSMs matching *target_protein_id*.

    Proteins are resolved through PeptideEvidenceRef -> PeptideEvidence ->
    DBSequence accessions (the ProteinAccession/DBSequence elements are NOT
    part of the SpectrumIdentificationItem subtree). Decoy PeptideEvidence
    entries (isDecoy="true") are skipped, so target-decoy FDR-style hits
    never match the target protein.

    Parameters
    ----------
    mzid_path : str
        Path to the .mzid file.
    target_protein_id : str
        Protein accession to search for.

    Returns
    -------
    list[dict]
        Each dict contains keys: scan, spectrum_id, sequence, charge, mods.
    """
    tree = ET.parse(mzid_path)
    root = tree.getroot()

    clean_target = _primary_accession(target_protein_id)
    logger.debug("Cleaned target protein ID: %s", clean_target)

    # DBSequence id -> accession
    dbseq_acc: Dict[str, str] = {}
    for db in _iter_elems(root, "DBSequence"):
        dbseq_acc[db.get("id", "")] = db.get("accession", "")

    # PeptideEvidence id -> list of protein accessions (decoys skipped)
    pe_proteins: Dict[str, List[str]] = {}
    for pe in _iter_elems(root, "PeptideEvidence"):
        if (pe.get("isDecoy") or "").lower() == "true":
            continue
        acc = dbseq_acc.get(pe.get("dBSequence_ref", ""), "")
        if acc:
            pe_id = pe.get("id", "")
            proteins = pe_proteins.setdefault(pe_id, [])
            if acc not in proteins:
                proteins.append(acc)

    # Peptide id -> (sequence, mods). Modifications live under the Peptide
    # element in the mzIdentML schema, not under SpectrumIdentificationItem.
    peptide_map: Dict[str, Dict[str, Any]] = {}
    for pep in _iter_elems(root, "Peptide"):
        seq_elem = _first_child(pep, "PeptideSequence")
        sequence = (seq_elem.text or "").strip() if seq_elem is not None else ""
        mods: List[Tuple[int, float]] = []
        for mod in _iter_elems(pep, "Modification"):
            try:
                location = int(mod.get("location") or 0)
                mono = (mod.get("monoisotopicMassDelta")
                        or mod.get("avgMassDelta") or "0")
                mods.append((location, float(mono)))
            except ValueError as exc:
                logger.warning("Skipping unparseable modification: %s", exc)
        peptide_map[pep.get("id", "")] = {"sequence": sequence, "mods": mods}

    psm_list: List[Dict[str, Any]] = []

    for sir in _iter_elems(root, "SpectrumIdentificationResult"):
        spectrum_id = sir.get("spectrumID", "")
        scan_num = _extract_scan(spectrum_id)

        for sii in _iter_elems(sir, "SpectrumIdentificationItem"):
            protein_refs: List[str] = []
            for per in _iter_elems(sii, "PeptideEvidenceRef"):
                protein_refs.extend(
                    pe_proteins.get(per.get("peptideEvidence_ref", ""), [])
                )
            # Fallback for tools that nest DBSequence elements in the item
            for db in _iter_elems(sii, "DBSequence"):
                acc = db.get("accession", "")
                if acc:
                    protein_refs.append(acc)

            matched = any(
                _protein_matches(ref, clean_target)
                for ref in protein_refs if ref
            )
            if not matched:
                continue

            pep_ref = sii.get("peptide_ref", "")
            pep_entry = peptide_map.get(pep_ref, {})
            sequence = pep_entry.get("sequence", "")
            if not sequence:
                seq_elem = _first_child(sii, "PeptideSequence")
                sequence = (seq_elem.text or "").strip() if seq_elem is not None else ""
            sequence = re.sub(r"\[.*?\]", "", sequence).upper()

            try:
                charge = int(sii.get("chargeState") or 1)
            except ValueError:
                charge = 1

            psm_list.append({
                "scan": scan_num,
                "spectrum_id": spectrum_id,
                "sequence": sequence,
                "charge": charge,
                "mods": list(pep_entry.get("mods", [])),
            })

    logger.debug("Found %d PSMs", len(psm_list))
    return psm_list


# ---------------------------------------------------------------------------
# Fragment ion calculation
# ---------------------------------------------------------------------------

def calculate_fragments(
    sequence: str,
    charge: int,
    modifications: Optional[List[Tuple[int, float]]] = None,
) -> Dict[str, List[float]]:
    """Calculate theoretical singly-charged b/y ion m/z values.

    b_i = sum(res[0:i]) + proton, y_i = sum(res[n-i:]) + H2O + proton.
    Modifications are (location, delta) with 1-based residue locations;
    location 0 denotes the N-terminus (included in every b ion) and
    location len+1 the C-terminus (included in every y ion).

    Parameters
    ----------
    sequence : str
        Peptide amino-acid sequence (one-letter codes).
    charge : int
        Precursor charge state (fragment ions are reported at z=1).
    modifications : list[tuple[int, float]], optional
        Each tuple is (position, mass_delta).

    Returns
    -------
    dict[str, list[float]]
        Keys ``'b'`` and ``'y'`` with m/z lists (ascending fragment index).

    Raises
    ------
    ValueError
        If the sequence contains an unknown residue letter.
    """
    if modifications is None:
        modifications = []

    unknown = sorted(set(sequence) - set(STD_AA_MASS))
    if unknown:
        raise ValueError(
            f"Unknown residue(s) {''.join(unknown)} in sequence '{sequence}'"
        )

    n = len(sequence)
    frags: Dict[str, List[float]] = {"b": [], "y": []}

    for i in range(1, n):
        # b-ion: prefix of length i (residues 1..i)
        b_mass = sum(STD_AA_MASS[aa] for aa in sequence[:i]) + PROTON_MASS
        for pos, delta in modifications:
            if pos <= i:  # location 0 (N-term) included everywhere
                b_mass += delta
        frags["b"].append(b_mass)

        # y-ion: suffix of length i (residues n-i+1..n)
        y_mass = sum(STD_AA_MASS[aa] for aa in sequence[n - i:]) + H2O_MASS + PROTON_MASS
        for pos, delta in modifications:
            if pos >= n - i + 1:  # location n+1 (C-term) included everywhere
                y_mass += delta
        frags["y"].append(y_mass)

    return frags


# ---------------------------------------------------------------------------
# Spectrum loading
# ---------------------------------------------------------------------------

def _mgf_scan_number(params: Dict[str, str]) -> Optional[int]:
    """Derive a scan number from MGF params (scans/scan param, then title)."""
    raw = params.get("scans") or params.get("scan") or ""
    token = str(raw).split(",")[0].strip()
    match = re.match(r"(\d+)", token)
    if match:
        return int(match.group(1))
    title = params.get("title", "")
    match = re.search(r"scan[ =](\d+)", title, re.IGNORECASE)
    return int(match.group(1)) if match else None


def iter_mgf_spectra(file_path: str) -> Iterator[Dict[str, Any]]:
    """Yield MGF spectra as dicts with keys scan, mz, inten."""
    params: Dict[str, str] = {}
    mz_list: List[float] = []
    inten_list: List[float] = []
    with open(file_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            upper = line.upper()
            if upper == "BEGIN IONS":
                params, mz_list, inten_list = {}, [], []
            elif upper == "END IONS":
                if mz_list:
                    yield {
                        "scan": _mgf_scan_number(params),
                        "mz": np.array(mz_list),
                        "inten": np.array(inten_list),
                    }
                params, mz_list, inten_list = {}, [], []
            elif "=" in line and not mz_list:
                key, _, value = line.partition("=")
                params[key.strip().lower()] = value.strip()
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz_list.append(float(parts[0]))
                        inten_list.append(float(parts[1]))
                    except ValueError:
                        logger.warning("Skipping unparseable MGF peak line: %s", line)


def load_spectra_mgf(
    file_path: str, scan_nums: set
) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
    """Single pass over an MGF file; maps scan number -> (mz, intensity).

    Handles 'SCANS=n' and 'SCANS=lo-hi' range forms (attributed to lo).
    """
    found: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    for spec in iter_mgf_spectra(file_path):
        scan = spec["scan"]
        if scan is not None and scan in scan_nums and scan not in found:
            found[scan] = (spec["mz"], spec["inten"])
    return found


def load_spectra_mzml(
    file_path: str, scan_nums: set
) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
    """Single pass over an mzML file (requires pymzml); scan number -> arrays.

    Spectra are matched by the scan= field of the spectrum ID (MS:1000016
    is retention time and was wrongly used before v1.1.0).
    """
    if pymzml is None:
        logger.error("pymzml is required for mzML files. Install with: pip install pymzml")
        return {}
    found: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    for spec in pymzml.run.Reader(file_path):
        if spec.ms_level != 2:
            continue
        match = re.search(r"scan=(\d+)", str(spec.ID))
        if not match:
            continue
        scan = int(match.group(1))
        if scan in scan_nums and scan not in found:
            found[scan] = (np.asarray(spec.mz), np.asarray(spec.i))
    return found


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_ms2(
    mz: np.ndarray,
    intensity: np.ndarray,
    frags: Dict[str, List[float]],
    scan_num: int,
    sequence: str,
    protein_id: str,
    output_path: str,
    tolerance: float = 0.1,
) -> None:
    """Draw a single annotated MS2 spectrum and save to *output_path*."""
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.stem(mz, intensity, linefmt="gray", markerfmt=" ", basefmt=" ")

    for ion_type in ("b", "y"):
        for i, theo_mz in enumerate(frags[ion_type]):
            nearest_idx = int(np.argmin(np.abs(mz - theo_mz)))
            if abs(mz[nearest_idx] - theo_mz) < tolerance:
                label = f"{ion_type}{i + 1}"
                ax.vlines(
                    theo_mz, 0, intensity[nearest_idx],
                    colors=COLORS[ion_type], linestyles="dashed", alpha=0.7,
                )
                ax.text(
                    theo_mz + ION_OFFSET,
                    intensity[nearest_idx],
                    label,
                    color=COLORS[ion_type],
                    fontsize=8,
                )

    ax.set_title(f"Scan {scan_num} | {sequence}\nProtein: {protein_id}", fontsize=12)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")

    seq_text = " ".join(list(sequence))
    fig.text(
        0.5, 0.01, seq_text,
        ha="center", fontsize=10,
        bbox=dict(facecolor="lightgray", alpha=0.5),
    )

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved: %s", output_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plot_ms2_spectrum.py",
        description="Visualize MS2 spectra with theoretical b/y ion annotations "
                    "from mzIdentML + mzML/MGF files.",
        epilog=(
            "Examples:\n"
            "  %(prog)s --mzid result.mzid --spectra data.mgf  --protein P12345 -o spec\n"
            "  %(prog)s --mzid result.mzid --spectra data.mzML --protein P12345 -o spec\n"
            "\n"
            "mzML support requires the optional pymzml package.\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument("--mzid", required=True, help="Path to the mzIdentML (.mzid) file.")
    parser.add_argument("--spectra", required=True, help="Path to the MS/MS file (.mzML or .mgf).")
    parser.add_argument("--protein", required=True, help="Target protein accession ID.")
    parser.add_argument(
        "-o", "--output", default="ms2_spectrum",
        help="Output prefix (figures saved as <prefix>_<N>.png; parent "
             "directories are created on demand).",
    )
    parser.add_argument(
        "--tolerance", type=float, default=0.1,
        help="Mass tolerance (Da) for matching theoretical fragments.",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    # 1. Parse mzIdentML
    psms = parse_mzid(args.mzid, args.protein)
    if not psms:
        logger.error("No PSMs found for protein '%s' in %s", args.protein, args.mzid)
        sys.exit(1)

    # 2. Load spectra (single pass over the file for all PSMs)
    spectra_path = Path(args.spectra)
    ext = spectra_path.suffix.lower()
    scan_nums = {psm["scan"] for psm in psms if psm["scan"] is not None}
    for psm in psms:
        if psm["scan"] is None:
            logger.warning("PSM with spectrumID '%s' has no scan number; cannot look it up.",
                           psm["spectrum_id"])

    if ext == ".mgf":
        spectra = load_spectra_mgf(str(spectra_path), scan_nums)
    elif ext in (".mzml", ".mzxml"):
        spectra = load_spectra_mzml(str(spectra_path), scan_nums)
    else:
        logger.error("Unsupported spectrum format: %s", ext)
        sys.exit(1)

    # 3. Plot
    n_plotted = 0
    for idx, psm in enumerate(psms, 1):
        if psm["scan"] is None or psm["scan"] not in spectra:
            logger.warning("Scan %s not found in %s", psm["scan"], spectra_path)
            continue

        mz, intensity = spectra[psm["scan"]]
        try:
            frags = calculate_fragments(psm["sequence"], psm["charge"], psm["mods"])
        except ValueError as exc:
            logger.warning("Skipping scan %s: %s", psm["scan"], exc)
            continue

        out_path = f"{args.output}_{idx}.png"
        out = Path(out_path)
        if out.parent != Path("."):
            out.parent.mkdir(parents=True, exist_ok=True)
        plot_ms2(
            mz, intensity, frags,
            psm["scan"], psm["sequence"], args.protein,
            out_path, tolerance=args.tolerance,
        )
        n_plotted += 1

    if n_plotted == 0:
        logger.error("No spectra were plotted (0 of %d PSMs matched a scan).", len(psms))
        sys.exit(1)
    logger.info("Done — %d/%d spectra plotted.", n_plotted, len(psms))


if __name__ == "__main__":
    main()
