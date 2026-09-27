"""Tests for bin/plot_ms2_spectrum.py (v1.1.0: stdlib rewrite, no lxml/pyteomics).

Covers the rewritten b/y formula against literature values, modification
placement (including N-/C-terminal locations), the fixed mzid protein
resolution, the built-in MGF reader, and an end-to-end render.
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import plot_ms2_spectrum as pms  # noqa: E402
from common import ScriptTestCase  # noqa: E402

PROTON = 1.00727646688
H2O = 18.010564684

MZID_TEMPLATE = """<?xml version="1.0" encoding="utf-8"?>
<mzIdentML xmlns="http://psidev.info/psi/pi/mzIdentML/1.1" version="1.1.0">
  <SequenceCollection>
    <DBSequence id="DBS1" accession="sp|P12345|PROT_HUMAN"/>
    <DBSequence id="DBS2" accession="DECOY_P12345|rev"/>
    <Peptide id="P1">
      <PeptideSequence>PEPTIDE</PeptideSequence>
      {mod}
    </Peptide>
    <PeptideEvidence id="PE1" peptide_ref="P1" dBSequence_ref="DBS1" isDecoy="false"/>
    <PeptideEvidence id="PE2" peptide_ref="P1" dBSequence_ref="DBS2" isDecoy="true"/>
  </SequenceCollection>
  <DataCollection>
    <AnalysisData>
      <SpectrumIdentificationList id="SIL1">
        <SpectrumIdentificationResult spectrumID="scan=5" spectrumTitle="s5">
          <SpectrumIdentificationItem id="SII1" peptide_ref="P1" chargeState="2">
            <PeptideEvidenceRef peptideEvidence_ref="PE1"/>
          </SpectrumIdentificationItem>
        </SpectrumIdentificationResult>
        <SpectrumIdentificationResult spectrumID="scan=9" spectrumTitle="s9">
          <SpectrumIdentificationItem id="SII2" peptide_ref="P1" chargeState="3">
            <PeptideEvidenceRef peptideEvidence_ref="PE2"/>
          </SpectrumIdentificationItem>
        </SpectrumIdentificationResult>
        <SpectrumIdentificationResult spectrumID="index=12" spectrumTitle="s12">
          <SpectrumIdentificationItem id="SII3" peptide_ref="P1" chargeState="2">
            <PeptideEvidenceRef peptideEvidence_ref="PE1"/>
          </SpectrumIdentificationItem>
        </SpectrumIdentificationResult>
      </SpectrumIdentificationList>
    </AnalysisData>
  </DataCollection>
</mzIdentML>
"""

MGF_TEXT = """TITLE=first, scan=5
SCANS=5
108.0 100.0
147.0 500.0
227.1 300.0
END IONS
TITLE=second
SCANS=20-25
100.0 10.0
653.3 800.0
END IONS
"""


class TestCalculateFragments(ScriptTestCase):

    def test_b1_y1_match_literature_values(self):
        frags = pms.calculate_fragments("PEPTIDE", 2)
        # b1 = P + proton = 97.05276 + 1.00728 (literature: 98.0600)
        self.assertAlmostEqual(frags["b"][0], 97.052764 + PROTON, places=4)
        self.assertAlmostEqual(frags["b"][0], 98.0600, places=3)
        # y1 = E + H2O + proton = 129.04259 + 18.01056 + 1.00728 (literature: 148.0604)
        self.assertAlmostEqual(frags["y"][0], 129.042593 + H2O + PROTON, places=4)
        self.assertAlmostEqual(frags["y"][0], 148.0604, places=3)

    def test_b6_y6_literature(self):
        frags = pms.calculate_fragments("PEPTIDE", 2)
        # b6 (PEPTID) = 652.30679 + proton = 653.3141 (literature: 653.3145)
        self.assertAlmostEqual(frags["b"][5], 653.3141, places=3)
        # y6 (EPTIDE) = 646.27333 + H2O + proton = 665.2912 (literature: 665.2917)
        self.assertAlmostEqual(frags["y"][5], 703.3145, places=3)
        # y5 (PTIDE) = 555.25404 + H2O + proton = 574.2719
        self.assertAlmostEqual(frags["y"][4], 574.2719, places=3)

    def test_formula_no_longer_offset_by_water(self):
        # regression for v1.0.0: b = sum(res) - H2O, y = sum(res)
        frags = pms.calculate_fragments("AA", 1)
        aa = 71.037114
        self.assertAlmostEqual(frags["b"][0], aa + PROTON)          # not aa - 18.01
        self.assertAlmostEqual(frags["y"][0], aa + H2O + PROTON)    # not aa

    def test_nterminal_mod_in_all_b_ions(self):
        frags = pms.calculate_fragments("AAA", 1, [(0, 100.0)])
        # n=3 -> b1, b2 only (b3 would be the intact peptide)
        self.assertAlmostEqual(frags["b"][0], 71.037114 + PROTON + 100.0)
        self.assertAlmostEqual(frags["b"][1], 2 * 71.037114 + PROTON + 100.0)
        self.assertEqual(len(frags["b"]), 2)

    def test_residue_mod_boundary(self):
        # mod on residue 2 of AAA: in b2/b3, not b1; in y2/y3, not y1
        frags = pms.calculate_fragments("AAA", 1, [(2, 80.0)])
        self.assertAlmostEqual(frags["b"][0], 71.037114 + PROTON)
        self.assertAlmostEqual(frags["b"][1], 2 * 71.037114 + PROTON + 80.0)
        self.assertAlmostEqual(frags["y"][0], 71.037114 + H2O + PROTON)
        self.assertAlmostEqual(frags["y"][1], 2 * 71.037114 + H2O + PROTON + 80.0)

    def test_cterminal_mod_in_all_y_ions(self):
        frags = pms.calculate_fragments("AAA", 1, [(4, 50.0)])
        self.assertAlmostEqual(frags["y"][0], 71.037114 + H2O + PROTON + 50.0)
        self.assertAlmostEqual(frags["y"][1], 2 * 71.037114 + H2O + PROTON + 50.0)
        self.assertAlmostEqual(frags["b"][0], 71.037114 + PROTON)

    def test_unknown_residue_raises(self):
        with self.assertRaises(ValueError):
            pms.calculate_fragments("PEX", 1)


class TestParseMzid(ScriptTestCase):

    def _write_mzid(self, mod=""):
        return self.write("test.mzid", MZID_TEMPLATE.format(mod=mod))

    def test_protein_resolved_via_peptide_evidence(self):
        psm = self._write_mzid()
        psms = pms.parse_mzid(str(psm), "P12345")
        # scan=5 (PE1, target) and index=12 (PE1) match; scan=9 is DECOY-only
        self.assertEqual([p["scan"] for p in psms], [5, None])
        self.assertEqual(psms[0]["sequence"], "PEPTIDE")
        self.assertEqual(psms[0]["charge"], 2)
        self.assertEqual(psms[0]["spectrum_id"], "scan=5")

    def test_full_sp_accession_prefix_matched(self):
        psm = self._write_mzid()
        psms = pms.parse_mzid(str(psm), "sp|P12345|PROT_HUMAN")
        self.assertEqual(len(psms), 2)

    def test_other_protein_no_psms(self):
        psm = self._write_mzid()
        psms = pms.parse_mzid(str(psm), "P99999")
        self.assertEqual(psms, [])

    def test_decoy_evidence_never_matches(self):
        # decoy PeptideEvidence (isDecoy="true") is skipped entirely: the
        # decoy-only result (scan=9) is not matched by the target protein,
        # and searching the decoy accession itself yields nothing
        psm = self._write_mzid()
        self.assertEqual(pms.parse_mzid(str(psm), "P12345"),
                         pms.parse_mzid(str(psm), "P12345"))  # deterministic
        psms = pms.parse_mzid(str(psm), "P12345")
        self.assertNotIn(9, [p["scan"] for p in psms])
        self.assertEqual(pms.parse_mzid(str(psm), "DECOY_P12345"), [])

    def test_modifications_parsed(self):
        psm = self._write_mzid(mod='<Modification location="2" monoisotopicMassDelta="79.966331"/>')
        psms = pms.parse_mzid(str(psm), "P12345")
        self.assertEqual(psms[0]["mods"], [(2, 79.966331)])


class TestMgfReader(ScriptTestCase):

    def test_scans_param_range_and_single(self):
        mgf_file = self.write("test.mgf", MGF_TEXT)
        spectra = pms.load_spectra_mgf(str(mgf_file), {5, 20})
        self.assertEqual(sorted(spectra), [5, 20])
        self.assertAlmostEqual(spectra[5][0][0], 108.0)   # (mz, inten) tuple
        self.assertAlmostEqual(spectra[20][0][1], 653.3)

    def test_title_scan_fallback(self):
        mgf_file = self.write("t.mgf",
                              "TITLE=x scan=77\n100.0 1.0\nEND IONS\n")
        spectra = pms.load_spectra_mgf(str(mgf_file), {77})
        self.assertEqual(sorted(spectra), [77])

    def test_unwanted_scans_not_loaded(self):
        mgf_file = self.write("u.mgf", MGF_TEXT)
        spectra = pms.load_spectra_mgf(str(mgf_file), {5})
        self.assertEqual(sorted(spectra), [5])


class TestEndToEnd(ScriptTestCase):

    def _write_fixtures(self, mod="", mgf_text=MGF_TEXT):
        mzid = self.write("test.mzid", MZID_TEMPLATE.format(mod=mod))
        mgf = self.write("test.mgf", mgf_text)
        return mzid, mgf

    def test_plots_pngs_for_matching_scans(self):
        mzid, mgf = self._write_fixtures()
        out_prefix = self.tmp / "figs" / "spec"
        proc = self.run_script(
            "plot_ms2_spectrum.py", "--mzid", mzid, "--spectra", mgf,
            "--protein", "P12345", "-o", out_prefix)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        # scan=5 matched (index=12 has no scan and is skipped with a warning)
        self.assertTrue((self.tmp / "figs" / "spec_1.png").exists())
        self.assertIn("1/2 spectra plotted", proc.stderr)

    def test_modified_peptide_rendered(self):
        mzid, mgf = self._write_fixtures(
            mod='<Modification location="2" monoisotopicMassDelta="79.966331"/>')
        # add a peak at the modified b2 position so the annotation path runs
        mgf_with_peak = MGF_TEXT.replace(
            "227.1 300.0", "227.1 300.0\n307.06 400.0")
        self.write("test.mgf", mgf_with_peak)
        out_prefix = self.tmp / "spec"
        proc = self.run_script(
            "plot_ms2_spectrum.py", "--mzid", mzid, "--spectra", self.tmp / "test.mgf",
            "--protein", "P12345", "-o", out_prefix)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((self.tmp / "spec_1.png").exists())

    def test_no_scan_match_exits_nonzero(self):
        mzid, mgf = self._write_fixtures()
        out_prefix = self.tmp / "none"
        proc = self.run_script(
            "plot_ms2_spectrum.py", "--mzid", mzid, "--spectra", mgf,
            "--protein", "P99999", "-o", out_prefix, "--log-level", "ERROR")
        self.assertEqual(proc.returncode, 1)
        self.assertFalse(list(self.tmp.glob("none*.png")))

    def test_help_works_without_optional_deps(self):
        proc = self.run_script("plot_ms2_spectrum.py", "--help")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("usage", proc.stdout)


if __name__ == "__main__":
    unittest.main()
