"""End-to-end tests for bin/deseq2_multigroup.R.

Runs the REAL DESeq2 pipeline (R_LIBS resolved from BUC_R_LIBS or ~/R/Rlib_*)
on a synthetic 50-gene x 6-sample matrix with planted DE genes. Validation
tests exit before the DESeq run and are fast; the end-to-end tests take
~30 s each.
"""

import os
import random
import subprocess
import unittest
from pathlib import Path

from common import BIN, ScriptTestCase, find_r_pair  # noqa: E402

SCRIPT = "deseq2_multigroup.R"
NEEDED_PKGS = ("DESeq2", "ashr", "pheatmap", "BiocParallel")

RSCRIPT, R_LIBS = find_r_pair(*NEEDED_PKGS)


def _make_counts(n_genes=50, n_per_group=2):
    """Deterministic synthetic matrix: g1-g10 up in treatA, g11-g20 down
    in treatB, g21+ background; gNA has a single read (NA padj candidate)."""
    rng = random.Random(42)
    samples = ([("control", i) for i in range(n_per_group)] +
               [("treatA", i) for i in range(n_per_group)] +
               [("treatB", i) for i in range(n_per_group)])
    header = ["gene"] + [f"S{g}{i}" for g, i in samples]
    rows = ["\t".join(header)]
    for g in range(1, n_genes + 1):
        base = rng.randint(10, 60)
        vals = []
        for grp, _ in samples:
            v = base + rng.randint(-5, 5)
            if g <= 10 and grp == "treatA":
                v = base + 60 + rng.randint(0, 20)
            elif 10 < g <= 20 and grp == "treatB":
                v = max(1, base // 4 + rng.randint(0, 3))
            vals.append(v)
        rows.append(f"g{g}\t" + "\t".join(map(str, vals)))
    rows.append("gNA\t" + "\t".join(["0"] * 5 + ["1"]))
    return "\n".join(rows) + "\n", [s for s, _ in samples]


SAMPLES_CSV = "id,group\nScontrol0,control\nScontrol1,control\n" \
              "StreatA0,treatA\nStreatA1,treatA\nStreatB0,treatB\nStreatB1,treatB\n"


@unittest.skipIf(R_LIBS is None,
                 f"R with {', '.join(NEEDED_PKGS)} not available (set BUC_R_LIBS)")
class TestDeseq2Multigroup(ScriptTestCase):

    def _env(self):
        env = dict(os.environ)
        if R_LIBS:
            env["R_LIBS"] = R_LIBS
        return env

    def _write_fixtures(self, counts_text=None, samples=SAMPLES_CSV):
        counts, _ = _make_counts()
        cpath = self.write("counts.tsv", counts_text or counts)
        spath = self.write("samples.csv", samples)
        return cpath, spath

    def _run(self, *args):
        return subprocess.run(
            [RSCRIPT, str(BIN / SCRIPT)] + [str(a) for a in args],
            capture_output=True, text=True, timeout=600, env=self._env())

    def _base(self, cpath, spath, out):
        return ("-c", str(cpath), "-s", str(spath), "-o", str(out), "-t", "2")

    def test_end_to_end_outputs(self):
        cpath, spath = self._write_fixtures()
        out = self.tmp / "res"
        proc = self._run(*self._base(cpath, spath, out))
        self.assertEqual(proc.returncode, 0, proc.stderr)

        # all logs on stderr, nothing on stdout
        self.assertEqual(proc.stdout, "")
        self.assertIn("[INFO] Running DESeq2...", proc.stderr)

        res_dir = out / "DEG_tables"
        deg = res_dir / "treatA_vs_control_DESeq2.tsv"
        self.assertTrue(deg.exists(), sorted(p.name for p in res_dir.iterdir()))
        header = deg.read_text().splitlines()[0]
        self.assertIn("regulation", header.split("\t"))
        self.assertIn("log2FoldChange", header.split("\t"))

        regs = {ln.split("\t")[-2] for ln in deg.read_text().splitlines()[1:]}
        self.assertIn("Up", regs)  # planted up-genes must be recovered

        for name in ("vst_Pearson_heatmap.pdf", "vst_PCA_plot.pdf",
                     "normalized_counts.tsv", "DEG_tables/treatB_vs_control_volcano.pdf"):
            self.assertTrue((out / name).exists(), name)

    def test_na_padj_reported_not_silent(self):
        cpath, spath = self._write_fixtures()
        proc = self._run(*self._base(cpath, spath, self.tmp / "res"))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("NA padj", proc.stderr)
        self.assertLessEqual(proc.stderr.count("[WARN]"), 4)

    def test_batch_flag_without_column_warns(self):
        counts, _ = _make_counts()
        cpath = self.write("counts.tsv", counts)
        spath = self.write("samples.csv", SAMPLES_CSV)
        proc = self._run("-c", str(cpath), "-s", str(spath),
                         "-o", str(self.tmp / "res"), "-t", "2", "-b")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("no 'batch' column", proc.stderr)

    def test_missing_group_column_fails(self):
        counts, _ = _make_counts()
        cpath = self.write("counts.tsv", counts)
        spath = self.write("samples.csv", "id\nScontrol0\nStreatA0\nStreatB0\n")
        proc = self._run("-c", str(cpath), "-s", str(spath),
                         "-o", str(self.tmp / "res"))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("missing required column(s): group", proc.stderr)

    def test_duplicate_sample_ids_fail(self):
        counts, _ = _make_counts()
        cpath = self.write("counts.tsv", counts)
        spath = self.write("samples.csv",
                           "id,group\nScontrol0,control\nScontrol0,control\n"
                           "StreatA0,treatA\nStreatA1,treatA\nStreatB0,treatB\nStreatB1,treatB\n")
        proc = self._run("-c", str(cpath), "-s", str(spath),
                         "-o", str(self.tmp / "res"))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Duplicate sample ids", proc.stderr)

    def test_non_integer_counts_fail(self):
        counts, samples = _make_counts()
        lines = counts.splitlines()
        cells = lines[1].split("\t")
        cells[1] = "10.5"
        lines[1] = "\t".join(cells)
        cpath = self.write("counts.tsv", "\n".join(lines) + "\n")
        spath = self.write("samples.csv", SAMPLES_CSV)
        proc = self._run("-c", str(cpath), "-s", str(spath),
                         "-o", str(self.tmp / "res"))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("non-integer", proc.stderr)

    def test_na_counts_fail(self):
        counts, _ = _make_counts()
        lines = counts.splitlines()
        cells = lines[1].split("\t")
        cells[1] = "NA"
        lines[1] = "\t".join(cells)
        cpath = self.write("counts.tsv", "\n".join(lines) + "\n")
        spath = self.write("samples.csv", SAMPLES_CSV)
        proc = self._run("-c", str(cpath), "-s", str(spath),
                         "-o", str(self.tmp / "res"))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("NA value(s)", proc.stderr)

    def test_unknown_control_fails_with_available_groups(self):
        counts, _ = _make_counts()
        cpath = self.write("counts.tsv", counts)
        spath = self.write("samples.csv", SAMPLES_CSV)
        proc = self._run("-c", str(cpath), "-s", str(spath),
                         "-o", str(self.tmp / "res"), "-r", "nosuch")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Control group 'nosuch' not found", proc.stderr)
        self.assertIn("control", proc.stderr)

    def test_threads_must_be_positive(self):
        counts, _ = _make_counts()
        cpath = self.write("counts.tsv", counts)
        spath = self.write("samples.csv", SAMPLES_CSV)
        proc = self._run("-c", str(cpath), "-s", str(spath),
                         "-o", str(self.tmp / "res"), "-t", "0")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("--threads must be >= 1", proc.stderr)

    def test_help_goes_to_stdout(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Usage: Rscript deseq2_multigroup.R", proc.stdout)
        self.assertIn("--fdr", proc.stdout)


if __name__ == "__main__":
    unittest.main()
