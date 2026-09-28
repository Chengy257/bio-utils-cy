"""Functional tests for bin/batch_parse_sanger.R.

Builds a minimal but valid ABIF (.ab1) trace in pure Python: big-endian
header + directory entries (sangerseqR reads element type 4 as 2-byte
ints and stores <=4-byte payloads inline in the entry), then runs the
real sangerseqR pipeline end to end.
"""

import math
import os
import struct
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "sangerseqR", "Biostrings")


def build_ab1(seq, points_per_base=10, peak_amp=3000, baseline=200):
    n = len(seq)
    total = n * points_per_base
    ploc = [i * points_per_base + points_per_base // 2 for i in range(n)]
    chan_of = {"A": 9, "C": 10, "G": 11, "T": 12}
    channels = {v: [baseline] * total for v in range(9, 13)}
    for i, base in enumerate(seq):
        center = ploc[i]
        for off in range(-4, 5):
            p = center + off
            if 0 <= p < total:
                channels[chan_of[base]][p] += int(
                    peak_amp * math.exp(-(off ** 2) / 4.0))

    entries, blobs = [], []

    def add(name, tagnum, elemtype, elemsize, payload, inline=False):
        nelem = len(payload) // elemsize
        if inline:
            datafield = payload + b"\x00" * (8 - len(payload))
            entries.append(struct.pack(">4sihhiiii", name.encode(), tagnum,
                                       elemtype, elemsize, nelem,
                                       len(payload), 0, 0)[:20] + datafield)
            blobs.append(None)
        else:
            entries.append(struct.pack(">4sihhiiii", name.encode(), tagnum,
                                       elemtype, elemsize, nelem,
                                       len(payload), 0, 0))
            blobs.append(payload)

    add("FWO_", 1, 2, 1, b"ACGT", inline=True)
    add("PBAS", 2, 2, 1, seq.encode())
    add("PLOC", 2, 4, 2, struct.pack(">%dh" % n, *ploc))
    for tag in (9, 10, 11, 12):
        add("DATA", tag, 4, 2, struct.pack(">%dh" % total, *channels[tag]))

    header = struct.pack(">4s h 4s i h h i i i i",
                         b"ABIF", 1, b"sang", 0, 0, 28, len(entries), 28,
                         128, 0)
    header += b"\x00" * (128 - len(header))

    offset = 128 + 28 * len(entries)
    fixed = []
    for e, blob in zip(entries, blobs):
        if blob is None:
            fixed.append(e)
            continue
        e = bytearray(e)
        e[20:24] = struct.pack(">i", offset)
        fixed.append(bytes(e))
        offset += len(blob) + (len(blob) % 2)

    body = b"".join(
        blob + (b"\x00" if len(blob) % 2 else b"")
        for blob in blobs if blob is not None)
    return header + b"".join(fixed) + body


@unittest.skipUnless(RS, "R with getopt+sangerseqR+Biostrings not available")
class TestBatchParseSanger(ScriptTestCase):
    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        old_cwd = os.getcwd()
        os.chdir(self.tmp)  # default input/output dirs are the CWD
        try:
            return subprocess.run([RS, str(BIN / "batch_parse_sanger.R")]
                                  + [str(a) for a in args],
                                  capture_output=True, text=True, timeout=600,
                                  env=env)
        finally:
            os.chdir(old_cwd)

    def _ref(self, seq="ACGT" * 20):
        return self.write("ref.txt", seq + "\n")

    def test_full_success_path(self):
        ab1 = self.tmp / "good.ab1"
        ab1.write_bytes(build_ab1("ACGT" * 20))
        ref = self._ref()
        out = self.tmp / "res"
        proc = self._run("-r", ref, "-i", self.tmp, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Processed: 1 / 1 files successfully", proc.stderr)
        self.assertTrue((out / "good_chromatogram.pdf").exists())
        self.assertTrue((out / "good_pairwiseAlignment_ref.txt").exists())
        self.assertTrue((out / "good_pairwiseAlignment_allele.txt").exists())

    def test_all_fail_exits_1(self):
        # regression: the batch always exited 0 even when nothing parsed
        (self.tmp / "broken.ab1").write_bytes(b"NOTABIF" * 20)
        ref = self._ref()
        proc = self._run("-r", ref, "-i", self.tmp, "-o", "res")
        self.assertEqual(proc.returncode, 1, proc.stderr)
        self.assertIn("No .ab1 file was processed successfully", proc.stderr)

    def test_partial_failure_still_reports(self):
        (self.tmp / "good.ab1").write_bytes(build_ab1("ACGT" * 20))
        (self.tmp / "bad.ab1").write_bytes(b"garbage" * 10)
        ref = self._ref()
        proc = self._run("-r", ref, "-i", self.tmp, "-o", "res")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Processed: 1 / 2 files successfully", proc.stderr)
        self.assertIn("Error processing bad.ab1", proc.stderr)

    def test_signal_cutoff_validation(self):
        ref = self._ref()
        for bad in ("0", "1", "1.5", "-0.1"):
            proc = self._run("-r", ref, "-i", self.tmp, "-o", "res",
                             "-s", bad)
            self.assertEqual(proc.returncode, 1, bad)
            self.assertIn("strictly between 0 and 1", proc.stderr)

    def test_trim5_validation(self):
        ref = self._ref()
        proc = self._run("-r", ref, "-i", self.tmp, "-o", "res", "-t", "-3")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("non-negative integer", proc.stderr)

    def test_no_ab1_files(self):
        ref = self._ref()
        proc = self._run("-r", ref, "-i", self.tmp, "-o", "res")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No .ab1 files found", proc.stderr)

    def test_missing_reference(self):
        proc = self._run("-r", self.tmp / "nope.txt", "-i", self.tmp)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Reference file not found", proc.stderr)

    def test_empty_reference(self):
        ref = self.write("empty_ref.txt", ">only_a_header\n\n")
        proc = self._run("-r", ref, "-i", self.tmp, "-o", "res")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("no sequence data", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
