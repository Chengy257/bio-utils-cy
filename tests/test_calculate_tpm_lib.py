"""Tests for lib/calculate_tpm.R — the shared expression-unit library.

An R driver sources the library and asserts hand-computed values plus the
guard rails. Pure base R; any working Rscript qualifies.
"""

import subprocess
import unittest
from pathlib import Path

from common import REPO_ROOT, RSCRIPT, ScriptTestCase  # noqa: E402

LIB = REPO_ROOT / "lib" / "calculate_tpm.R"

DRIVER = r'''
source("{lib}")

ok <- function(name, cond) cat(if (isTRUE(cond)) "OK" else "FAIL", name, "\n")

counts <- c(10, 20, 30); eff <- c(100, 200, 300)

# TPM: rate = c(0.1, 0.1, 0.1) -> equal shares of 1e6
tpm <- countToTpm(counts, eff)
ok("tpm_hand", isTRUE(all.equal(tpm, rep(1e6 / 3, 3), tolerance = 1e-8)))
ok("tpm_sums_1e6", isTRUE(all.equal(sum(tpm), 1e6, tolerance = 1e-8)))

# FPKM: N=60 -> count/(eff/1e3)/(6e-5) = 1666666.67 for all three
fpkm <- countToFpkm(counts, eff)
ok("fpkm_hand", isTRUE(all.equal(fpkm, rep(10 / 0.1 / (60 / 1e6), 3), tolerance = 1e-8)))

# fpkmToTpm(countToFpkm) == countToTpm (unit identity)
ok("fpkm_tpm_identity", isTRUE(all.equal(fpkmToTpm(fpkm), tpm, tolerance = 1e-8)))

# RPM: 10/60, 20/60, 30/60 of 1e6
ok("rpm_hand", isTRUE(all.equal(countToRpm(counts),
                                c(10, 20, 30) / 60 * 1e6, tolerance = 1e-8)))

# effective counts: counts * len/effLen
ok("effcounts_hand", isTRUE(all.equal(countToEffCounts(counts, c(50, 100, 150), eff),
                                      c(5, 10, 15), tolerance = 1e-8)))

# --- guards: each must stop ---
guard <- function(name, expr, pattern) {{
    caught <- tryCatch({{ expr; "NO-ERROR" }}, error = function(e) conditionMessage(e))
    cat(if (grepl(pattern, caught, fixed = TRUE)) "OK" else "FAIL", name, "|",
        substr(caught, 1, 60), "\n")
}}
guard("len_mismatch", countToTpm(c(1, 2), c(100)), "equal length")
guard("zero_efflen", countToTpm(c(1, 2), c(100, 0)), "positive")
guard("allzero_tpm", countToTpm(c(0, 0, 0), eff), "undefined")
guard("allzero_fpkm", countToFpkm(c(0, 0, 0), eff), "undefined")
guard("allzero_rpm", countToRpm(c(0, 0, 0)), "undefined")
guard("allzero_fpkm2tpm", fpkmToTpm(c(0, 0, 0)), "undefined")
guard("empty_fpkm2tpm", fpkmToTpm(numeric(0)), "non-empty")
guard("effcounts_mismatch", countToEffCounts(c(1, 2), c(1), c(10, 10)), "equal length")
guard("effcounts_zero_eff", countToEffCounts(c(1, 2), c(1, 2), c(10, 0)), "positive")
'''


@unittest.skipUnless(RSCRIPT and Path(RSCRIPT).exists(), "no working Rscript")
class TestCalculateTpmLib(ScriptTestCase):

    def test_library_functions_and_guards(self):
        driver = self.write("driver.R", DRIVER.format(lib=str(LIB)))
        proc = subprocess.run([RSCRIPT, str(driver)], capture_output=True,
                              text=True, timeout=300)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        lines = {ln.split()[1]: ln for ln in proc.stdout.splitlines() if ln}
        fails = [k for k, ln in lines.items() if ln.startswith("FAIL")]
        self.assertEqual(fails, [], proc.stdout)
        for name in ("tpm_hand", "tpm_sums_1e6", "fpkm_hand", "fpkm_tpm_identity",
                     "rpm_hand", "effcounts_hand", "len_mismatch", "zero_efflen",
                     "allzero_tpm", "allzero_fpkm", "allzero_rpm", "allzero_fpkm2tpm",
                     "empty_fpkm2tpm", "effcounts_mismatch", "effcounts_zero_eff"):
            self.assertIn(name, lines, proc.stdout)

    def test_no_stdout_pollution_on_source(self):
        driver = self.write("quiet.R", f'source("{LIB}")\ncat("SOURCED-QUIETLY\\n")\n')
        proc = subprocess.run([RSCRIPT, str(driver)], capture_output=True,
                              text=True, timeout=300)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout, "SOURCED-QUIETLY\n")  # no banner


if __name__ == "__main__":
    unittest.main()
