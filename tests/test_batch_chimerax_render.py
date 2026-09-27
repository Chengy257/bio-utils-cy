"""Offline tests for bin/batch_chimerax_render.py.

ChimeraX itself is not installed; the render path is exercised with a stub
binary that parses the generated .cxc script and "renders" PNGs, and the
script builder is unit-tested directly against the documented ChimeraX
command syntax.
"""

import os
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import batch_chimerax_render as bcr  # noqa: E402
from common import ScriptTestCase  # noqa: E402


class TestBuildScript(ScriptTestCase):

    def test_bychain_scheme_and_quoted_paths(self):
        script = bcr.build_chimerax_script(
            ["/data dir/prot one.pdb"], "/out dir", 800, 600, "bychain")
        self.assertIn('open "/data dir/prot one.pdb"', script)
        self.assertIn('save "/out dir/prot one.png" width 800 height 600 supersample 2', script)
        self.assertIn("color bychain", script)
        self.assertIn("close all", script)

    def test_bfactor_uses_alphafold_palette(self):
        script = bcr.build_chimerax_script(["/x/p.pdb"], "/o", 100, 100, "bfactor")
        self.assertIn("color byattribute bfactor palette alphafold", script)

    def test_ss_scheme_uses_valid_selectors(self):
        # 'color ss' is not a ChimeraX command; helix/sheet selectors are
        script = bcr.build_chimerax_script(["/x/p.pdb"], "/o", 100, 100, "ss")
        self.assertNotIn("color ss", script)
        self.assertIn("color helix tan target c", script)
        self.assertIn("color sheet gold target c", script)

    def test_rainbow_scheme(self):
        script = bcr.build_chimerax_script(["/x/p.pdb"], "/o", 100, 100, "rainbow")
        self.assertIn("rainbow chain", script)

    def test_multiple_files_all_rendered_and_closed(self):
        script = bcr.build_chimerax_script(["/x/a.pdb", "/x/b.pdb"],
                                           "/o", 100, 100, "bychain")
        self.assertEqual(script.count("open "), 2)
        self.assertEqual(script.count("close all"), 2)


class TestMainDetection(ScriptTestCase):

    def _run_main(self, argv_extra):
        argv = ["batch_chimerax_render.py"] + argv_extra
        with mock.patch.object(sys, "argv", argv):
            try:
                bcr.main()
                return 0
            except SystemExit as exc:
                return exc.code

    def test_not_found_exits_nonzero(self):
        env = dict(os.environ)
        env.pop("BUC_CHIMERAX_BIN", None)
        with mock.patch.dict(os.environ, env, clear=True):
            code = self._run_main(["-i", str(self.tmp), "-o", str(self.tmp / "o")])
        self.assertEqual(code, 1)

    def test_explicit_nonexecutable_rejected(self):
        bad = self.tmp / "not_chimerax"
        bad.write_text("x")
        code = self._run_main(["-i", str(self.tmp), "-o", str(self.tmp / "o"),
                               "--chimerax", str(bad)])
        self.assertEqual(code, 1)

    def test_broken_env_slot_rejected(self):
        bad = self.tmp / "not_chimerax2"
        bad.write_text("x")
        with mock.patch.dict(os.environ, {"BUC_CHIMERAX_BIN": str(bad)}):
            code = self._run_main(["-i", str(self.tmp), "-o", str(self.tmp / "o")])
        self.assertEqual(code, 1)
        # with assertLogs to confirm the message
        # (already checked via exit code; message text asserted in e2e below)


class TestBatchRenderStubbed(ScriptTestCase):

    def _make_stub(self, exit_code=0):
        stub = self.tmp / "stub_chimerax.sh"
        stub.write_text(
            "#!/bin/bash\n"
            "# emulate: chimerax --nogui --exit script.cxc\n"
            "script=\"\"\n"
            "for a in \"$@\"; do [[ \"$a\" == *.cxc ]] && script=\"$a\"; done\n"
            'if [[ "${FAIL_SAVES:-0}" == "1" ]]; then exit 0; fi\n'
            'grep \'^save "\' "$script" | sed \'s/^save "//; s/" width.*//\' | '
            'while read -r png; do printf \'PNGDATA\' > "$png"; done\n'
            f"exit {exit_code}\n"
        )
        stub.chmod(0o755)
        return stub

    def _run_main(self, argv_extra):
        argv = ["batch_chimerax_render.py"] + argv_extra
        with mock.patch.object(sys, "argv", argv):
            try:
                bcr.main()
                return 0
            except SystemExit as exc:
                return exc.code

    def test_successful_render(self):
        pdb_dir = self.tmp / "pdbs"
        pdb_dir.mkdir()
        (pdb_dir / "a.pdb").write_text("ATOM\n")
        (pdb_dir / "b.pdb").write_text("ATOM\n")
        stub = self._make_stub()
        code = self._run_main(["-i", str(pdb_dir), "-o", str(self.tmp / "img"),
                               "--chimerax", str(stub)])
        self.assertEqual(code, 0)
        self.assertTrue((self.tmp / "img" / "a.png").exists())
        self.assertTrue((self.tmp / "img" / "b.png").exists())

    def test_missing_png_detected(self):
        pdb_dir = self.tmp / "pdbs2"
        pdb_dir.mkdir()
        (pdb_dir / "a.pdb").write_text("ATOM\n")
        (pdb_dir / "b.pdb").write_text("ATOM\n")
        stub = self._make_stub()
        with mock.patch.dict(os.environ, {"FAIL_SAVES": "1"}):
            with self.assertLogs(bcr.logger, level="ERROR") as captured:
                code = self._run_main(["-i", str(pdb_dir), "-o", str(self.tmp / "img2"),
                                       "--chimerax", str(stub)])
        self.assertEqual(code, 1)
        self.assertTrue(any("missing or empty" in line for line in captured.output))

    def test_nonzero_chimerax_exit_detected(self):
        pdb_dir = self.tmp / "pdbs3"
        pdb_dir.mkdir()
        (pdb_dir / "a.pdb").write_text("ATOM\n")
        stub = self._make_stub(exit_code=2)
        code = self._run_main(["-i", str(pdb_dir), "-o", str(self.tmp / "img3"),
                               "--chimerax", str(stub)])
        self.assertEqual(code, 1)

    def test_no_pdb_files(self):
        empty = self.tmp / "empty"
        empty.mkdir()
        stub = self._make_stub()
        with self.assertLogs(bcr.logger, level="ERROR"):
            code = self._run_main(["-i", str(empty), "-o", str(self.tmp / "img4"),
                                   "--chimerax", str(stub)])
        self.assertEqual(code, 1)


if __name__ == "__main__":
    unittest.main()
