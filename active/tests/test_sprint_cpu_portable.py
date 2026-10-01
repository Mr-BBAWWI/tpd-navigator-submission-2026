"""Unit tests for the isolated v6 CPU launcher and its packaging hook."""
from __future__ import annotations

import errno
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
DELIVERY = ROOT / "reviewer_delivery"
if str(DELIVERY) not in sys.path:
    sys.path.insert(0, str(DELIVERY))

import build_sprint_package as sprint
import sprint_cpu_portable as cpu
from packages.science import dual_e3
from scripts import run_design_panel


def _success_job(attempted=8, completed=8):
    return {
        "state": "completed",
        "error_code": None,
        "result": {
            "summary": {
                "docking_attempted": attempted,
                "stage_counts": {"dock_all_eligible": {"attempted": attempted, "completed": completed}},
                "qualified_analogs": 0,
                "selected_analogs": 0,
            },
            "parent_redocking": {"real_docking": True, "status": "completed_with_limits"},
        },
    }


class CpuPortableTests(unittest.TestCase):
    def _run(self, root: Path, job: dict, returncode=0, side_effect=None, check_manifest=None, app=ROOT):
        calls = []

        def runner(argv, **kwargs):
            calls.append((argv, kwargs))
            if side_effect is not None:
                raise side_effect
            receipt = Path(argv[argv.index("--receipt") + 1])
            receipt.write_text(json.dumps(job), encoding="utf-8")
            kwargs["stdout"].write("mock child output\n")
            return subprocess.CompletedProcess(argv, returncode)

        with mock.patch.object(cpu, "ROOT", root), \
                mock.patch.object(cpu, "APP", app), \
                mock.patch.object(cpu.sys, "executable", sys.executable):
            manifest = check_manifest if check_manifest is not None else mock.Mock()
            code = cpu.run_cpu(check_manifest=manifest, runner=runner)
        return code, calls

    def test_exact_strict_argv_and_actual_docking_success(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            code, calls = self._run(root, _success_job())
            self.assertEqual(code, 0)
            runs = list((root / "reviewer-data/cpu-runs").iterdir())
            self.assertEqual(len(runs), 1)
            run = runs[0]
            expected = [
                sys.executable, "-X", "utf8", "-I", "-B",
                str(ROOT / "scripts/run_design_panel.py"),
                "--data-dir", str(run / "state"),
                "--panel-size", "10", "--linker-id", "alkyl_c6",
                "--max-elapsed-seconds", "900",
                "--receipt", str(run / "job-receipt.json"),
            ]
            self.assertEqual(calls[0][0], expected)
            self.assertNotIn("--no-dock", calls[0][0])
            self.assertNotIn("--exploratory", calls[0][0])
            self.assertNotIn("--api", calls[0][0])
            self.assertIs(calls[0][1]["stderr"], subprocess.STDOUT)
            self.assertFalse(calls[0][1]["check"])
            self.assertEqual(calls[0][1]["cwd"], root)
            self.assertNotIn("DACON_API_KEY", calls[0][1]["env"])
            receipt = json.loads((run / "launcher-receipt.json").read_text(encoding="utf-8"))
            self.assertEqual(receipt["status"], "completed")
            self.assertEqual(receipt["mode"], "strict_cpu_docking")
            self.assertEqual(receipt["panel_size_requested"], 10)
            self.assertTrue(receipt["parent_real_docking"])
            self.assertEqual(receipt["actual_summary"]["qualified_analogs"], 0)
            self.assertFalse(receipt["scientific_approval"])

    def test_helper_argv_builds_payload_accepted_by_real_validator(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            code, calls = self._run(root, _success_job())
            self.assertEqual(code, 0)
            argv = calls[0][0]
            script_index = argv.index(str(ROOT / "scripts/run_design_panel.py"))
            args = run_design_panel.build_parser().parse_args(argv[script_index + 1:])
            # Match the exact CLI payload construction without executing science.
            payload = {
                "exploratory": args.exploratory,
                "dock": not args.no_dock,
                "use_api": args.api,
                "panel_size": args.panel_size,
            }
            if args.linker_id:
                payload["linker_ids"] = args.linker_id
            validated = dual_e3.validate(payload)
            self.assertEqual(validated["panel_size"], 10)
            self.assertTrue(validated["dock"])
            self.assertFalse(validated["exploratory"])
            self.assertFalse(validated["use_api"])
            self.assertEqual(validated["linker_ids"], ["alkyl_c6"])
            with self.assertRaises(ValueError):
                dual_e3.validate({"panel_size": 6})

    def test_failed_receipt_and_zero_docking_are_nonzero(self):
        fixtures = [
            ({"state": "failed", "error_code": "DOCK_FAILED", "result": None}, 1),
            (_success_job(attempted=0, completed=0), 0),
        ]
        for job, child_code in fixtures:
            with self.subTest(job=job), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                code, _ = self._run(root, job, child_code)
                self.assertNotEqual(code, 0)
                run = next((root / "reviewer-data/cpu-runs").iterdir())
                receipt = json.loads((run / "launcher-receipt.json").read_text(encoding="utf-8"))
                self.assertEqual(receipt["status"], "failed")
                self.assertTrue((run / "failure-note.txt").is_file())
                self.assertIn("STRICT_CPU_DOCKING_FAILED", (run / "run.log").read_text(encoding="utf-8"))

    def test_runner_exception_has_safe_log_and_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            code, _ = self._run(root, {}, side_effect=RuntimeError("secret fixture text"))
            self.assertNotEqual(code, 0)
            run = next((root / "reviewer-data/cpu-runs").iterdir())
            log = (run / "run.log").read_text(encoding="utf-8")
            receipt = json.loads((run / "launcher-receipt.json").read_text(encoding="utf-8"))
            self.assertNotIn("secret fixture text", log)
            self.assertEqual(receipt["error_code"], "CPU_LAUNCH_EXCEPTION")
            self.assertTrue((run / "failure-note.txt").is_file())

    def test_manifest_failure_has_safe_stage_specific_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            failed_manifest = mock.Mock(side_effect=RuntimeError("PACKAGE_HASH_MISMATCH: secret detail"))
            code, calls = self._run(root, {}, check_manifest=failed_manifest)
            self.assertNotEqual(code, 0)
            self.assertEqual(calls, [])
            run = next((root / "reviewer-data/cpu-runs").iterdir())
            receipt_path = run / "launcher-receipt.json"
            self.assertTrue((run / "state").is_dir())
            self.assertTrue(receipt_path.is_file())
            self.assertTrue((run / "failure-note.txt").is_file())
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            self.assertEqual(receipt["error_code"], "CPU_MANIFEST_CHECK_FAILED")
            artifacts = "\n".join(
                path.read_text(encoding="utf-8")
                for path in (run / "run.log", run / "failure-note.txt", receipt_path)
            )
            self.assertIn("CPU_MANIFEST_CHECK_FAILED", artifacts)
            self.assertNotIn("PACKAGE_HASH_MISMATCH", artifacts)
            self.assertNotIn("secret detail", artifacts)

    def test_missing_script_has_specific_error_and_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            missing_app = root / "missing-app"
            code, calls = self._run(root, {}, app=missing_app)
            self.assertNotEqual(code, 0)
            self.assertEqual(calls, [])
            run = next((root / "reviewer-data/cpu-runs").iterdir())
            receipt_path = run / "launcher-receipt.json"
            self.assertTrue(receipt_path.is_file())
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            self.assertEqual(receipt["error_code"], "CPU_SCRIPT_MISSING")
            self.assertTrue((run / "failure-note.txt").is_file())

    def test_existing_user_state_preserved_and_runs_unique(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "design-data").mkdir()
            (root / "reviewer-data").mkdir()
            design = root / "design-data/user.bin"
            reviewer = root / "reviewer-data/user.json"
            design.write_bytes(b"design-state")
            reviewer.write_bytes(b"review-state")
            self._run(root, _success_job())
            self._run(root, _success_job())
            self.assertEqual(design.read_bytes(), b"design-state")
            self.assertEqual(reviewer.read_bytes(), b"review-state")
            self.assertEqual(len(list((root / "reviewer-data/cpu-runs").iterdir())), 2)

    def test_linked_reviewer_data_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            outside = root / "outside"
            outside.mkdir()
            link = root / "reviewer-data"
            try:
                link.symlink_to(outside, target_is_directory=True)
            except OSError as exc:
                denied = exc.errno in {errno.EACCES, errno.EPERM} or getattr(exc, "winerror", None) == 1314
                if os.name == "nt" and denied:
                    self.skipTest("Windows account lacks symlink permission")
                raise
            checked = mock.Mock()
            with mock.patch.object(cpu, "ROOT", root):
                with self.assertRaisesRegex(RuntimeError, "PARENT_INVALID"):
                    cpu.run_cpu(check_manifest=checked, runner=mock.Mock())
            checked.assert_not_called()
            self.assertEqual(list(outside.iterdir()), [])

    def test_isolated_help_from_unrelated_cwd(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = subprocess.run(
                [sys.executable, "-X", "utf8", "-I", str(DELIVERY / "sprint_cpu_portable.py"), "--help"],
                cwd=temporary,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                check=False,
            )
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("cpu", result.stdout)


class SprintCpuPackagingTests(unittest.TestCase):
    def test_release_payload_requests_validated_minimum_panel(self):
        payload = sprint._release_payload()
        self.assertEqual(payload["cpu_trial"]["panel_size_requested"], 10)
        self.assertEqual(payload["cpu_trial"]["linker_id"], "alkyl_c6")

    def test_stage_hook_copies_helper_and_replaces_only_v6_launcher(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            active = root / "active"
            source = active / "reviewer_delivery/sprint_cpu_portable.py"
            source.parent.mkdir(parents=True)
            source.write_text("# helper fixture\n", encoding="utf-8")
            stage = root / "stage"
            (stage / "program").mkdir(parents=True)
            inherited = b'"%~dp0runtime\\python.exe" -X utf8 -I -B "%~dp0program\\scientific_portable.py" cpu\r\n'
            (stage / "02_CPU_직접시험.cmd").write_bytes(inherited)
            original_callback = sprint._ORIGINAL_STAGE_RELEASE_FILES
            with mock.patch.object(sprint, "ACTIVE", active), \
                    mock.patch.object(sprint, "_ORIGINAL_STAGE_RELEASE_FILES", mock.Mock()):
                sprint._stage_release_files(stage)
            launcher = (stage / "02_CPU_직접시험.cmd").read_text(encoding="utf-8")
            self.assertIn("program\\sprint_cpu_portable.py", launcher)
            self.assertNotIn("scientific_portable.py", launcher)
            self.assertIn("reviewer-data\\cpu-runs", launcher)
            self.assertIn("pause", launcher.lower())
            self.assertEqual((stage / "program/sprint_cpu_portable.py").read_text(), "# helper fixture\n")
            self.assertIs(sprint._ORIGINAL_STAGE_RELEASE_FILES, original_callback)

    def test_stage_hook_refuses_existing_helper(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            active = root / "active"
            source = active / "reviewer_delivery/sprint_cpu_portable.py"
            source.parent.mkdir(parents=True)
            source.write_text("source", encoding="utf-8")
            stage = root / "stage"
            (stage / "program").mkdir(parents=True)
            (stage / "program/sprint_cpu_portable.py").write_text("existing", encoding="utf-8")
            (stage / "02_CPU_직접시험.cmd").write_text("scientific_portable.py cpu", encoding="utf-8")
            with mock.patch.object(sprint, "ACTIVE", active), \
                    mock.patch.object(sprint, "_ORIGINAL_STAGE_RELEASE_FILES", mock.Mock()):
                with self.assertRaisesRegex(RuntimeError, "HELPER_EXISTS"):
                    sprint._stage_release_files(stage)


if __name__ == "__main__":
    unittest.main()
