from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import MappingProxyType
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "verify_parent_funnel.py"
VALID_SHA256 = "0" * 64


class VerifyParentFunnelCliTests(unittest.TestCase):
    def run_isolated(self, *arguments: str, cwd: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, "-I", str(SCRIPT), *arguments],
            cwd=cwd,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_help_works_under_isolated_python_from_temporary_cwd(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            completed = self.run_isolated("--help", cwd=Path(directory))
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--bundle", completed.stdout)
        self.assertIn("--expected-manifest-sha256", completed.stdout)

    def test_malformed_sha_is_rejected_by_actual_verifier(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bundle = root / "bundle"
            bundle.mkdir()
            before = tuple(bundle.iterdir())
            completed = self.run_isolated(
                "--bundle",
                str(bundle),
                "--expected-manifest-sha256",
                "not-a-sha256",
                cwd=root,
            )
            after = tuple(bundle.iterdir())
        self.assertEqual(completed.returncode, 2)
        self.assertEqual(completed.stdout, "")
        self.assertIn("VERIFY_FAILED", completed.stderr)
        self.assertIn("EXPECTED_MANIFEST_SHA256_INVALID", completed.stderr)
        self.assertEqual(before, after)

    def test_missing_bundle_is_rejected_by_actual_verifier(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            missing = root / "does-not-exist"
            completed = self.run_isolated(
                "--bundle",
                str(missing),
                "--expected-manifest-sha256",
                VALID_SHA256,
                cwd=root,
            )
            self.assertFalse(missing.exists())
        self.assertEqual(completed.returncode, 2)
        self.assertEqual(completed.stdout, "")
        self.assertIn("VERIFY_FAILED", completed.stderr)

    def test_success_interface_uses_only_mocked_verifier_result(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "verify_parent_funnel_cli_test_module",
            SCRIPT,
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        parent = MappingProxyType({
            "parent_id": "PARENT-1",
            "parent_redocking": "pass",
            "actual_counts": MappingProxyType({
                "cheap_passed": 3,
                "pose_qualified": 1,
                "selected": 0,
                "strict_qualified": 0,
            }),
            "evidence_tiers": MappingProxyType({
                "exact_source_join": False,
                "missing_evidence": True,
            }),
            "actual_families": MappingProxyType({
                "cheap_passed": ("ring_modification",),
                "pose_qualified": ("ring_modification",),
                "strict_qualified": (),
            }),
            "unbounded_internal_detail": "must not be printed",
        })
        result = MappingProxyType({
            "format": "parent-funnel-dossier/1.0",
            "status": "verified_computed_diagnostic",
            "manifest_sha256": VALID_SHA256,
            "parent_count": 1,
            "human_selection_count": 0,
            "strict_trusted_selected_count": 0,
            "scientific_approval": False,
            "m2_registered": False,
            "reference_verification": "build_time_external_check_only_not_independently_reverified",
            "parents": (parent,),
            "compact_summary": {"large": "not part of the CLI projection"},
        })

        stdout = io.StringIO()
        stderr = io.StringIO()
        with mock.patch.object(module, "verify_funnel_bundle", return_value=result) as verifier:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                status = module.main([
                    "--bundle",
                    "read-only-bundle",
                    "--expected-manifest-sha256",
                    VALID_SHA256,
                ])

        self.assertEqual(status, 0)
        self.assertEqual(stderr.getvalue(), "")
        verifier.assert_called_once_with("read-only-bundle", VALID_SHA256)
        payload = json.loads(stdout.getvalue())
        self.assertEqual(payload["top_counts"]["parent_count"], 1)
        self.assertEqual(payload["parents"][0]["actual_counts"]["pose_qualified"], 1)
        self.assertEqual(
            payload["parents"][0]["actual_families"]["cheap_passed"],
            ["ring_modification"],
        )
        self.assertNotIn("compact_summary", payload)
        self.assertNotIn("unbounded_internal_detail", payload["parents"][0])


if __name__ == "__main__":
    unittest.main()
