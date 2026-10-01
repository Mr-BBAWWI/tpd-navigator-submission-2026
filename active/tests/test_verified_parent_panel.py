"""Strict verified-parent CLI and real Workbench policy-boundary tests."""
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest import mock
import io
import json
import tempfile
import unittest

from packages.contracts import ContractError
from packages.platform.store import Store
from packages.platform.workbench import WorkbenchService
from scripts import run_verified_parent_panel as runner


class RealPolicyBoundaryTests(unittest.TestCase):
    def test_missing_policy_rejects_through_workbench_without_job(self):
        with tempfile.TemporaryDirectory(prefix="verified-policy-") as directory:
            store = Store(Path(directory) / "store")
            service = WorkbenchService(store, "verified-policy-test")
            with self.assertRaisesRegex(
                ContractError,
                "JOB_SCIENTIFIC_POLICY_INVALID",
            ):
                service._resolve_design_policy(
                    "missing-policy",
                    "SMARCA2-FX5",
                )
            with store.db() as db:
                count = db.execute(
                    "SELECT COUNT(*) FROM workbench_jobs",
                ).fetchone()[0]
            self.assertEqual(count, 0)

    def test_built_payloads_validate_for_fx5_and_reference_parent(self):
        import packages.science.dual_e3 as dual_e3

        base = {
            "scientific_policy_id": "stored-m2-policy",
            "panel_size": 20,
            "linker_id": None,
        }
        fx5 = runner.build_payload(SimpleNamespace(
            parent_id="SMARCA2-FX5",
            **base,
        ))
        reference = runner.build_payload(SimpleNamespace(
            parent_id="SMARCA2-9D12-A1A1P",
            **base,
        ))

        expected_fx5 = {
            "exploratory": False,
            "dock": True,
            "use_api": False,
            "target": "SMARCA2",
            "panel_size": 20,
            "linker_ids": ["alkyl_c6"],
        }
        expected_reference = {
            **expected_fx5,
            "parent_id": "SMARCA2-9D12-A1A1P",
        }
        fx5_for_validation = {
            key: value
            for key, value in fx5.items()
            if key != "scientific_policy_id"
        }
        reference_for_validation = {
            key: value
            for key, value in reference.items()
            if key != "scientific_policy_id"
        }

        self.assertEqual(
            set(fx5) - set(fx5_for_validation),
            {"scientific_policy_id"},
        )
        self.assertEqual(
            set(reference) - set(reference_for_validation),
            {"scientific_policy_id"},
        )
        self.assertEqual(fx5_for_validation, expected_fx5)
        self.assertEqual(reference_for_validation, expected_reference)
        self.assertEqual(dual_e3.validate(fx5_for_validation), expected_fx5)
        self.assertEqual(
            dual_e3.validate(reference_for_validation),
            expected_reference,
        )
        self.assertNotIn("parent_id", fx5)
        self.assertEqual(fx5["panel_size"], 20)
        self.assertEqual(reference["parent_id"], "SMARCA2-9D12-A1A1P")


class CliTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="verified-parent-panel-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.data = self.root / "store"
        self.data.mkdir()
        self.index = self.data / "index.sqlite3"
        self.index.write_bytes(b"existing-store")
        self.receipt = self.root / "receipts" / "receipt.json"
        self.receipt.parent.mkdir()

    def arguments(self, receipt=None):
        return [
            "--data-dir", str(self.data),
            "--project", "project-a",
            "--parent-id", "SMARCA2-FX5",
            "--scientific-policy-id", "policy-17",
            "--receipt", str(receipt or self.receipt),
        ]

    def test_valid_orchestration_forwards_only_strict_payload(self):
        events = []
        expected_payload = {
            "exploratory": False,
            "dock": True,
            "use_api": False,
            "target": "SMARCA2",
            "scientific_policy_id": "policy-17",
            "panel_size": 20,
            "linker_ids": ["alkyl_c6"],
        }

        class FakeStore:
            def __init__(self, path):
                events.append(("store", path))

        class FakeService:
            def __init__(self, store, project):
                events.append(("service", project))

            def _resolve_design_policy(self, policy_id, parent_id):
                events.append(("resolve", policy_id, parent_id))
                return {"validated": True}

            def create(self, result_id, operation, request_key, payload):
                events.append(("create", result_id, operation, request_key, payload))
                return ({"id": "job-real", "state": "queued"}, True)

            def execute(self, job_id):
                events.append(("execute", job_id))

            def view(self, job_id):
                events.append(("view", job_id))
                return {
                    "id": job_id,
                    "state": "completed",
                    "error_code": None,
                    "result": {
                        "summary": {"selected_analogs": 13},
                    },
                    "parameters": expected_payload,
                }

        runtime = (FakeStore, FakeService, ContractError)
        output = io.StringIO()
        with mock.patch.object(runner, "_runtime", return_value=runtime):
            with redirect_stdout(output):
                code = runner.main(self.arguments())

        self.assertEqual(code, 0)
        create = next(event for event in events if event[0] == "create")
        self.assertEqual(create[1:3], ("design:SMARCA2", "design_panel"))
        self.assertTrue(create[3].startswith("verified-parent-"))
        self.assertEqual(create[4], expected_payload)
        self.assertEqual(
            events[2],
            ("resolve", "policy-17", "SMARCA2-FX5"),
        )
        self.assertFalse(any("approve" in event[0] for event in events))
        receipt = json.loads(self.receipt.read_text(encoding="utf-8"))
        self.assertEqual(receipt["id"], "job-real")
        self.assertEqual(receipt["parameters"], expected_payload)
        printed = json.loads(output.getvalue().splitlines()[-1])
        self.assertEqual(printed["selected_actual_count"], 13)

    def test_execution_failure_saves_actual_view_and_returns_failure(self):
        actual = {
            "id": "job-failed",
            "state": "completed",
            "error_code": None,
            "result": {"summary": {"selected_analogs": 11}},
        }

        class FakeStore:
            def __init__(self, path):
                pass

        class FakeService:
            def __init__(self, store, project):
                pass

            def _resolve_design_policy(self, policy_id, parent_id):
                return {"validated": True}

            def create(self, *args):
                return ({"id": "job-failed"}, True)

            def execute(self, job_id):
                raise RuntimeError("worker disconnected")

            def view(self, job_id):
                return actual

        error = io.StringIO()
        runtime = (FakeStore, FakeService, ContractError)
        with mock.patch.object(runner, "_runtime", return_value=runtime):
            with redirect_stdout(io.StringIO()), redirect_stderr(error):
                code = runner.main(self.arguments())
        self.assertEqual(code, 1)
        self.assertEqual(
            json.loads(self.receipt.read_text(encoding="utf-8")),
            actual,
        )
        self.assertIn("worker disconnected", error.getvalue())

    def test_policy_preflight_failure_writes_no_receipt_or_job(self):
        original_index = self.index.read_bytes()

        class FakeStore:
            def __init__(self, path):
                self.path = path

        class RejectingService:
            created = False
            executed = False

            def __init__(self, store, project):
                pass

            def _resolve_design_policy(self, policy_id, parent_id):
                raise ContractError("JOB_SCIENTIFIC_POLICY_INVALID")

            def create(self, *args):
                type(self).created = True

            def execute(self, *args):
                type(self).executed = True

        runtime = (FakeStore, RejectingService, ContractError)
        error = io.StringIO()
        with mock.patch.object(runner, "_runtime", return_value=runtime):
            with redirect_stderr(error):
                code = runner.main(self.arguments())

        self.assertEqual(code, 2)
        self.assertFalse(self.receipt.exists())
        self.assertFalse(RejectingService.created)
        self.assertFalse(RejectingService.executed)
        self.assertEqual(self.index.read_bytes(), original_index)
        self.assertIn("JOB_SCIENTIFIC_POLICY_INVALID", error.getvalue())
        self.assertIn("추가 개발 필요", error.getvalue())

    def test_existing_receipt_is_never_overwritten_before_runtime(self):
        self.receipt.write_text("keep-me", encoding="utf-8")
        with mock.patch.object(runner, "_runtime") as runtime:
            with redirect_stderr(io.StringIO()):
                code = runner.main(self.arguments())
        self.assertEqual(code, 2)
        self.assertEqual(self.receipt.read_text(encoding="utf-8"), "keep-me")
        runtime.assert_not_called()

    def test_missing_existing_store_index_fails_cleanly_before_runtime(self):
        self.index.unlink()
        error = io.StringIO()
        with mock.patch.object(runner, "_runtime") as runtime:
            with redirect_stderr(error):
                code = runner.main(self.arguments())
        self.assertEqual(code, 2)
        self.assertFalse(self.receipt.exists())
        self.assertIn("EXISTING_STORE_INDEX_REQUIRED", error.getvalue())
        self.assertNotIn("Traceback", error.getvalue())
        runtime.assert_not_called()

    def test_selected_actual_count_uses_only_selected_analogs(self):
        self.assertEqual(
            runner._selected_actual_count({
                "selected_analogs": 7,
                "requested_panel": 20,
            }),
            7,
        )
        self.assertIsNone(runner._selected_actual_count({
            "requested_panel": 20,
            "selected_count": 19,
        }))
        self.assertIsNone(runner._selected_actual_count({
            "selected_analogs": True,
        }))


if __name__ == "__main__":
    unittest.main()
