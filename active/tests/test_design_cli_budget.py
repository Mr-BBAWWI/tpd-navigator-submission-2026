"""Unit tests for the design-panel CLI execution budget."""
import copy
import io
import unittest
from contextlib import redirect_stderr

from packages.platform import workbench
from scripts import run_design_panel


class CapturingService:
    def __init__(self, error=None):
        self.error = error
        self.limits = None
        self.arguments = None

    def create(self, result_id, kind, request_key, payload):
        self.limits = copy.deepcopy(workbench.LIMITS)
        self.arguments = (result_id, kind, request_key, copy.deepcopy(payload))
        if self.error is not None:
            raise self.error
        return {"id": "synthetic-job"}, False


class DesignCliBudgetTests(unittest.TestCase):
    def _parse(self, *arguments):
        return run_design_panel.build_parser().parse_args([
            *arguments,
            "--receipt",
            "synthetic-receipt.json",
        ])

    def test_argparse_accepts_inclusive_budget_bounds(self):
        self.assertEqual(self._parse("--max-elapsed-seconds", "1").max_elapsed_seconds, 1)
        self.assertEqual(
            self._parse("--max-elapsed-seconds", "7200").max_elapsed_seconds,
            7200,
        )

    def test_argparse_rejects_invalid_budgets(self):
        for value in ("0", "-1", "7201", "not-an-integer"):
            with self.subTest(value=value), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    self._parse("--max-elapsed-seconds", value)

    def test_omitted_override_leaves_default_900_unchanged(self):
        original_limits = workbench.LIMITS
        self.assertEqual(original_limits["max_elapsed_seconds"], 900)
        arguments = self._parse()
        self.assertIsNone(arguments.max_elapsed_seconds)

        service = CapturingService()
        run_design_panel.create_job(
            service,
            "design:synthetic",
            "synthetic-request",
            {"panel_size": 1},
            arguments.max_elapsed_seconds,
        )

        self.assertEqual(service.limits["max_elapsed_seconds"], 900)
        self.assertIs(workbench.LIMITS, original_limits)

    def test_scoped_override_is_restored_after_success_and_failure(self):
        original_limits = workbench.LIMITS

        successful = CapturingService()
        run_design_panel.create_job(
            successful,
            "design:synthetic",
            "successful-request",
            {},
            123,
        )
        self.assertEqual(successful.limits["max_elapsed_seconds"], 123)
        self.assertIs(workbench.LIMITS, original_limits)

        failing = CapturingService(RuntimeError("synthetic creation failure"))
        with self.assertRaisesRegex(RuntimeError, "synthetic creation failure"):
            run_design_panel.create_job(
                failing,
                "design:synthetic",
                "failing-request",
                {},
                456,
            )
        self.assertEqual(failing.limits["max_elapsed_seconds"], 456)
        self.assertIs(workbench.LIMITS, original_limits)

    def test_creation_receives_request_budget_and_preserves_all_other_limits(self):
        original_limits = workbench.LIMITS
        original_values = copy.deepcopy(original_limits)
        payload = {"exploratory": False, "dock": True, "panel_size": 12}
        service = CapturingService()

        run_design_panel.create_job(
            service,
            "design:synthetic",
            "budget-request",
            payload,
            321,
        )

        expected = copy.deepcopy(original_values)
        expected["max_elapsed_seconds"] = 321
        self.assertEqual(service.limits, expected)
        for name, value in original_values.items():
            if name != "max_elapsed_seconds":
                self.assertEqual(service.limits[name], value)
        self.assertEqual(
            service.arguments,
            (
                "design:synthetic",
                "design_panel",
                "budget-request",
                payload,
            ),
        )
        self.assertIs(workbench.LIMITS, original_limits)
        self.assertEqual(workbench.LIMITS, original_values)


if __name__ == "__main__":
    unittest.main()
