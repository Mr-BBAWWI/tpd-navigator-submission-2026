import os
from pathlib import Path
import tempfile
import unittest

from packages.platform.scientific_acceptance import ScientificAcceptanceService


def _service():
    return object.__new__(ScientificAcceptanceService)


def _files(run_root: Path):
    prediction = run_root / "seed-23" / "boltz_output" / "model.cif"
    confidence = run_root / "seed-23" / "boltz_output" / "confidence_model.json"
    prediction.parent.mkdir(parents=True)
    prediction.write_text("prediction", encoding="utf-8")
    confidence.write_text("{}", encoding="utf-8")
    return prediction, confidence


def _inspection(prediction, confidence, metric=0.75):
    return {
        "status": "computed_hypothesis",
        "actual_computation": True,
        "prediction": str(prediction),
        "confidence_file": str(confidence),
        "descriptive_metrics": {"interface_confidence": metric},
        "mapping": {"chain_a": "A", "chain_b": "B"},
        "confidence": "high",
    }


class TernaryInspectionPathTests(unittest.TestCase):
    def setUp(self):
        self._temporary_directory = tempfile.TemporaryDirectory(prefix="삼원-검사-")
        self.addCleanup(self._temporary_directory.cleanup)
        self.root = Path(self._temporary_directory.name)

    def change_directory(self, path):
        previous = Path.cwd()
        os.chdir(path)
        self.addCleanup(os.chdir, previous)

    def assert_not_matching(self, run_root, rerun, recorded):
        self.assertFalse(
            _service()._matching_ternary_inspection(run_root, rerun, recorded)
        )

    def test_matching_accepts_same_bounded_files_as_relative_and_absolute(self):
        self.change_directory(self.root)
        run_root = (
            self.root
            / ".localdata"
            / "expert-closure-20260930"
            / "priority-batch"
            / "한국-W1"
        )
        prediction, confidence = _files(run_root)
        relative_prediction = prediction.relative_to(self.root)
        relative_confidence = confidence.relative_to(self.root)

        for rerun_relative in (False, True):
            with self.subTest(rerun_relative=rerun_relative):
                if rerun_relative:
                    rerun = _inspection(relative_prediction, relative_confidence)
                    recorded = _inspection(prediction, confidence)
                else:
                    rerun = _inspection(prediction, confidence)
                    recorded = _inspection(relative_prediction, relative_confidence)

                rerun_before = dict(rerun)
                recorded_before = dict(recorded)
                self.assertTrue(
                    _service()._matching_ternary_inspection(
                        run_root,
                        rerun,
                        recorded,
                    )
                )
                self.assertEqual(rerun, rerun_before)
                self.assertEqual(recorded, recorded_before)

    def test_rejects_different_prediction_file_with_same_filename(self):
        run_root = self.root / "run"
        prediction, confidence = _files(run_root)
        other_prediction = run_root / "other" / prediction.name
        other_prediction.parent.mkdir()
        other_prediction.write_text("different", encoding="utf-8")

        rerun = _inspection(prediction, confidence)
        recorded = _inspection(other_prediction, confidence)

        self.assert_not_matching(run_root, rerun, recorded)

    def test_rejects_different_confidence_file_with_same_filename(self):
        run_root = self.root / "run"
        prediction, confidence = _files(run_root)
        other_confidence = run_root / "other" / confidence.name
        other_confidence.parent.mkdir()
        other_confidence.write_text("{}", encoding="utf-8")

        rerun = _inspection(prediction, confidence)
        recorded = _inspection(prediction, other_confidence)

        self.assert_not_matching(run_root, rerun, recorded)

    def test_rejects_traversal_even_when_metrics_match(self):
        self.change_directory(self.root)
        run_root = self.root / "run"
        prediction, confidence = _files(run_root)
        outside = self.root / "outside" / prediction.name
        outside.parent.mkdir()
        outside.write_text("outside", encoding="utf-8")

        rerun = _inspection(prediction, confidence)
        recorded = _inspection(
            run_root / ".." / "outside" / prediction.name,
            confidence,
        )

        self.assert_not_matching(run_root, rerun, recorded)

    def test_rejects_absolute_file_outside_run_root(self):
        run_root = self.root / "run"
        prediction, confidence = _files(run_root)
        outside_directory = self.root / "외부"
        outside_directory.mkdir()
        outside_prediction = outside_directory / "model.cif"
        outside_prediction.write_text("outside", encoding="utf-8")

        rerun = _inspection(prediction, confidence)
        recorded = _inspection(outside_prediction, confidence)

        self.assert_not_matching(run_root, rerun, recorded)

    def test_rejects_changed_metric_mapping_or_confidence(self):
        run_root = self.root / "run"
        prediction, confidence_file = _files(run_root)
        baseline = _inspection(prediction, confidence_file)

        changed_metric = _inspection(prediction, confidence_file, metric=0.76)
        changed_mapping = _inspection(prediction, confidence_file)
        changed_mapping["mapping"] = {"chain_a": "A", "chain_b": "C"}
        changed_confidence = _inspection(prediction, confidence_file)
        changed_confidence["confidence"] = "medium"

        for label, recorded in (
            ("metric", changed_metric),
            ("mapping", changed_mapping),
            ("confidence", changed_confidence),
        ):
            with self.subTest(field=label):
                self.assert_not_matching(run_root, baseline, recorded)

    def test_rejects_arbitrary_extensions(self):
        for prediction_name, confidence_name in (
            ("model.exe", "confidence_model.json"),
            ("model.cif", "confidence_model.txt"),
        ):
            with self.subTest(
                prediction_name=prediction_name,
                confidence_name=confidence_name,
            ):
                run_root = self.root / f"run-{prediction_name}-{confidence_name}"
                output = run_root / "seed-23" / "boltz_output"
                output.mkdir(parents=True)
                prediction = output / prediction_name
                confidence = output / confidence_name
                prediction.write_text("prediction", encoding="utf-8")
                confidence.write_text("{}", encoding="utf-8")
                rerun = _inspection(prediction, confidence)
                recorded = _inspection(prediction, confidence)

                self.assert_not_matching(run_root, rerun, recorded)


if __name__ == "__main__":
    unittest.main()
