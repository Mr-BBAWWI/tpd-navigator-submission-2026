import copy
import unittest

from packages.science.e3_benchmark import compare_paired


class E3BenchmarkTests(unittest.TestCase):
    def _coordinates(self):
        xyz = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        reference = {
            role: {"ids": ["1", "2", "3"], "xyz": copy.deepcopy(xyz)}
            for role in ("target", "e3", "ligand")
        }
        return reference, copy.deepcopy(reference)

    def _provenance(self):
        return {
            "model": "actual-test-checkpoint",
            "seed": 0,
            "input_sha256": "a" * 64,
            "prediction_sha256": "b" * 64,
            "e3_type": "CRBN",
        }

    def test_seed_zero_is_valid_and_actual_metadata_is_preserved(self):
        reference, prediction = self._coordinates()
        prediction["ligand"]["xyz"] = [[1, 0, 0], [2, 0, 0], [1, 1, 0]]
        provenance = self._provenance()
        provenance["applied_model_metadata"] = {
            "checkpoint": "actual-test-checkpoint",
            "seed": 0,
        }

        result = compare_paired(reference, prediction, provenance)

        self.assertEqual(result["provenance"], provenance)
        self.assertEqual(result["provenance"]["seed"], 0)
        self.assertEqual(result["provenance"]["applied_model_metadata"]["seed"], 0)
        self.assertEqual(result["evaluation_kind"], "prediction_reference_comparison")
        self.assertAlmostEqual(
            result["metrics"]["ligand_RMSD_after_target_alignment_A"], 1.0
        )

    def test_empty_null_and_invalid_provenance_are_rejected(self):
        reference, prediction = self._coordinates()
        invalid_updates = [
            {"model": ""},
            {"model": None},
            {"input_sha256": ""},
            {"input_sha256": None},
            {"prediction_sha256": ""},
            {"prediction_sha256": None},
            {"e3_type": ""},
            {"e3_type": None},
            {"seed": None},
            {"seed": True},
            {"seed": -1},
        ]
        for update in invalid_updates:
            with self.subTest(update=update):
                provenance = self._provenance()
                provenance.update(update)
                with self.assertRaisesRegex(ValueError, "BENCHMARK_PROVENANCE_REQUIRED"):
                    compare_paired(reference, prediction, provenance)

    def test_only_supported_e3_types_are_accepted(self):
        reference, prediction = self._coordinates()
        provenance = self._provenance()
        provenance["e3_type"] = "DDB1"
        with self.assertRaisesRegex(ValueError, "BENCHMARK_E3"):
            compare_paired(reference, prediction, provenance)

        provenance["e3_type"] = "VHL"
        self.assertEqual(
            compare_paired(reference, prediction, provenance)["e3_type"], "VHL"
        )


if __name__ == "__main__":
    unittest.main()
