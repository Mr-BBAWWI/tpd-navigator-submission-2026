import math
import unittest

from packages.science import boltz_worker as worker
from scripts import run_calibration_seed_holdout as protocol


class CalibrationSeedHoldoutTests(unittest.TestCase):
    @staticmethod
    def models(iplddt, confidence=None):
        if confidence is None:
            confidence = [0.1 + index / 10 for index in range(5)]
        return [
            {
                "model_index": index,
                "confidence_raw": {
                    "complex_iplddt": iplddt[index],
                    "confidence_score": confidence[index],
                    "iplddt": 0.99 - index / 100,
                },
            }
            for index in range(5)
        ]

    def test_selector_uses_complex_iplddt_not_iplddt_or_confidence_score(self):
        models = self.models(
            [0.20, 0.30, 0.95, 0.40, 0.50],
            [0.99, 0.98, 0.10, 0.97, 0.96],
        )
        self.assertEqual(protocol.select_model(models), 2)
        self.assertEqual(protocol.select_confidence_baseline(models), 0)

    def test_selector_rejects_nan(self):
        models = self.models([0.2, 0.3, math.nan, 0.4, 0.5])
        with self.assertRaisesRegex(worker.BoltzWorkerError,
                                    "FINITE_RAW_COMPLEX_IPLDDT_REQUIRED"):
            protocol.select_model(models)

    def test_selector_tie_uses_lower_model_index(self):
        models = self.models([0.2, 0.91, 0.4, 0.91, 0.5])
        self.assertEqual(protocol.select_model(models), 1)

    def test_selector_rejects_missing_model_index(self):
        models = self.models([0.2, 0.3, 0.4, 0.5, 0.6])
        del models[3]["model_index"]
        with self.assertRaisesRegex(worker.BoltzWorkerError, "MODEL_INDEX_REQUIRED"):
            protocol.select_model(models)

    def test_missing_complex_iplddt_does_not_fall_back(self):
        models = self.models([0.2, 0.3, 0.4, 0.5, 0.6])
        del models[4]["confidence_raw"]["complex_iplddt"]
        self.assertIn("confidence_score", models[4]["confidence_raw"])
        self.assertIn("iplddt", models[4]["confidence_raw"])
        with self.assertRaisesRegex(worker.BoltzWorkerError,
                                    "RAW_COMPLEX_IPLDDT_REQUIRED"):
            protocol.select_model(models)

    def test_heldout_seed_set_is_exact_and_disjoint(self):
        self.assertEqual(protocol.SEEDS, [101, 127, 149, 173, 197])
        self.assertEqual(protocol.PRIOR_SEEDS, [23, 41, 61, 79, 97])
        self.assertTrue(set(protocol.SEEDS).isdisjoint(protocol.PRIOR_SEEDS))
        self.assertEqual(
            {(seed, index) for seed in protocol.SEEDS
             for index in protocol.MODEL_INDICES},
            {(seed, index) for seed in [101, 127, 149, 173, 197]
             for index in range(5)},
        )

    def test_failed_seed_cannot_pass_even_with_three_marked_passes(self):
        receipts = [
            {"seed": seed, "status": "success", "selected_model_pass": index < 3}
            for index, seed in enumerate(protocol.SEEDS)
        ]
        self.assertTrue(protocol.protocol_pass(receipts))
        receipts[4]["status"] = "failure"
        receipts[4]["selected_model_pass"] = True
        self.assertFalse(protocol.protocol_pass(receipts))

    def test_missing_seed_cannot_be_replaced_by_duplicate(self):
        receipts = [
            {"seed": seed, "status": "success", "selected_model_pass": True}
            for seed in protocol.SEEDS
        ]
        receipts[-1]["seed"] = protocol.SEEDS[0]
        self.assertFalse(protocol.protocol_pass(receipts))


if __name__ == "__main__":
    unittest.main()
