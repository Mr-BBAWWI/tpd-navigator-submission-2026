import copy
import unittest

from packages.science import scientific_assessment


class CalibrationSupplementAssessmentTests(unittest.TestCase):
    def _policy(self):
        return {
            "_validated_by_platform": True,
            "numerical_criteria": {},
            "decisions": {},
            "scope": {"job_id": "job-1"},
        }

    def _bound_report(self, policy):
        digest = scientific_assessment._digest(policy)
        rows = []
        values = {
            23: 1.5420394953522754,
            41: 12.38776540052114,
            61: 1.1780957439827384,
            79: 7.4289082360252525,
            97: 2.505346661961516,
        }
        for seed, value in values.items():
            rows.append({
                "seed": seed,
                "technical_execution_success": True,
                "execution_status": "completed_exit0",
                "inspection_status": "success",
                "original_status": "failure" if seed == 23 else "success",
                "original_failure_reason": (
                    "PREDICTION_CHAIN_MAPPING_AMBIGUOUS" if seed == 23 else None),
                "metrics": {
                    "comparison.metrics.e3_CA_RMSD_after_target_alignment_A": value
                },
                "model_confidence": {"confidence_score": 0.5},
            })
        return {
            "_validated_by_platform": True,
            "job_id": "job-1",
            "current_job_id": "job-1",
            "policy_digest": digest,
            "id": "CRBN_6BOY",
            "kind": "source_independent_benchmark",
            "measurement_binding": {
                "source_id": "CRBN_6BOY",
                "module_version": "crbn-calibration-import/1",
                "module_fingerprint": "platform-provided-fingerprint",
            },
            "source": {"case": "6BOY", "e3_type": "CRBN"},
            "registered_artifacts": {"frozen-base-manifest": {"sha256": "base"}},
            "distribution": {
                "seeds": rows,
                "technical_success_count": 5,
                "scientific_model_quality_success_count": None,
            },
        }

    def test_five_real_rows_remain_pending_without_expert_threshold(self):
        policy = self._policy()
        report = self._bound_report(policy)
        status, observed = scientific_assessment._calibration_covered(
            {}, {}, 3, [report], scientific_assessment._digest(policy), policy["scope"])
        self.assertEqual(status, "pending")
        self.assertEqual(observed["technical_success_count"], 5)
        self.assertIsNone(observed["scientific_model_quality_success_count"])
        self.assertTrue(observed["known_seed_41_retained"])
        metric = observed["metric_distributions"][
            "e3_CA_RMSD_after_target_alignment_A"]
        self.assertEqual(metric["median"], 2.505346661961516)
        self.assertEqual(metric["range"], [1.1780957439827384, 12.38776540052114])
        self.assertAlmostEqual(metric["iqr"], 5.886868740672977)
        seed23 = next(row for row in observed["all_outcomes"] if row["seed"] == 23)
        self.assertEqual(seed23["original_status"], "failure")
        self.assertEqual(seed23["original_failure_reason"],
                         "PREDICTION_CHAIN_MAPPING_AMBIGUOUS")
        self.assertEqual(observed["model_confidence_by_seed"]["23"],
                         [{"confidence_score": 0.5}])
        self.assertNotIn("confidence_score", observed["metric_distributions"])

    def test_human_approved_bare_metric_covers_prefixed_sidecar_rows(self):
        policy = self._policy()
        report = self._bound_report(policy)
        decisions = {"calibration_criterion": {
            "metric": "e3_CA_RMSD_after_target_alignment_A",
            "threshold": 20.0,
            "comparison": "max_lte",
            "expected_seeds": [23, 41, 61, 79, 97],
            "scope": "synthetic unittest human-approved bound criterion",
        }}
        status, observed = scientific_assessment._calibration_covered(
            {}, decisions, 3, [report], scientific_assessment._digest(policy),
            policy["scope"])
        self.assertEqual(status, "pass")
        self.assertEqual(set(observed["metric_by_seed"]), {23, 41, 61, 79, 97})

        decisions["calibration_criterion"]["threshold"] = 10.0
        status, _ = scientific_assessment._calibration_covered(
            {}, decisions, 3, [report], scientific_assessment._digest(policy),
            policy["scope"])
        self.assertEqual(status, "failed")

    def test_unbound_policy_or_job_cannot_expand_frozen_receipts(self):
        policy = self._policy()
        report = self._bound_report(policy)
        for mutation in (
            lambda value: value.update(policy_digest="tampered"),
            lambda value: value.update(job_id="other-job"),
            lambda value: value.update(current_job_id="other-job"),
            lambda value: value.update(_validated_by_platform=False),
        ):
            candidate = copy.deepcopy(report)
            mutation(candidate)
            result = {"calibration": {"CRBN": {"seed_receipts": [
                {"seed": 23, "execution_success": True, "metrics": {"x": 1}},
                {"seed": 41, "execution_success": True, "metrics": {"x": 9}},
                {"seed": 61, "execution_success": True, "metrics": {"x": 2}},
            ]}}}
            status, observed = scientific_assessment._calibration_covered(
                result, {}, 3, [candidate], scientific_assessment._digest(policy),
                policy["scope"])
            self.assertEqual(status, "pending")
            self.assertEqual(observed["receipt_seeds"], [23, 41, 61])
            self.assertIsNone(observed["supplemental_source"])

    def test_failed_seed_is_retained_and_cannot_become_success_scalar(self):
        policy = self._policy()
        report = self._bound_report(policy)
        row = next(row for row in report["distribution"]["seeds"] if row["seed"] == 97)
        row.update({
            "technical_execution_success": False,
            "execution_success": True,
            "execution_status": "process_failure",
            "inspection_status": "failed_or_missing",
            "original_status": "failure",
            "original_failure_reason": "BOLTZ_NONZERO_EXIT",
            "metrics": None,
        })
        status, observed = scientific_assessment._calibration_covered(
            {}, {}, 3, [report], scientific_assessment._digest(policy), policy["scope"])
        self.assertEqual(status, "pending")
        self.assertEqual(observed["technical_success_count"], 4)
        failed = next(row for row in observed["all_outcomes"] if row["seed"] == 97)
        self.assertFalse(failed["execution_success"])
        self.assertEqual(failed["failure"]["original_failure_reason"],
                         "BOLTZ_NONZERO_EXIT")


if __name__ == "__main__":
    unittest.main()
