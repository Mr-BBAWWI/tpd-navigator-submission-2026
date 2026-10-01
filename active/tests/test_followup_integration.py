from __future__ import annotations

import copy
import unittest

from packages.science import scientific_assessment
from packages.science.expert_followup import EXPECTED_SEEDS, load_followup


LIMITS = {
    "target_CA_RMSD_A": 1.5,
    "ligand_heavy_atom_RMSD_after_target_alignment_A": 3.0,
    "e3_CA_RMSD_after_target_alignment_A": 10.0,
    "contact_jaccard": .4,
}


def result():
    return {
        "parent_scope": {"actual_design_parent_id": "SMARCA2-FX5"},
        "sites": {"atoms": [{"atom_map": 19, "state": "MODIFIABLE"}]},
        "analogs": [],
        "protac_candidates": [],
        "summary": {},
    }


def policy(with_followup=True):
    value = {
        "_validated_by_platform": True,
        "numerical_criteria": {},
        "scope": {"project_id": "p", "job_id": "job-1"},
        "decisions": {},
        "opinion": {"scope": {"fx5_bound_scope": True}},
    }
    if with_followup:
        value["followup"] = load_followup()
    return value


def calibration_report(policy_value, passing):
    digest = scientific_assessment._digest(policy_value)
    rows = []
    for seed in EXPECTED_SEEDS:
        metrics = dict(LIMITS)
        if seed not in passing:
            metrics["target_CA_RMSD_A"] = 2.0
        rows.append({
            "seed": seed,
            "technical_execution_success": True,
            "metrics": metrics,
            "original_status": "failed" if seed == 41 else "completed",
            "original_failure_reason": "retained" if seed == 41 else None,
        })
    return {
        "id": "CRBN_6BOY",
        "kind": "source_independent_benchmark",
        "distribution": {"seeds": rows},
        "measurement_binding": {
            "source_id": "CRBN_6BOY",
            "module_version": "test",
            "module_fingerprint": "test",
        },
        "_validated_by_platform": True,
        "job_id": "job-1",
        "current_job_id": "job-1",
        "policy_digest": digest,
    }


def criterion(assessment, identifier):
    return next(item for item in assessment["criteria"] if item["id"] == identifier)


class FollowupIntegrationTests(unittest.TestCase):
    def test_untrusted_or_absent_followup_is_ignored(self):
        trusted_without = scientific_assessment.evaluate_design(result(), policy=policy(False))
        self.assertIsNone(trusted_without["expert_followup"])
        forged = policy(True)
        forged["_validated_by_platform"] = False
        untrusted = scientific_assessment.evaluate_design(result(), policy=forged)
        self.assertIsNone(untrusted["expert_followup"])

    def test_two_of_five_fails_and_three_of_five_passes_calibration_only(self):
        for passing, expected in (({23, 61}, "failed"), ({23, 61, 79}, "pass")):
            bound_policy = policy()
            report = calibration_report(bound_policy, passing)
            assessment = scientific_assessment.evaluate_design(
                result(), policy=bound_policy,
                supplements={"calibration": [report]},
            )
            observed = criterion(assessment, "known_crbn_calibration")
            self.assertEqual(observed["status"], expected)
            self.assertEqual(observed["observed"]["expert_followup"]["passing_seed_count"],
                             len(passing))
            self.assertEqual(criterion(assessment, "novel_ternary_geometry")["status"],
                             "pending")
            self.assertEqual(criterion(assessment, "formal_expert_decision")["status"],
                             "pending")
            self.assertEqual(assessment["demo_scope"], "DEMO/INCOMPLETE")
            self.assertFalse(assessment["scientific_accepted"])

    def test_old_policy_digest_report_requires_reimport(self):
        bound_policy = policy()
        report = calibration_report(bound_policy, {23, 61, 79})
        report["policy_digest"] = "0" * 64
        assessment = scientific_assessment.evaluate_design(
            result(), policy=bound_policy,
            supplements={"calibration": [report]},
        )
        self.assertNotEqual(criterion(assessment, "known_crbn_calibration")["status"],
                            "pass")
        self.assertIsNone(assessment["expert_followup"]["quantitative_subchecks"]
                          ["calibration"])

    def test_direct_authenticated_calibration_decision_keeps_precedence(self):
        bound_policy = policy()
        bound_policy["decisions"]["calibration_criterion"] = {
            "expected_seeds": list(EXPECTED_SEEDS),
            "metric": "target_CA_RMSD_A",
            "threshold": 0.1,
            "comparison": "max_lte",
            "scope": "authenticated direct test",
        }
        report = calibration_report(bound_policy, set(EXPECTED_SEEDS))
        assessment = scientific_assessment.evaluate_design(
            result(), policy=bound_policy,
            supplements={"calibration": [report]},
        )
        self.assertEqual(criterion(assessment, "known_crbn_calibration")["status"],
                         "failed")
        self.assertIsNone(assessment["expert_followup"]["quantitative_subchecks"]
                          ["calibration"])

    def test_source_provenance_and_engine_fingerprint_change_inputs(self):
        bound_policy = policy()
        assessment = scientific_assessment.evaluate_design(result(), policy=bound_policy)
        followup = assessment["expert_followup"]
        self.assertEqual(followup["source"]["sha256"],
                         bound_policy["followup"]["source"]["sha256"])
        self.assertEqual(len(assessment["engine_fingerprint"]), 64)
        mutated = copy.deepcopy(bound_policy)
        mutated["followup"]["source"]["sha256"] = "f" * 64
        ignored = scientific_assessment.evaluate_design(result(), policy=mutated)
        self.assertIsNone(ignored["expert_followup"])


if __name__ == "__main__":
    unittest.main()
