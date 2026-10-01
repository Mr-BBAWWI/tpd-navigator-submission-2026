"""Synthetic regression fixtures for the bundled-followup population guard."""
from __future__ import annotations

import copy
import unittest
from unittest import mock

from packages.science import scientific_assessment


def criterion(value, identifier):
    return next(row for row in value["criteria"] if row["id"] == identifier)


def fixtures(count=1):
    smiles = ["CC", "CCC", "CCCC"]
    ids = [f"synthetic-A-{index}" for index in range(count)]
    result = {
        "sites": {"atoms": [{"atom_map": 19, "state": "MODIFIABLE"}]},
        "analogs": [{
            "id": identifier, "mapped_smiles": smiles[index],
            "canonical_smiles": smiles[index], "attachment_site_atom_maps": [19],
            "transformation_class": "ring_expansion", "cheap_filter": {"valid": True},
            "selected": True, "pipeline_status": "qualified",
            "assembly_eligible": True, "parent_redocking_supported": True,
            "docking": {"status": "completed_with_limits", "pose_preserved": True,
                        "passing_pose_count": 1},
        } for index, identifier in enumerate(ids)],
        "protac_candidates": [], "calibration": {"CRBN": {"seed_receipts": []}},
    }
    policy = {
        "_validated_by_platform": True,
        "numerical_criteria": {"panel_min": 1, "panel_max": 20},
        "scope": {"project_id": "synthetic-project", "job_id": "synthetic-job"},
        "decisions": {"microstate_expert_decision": {
            "accepted": True, "state_choices": {identifier: 0 for identifier in ids},
            "pH_conditions": {"pH": 7.4, "description": "synthetic fixture"},
        }},
    }
    digest = scientific_assessment.evaluate_design(result, policy=policy)["policy_digest"]
    reports = [{
        "_validated_by_platform": True, "job_id": "synthetic-job",
        "current_job_id": "synthetic-job", "policy_digest": digest,
        "bound_result_sha256": "1" * 64, "analog_id": identifier,
        "report_id": "interaction-" + identifier,
        "requirements": [{"id": "contact", "required": True,
                          "status": "preserved", "computed_pass": True}],
        "failures": [],
        "state_alternatives": [{
            "status": "computed", "mapped_smiles": smiles[index], "pH_context": 7.4,
            "requirements": [{"id": "contact", "required": True, "observed": True,
                              "unambiguous": True,
                              "pending_missing_protein_hydrogen": False}],
        }],
    } for index, identifier in enumerate(ids)]
    return result, policy, reports


def population_row(identifier, smiles):
    return {
        "_validated_by_platform": True, "project_id": "synthetic-project",
        "job_id": "synthetic-job", "analog_id": identifier,
        "report_id": "interaction-" + identifier, "state_index": 0, "pH": 7.4,
        "analog_canonical_isomeric_graph": smiles,
        "state_canonical_isomeric_graph": smiles,
        "bound_result_sha256": "1" * 64,
        "source_interaction_ref": {"artifact_id": "synthetic-interaction"},
        "source_interaction_sha256": "2" * 64,
        "quantity_kind": "state_population", "value": 0.0,
        "source": "synthetic:quantitative-record", "source_kind": "computed",
        "method": "synthetic method", "protocol": "synthetic protocol",
        "conditions": "synthetic pH 7.4 conditions",
        "uncertainty_description": "synthetic uncertainty",
        "scientist_interpretation": "Synthetic interpretation selecting state zero",
        "evidence_ref": {"artifact_id": "synthetic-source"},
        "source_locator": "/records/0",
        "review_ref": {"artifact_id": "synthetic-review"},
    }


class MicrostatePopulationGuardTests(unittest.TestCase):
    def evaluate(self, result, policy, reports, active=True):
        followup = {"source": {"synthetic": True}, "source_policy_digest": "fixture",
                    "manifest": {}, "rules": {}}
        with mock.patch.object(scientific_assessment, "_trusted_followup",
                               return_value=followup if active else None):
            return scientific_assessment.evaluate_design(
                result, policy=policy, supplements={"interactions": reports})

    def rebind(self, result, policy, reports):
        digest = scientific_assessment.evaluate_design(result, policy=policy)["policy_digest"]
        for report in reports:
            report["policy_digest"] = digest

    def test_missing_quantitative_evidence_downgrades_both_but_nofollowup_is_unchanged(self):
        result, policy, reports = fixtures()
        guarded = self.evaluate(result, policy, reports)
        self.assertEqual(criterion(guarded, "core_interaction_preservation")["status"], "pending")
        self.assertEqual(criterion(guarded, "microstates_h_direction")["status"], "pending")
        unguarded = self.evaluate(result, policy, reports, active=False)
        self.assertEqual(criterion(unguarded, "core_interaction_preservation")["status"], "pass")
        self.assertEqual(criterion(unguarded, "microstates_h_direction")["status"], "pass")

    def test_fully_bound_zero_population_record_allows_existing_checks(self):
        result, policy, reports = fixtures()
        policy["decisions"]["microstate_expert_decision"]["population_evidence"] = [
            population_row("synthetic-A-0", "CC")]
        self.rebind(result, policy, reports)
        value = self.evaluate(result, policy, reports)
        self.assertEqual(criterion(value, "core_interaction_preservation")["status"], "pass")
        self.assertEqual(criterion(value, "microstates_h_direction")["status"], "pass")
        coverage = criterion(value, "microstates_h_direction")["observed"]["population_evidence_coverage"]
        self.assertEqual(coverage[0]["value"], 0.0)
        self.assertIs(value["scientific_accepted"], False)

    def test_missing_choice_bad_index_empty_ph_and_incomplete_coverage_remain_pending(self):
        mutations = (
            lambda p, r: p["decisions"]["microstate_expert_decision"].update(state_choices={}),
            lambda p, r: p["decisions"]["microstate_expert_decision"].update(
                state_choices={"synthetic-A-0": 9}),
            lambda p, r: r[0]["state_alternatives"][0].pop("pH_context"),
        )
        for mutate in mutations:
            result, policy, reports = fixtures()
            mutate(policy, reports); self.rebind(result, policy, reports)
            value = self.evaluate(result, policy, reports)
            self.assertEqual(criterion(value, "core_interaction_preservation")["status"], "pending")
        result, policy, reports = fixtures(3)
        policy["decisions"]["microstate_expert_decision"]["population_evidence"] = [
            population_row("synthetic-A-0", "CC"), population_row("synthetic-A-1", "CCC")]
        self.rebind(result, policy, reports)
        value = self.evaluate(result, policy, reports)
        self.assertEqual(criterion(value, "core_interaction_preservation")["status"], "pending")

    def test_malformed_reused_and_failed_records_do_not_pass(self):
        changes = (("state_index", 1), ("pH", 6.5),
                   ("state_canonical_isomeric_graph", "CCC"),
                   ("bound_result_sha256", "9" * 64),
                   ("project_id", "other"), ("value", float("nan")), ("value", 1.1))
        for field, changed in changes:
            result, policy, reports = fixtures()
            row = population_row("synthetic-A-0", "CC"); row[field] = changed
            policy["decisions"]["microstate_expert_decision"]["population_evidence"] = [row]
            self.rebind(result, policy, reports)
            self.assertEqual(criterion(self.evaluate(result, policy, reports),
                                       "microstates_h_direction")["status"], "pending")
        result, policy, reports = fixtures()
        row = population_row("synthetic-A-0", "CC")
        policy["decisions"]["microstate_expert_decision"]["population_evidence"] = [row, copy.deepcopy(row)]
        self.rebind(result, policy, reports)
        self.assertEqual(criterion(self.evaluate(result, policy, reports),
                                   "microstates_h_direction")["status"], "pending")
        reports[0]["failures"] = [{"reason": "synthetic geometry failure"}]
        self.assertEqual(criterion(self.evaluate(result, policy, reports),
                                   "core_interaction_preservation")["status"], "failed")


if __name__ == "__main__":
    unittest.main()
