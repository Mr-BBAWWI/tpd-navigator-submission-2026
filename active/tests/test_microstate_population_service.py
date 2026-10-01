"""Service integration and absent-protein-H tests using synthetic fixtures only."""
from __future__ import annotations

import copy
import json
import unittest

from packages.contracts import ContractError
from packages.platform.review_identity import ReviewAuthError
from packages.science import interaction_review

try:
    import test_scientific_acceptance_service as scientific_fixtures
except ImportError:
    from tests import test_scientific_acceptance_service as scientific_fixtures


class MicrostatePopulationServiceTests(unittest.TestCase):
    project = scientific_fixtures.ScientificAcceptanceServiceTests.project
    setUp = scientific_fixtures.ScientificAcceptanceServiceTests.setUp
    _completed_job = scientific_fixtures.ScientificAcceptanceServiceTests._completed_job
    _replace_job = scientific_fixtures.ScientificAcceptanceServiceTests._replace_job

    def snapshot(self):
        with self.store.db() as db:
            _, result, binding, _ = self.service._verified_result(db, self.job_id)
            policy = self.service._policy(db, self.job_id, result)
        return result, binding, policy

    def report(self, marker="source"):
        return {
            "format": "synthetic-interaction-fixture", "fixture_marker": marker,
            "report_id": "interaction-A-1", "analog_id": "A-1",
            "requirements": [{"id": "synthetic-contact", "required": True,
                              "status": "preserved", "computed_pass": True}],
            "failures": [], "pH_conditions": {"pH": 7.4},
            "state_alternatives": [{
                "index": 0, "status": "computed", "mapped_smiles": "[CH3:1][CH3:2]",
                "pH_context": 7.4,
                "requirements": [{"id": "synthetic-contact", "required": True,
                                  "observed": True, "unambiguous": True,
                                  "pending_missing_protein_hydrogen": False}],
            }], "scientific_approved": False,
        }

    def register_interaction(self, report, policy):
        result, binding, _ = self.snapshot(); paths = []
        evidence_binding = {
            "project": self.project, "job_id": self.job_id,
            "result_sha256": binding["result_sha256"],
            "evaluation_policy_digest": self.service._policy_digest(policy),
            "evaluation_policy": {"synthetic_fixture": True},
        }
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            made = self.service._register_evidence(
                db, paths, self.job_id, "interaction", "A-1",
                evidence_binding, report)
        return made["ref"]

    def source_document(self, result, interaction_ref, **changes):
        record = {
            "project_id": self.project, "job_id": self.job_id,
            "result_sha256": result["files"]["json"]["sha256"],
            "interaction_ref": interaction_ref, "analog_id": "A-1",
            "selected_state_index": 0, "pH": 7.4,
            "analog_canonical_isomeric_graph": "CC",
            "state_canonical_isomeric_graph": "CC",
            "source_identifier": "synthetic:population", "source_kind": "computed",
            "method": "synthetic calculation", "protocol": "synthetic protocol",
            "conditions": "synthetic pH 7.4 conditions",
            "uncertainty_description": "synthetic uncertainty",
            "interpretation": "Synthetic record identifies selected state zero",
            "state_population": 0.0,
        }
        record.update(changes)
        return {"records": [record], "fixture_notice": "synthetic, not scientific evidence"}

    def register_source(self, result, interaction_ref, **changes):
        text = json.dumps(
            self.source_document(result, interaction_ref, **changes),
            allow_nan=False)
        return self.service.register_quantitative_source(
            self.job_id, text, self.token)

    def submission(self, assessment, source_ref, review_ref, rows=None):
        row = {"analog_id": "A-1", "report_id": "interaction-A-1",
               "state_index": 0, "pH": 7.4, "evidence_ref": source_ref,
               "locator": "/records/0", "review_ref": review_ref,
               "interpretation": "Synthetic reviewer interpretation selecting state zero"}
        return {
            "kind": "microstate_decision", "action": "accept",
            "intent": "Synthetic diagnostic state selection",
            "reason": "Synthetic fixture; not scientific evidence",
            "expected_policy_revision": assessment["policy_revision"],
            "data": {"state_choices": {"A-1": 0},
                     "pH_conditions": {"pH": 7.4, "description": "synthetic"},
                     "population_evidence": rows if rows is not None else [row]},
        }

    def test_old_policy_ignored_then_current_recompute_passes_without_policy_digest_churn(self):
        result, _, policy = self.snapshot()
        source_interaction = self.register_interaction(self.report(), policy)
        source_ref = self.register_source(
            result, source_interaction)["ref"]
        review = self.service.create_statement(
            "Synthetic reviewer statement; not quantitative evidence.",
            "fixture:microstate-review", self.token)
        assessment, _ = self.service.create(self.job_id)
        self.service.record_decision(
            assessment["id"], self.submission(assessment, source_ref, review["ref"]),
            self.token)
        pending, _ = self.service.create(self.job_id)
        rows = {row["id"]: row for row in pending["criteria"]}
        self.assertEqual(rows["core_interaction_preservation"]["status"], "pending")
        _, _, current_policy = self.snapshot()
        digest_before = self.service._policy_digest(current_policy)
        self.register_interaction(self.report("current-policy-recompute"), current_policy)
        _, _, reconstructed = self.snapshot()
        self.assertEqual(self.service._policy_digest(reconstructed), digest_before)
        current, _ = self.service.create(self.job_id)
        rows = {row["id"]: row for row in current["criteria"]}
        self.assertEqual(rows["core_interaction_preservation"]["status"], "pass")
        self.assertEqual(rows["microstates_h_direction"]["status"], "pass")
        self.assertIs(current["scientific_accepted"], False)

    def test_canonical_source_reupload_keeps_active_policy_and_original_provenance(self):
        result, _, policy = self.snapshot()
        source_interaction = self.register_interaction(self.report(), policy)
        document = self.source_document(result, source_interaction)
        pretty = json.dumps(document, allow_nan=False, indent=2)
        compact = json.dumps(
            document, allow_nan=False, separators=(",", ":"))
        first = self.service.register_quantitative_source(
            self.job_id, pretty, self.token)
        review = self.service.create_statement(
            "Synthetic reviewer statement; not quantitative evidence.",
            "fixture:canonical-source-review", self.token)
        assessment, _ = self.service.create(self.job_id)
        self.service.record_decision(
            assessment["id"],
            self.submission(assessment, first["ref"], review["ref"]),
            self.token)

        _, _, active_policy = self.snapshot()
        digest_before = self.service._policy_digest(active_policy)
        second = self.service.register_quantitative_source(
            self.job_id, compact, self.token)
        self.assertFalse(second["created"])
        self.assertEqual(second["id"], first["id"])
        self.assertEqual(second["ref"], first["ref"])
        self.assertEqual(
            second["source_upload_sha256"], first["source_upload_sha256"])
        self.assertEqual(
            second["source_document_sha256"],
            first["source_document_sha256"])
        with self.store.db() as db:
            source_count = db.execute("""SELECT COUNT(*) AS n
                FROM scientific_evidence WHERE project=? AND job_id=?
                AND kind='microstate_population_source'""",
                (self.project, self.job_id)).fetchone()["n"]
        self.assertEqual(source_count, 1)
        _, _, reconstructed = self.snapshot()
        self.assertEqual(
            self.service._policy_digest(reconstructed), digest_before)

        self.register_interaction(
            self.report("canonical-source-current-policy"), reconstructed)
        _, _, recomputed = self.snapshot()
        self.assertEqual(self.service._policy_digest(recomputed), digest_before)
        current, _ = self.service.create(self.job_id)
        rows = {row["id"]: row for row in current["criteria"]}
        self.assertEqual(rows["core_interaction_preservation"]["status"], "pass")
        self.assertEqual(rows["microstates_h_direction"]["status"], "pass")
        self.assertIs(current["scientific_accepted"], False)

    def test_bad_source_scope_state_quantity_duplicates_and_review_separation_rejected(self):
        result, _, policy = self.snapshot()
        interaction_ref = self.register_interaction(self.report(), policy)
        review = self.service.create_statement(
            "Synthetic reviewer statement.", "fixture:review", self.token)
        assessment, _ = self.service.create(self.job_id)
        mutations = ({"project_id": "other"}, {"job_id": "other"},
                     {"result_sha256": "0" * 64}, {"pH": 6.5},
                     {"state_canonical_isomeric_graph": "CCC"},
                     {"selected_state_index": True}, {"state_population": 1.1})
        for mutation in mutations:
            raw = json.dumps(self.source_document(result, interaction_ref, **mutation),
                             allow_nan=False).encode()
            ref = self.service.port.put_raw(raw, "application/json", "source")
            with self.assertRaises(ContractError):
                self.service.record_decision(
                    assessment["id"], self.submission(assessment, ref, review["ref"]),
                    self.token)
        nan_raw = json.dumps(self.source_document(
            result, interaction_ref, state_population=float("nan")),
            allow_nan=True).encode()
        nan_ref = self.service.port.put_raw(nan_raw, "application/json", "source")
        with self.assertRaises(ContractError):
            self.service.record_decision(
                assessment["id"], self.submission(assessment, nan_ref, review["ref"]),
                self.token)
        good_ref = self.register_source(result, interaction_ref)["ref"]
        same = self.submission(assessment, good_ref, good_ref)
        with self.assertRaises(ContractError):
            self.service.record_decision(assessment["id"], same, self.token)
        duplicate = self.submission(assessment, good_ref, review["ref"])
        duplicate["data"]["population_evidence"].append(
            copy.deepcopy(duplicate["data"]["population_evidence"][0]))
        with self.assertRaises(ContractError):
            self.service.record_decision(assessment["id"], duplicate, self.token)

    def test_authenticated_registration_lists_source_and_only_own_statement(self):
        result, _, policy = self.snapshot()
        interaction_ref = self.register_interaction(self.report(), policy)
        statement = self.service.create_statement(
            "Synthetic expert review statement, separate from raw quantity.",
            "fixture:population-review", self.token)
        made = self.register_source(result, interaction_ref)
        listed = self.service.list_quantitative_sources(
            self.job_id, self.token)
        self.assertEqual(listed["sources"][0]["ref"], made["ref"])
        self.assertEqual(listed["review_statements"][0]["ref"], statement["ref"])
        self.assertEqual(listed["interaction_reports"][0]["ref"], interaction_ref)
        self.assertEqual(
            listed["sources"][0]["state"],
            "unverified_quantity_pending_expert_review")
        assessment, _ = self.service.create(self.job_id)
        self.assertFalse(assessment["scientific_accepted"])

    def test_registration_rejects_wrong_scope_human_source_and_malformed_json(self):
        result, _, policy = self.snapshot()
        interaction_ref = self.register_interaction(self.report(), policy)
        for changes in ({"job_id": "other-job"},
                        {"project_id": "other-project"},
                        {"selected_state_index": True}):
            with self.subTest(changes=changes), self.assertRaises(ContractError):
                self.register_source(result, interaction_ref, **changes)
        with self.assertRaises(ContractError):
            self.service.register_quantitative_source(
                self.job_id, '{"type":"human_statement","records":[]}',
                self.token)
        with self.assertRaises(ContractError):
            self.service.register_quantitative_source(
                self.job_id, '{"records":[],"records":[]}', self.token)
        with self.assertRaises(ReviewAuthError):
            self.service.list_quantitative_sources(self.job_id, "invalid-session")

    def test_unregistered_quantitative_artifact_cannot_support_decision(self):
        result, _, policy = self.snapshot()
        interaction_ref = self.register_interaction(self.report(), policy)
        raw_ref = self.service.port.put_raw(
            json.dumps(self.source_document(result, interaction_ref)).encode(),
            "application/json", "source")
        review = self.service.create_statement(
            "Synthetic reviewer statement.", "fixture:review", self.token)
        assessment, _ = self.service.create(self.job_id)
        with self.assertRaises(ContractError):
            self.service.record_decision(
                assessment["id"],
                self.submission(assessment, raw_ref, review["ref"]), self.token)

    def test_strict_parser_rejects_noncanonical_or_unsafe_json_numbers(self):
        invalid = (
            b'{"records":[],"records":[]}',
            b'{"records":[{"value":1,"value":2}]}',
            b'{"value":NaN}', b'{"value":Infinity}', b'{"value":-Infinity}',
            b'{"value":1e10000}', b'[]', b'{"value":1} trailing',
            b'\xff{"value":1}',
            b'{"value":' + (b'9' * 10000) + b'}',
        )
        for raw in invalid:
            with self.subTest(raw=raw[:80]):
                with self.assertRaises(ContractError):
                    self.service._strict_quantitative_document(raw)

    def test_direction_specific_absent_protein_h_receipt_uses_engine_record_format(self):
        requirement = {
            "id": "hbond", "kind": "directional_hbond", "required": True,
            "ligand_maps": [19],
            "protein_atom_ids": ["A:42:ASP:OD1"],
        }
        actual_match = {
            "id": "identifier-is-sorted-separately",
            "kind": "directional_hbond",
            "participants": ["L:19", "A:42:ASP:OD1"],
            "distance_DA_A": 2.8,
            "distance_HA_A": 1.8,
            "angle_DHA_deg": 171.0,
            "requires_review": False,
            "description": "engine-format directional geometry",
        }
        self.assertTrue(interaction_review._participant_match(
            actual_match, requirement))
        row = {
            **requirement,
            "status": "preserved", "computed_pass": True,
            "pending_missing_protein_hydrogen": False,
            "needs_expert": False,
            "before_matches": [copy.deepcopy(actual_match)],
            "after_matches": [copy.deepcopy(actual_match)],
        }
        base = {
            "protein_hydrogen_status": "not_required_or_not_implicated",
            "failures": [], "requirements": [row],
        }
        receipt = self.service._absent_protein_hydrogen_receipt(base)
        self.assertFalse(receipt["requires_review"])

        for field in ("before_matches", "after_matches"):
            changed = copy.deepcopy(base)
            changed["requirements"][0][field][0]["participants"] = [
                "A:42:ASP:OD1", "L:19"]
            self.assertTrue(self.service._absent_protein_hydrogen_receipt(
                changed)["requires_review"])

        wrong_ligand = copy.deepcopy(base)
        wrong_ligand["requirements"][0]["after_matches"][0]["participants"][0] = "L:20"
        self.assertTrue(self.service._absent_protein_hydrogen_receipt(
            wrong_ligand)["requires_review"])
        wrong_protein = copy.deepcopy(base)
        wrong_protein["requirements"][0]["after_matches"][0]["participants"][1] = "A:43:ASP:OD1"
        self.assertTrue(self.service._absent_protein_hydrogen_receipt(
            wrong_protein)["requires_review"])
        review_match = copy.deepcopy(base)
        review_match["requirements"][0]["after_matches"][0]["requires_review"] = True
        self.assertTrue(self.service._absent_protein_hydrogen_receipt(
            review_match)["requires_review"])
        missing_measurement = copy.deepcopy(base)
        missing_measurement["requirements"][0]["after_matches"][0].pop("distance_HA_A")
        self.assertTrue(self.service._absent_protein_hydrogen_receipt(
            missing_measurement)["requires_review"])
        expert = copy.deepcopy(base)
        expert["requirements"][0]["needs_expert"] = True
        self.assertTrue(self.service._absent_protein_hydrogen_receipt(
            expert)["requires_review"])
        for requirements in (None, {}, [], ["invalid"]):
            changed = copy.deepcopy(base)
            if requirements is None:
                changed.pop("requirements")
            else:
                changed["requirements"] = requirements
            self.assertTrue(self.service._absent_protein_hydrogen_receipt(
                changed)["requires_review"])


if __name__ == "__main__":
    unittest.main()
