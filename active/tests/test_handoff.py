"""Read-only consumption of the real B example; no raw cache or science imports."""
import copy
import json
from pathlib import Path
import unittest

from packages.science.handoff import (HandoffError, artifact_ref, compare_observation, compare_versions,
    encoded, match_compound, parse, read_material, sha, validate_material)

ACTIVE = Path(__file__).resolve().parents[1]
EXAMPLE = ACTIVE / "outputs/b_handoff_20260923"


def reseal(material):
    material["material_digest"] = sha(encoded({k:v for k,v in material.items() if k != "material_digest"}))
    return material


class HandoffConsumerTests(unittest.TestCase):
    def setUp(self):
        self.material = parse((EXAMPLE / "material.json").read_bytes())
        self.query = {"doi": "10.1038/s41589-019-0294-6", "paper_name": "PROTAC 1", "compound_number": "2"}

    def test_real_aliases_do_not_confuse_protac_and_compound_numbers(self):
        self.assertEqual(match_compound(self.material, **self.query)["compound_id"], "C01")
        self.query.update(paper_name="PROTAC 2", compound_number="3")
        self.assertEqual(match_compound(self.material, **self.query)["compound_id"], "C02")
        self.query.update(paper_name="PROTAC 1", compound_number="1")
        self.assertEqual(match_compound(self.material, **self.query)["status"], "conflict")

    def test_another_article_and_acbi1_are_not_matched(self):
        for update in ({"doi": "10.0000/another-paper"}, {"paper_name": "ACBI1", "compound_number": "5"}):
            self.assertEqual(match_compound(self.material, **{**self.query, **update})["status"], "not_found")

    def test_missing_identity_and_wrong_number_type_rejected(self):
        with self.assertRaisesRegex(HandoffError, "IDENTITY_REQUIRED"):
            match_compound(self.material, doi=self.query["doi"])
        with self.assertRaisesRegex(HandoffError, "NUMBER_INVALID"):
            match_compound(self.material, doi=self.query["doi"], compound_number=2)

    def test_stale_material_cannot_be_silently_consumed(self):
        with self.assertRaisesRegex(HandoffError, "STALE"):
            match_compound(self.material, **self.query, expected_digest="0" * 64)

    def test_human_and_dispatch_flags_cannot_be_promoted(self):
        for key,value in (("human_review", "approved"), ("dispatch_authorized", True), ("approval_record_created", True)):
            material = copy.deepcopy(self.material)
            material["authority"][key] = value
            with self.assertRaisesRegex(HandoffError, "SCHEMA"):
                validate_material(reseal(material))

    def test_experiment_cannot_be_reassigned_to_another_compound(self):
        self.material["h1"]["observations"][0]["original"]["candidate_id"] = "C02"
        with self.assertRaisesRegex(HandoffError, "OBSERVATION_COMPOUND"):
            validate_material(reseal(self.material))

    def test_missing_experiment_cannot_gain_parent_value(self):
        self.material["h1"]["observations"][1]["original"]["value"] = 300
        with self.assertRaisesRegex(HandoffError, "MISSING_OBSERVATION"):
            validate_material(reseal(self.material))

    def test_molecule_version_and_hypothesis_subject_are_bound(self):
        material = copy.deepcopy(self.material)
        material["h1"]["compounds"][1]["molecule_version"] = "0" * 64
        with self.assertRaisesRegex(HandoffError, "MOLECULE_VERSION"):
            validate_material(reseal(material))
        self.material["h2"]["subject_molecule_id"] = self.material["h1"]["compounds"][1]["molecule_id"]
        with self.assertRaisesRegex(HandoffError, "HYPOTHESIS_SUBJECT"):
            validate_material(reseal(self.material))

    def test_changed_units_relations_and_missing_conditions_require_review(self):
        original = self.material["h1"]["observations"][0]["original"]
        for field,value in (("unit", "uM"), ("relation", "<"), ("exposure_time_h", 24),
                            ("kind", "computed_prediction"), ("locator", ".//fig[@id='F3']/caption")):
            incoming = {**original, field:value}
            result = compare_observation(self.material, self.query, incoming)
            self.assertEqual(result["status"], "difference_requires_review")
            self.assertIn(field, result["differences"])
            self.assertEqual(result["incoming"], incoming)
        result = compare_observation(self.material, self.query, original)
        self.assertEqual(result["status"], "same_record_pending_review")
        self.assertIsNone(result["incoming"]["exposure_time_h"])

    def test_new_value_for_c02_is_not_automatically_accepted(self):
        incoming = {"endpoint": "DC50", "value": 300, "unit": "nM"}
        result = compare_observation(self.material, {**self.query, "paper_name": "PROTAC 2", "compound_number": "3"}, incoming)
        self.assertEqual(result["status"], "reference_missing_requires_review")
        self.assertIsNone(self.material["h1"]["observations"][1]["original"]["value"])

    def test_observation_name_conflicts_even_when_query_only_supplies_number(self):
        incoming = {**self.material["h1"]["observations"][0]["original"], "paper_name": "ACBI1"}
        result = compare_observation(self.material, {"doi": self.query["doi"], "compound_number": "2"}, incoming)
        self.assertEqual(result["status"], "identity_conflict")

    def test_existing_artifact_ref_with_injected_reader(self):
        manifest = parse((EXAMPLE / "manifest.json").read_bytes())
        data = (EXAMPLE / "material.json").read_bytes()
        ref = manifest["material_ref"]
        self.assertEqual(read_material(ref, lambda _: data), self.material)
        with self.assertRaisesRegex(HandoffError, "HASH_MISMATCH"):
            read_material(ref, lambda _: data + b" ")
        with self.assertRaisesRegex(HandoffError, "REF_KIND"):
            read_material({**ref, "schema_id": "urn:approval"}, lambda _: data)

    def test_input_change_requires_reassessment_without_changing_approval(self):
        modified = copy.deepcopy(self.material)
        modified["h1"]["coverage"] += "; new coverage annotation"
        result = compare_versions(self.material, reseal(modified))
        self.assertTrue(result["requires_M2_reassessment"])
        self.assertEqual(result["changed_sections"], ["h1"])
        self.assertFalse(result["approval_state_modified"])
        self.assertFalse(compare_versions(self.material, self.material)["requires_M2_reassessment"])

    def test_duplicate_and_nonfinite_json_rejected(self):
        for data in (b'{"x":1,"x":2}', b'{"x":NaN}'):
            with self.assertRaises(HandoffError): parse(data)


if __name__ == "__main__":
    unittest.main()
