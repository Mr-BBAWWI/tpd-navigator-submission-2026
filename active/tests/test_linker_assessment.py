import json
import tempfile
import unittest
from pathlib import Path

from packages.science import linker_assessment as la
from packages.science import linker_design


class ExpandedLibraryTests(unittest.TestCase):
    def test_original_ids_and_structures_are_preserved_and_families_are_covered(self):
        base = linker_design.library()
        expanded = la.expanded_library()
        original = {item["id"]: item["smiles"] for item in base["templates"]}
        resulting = {item["id"]: item["smiles"] for item in expanded["templates"]}
        self.assertGreaterEqual(len(resulting), 18)
        for template_id, smiles in original.items():
            self.assertEqual(resulting[template_id], smiles)
        report = la.validate_library(expanded)
        self.assertTrue(report["valid"])
        self.assertEqual(set(report["families"]), {
            "PEG", "alkyl", "PEG-alkyl", "piperazine", "triazole", "rigid-aromatic"
        })
        for lengths in report["families"].values():
            self.assertEqual(set(lengths), {"long", "medium", "short"})

    def test_generated_templates_are_hypotheses_without_exact_precedence_claim(self):
        library = la.expanded_library()
        generated = library["precedence"]["peg_short"]
        self.assertIn("가설", generated["design_status"])
        self.assertEqual(generated["exact_structure_precedence"], "확인된 근거 없음")
        json.dumps(library, allow_nan=False)

    def test_duplicate_id_and_bad_dummy_maps_are_rejected(self):
        library = la.expanded_library()
        library["templates"][1]["id"] = library["templates"][0]["id"]
        with self.assertRaisesRegex(ValueError, "DUPLICATE"):
            la.validate_library(library)

        library = la.expanded_library()
        library["templates"][0]["smiles"] = "[*:1001]CC[*:1001]"
        with self.assertRaisesRegex(ValueError, "DUMMY_MAPS"):
            la.validate_library(library)


class LinkerAssessmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.templates = {item["id"]: item for item in la.expanded_library()["templates"]}

    def test_missing_ternary_geometry_remains_unknown_and_descriptors_are_finite(self):
        result = la.assess_linker(self.templates["peg3"])
        self.assertFalse(result["hard_rejected"])
        self.assertEqual(result["ternary_geometry"]["status"], "미지")
        self.assertIn("warhead-only", result["ternary_geometry"]["reason"])
        self.assertGreater(result["graph_span"]["terminal_heavy_atom_shortest_path_bonds"], 0)
        self.assertGreater(result["contour_bound"]["conservative_total_upper_A"], 0)
        self.assertIn(result["sampled_3d"]["status"], {"계산됨", "계산되지 않음"})
        if result["sampled_3d"]["status"] == "계산됨":
            observed = result["sampled_3d"]["end_to_end_range_A"]
            self.assertLessEqual(observed["min_observed"], observed["max_observed"])
            self.assertIn("이론적 최대", result["sampled_3d"]["interpretation"])
        json.dumps(result, allow_nan=False)

    def test_only_distance_beyond_conservative_bound_hard_rejects(self):
        baseline = la.assess_linker(self.templates["alkyl_short"])
        bound = baseline["contour_bound"]["conservative_total_upper_A"]
        accepted = la.assess_linker(
            self.templates["alkyl_short"], required_attachment_distance_A=bound
        )
        rejected = la.assess_linker(
            self.templates["alkyl_short"], required_attachment_distance_A=bound + 0.01
        )
        self.assertFalse(accepted["hard_rejected"])
        self.assertTrue(rejected["hard_rejected"])
        self.assertIsNotNone(rejected["hard_rejection_reason"])

    def test_endpoint_reactions_are_labeled_as_hypotheses(self):
        result = la.assess_linker(
            self.templates["triazole"],
            endpoint_chemistry={"left": "azide", "right": "terminal alkyne"},
        )
        proposals = result["endpoint_chemistry"]["proposals"]
        self.assertEqual(len(proposals), 1)
        self.assertIn("click", proposals[0]["type"])
        self.assertIn("가설", proposals[0]["status"])
        self.assertIn("검증하지", result["endpoint_chemistry"]["interpretation"])

    def test_invalid_distance_is_not_silently_coerced(self):
        with self.assertRaises(TypeError):
            la.assess_linker(self.templates["peg_short"], True)
        with self.assertRaises(ValueError):
            la.assess_linker(self.templates["peg_short"], float("nan"))


class RouteEvidenceTests(unittest.TestCase):
    def complete_record(self, compound_id="CMP-1", record_id="R1"):
        return {
            "record_id": record_id,
            "compound_id": compound_id,
            "source_url": "https://example.org/supporting-information",
            "locator": "Supporting information, page 12, compound CMP-1",
            "steps": ["중간체 A와 B를 반응시켰다고 source에 기재됨"],
            "materials": ["중간체 A", "중간체 B"],
            "yield": "source-reported 42%",
            "purification": "source-reported chromatography",
            "characterization": ["source-reported MS entry"],
        }

    def test_different_compound_never_satisfies_exact_route(self):
        precedent = self.complete_record("OTHER-9")
        precedent["reusable_precedent"] = True
        result = la.route_evidence([precedent], "CMP-1")
        self.assertEqual(result["complete_exact_documentation_count"], 0)
        self.assertEqual(len(result["exact_compound_records"]), 0)
        self.assertEqual(len(result["reusable_precedent_records"]), 1)
        self.assertIn("정확한 compound ID route 근거 없음", result["status"])

    def test_missing_route_fields_are_reported_not_fabricated(self):
        incomplete = {
            "record_id": "R2",
            "compound_id": "CMP-1",
            "source_url": "not a URL",
            "locator": "",
            "steps": [],
        }
        result = la.route_evidence([incomplete], "CMP-1")
        report = result["exact_compound_records"][0]
        self.assertFalse(report["documentation_complete"])
        self.assertIn("materials", report["missing_or_invalid_fields"])
        self.assertIn("source_url", report["missing_or_invalid_fields"])
        self.assertNotIn("materials", report["source_record"])

    def test_complete_exact_documentation_is_not_called_synthetic_validation(self):
        result = la.route_evidence(
            [self.complete_record()],
            "CMP-1",
            generated_candidate={"novel_transformations": [{"type": "amidation"}]},
        )
        self.assertEqual(result["complete_exact_documentation_count"], 1)
        self.assertIn("독립적 재현성", result["conclusion"])
        self.assertIn("별도 검토", result["novel_transformations"]["evidence_status"])

    def test_local_json_ingest_and_duplicate_record_validation(self):
        records = [self.complete_record(), self.complete_record("CMP-2", "R2")]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "evidence.json"
            path.write_text(json.dumps({"records": records}), encoding="utf-8")
            loaded = la.ingest_source_json(path)
            self.assertEqual(len(loaded), 2)
            report = la.validate_source_json(path)
            self.assertEqual(report["record_count"], 2)

        duplicate = [self.complete_record(), self.complete_record("CMP-2", "R1")]
        with self.assertRaisesRegex(ValueError, "DUPLICATE"):
            la.validate_source_json(duplicate)

    def test_malformed_source_json_fails(self):
        with self.assertRaises(json.JSONDecodeError):
            la.ingest_source_json("[{bad json]")


class SummaryTests(unittest.TestCase):
    def test_sa_descriptor_does_not_assert_synthesizability(self):
        template = la.expanded_library()["templates"][0]
        assessment = la.assess_linker(template)
        summary = la.assessment_summary(assessment, sa_score=3.2)
        self.assertEqual(summary["synthesizability_conclusion"], "미지/검토 필요")
        self.assertIn("synthesizability", summary["descriptor_only_sa"]["interpretation"])
        json.dumps(summary, allow_nan=False)

    def test_nonfinite_sa_is_rejected(self):
        with self.assertRaises(ValueError):
            la.assessment_summary({}, sa_score=float("inf"))


if __name__ == "__main__":
    unittest.main()
