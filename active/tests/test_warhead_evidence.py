import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from packages.science.warhead_evidence import (
    SCHEMA_VERSION,
    classify_construct_domain,
    deduplicate_records,
    evaluate_candidate,
    load_catalog,
)


class WarheadEvidenceTests(unittest.TestCase):
    def eligible(self, **changes):
        candidate = {
            "target_accessions": ["P51531"],
            "heavy_atom_count": 20,
            "molecular_weight_Da": 300.0,
            "contains_carbon": True,
            "ccd": "ABC",
            "name": "organic inhibitor",
            "multiligand_crosslinked": False,
            "complete_protac": False,
            "contacting_ligand_heavy_atoms": 4,
            "contacting_protein_heavy_atoms": 5,
        }
        candidate.update(changes)
        return candidate

    def test_non_target_is_excluded(self):
        reason = evaluate_candidate(self.eligible(target_accessions=["P51532"]))
        self.assertEqual(reason, "non_target_chain")

    def test_complete_protac_is_excluded(self):
        reason = evaluate_candidate(self.eligible(name="complete PROTAC degrader"))
        self.assertEqual(reason, "complete_PROTAC_not_a_warhead")

    def test_multiligand_crosslink_is_excluded(self):
        reason = evaluate_candidate(self.eligible(multiligand_crosslinked=True))
        self.assertEqual(reason, "multiligand_crosslinked")

    def test_adp_is_excluded_as_natural_cofactor(self):
        reason = evaluate_candidate(self.eligible(ccd="ADP", name="ADP"))
        self.assertEqual(reason, "natural_metabolite_or_cofactor")

    def test_different_construct_domain_is_not_qualified(self):
        result = classify_construct_domain(
            [{
                "pdb_chain_ids": ["A"],
                "uniprot_sequence_begin": 705,
                "uniprot_sequence_end": 955,
            }],
            ["AA"],
            {"AA": ["A"]},
        )
        self.assertFalse(result["qualified_for_cards"])
        self.assertEqual(result["domain_annotation"], "different_construct_range")

    def test_genuine_fx5_bromodomain_construct_is_qualified(self):
        result = classify_construct_domain(
            [{
                "pdb_chain_ids": ["A"],
                "uniprot_sequence_begin": 1373,
                "uniprot_sequence_end": 1493,
            }],
            ["AA"],
            {"AA": ["A"]},
        )
        self.assertTrue(result["qualified_for_cards"])
        self.assertEqual(result["domain_annotation"], "bromodomain")

    def test_duplicate_keeps_all_experimental_references(self):
        common = {
            "parent_connectivity_key": "C[NH2+]C",
            "identity": {"canonical_isomeric_smiles": "C[NH2+]C"},
            "ccd": "LIG",
        }
        records = [
            {**common, "record_id": "1AAA:C:LIG", "pdb": "1AAA", "ligand_label_asym_id": "C"},
            {**common, "record_id": "2BBB:D:LIG", "pdb": "2BBB", "ligand_label_asym_id": "D"},
        ]
        groups = deduplicate_records(records)
        self.assertEqual(len(groups), 1)
        self.assertEqual(
            [version["record_id"] for version in groups[0]["experimental_versions"]],
            ["1AAA:C:LIG", "2BBB:D:LIG"],
        )

    def test_changed_source_hash_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source_snapshots" / "abc" / "source.cif"
            source.parent.mkdir(parents=True)
            original = b"original source\n"
            source.write_bytes(original)
            catalog = {
                "schema_version": SCHEMA_VERSION,
                "sources_root": "source_snapshots",
                "sources": [{
                    "relative_path": "abc/source.cif",
                    "sha256": hashlib.sha256(original).hexdigest(),
                }],
                "artifacts": [],
            }
            path = root / "catalog.json"
            path.write_text(json.dumps(catalog), encoding="utf-8")
            source.write_bytes(b"changed source\n")
            with self.assertRaisesRegex(ValueError, "Source hash verification failed"):
                load_catalog(path)

    def test_no_measured_kd_is_inferred_by_candidate_filter(self):
        candidate = self.eligible()
        self.assertIsNone(evaluate_candidate(candidate))
        self.assertNotIn("measured_Kd", candidate)
        self.assertNotIn("Kd", candidate)

    def test_contact_requirement_applies_to_both_sides(self):
        self.assertEqual(
            evaluate_candidate(self.eligible(contacting_ligand_heavy_atoms=2)),
            "fewer_than_3_contacting_ligand_heavy_atoms",
        )
        self.assertEqual(
            evaluate_candidate(self.eligible(contacting_protein_heavy_atoms=2)),
            "fewer_than_3_contacting_protein_heavy_atoms",
        )


if __name__ == "__main__":
    unittest.main()
