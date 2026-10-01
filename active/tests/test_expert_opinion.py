import copy
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from packages.contracts import ContractError
from packages.science.expert_opinion import (
    FX5_SOURCE_FILES,
    MANIFEST_RELATIVE_PATH,
    ROOT,
    SOURCE_RELATIVE_PATH,
    apply_bundled_opinion,
    load_bundled_opinion,
)


def base_policy():
    return {
        "_validated_by_platform": True,
        "numerical_criteria": {
            "modifiable_sites_min": 2,
            "distinct_graphs_per_site_min": 30,
            "calibration_seeds_min": 3,
            "novel_ternary_repeats_min": 2,
        },
        "decisions": {},
    }


def source_binding():
    return {
        path: {"status": "configured_hash_verified", "sha256": digest}
        for path, digest in FX5_SOURCE_FILES.items()
    }


def mapping():
    return [
        {"atom_map": atom_map, "rdkit_index_zero_based": atom_map - 1}
        for atom_map in range(1, 21)
    ]


def site_atoms():
    values = []
    for atom_map in range(1, 21):
        if atom_map == 19:
            state = "MODIFIABLE"
        elif atom_map in (8, 11, 12):
            state = "UNKNOWN"
        else:
            state = "PROTECTED"
        values.append({"atom_map": atom_map, "state": state})
    return values


def actual_shape_result():
    common = list(range(1, 21))
    core = [1, 2, 3, 4, 5, 6, 7, 9, 10, 13, 14, 15, 16, 17, 18, 20]
    diagnostics = []
    for index in range(5):
        diagnostics.append({
            "pose_index_zero_based": index,
            "common_mapped_heavy_atom_count": 20,
            "core_mapped_atom_count": 16,
            "common_atom_maps": common,
            "core_atom_maps": core,
            "docking_pose_preserved": index == 0,
            "status": "pass" if index == 0 else "fail",
        })
    return {
        "input_binding": {"source_files": source_binding()},
        "parent_scope": {
            "actual_design_parent_id": "SMARCA2-FX5",
            "receptor_frame": "6HAZ chain A",
            "parent_redock_required_for_qualification": True,
        },
        "warhead": {
            "id": "SMARCA2-FX5",
            "pdb": "6HAZ",
            "target_chain": "A",
            "parent_identity": {
                "canonical_isomeric_smiles": "Nc1nnc(-c2ccccc2O)cc1N1CC[NH2+]CC1",
                "inchikey": "SZKHGLTYXIDOFH-UHFFFAOYSA-O",
                "formal_charge": 1,
                "heavy_atom_count": 20,
            },
            "design_identity": {
                "canonical_isomeric_smiles": "Nc1nnc(-c2ccccc2O)cc1N1CCNCC1",
                "inchikey": "SZKHGLTYXIDOFH-UHFFFAOYSA-N",
                "formal_charge": 0,
                "heavy_atom_count": 20,
            },
            "atom_mapping": mapping(),
        },
        "parent_redocking": {
            "status": "completed_with_limits",
            "real_docking": True,
            "parent_redocking": True,
            "results": {
                "pose_preservation": {
                    "docking_pose_preserved": True,
                    "status": "pass",
                    "all_poses_diagnostics": diagnostics,
                }
            },
        },
        "sites": {
            "atoms": site_atoms(),
            "method": {
                "SASA_tool": "RDKit rdFreeSASA",
                "SASA_algorithm": "LeeRichards",
            },
        },
    }


class ExpertOpinionTests(unittest.TestCase):
    def apply(self, result=None, decisions=None, policy=None, verified=True, root=ROOT):
        return apply_bundled_opinion(
            base_policy() if policy is None else policy,
            actual_shape_result() if result is None else result,
            [] if decisions is None else decisions,
            archived_result_verified=verified,
            root=root,
        )

    def test_actual_shape_applies_fx5_waiver_and_both_minimum_floors(self):
        policy = self.apply()
        numerical = policy["numerical_criteria"]
        opinion = policy["opinion"]
        self.assertTrue(opinion["scope"]["fx5_bound_scope"])
        self.assertTrue(opinion["site_waiver"]["applied"])
        self.assertEqual(numerical["modifiable_sites_min"], 1)
        self.assertEqual(numerical["distinct_graphs_per_site_min"], 30)
        self.assertEqual(numerical["calibration_seeds_min"], 5)
        self.assertEqual(numerical["novel_ternary_repeats_min"], 5)
        self.assertEqual(opinion["analysis_provenance"], {"status": "pending_metadata"})
        self.assertEqual(opinion["planning"]["execution_status"], "pending")
        self.assertFalse(opinion["formal_approval"])
        self.assertFalse(opinion["scientific_accepted"])

    def test_other_parent_does_not_change_any_minimum(self):
        result = actual_shape_result()
        result["parent_scope"]["actual_design_parent_id"] = "SMARCA2-OTHER"
        policy = self.apply(result)
        self.assertEqual(policy["numerical_criteria"]["modifiable_sites_min"], 2)
        self.assertEqual(policy["numerical_criteria"]["calibration_seeds_min"], 3)
        self.assertEqual(policy["numerical_criteria"]["novel_ternary_repeats_min"], 2)

    def test_additional_source_fields_and_files_do_not_hide_required_binding_errors(self):
        result = actual_shape_result()
        for record in result["input_binding"]["source_files"].values():
            record["bytes"] = 123
            record["verified"] = True
        result["input_binding"]["source_files"]["design_sources/additional-provenance.json"] = {
            "status": "configured_hash_verified",
            "sha256": "a" * 64,
        }
        policy = self.apply(result)
        self.assertTrue(policy["opinion"]["scope"]["fx5_bound_scope"])
        self.assertTrue(policy["opinion"]["site_waiver"]["applied"])

        result["input_binding"]["source_files"]["design_sources/6HAZ.cif"]["sha256"] = "f" * 64
        failures = self.apply(result)["opinion"]["scope"]["condition_failures"]
        self.assertIn("INPUT_BINDING_SOURCE_RECORD_MISMATCH", failures)

    def test_bad_frame_wrong_source_and_unverified_archive_are_rejected(self):
        bad_frame = actual_shape_result()
        bad_frame["parent_scope"]["receptor_frame"] = "6HAZ chain B"
        self.assertIn(
            "RECEPTOR_FRAME_NOT_EXACT_6HAZ_CHAIN_A",
            self.apply(bad_frame)["opinion"]["scope"]["condition_failures"],
        )
        bad_source = actual_shape_result()
        bad_source["input_binding"]["source_files"]["design_sources/6HAZ.cif"]["sha256"] = "f" * 64
        self.assertIn(
            "INPUT_BINDING_SOURCE_RECORD_MISMATCH",
            self.apply(bad_source)["opinion"]["scope"]["condition_failures"],
        )
        with self.assertRaises(ContractError):
            self.apply(verified=False)

    def test_true_fields_require_bool_not_numeric_or_text(self):
        for field in ("real_docking", "parent_redocking"):
            for value in (1, "true"):
                result = actual_shape_result()
                result["parent_redocking"][field] = value
                failures = self.apply(result)["opinion"]["site_waiver"]["condition_failures"]
                self.assertTrue(any("BOOLEAN_TRUE_REQUIRED" in item for item in failures))

    def test_duplicate_mapping_and_duplicate_site_map_are_rejected(self):
        graph_duplicate = actual_shape_result()
        graph_duplicate["warhead"]["atom_mapping"][1]["atom_map"] = 1
        self.assertIn(
            "ATOM_MAPPING_DUPLICATE_OR_INCOMPLETE",
            self.apply(graph_duplicate)["opinion"]["scope"]["condition_failures"],
        )
        site_duplicate = actual_shape_result()
        site_duplicate["sites"]["atoms"][1]["atom_map"] = 1
        self.assertIn(
            "SITE_MAP_DUPLICATE",
            self.apply(site_duplicate)["opinion"]["site_waiver"]["condition_failures"],
        )

    def test_unknown_sites_are_observed_but_never_changed(self):
        result = actual_shape_result()
        before = copy.deepcopy(result["sites"]["atoms"])
        policy = self.apply(result)
        self.assertEqual(result["sites"]["atoms"], before)
        observed = policy["opinion"]["site_waiver"]["site_states_observed_without_mutation"]
        self.assertEqual(observed["8"], "UNKNOWN")
        self.assertEqual(observed["11"], "UNKNOWN")
        self.assertEqual(observed["12"], "UNKNOWN")

    def test_direct_numerical_values_have_priority_for_each_key(self):
        policy = base_policy()
        policy["numerical_criteria"].update({
            "modifiable_sites_min": 4,
            "calibration_seeds_min": 2,
            "novel_ternary_repeats_min": 3,
        })
        decision = {
            "id": "direct-numerical",
            "action": "accept",
            "kind": "numerical_criteria",
            "data": {"values": {
                "modifiable_sites_min": 4,
                "calibration_seeds_min": 2,
                "novel_ternary_repeats_min": 3,
            }},
        }
        applied = self.apply(decisions=[decision], policy=policy)
        self.assertEqual(applied["numerical_criteria"]["modifiable_sites_min"], 4)
        self.assertEqual(applied["numerical_criteria"]["calibration_seeds_min"], 2)
        self.assertEqual(applied["numerical_criteria"]["novel_ternary_repeats_min"], 3)

    def test_wrong_diagnostic_counts_or_duplicate_pose_index_block_waiver(self):
        wrong = actual_shape_result()
        wrong["parent_redocking"]["results"]["pose_preservation"]["all_poses_diagnostics"][0]["core_mapped_atom_count"] = 15
        self.assertIn(
            "PARENT_POSE_MAPPING_COUNTS_MISMATCH",
            self.apply(wrong)["opinion"]["site_waiver"]["condition_failures"],
        )
        duplicate = actual_shape_result()
        duplicate["parent_redocking"]["results"]["pose_preservation"]["all_poses_diagnostics"][1]["pose_index_zero_based"] = 0
        self.assertIn(
            "PARENT_POSE_INDEX_DUPLICATE",
            self.apply(duplicate)["opinion"]["site_waiver"]["condition_failures"],
        )

    def _portable_root(self):
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)
        for relative in (MANIFEST_RELATIVE_PATH, SOURCE_RELATIVE_PATH):
            destination = root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / relative, destination)
        for relative in FX5_SOURCE_FILES:
            source = ROOT / "cases" / relative
            destination = root / "cases" / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
        return temporary, root

    def test_portable_case_clone_and_source_fingerprints(self):
        temporary, root = self._portable_root()
        try:
            policy = self.apply(root=root)
            self.assertTrue(policy["opinion"]["site_waiver"]["applied"])
            self.assertEqual(set(policy["opinion"]["verified_input_sources"]), set(FX5_SOURCE_FILES))
            self.assertEqual(len(policy["opinion"]["policy_digest"]), 64)
        finally:
            temporary.cleanup()

    def test_unicode_root_reads_actual_sdf_files(self):
        temporary, root = self._portable_root()
        unicode_root = Path(temporary.name) / "바이오공모전"
        try:
            shutil.copytree(root, unicode_root)
            policy = self.apply(root=unicode_root)
            self.assertTrue(policy["opinion"]["scope"]["fx5_bound_scope"])
            self.assertTrue(policy["opinion"]["site_waiver"]["applied"])
        finally:
            temporary.cleanup()

    def test_invalid_record_after_valid_sdf_molecule_is_rejected(self):
        temporary, root = self._portable_root()
        try:
            target = root / "cases/design_sources/SMARCA2-neutral-design.sdf"
            with target.open("ab") as stream:
                stream.write(b"invalid duplicate record\n$$$$\n")
            result = actual_shape_result()
            result["input_binding"]["source_files"]["design_sources/SMARCA2-neutral-design.sdf"]["sha256"] = __import__("hashlib").sha256(target.read_bytes()).hexdigest()
            failures = self.apply(result=result, root=root)["opinion"]["scope"]["condition_failures"]
            self.assertIn("PARENT_CANONICAL_GRAPH_UNVERIFIED", failures)
        finally:
            temporary.cleanup()

    def test_tampered_manifest_document_or_bound_source_is_rejected(self):
        for relative in (MANIFEST_RELATIVE_PATH, SOURCE_RELATIVE_PATH, "design_sources/6HAZ.cif"):
            temporary, root = self._portable_root()
            try:
                target = root / relative if relative in (MANIFEST_RELATIVE_PATH, SOURCE_RELATIVE_PATH) else root / "cases" / relative
                with target.open("ab") as stream:
                    stream.write(b"tamper")
                if relative in (MANIFEST_RELATIVE_PATH, SOURCE_RELATIVE_PATH):
                    with self.assertRaises(ContractError):
                        load_bundled_opinion(root)
                else:
                    failures = self.apply(root=root)["opinion"]["scope"]["condition_failures"]
                    self.assertIn("BOUND_SOURCE_CURRENT_HASH_MISMATCH", failures)
            finally:
                temporary.cleanup()

    def test_engine_fingerprint_contains_module_manifest_and_document(self):
        bundle = load_bundled_opinion()
        self.assertEqual(len(bundle["decision_extract_sha256"]), 64)
        self.assertEqual(len(bundle["conditions_extract_sha256"]), 64)
        from packages.platform.scientific_acceptance import ENGINE_FILES, engine_fingerprint
        self.assertIn("packages/science/expert_opinion.py", ENGINE_FILES)
        self.assertIn(MANIFEST_RELATIVE_PATH, ENGINE_FILES)
        self.assertIn(SOURCE_RELATIVE_PATH, ENGINE_FILES)
        self.assertEqual(len(engine_fingerprint()), 64)


if __name__ == "__main__":
    unittest.main()
