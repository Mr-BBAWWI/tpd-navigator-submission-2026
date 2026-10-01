import json
import unittest
from unittest import mock

from rdkit import Chem
from rdkit.Chem import AllChem

import packages.science.interaction_review as ir
from packages.science.interaction_review import assess_interactions
from packages.science.scientific_assessment import evaluate_design
from packages.science.synthesis_review import assess_synthesis


def mol3d(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise AssertionError(f"invalid fixture SMILES: {smiles}")
    for index, atom in enumerate(mol.GetAtoms(), 1):
        if atom.GetAtomicNum() > 1:
            atom.SetAtomMapNum(index)
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=23) != 0:
        raise AssertionError(f"could not embed fixture SMILES: {smiles}")
    mol = Chem.RemoveHs(mol)
    mol.GetConformer().Set3D(True)
    return mol


def protein(atom="OD1", residue="ASN", xyz=(2.8, 0.0, 0.0)):
    return [{
        "xyz": list(xyz),
        "type_symbol": "O",
        "label_comp_id": residue,
        "auth_seq_id": "1",
        "label_atom_id": atom,
        "label_asym_id": "A",
    }]


def candidate():
    return {
        "candidate_id": "D-1",
        "mapped_smiles": "[CH3:1][CH2:2][CH2:3][CH3:4]",
        "canonical_smiles": "CCCC",
        "atom_roles": {
            "warhead_maps": [1],
            "linker_maps": [2, 3],
            "recruiter_maps": [4],
        },
        "attachment_metadata": {
            "warhead_linker_bond": {
                "attachment_atom_map": 1,
                "partner_atom_map": 2,
                "bond_type": "SINGLE",
            },
            "recruiter_linker_bond": {
                "attachment_atom_map": 4,
                "partner_atom_map": 3,
                "bond_type": "SINGLE",
            },
        },
    }


def base_result():
    return {
        "summary": {"valid_analogs": 999},
        "sites": {
            "atoms": [
                {"atom_map": 19, "state": "MODIFIABLE"},
                {"atom_map": 8, "state": "UNKNOWN"},
                {"atom_map": 11, "state": "UNKNOWN"},
                {"atom_map": 12, "state": "UNKNOWN"},
            ]
        },
        "analogs": [],
        "protac_candidates": [],
        "calibration": {"CRBN": {"seed_receipts": []}},
    }


def valid_analog(identifier, site, smiles="CC", selected=False, preview=False):
    return {
        "id": identifier,
        "mapped_smiles": smiles,
        "canonical_smiles": smiles,
        "attachment_site_atom_maps": [site],
        "transformation_class": "ring_expansion",
        "cheap_filter": {"valid": True},
        "selected": selected,
        "pipeline_status": "qualified",
        "assembly_eligible": True,
        "parent_redocking_supported": True,
        "preview": preview,
        "docking": {
            "status": "completed_with_limits",
            "pose_preserved": True,
            "passing_pose_count": 1,
        },
    }


def criterion_by_id(assessment, identifier):
    return {item["id"]: item for item in assessment["criteria"]}[identifier]


class ScientificAssessmentTests(unittest.TestCase):
    def test_exact_participant_l1_does_not_match_l10_or_forged_protein_id(self):
        requirement = {
            "kind": "proximity",
            "ligand_maps": [1],
            "protein_atom_ids": ["A:1:ASN:OD1"],
        }
        self.assertFalse(ir._participant_match({
            "kind": "hydrophobic_contact",
            "participants": ["L:10", "A:1:ASN:OD1"],
        }, requirement))
        self.assertFalse(ir._participant_match({
            "kind": "hydrophobic_contact",
            "participants": ["L:1", "prefix-A:1:ASN:OD1-forged"],
        }, requirement))
        self.assertTrue(ir._participant_match({
            "kind": "hydrophobic_contact",
            "participants": ["L:1", "A:1:ASN:OD1"],
        }, requirement))

    def test_actual_aromatic_profile_format_and_exact_ring_map_set(self):
        interaction = {
            "kind": "aromatic_contact_geometry",
            "participants": [
                "L:ring:L:1,L:2,L:3,L:4,L:5,L:6",
                "A:1417:PHE",
            ],
        }
        requirement = {
            "kind": "aromatic",
            "ligand_maps": [1, 2, 3, 4, 5, 6],
            "protein_atom_ids": ["A:1417:PHE:CZ"],
        }
        self.assertTrue(ir._participant_match(interaction, requirement))

        wrong_maps = dict(requirement, ligand_maps=[1, 2, 3, 4, 5, 10])
        self.assertFalse(ir._participant_match(interaction, wrong_maps))

        wrong_atom = dict(
            requirement,
            protein_atom_ids=["A:1417:PHE:CA"],
        )
        self.assertFalse(ir._participant_match(interaction, wrong_atom))

    def test_gained_interaction_is_not_preserved(self):
        parent = mol3d("N")
        pose = Chem.Mol(parent)
        profiles = iter([
            ({"interactions": []}, {"requires_review": False}),
            ({
                "interactions": [{
                    "id": "gained-h1",
                    "kind": "directional_hbond",
                    "participants": ["L:1", "A:1:ASN:OD1"],
                }]
            }, {"requires_review": False}),
            ({"interactions": []}, {"requires_review": False}),
        ])
        requirements = [{
            "id": "h1",
            "kind": "directional_hbond",
            "ligand_maps": [1],
            "protein_atom_ids": ["A:1:ASN:OD1"],
            "required": True,
        }]

        with mock.patch.object(
            ir,
            "_profile_for",
            side_effect=lambda *args, **kwargs: next(profiles),
        ), mock.patch.object(
            ir,
            "enumerate_microstates",
            return_value=[Chem.Mol(pose)],
        ):
            report = assess_interactions(
                parent,
                pose,
                protein(),
                requirements=requirements,
                max_states=1,
            )

        requirement = report["requirements"][0]
        self.assertEqual(requirement["status"], "pending")
        self.assertEqual(
            requirement["reason"],
            "baseline_definition_required",
        )
        self.assertIs(requirement["computed_pass"], False)

    def test_missing_pose_map_is_graceful_loss(self):
        parent = mol3d("CN")
        pose = mol3d("C")
        requirements = [{
            "id": "x",
            "kind": "hydrophobic",
            "ligand_maps": [2],
            "protein_atom_ids": ["A:1:ASN:OD1"],
            "required": True,
        }]
        computed_profile = ({
            "interactions": [{
                "id": "x-contact",
                "kind": "hydrophobic_contact",
                "participants": ["L:2", "A:1:ASN:OD1"],
            }]
        }, {"requires_review": False})

        with mock.patch.object(
            ir,
            "_profile_for",
            return_value=computed_profile,
        ), mock.patch.object(
            ir,
            "enumerate_microstates",
            return_value=[],
        ):
            report = assess_interactions(
                parent,
                pose,
                protein(),
                requirements=requirements,
            )

        self.assertEqual(report["requirements"][0]["status"], "lost")

    def test_synthesized_false_and_foreign_route_remain_undocumented(self):
        design = candidate()
        design["synthesized"] = False
        route = [
            {
                "candidate_id": "D-1",
                "canonical_smiles": "CCO",
                "source": "paper",
                "locator": "page 1",
                "steps": [{
                    "reagents": ["x"],
                    "conditions": "room temperature",
                    "purification": "HPLC",
                    "characterization": "NMR",
                }],
            },
            "not-a-record",
        ]

        report = assess_synthesis(design, exact_route_records=route)

        self.assertIs(report["documented_route_for_graph"], False)
        self.assertEqual(len(report["route_evidence"]["records"]), 2)
        self.assertTrue(all(
            record["documentation_complete"] is False
            for record in report["route_evidence"]["records"]
        ))

    def test_unknown_site_graphs_cannot_satisfy_modifiable_site_criterion(self):
        result = base_result()
        result["analogs"] = [
            valid_analog(f"U-{site}-{number}", site, "C" * (number + 1))
            for site in (8, 11, 12)
            for number in range(30)
        ]

        criterion = criterion_by_id(
            evaluate_design(result),
            "distinct_constitutional_graphs",
        )

        self.assertEqual(criterion["status"], "failed")
        self.assertEqual(criterion["observed"]["qualifying_sites"], [])

    def test_docking_status_alone_does_not_qualify_and_assembly_preview_is_excluded(self):
        result = base_result()
        result["sites"]["atoms"].append({
            "atom_map": 20,
            "state": "MODIFIABLE",
        })

        preview = valid_analog(
            "preview",
            19,
            "CC",
            selected=True,
            preview=True,
        )
        preview["pipeline_status"] = "preview"
        preview["docking"] = {"status": "completed_with_limits"}

        status_only = valid_analog("status-only", 19, "CCC", selected=True)
        status_only["docking"] = {"status": "completed"}
        result["analogs"] = [preview, status_only]
        result["protac_candidates"] = [
            {
                "candidate_id": "preview-crbn",
                "e3_type": "CRBN",
                "mapped_smiles": "CC",
                "canonical_smiles": "CC",
                "preview": True,
            },
            {
                "candidate_id": "assembled-vhl",
                "e3_type": "VHL",
                "mapped_smiles": "CCC",
                "canonical_smiles": "CCC",
            },
        ]

        assessment = evaluate_design(result)
        panel = criterion_by_id(assessment, "qualified_panel")
        assembly = criterion_by_id(assessment, "both_e3_assembly")

        self.assertEqual(
            panel["observed"]["qualified_selected_unique"],
            0,
        )
        self.assertEqual(panel["status"], "failed")
        self.assertEqual(assembly["observed"], {"CRBN": 0, "VHL": 1})
        self.assertEqual(assembly["status"], "failed")

    def test_constitutional_graph_deduplication_ignores_stereochemistry(self):
        def stereoisomer(identifier, mapped_smiles):
            analog = valid_analog(
                identifier,
                19,
                mapped_smiles,
                selected=True,
            )
            mol = Chem.MolFromSmiles(mapped_smiles)
            self.assertIsNotNone(mol)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
            analog["canonical_smiles"] = Chem.MolToSmiles(
                mol,
                canonical=True,
                isomericSmiles=True,
            )
            analog["docking"] = {
                "status": "completed_with_limits",
                "pose_preserved": True,
                "passing_pose_count": 1,
            }
            return analog

        result = base_result()
        result["analogs"] = [
            stereoisomer("R", "[F:1][C@H:2]([Cl:3])[Br:4]"),
            stereoisomer("S", "[F:1][C@@H:2]([Cl:3])[Br:4]"),
        ]
        policy = {
            "_validated_by_platform": True,
            "numerical_criteria": {"panel_min": 2},
        }

        assessment = evaluate_design(result, policy=policy)
        graphs = criterion_by_id(
            assessment,
            "distinct_constitutional_graphs",
        )
        panel = criterion_by_id(assessment, "qualified_panel")

        self.assertEqual(
            graphs["observed"]["counts_by_site"],
            {"19": 1},
        )
        self.assertEqual(
            panel["observed"]["qualified_selected_unique"],
            2,
        )
        self.assertEqual(panel["status"], "pass")

    def test_invalid_or_constitutionally_disagreeing_graph_is_rejected(self):
        result = base_result()
        invalid = valid_analog("invalid", 19, "CC", selected=True)
        invalid["mapped_smiles"] = "not-smiles"
        invalid["canonical_smiles"] = "not-smiles"

        disagreeing = valid_analog("disagree", 19, "CC", selected=True)
        disagreeing["canonical_smiles"] = "CCC"
        result["analogs"] = [invalid, disagreeing]

        panel = criterion_by_id(
            evaluate_design(result),
            "qualified_panel",
        )

        self.assertEqual(
            panel["observed"]["qualified_selected_unique"],
            0,
        )

    def test_calibration_reads_inspection_comparison_metrics(self):
        result = base_result()
        result["calibration"]["CRBN"]["seed_receipts"] = [
            {
                "seed": seed,
                "execution_success": True,
                "inspection": {
                    "comparison": {"metrics": {"rmsd": value}}
                },
            }
            for seed, value in ((23, 1.0), (41, 2.0), (67, 1.5))
        ]
        policy = {
            "_validated_by_platform": True,
            "decisions": {
                "calibration_criterion": {
                    "metric": "rmsd",
                    "threshold": 2.0,
                    "expected_seeds": [23, 41, 67],
                    "scope": "CRBN",
                    "comparison": "max_lte",
                }
            },
        }

        criterion = criterion_by_id(
            evaluate_design(result, policy=policy),
            "known_crbn_calibration",
        )

        self.assertEqual(criterion["status"], "pass")
        self.assertEqual(
            criterion["observed"]["metric_distributions"]["rmsd"]["values"],
            [1.0, 2.0, 1.5],
        )

    def test_duplicated_seed_is_not_a_successful_repeat(self):
        result = base_result()
        result["protac_candidates"] = [
            {
                "candidate_id": "C",
                "e3_type": "CRBN",
                "mapped_smiles": "CC",
                "canonical_smiles": "CC",
            },
            {
                "candidate_id": "V",
                "e3_type": "VHL",
                "mapped_smiles": "CCC",
                "canonical_smiles": "CCC",
            },
        ]
        policy = {
            "_validated_by_platform": True,
            "decisions": {
                "ternary_criterion": {
                    "chosen_candidates": {"CRBN": ["C"], "VHL": ["V"]},
                    "expected_seeds": [23, 41],
                }
            },
        }
        digest = evaluate_design(result, policy=policy)["policy_digest"]
        receipts = [{
            "_validated_by_platform": True,
            "job_id": f"j{index}",
            "policy_digest": digest,
            "candidate_id": "C",
            "e3_type": "CRBN",
            "seed": 23,
            "actual_computation": True,
            "execution_success": True,
            "reference_free": True,
            "inspection_summary": {
                "compatibility": "compatible",
                "metrics": {"collision_score": 0.1},
            },
        } for index in range(2)]

        criterion = criterion_by_id(
            evaluate_design(
                result,
                policy=policy,
                supplements={"ternary": receipts},
            ),
            "novel_ternary_repeats",
        )

        self.assertEqual(
            criterion["observed"]["groups"]["C|CRBN"],
            1,
        )
        self.assertEqual(criterion["status"], "pending")

    def test_formal_reject_is_failed_and_engine_never_approves(self):
        result = base_result()
        policy = {
            "_validated_by_platform": True,
            "decisions": {"formal_expert_decision": "reject"},
        }

        assessment = evaluate_design(result, policy=policy)
        criterion = criterion_by_id(
            assessment,
            "formal_expert_decision",
        )

        self.assertEqual(criterion["status"], "failed")
        self.assertIs(assessment["scientific_accepted"], False)

    def test_formal_accept_still_cannot_make_pure_engine_approve(self):
        result = base_result()
        policy = {
            "_validated_by_platform": True,
            "decisions": {"formal_expert_decision": "accept"},
        }

        assessment = evaluate_design(result, policy=policy)

        self.assertEqual(
            criterion_by_id(
                assessment,
                "formal_expert_decision",
            )["status"],
            "pass",
        )
        self.assertIs(assessment["scientific_accepted"], False)

    def _selected_microstate_assessment(self, requirement):
        result = base_result()
        result["analogs"] = [
            valid_analog("A-1", 19, "CC", selected=True)
        ]
        policy = {
            "_validated_by_platform": True,
            "scope": {"job_id": "job-1"},
            "decisions": {
                "microstate_expert_decision": {
                    "accepted": True,
                    "state_choices": {"A-1": 0},
                    "pH_conditions": {
                        "pH": 7.4,
                        "description": "Explicit synthetic assay condition",
                    },
                }
            },
        }
        digest = evaluate_design(result, policy=policy)["policy_digest"]
        supplement = {
            "_validated_by_platform": True,
            "job_id": "job-1",
            "policy_digest": digest,
            "analog_id": "A-1",
            "report_id": "interaction-A-1",
            "pH_conditions": {
                "pH": 7.4,
                "description": "Explicit synthetic assay condition",
            },
            "requirements": [{
                "id": "h",
                "required": True,
                "status": "preserved",
                "computed_pass": True,
                "pending_missing_protein_hydrogen": False,
            }],
            "failures": [],
            "state_alternatives": [{
                "status": "computed",
                "hydrogen_receipt": {"requires_review": False},
                "pH_conditions": {"pH": 7.4},
                "requirements": [requirement],
            }],
        }
        return evaluate_design(
            result,
            policy=policy,
            supplements={"interactions": [supplement]},
        )

    def test_known_h_ambiguous_or_noncomputed_selected_state_is_failed(self):
        assessment = self._selected_microstate_assessment({
            "id": "h",
            "required": True,
            "computed_pass": False,
            "observed": False,
            "unambiguous": False,
            "pending_missing_protein_hydrogen": False,
        })
        criterion = criterion_by_id(
            assessment,
            "microstates_h_direction",
        )

        self.assertEqual(criterion["status"], "failed")
        self.assertTrue(any(
            item.get("analog_id") == "A-1"
            and item.get("requirement_id") == "h"
            and item.get("reason")
            == "required_interaction_not_unambiguous_in_selected_state"
            for item in criterion["observed"]["failures"]
        ))
        self.assertEqual(criterion["observed"]["pending"], [])

    def test_selected_state_missing_protein_hydrogen_is_pending(self):
        assessment = self._selected_microstate_assessment({
            "id": "h",
            "required": True,
            "computed_pass": False,
            "observed": False,
            "unambiguous": False,
            "pending_missing_protein_hydrogen": True,
        })
        criterion = criterion_by_id(
            assessment,
            "microstates_h_direction",
        )

        self.assertEqual(criterion["status"], "pending")
        self.assertEqual(criterion["observed"]["failures"], [])
        self.assertTrue(any(
            item.get("analog_id") == "A-1"
            and item.get("requirement_id") == "h"
            and item.get("reason") == "protein_hydrogen_pending"
            for item in criterion["observed"]["pending"]
        ))

    def test_minimizer_requires_review_remains_pending(self):
        result = base_result()
        result["analogs"] = [
            valid_analog("A-1", 19, "CC", selected=True)
        ]
        policy = {
            "_validated_by_platform": True,
            "scope": {"job_id": "job-1"},
            "decisions": {
                "microstate_expert_decision": {
                    "accepted": True,
                    "state_choices": {"interaction-A-1": 0},
                    "pH_conditions": {
                        "pH": 7.4,
                        "description": "Explicit synthetic assay condition",
                    },
                }
            },
        }
        digest = evaluate_design(result, policy=policy)["policy_digest"]
        supplement = {
            "_validated_by_platform": True,
            "job_id": "job-1",
            "policy_digest": digest,
            "analog_id": "A-1",
            "report_id": "interaction-A-1",
            "pH_conditions": {
                "pH": 7.4,
                "description": "Explicit synthetic assay condition",
            },
            "hydrogen_receipts": {
                "parent": {"requires_review": False},
                "pose": {"requires_review": True},
            },
            "requirements": [{
                "id": "h",
                "required": True,
                "status": "preserved",
                "computed_pass": True,
            }],
            "failures": [],
            "state_alternatives": [{
                "status": "computed",
                "hydrogen_receipts": {
                    "pose": {"requires_review": True},
                },
                "pH_conditions": {"pH": 7.4},
                "requirements": [{
                    "id": "h",
                    "required": True,
                    "computed_pass": True,
                    "unambiguous": True,
                    "pending_missing_protein_hydrogen": False,
                }],
            }],
        }

        assessment = evaluate_design(
            result,
            policy=policy,
            supplements={"interactions": [supplement]},
        )
        criterion = criterion_by_id(
            assessment,
            "microstates_h_direction",
        )

        self.assertEqual(criterion["status"], "pending")
        self.assertTrue(any(
            item.get("analog_id") == "A-1"
            and item.get("reason")
            == "selected_state_failed_or_minimizer_requires_review"
            for item in criterion["observed"]["pending"]
        ))

    def test_policy_digest_is_deterministic_and_old_public_boolean_is_untrusted(self):
        result = base_result()
        first = evaluate_design(result, policy={
            "_validated_by_platform": True,
            "numerical_criteria": {"panel_min": 1},
        })
        second = evaluate_design(result, policy={
            "numerical_criteria": {"panel_min": 1},
            "_validated_by_platform": True,
        })
        untrusted = evaluate_design(result, policy={
            "server_validated": True,
            "decisions": {"formal_expert_decision": "accept"},
        })

        self.assertEqual(first["policy_digest"], second["policy_digest"])
        self.assertIs(untrusted["policy_server_validated"], False)
        self.assertEqual(
            criterion_by_id(
                untrusted,
                "formal_expert_decision",
            )["status"],
            "pending",
        )
        json.dumps(first, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
