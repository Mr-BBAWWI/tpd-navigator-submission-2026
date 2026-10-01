from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from packages.science import novel_ternary as novel


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _ref(media_type="chemical/x-mdl-sdfile"):
    return {"artifact_id": "a-test", "version": 1, "sha256": "a" * 64,
            "media_type": media_type, "schema_id": "urn:tpd-navigator:raw:1",
            "provenance": "computed"}


def _analog():
    return {"mapped_smiles": "[CH3:1][CH2:2][OH:3]", "canonical_smiles": "CCO",
            "id": "WA-1", "pipeline_status": "qualified", "selected": True,
            "assembly_eligible": True, "parent_redocking_supported": True,
            "docking": {"status": "completed_with_limits", "pose_preserved": True,
                        "passing_pose_count": 1, "files": {"poses_sdf": _ref()}}}


def _qualification():
    analog = _analog()
    return {"id": analog["id"], "selected": analog["selected"],
            "pipeline_status": analog["pipeline_status"],
            "assembly_eligible": analog["assembly_eligible"],
            "parent_redocking_supported": analog["parent_redocking_supported"],
            "docking": {"status": analog["docking"]["status"],
                        "pose_preserved": analog["docking"]["pose_preserved"],
                        "passing_pose_count": analog["docking"]["passing_pose_count"],
                        "files": {"poses_sdf": analog["docking"]["files"]["poses_sdf"]}}}


def _source_cif(path: Path, deposit: str, sequence: str, description: str, chain: str,
                extras=()) -> None:
    three = {value: key for key, value in novel.worker.THREE_TO_ONE.items()}
    entities = [("1", "polymer", description)]
    asym = [(chain, "1")]
    polymers = [("1", "polypeptide(L)", sequence)]
    for index, (extra_chain, extra_description) in enumerate(extras, 2):
        entities.append((str(index), "polymer", extra_description))
        asym.append((extra_chain, str(index)))
        polymers.append((str(index), "polypeptide(L)", "AAA"))
    lines = ["data_test", f"_entry.id {deposit}", "loop_", "_entity.id", "_entity.type",
             "_entity.pdbx_description",
             *[f"{eid} {kind} '{desc}'" for eid, kind, desc in entities],
             "loop_", "_entity_poly.entity_id", "_entity_poly.type",
             "_entity_poly.pdbx_seq_one_letter_code_can",
             *[f"{eid} '{kind}' {seq}" for eid, kind, seq in polymers],
             "loop_", "_struct_asym.id", "_struct_asym.entity_id",
             *[f"{item_chain} {eid}" for item_chain, eid in asym],
             "loop_", "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
             "_atom_site.label_atom_id", "_atom_site.label_alt_id", "_atom_site.label_comp_id",
             "_atom_site.label_asym_id", "_atom_site.label_entity_id", "_atom_site.label_seq_id",
             "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
             "_atom_site.occupancy", "_atom_site.auth_asym_id", "_atom_site.pdbx_PDB_model_num"]
    for atom_id, residue in enumerate(sequence, 1):
        lines.append(f"ATOM {atom_id} C CA . {three[residue]} {chain} 1 {atom_id} "
                     f"{atom_id}.0 0.0 0.0 1.0 {chain} 1")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _candidate():
    return {"candidate_id": "D-test", "e3_type": "VHL",
            "assembly_mode": "pose_supported_hypothesis", "warhead_analog_id": "WA-1",
            "analog_qualification": _qualification(), "canonical_smiles": "CCCO",
            "mapped_smiles": "[CH3:1][CH2:2001][CH2:2002][OH:4001]",
            "atom_roles": {"warhead_maps": [1], "linker_maps": [2001, 2002],
                           "recruiter_maps": [4001]},
            "attachment_metadata": {
                "warhead_linker_bond": {"attachment_atom_map": 1,
                                         "partner_atom_map": 2001, "bond_type": "SINGLE"},
                "recruiter_linker_bond": {"attachment_atom_map": 4001,
                                           "partner_atom_map": 2002, "bond_type": "SINGLE"}}}


def _chem_block(path: Path, names, orders):
    lines = ["data_prediction", "loop_", "_chem_comp_bond.comp_id",
             "_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2",
             "_chem_comp_bond.value_order"]
    for left, right, order in orders:
        lines.append(f"LIG {names[left]} {names[right]} {order}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return novel._source_block(path)


class NovelTernaryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def plan(self):
        target, e3 = self.root / "6haz.cif", self.root / "6hay.cif"
        _source_cif(target, "6HAZ", "A" * 123,
                    "Probable global transcription activator SNF2L2 SMARCA2", "A")
        _source_cif(e3, "6HAY", "C" * 162,
                    "Von Hippel-Lindau disease tumor suppressor VHL", "B",
                    (("C", "Elongin-C"), ("D", "Elongin-B")))
        graph = novel._mapped_graph(_candidate())
        return novel.prepare_plan(
            _candidate(),
            target_source={"path": str(target), "label_asym_id": "A",
                           "expected_sha256": _sha(target), "role": "target",
                           "description_terms": ["SMARCA2"]},
            e3_source={"path": str(e3), "label_asym_id": "B",
                       "expected_sha256": _sha(e3), "role": "e3",
                       "description_terms": ["VHL"]},
            job_binding={"project": "p", "job_id": "j", "input_sha256": "1" * 64,
                         "result_sha256": "3" * 64,
                         "candidate_graph_sha256": graph["graph_sha256"]},
            policy_binding={"module": "policy", "revision": "r1", "digest": "2" * 64},
            seeds=[3, 7])

    def test_plan_qualification_binding_and_cofactor_normalization(self):
        plan = self.plan()
        novel.verify_plan(plan)
        self.assertTrue(plan["candidate_graph"]["canonical_full_graph"]["derived"])
        self.assertEqual("qualified",
                         plan["candidate_graph"]["analog_qualification"]["pipeline_status"])
        self.assertNotEqual(plan["bindings"]["job"]["result_sha256"],
                            plan["bindings"]["job"]["candidate_graph_sha256"])
        self.assertEqual({"elongin_C", "elongin_B"},
                         {item["category"] for item in
                          plan["sources"]["e3"]["excluded_proteins"]})

    def test_qualification_exact_snapshot_revalidated(self):
        plan = self.plan()
        tampered = copy.deepcopy(plan)
        tampered["candidate_graph"]["analog_qualification"]["assembly_eligible"] = False
        unsigned = dict(tampered)
        unsigned.pop("plan_digest")
        tampered["plan_digest"] = novel._digest(unsigned)
        with self.assertRaisesRegex(novel.NovelTernaryError, "ANALOG_NOT_ASSEMBLY_ELIGIBLE"):
            novel.verify_plan(tampered)
        extra = _candidate()
        extra["analog_qualification"]["invented"] = True
        with self.assertRaisesRegex(novel.NovelTernaryError, "FIELDS_INVALID"):
            novel._mapped_graph(extra)

    def test_graph_is_derived_and_attachment_checks_remain(self):
        candidate = _candidate()
        candidate["mapped_smiles"] = "[CH2:1]=[CH:2001][CH2:2002][OH:4001]"
        candidate["canonical_smiles"] = "C=CCO"
        with self.assertRaisesRegex(novel.NovelTernaryError, "BOND_TYPE_GRAPH_MISMATCH"):
            novel._mapped_graph(candidate)
        candidate = _candidate()
        candidate["atom_roles"]["warhead_maps"] = [True]
        with self.assertRaisesRegex(novel.NovelTernaryError, "ATOM_ROLE_MAPS_INVALID"):
            novel._mapped_graph(candidate)

    def test_source_deposit_and_binding_schema_fail(self):
        plan = self.plan()
        bad = copy.deepcopy(plan)
        bad["bindings"]["job"].pop("candidate_graph_sha256")
        unsigned = dict(bad)
        unsigned.pop("plan_digest")
        bad["plan_digest"] = novel._digest(unsigned)
        with self.assertRaisesRegex(novel.NovelTernaryError, "JOB_BINDING_FIELDS_INVALID"):
            novel.verify_plan(bad)
        source = Path(plan["sources"]["target"]["path"])
        source.write_text(source.read_text().replace("_entry.id 6HAZ", "_entry.id 6HAY"),
                          encoding="utf-8")
        plan["sources"]["target"]["expected_sha256"] = _sha(source)
        with self.assertRaisesRegex(novel.NovelTernaryError, "SOURCE_DEPOSIT_ID_MISMATCH"):
            novel._validate_source(plan["sources"]["target"], "target")

    def test_cli_actual_shape_and_exact_frozen_bytes_digest(self):
        target, e3 = self.root / "6haz.cif", self.root / "6hay.cif"
        _source_cif(target, "6HAZ", "A" * 123, "SMARCA2 SNF2L2", "A")
        _source_cif(e3, "6HAY", "C" * 162, "Von Hippel-Lindau VHL", "B")
        candidate = _candidate()
        candidate.pop("analog_qualification")
        design = {"format": "design-panel/20260930.4",
                  "input_binding": {"format": "design-panel/20260930.4",
                                    "source_files": [], "digest": "a" * 64},
                  "analogs": [_analog()], "protac_candidates": [candidate]}
        design_path, output = self.root / "design.json", self.root / "plan.json"
        raw = json.dumps(design, indent=2).encode("utf-8") + b"\n"
        design_path.write_bytes(raw)
        script = Path(__file__).resolve().parents[1] / "scripts" / "run_novel_ternary.py"
        result = subprocess.run([
            sys.executable, str(script), "--design-json", str(design_path),
            "--candidate-id", "D-test", "--target-cif", str(target), "--target-chain", "A",
            "--e3-cif", str(e3), "--e3-chain", "B", "--job-id", "declared-job",
            "--project", "declared-project", "--output", str(output),
            "--seeds", " 3, 7 ", "--prepare-only"], capture_output=True, text=True,
            cwd=self.root)
        self.assertEqual(0, result.returncode, result.stderr)
        plan = json.loads(output.read_text())
        self.assertEqual(hashlib.sha256(raw).hexdigest(),
                         plan["bindings"]["job"]["result_sha256"])
        self.assertEqual("a" * 64, plan["bindings"]["job"]["input_sha256"])
        self.assertEqual("WA-1", plan["candidate_graph"]["analog_qualification"]["id"])

    def test_cli_rejects_foreign_and_unqualified_analog(self):
        from scripts import run_novel_ternary as cli
        document = {"protac_candidates": [_candidate()], "analogs": []}
        with self.assertRaisesRegex(ValueError, "WARHEAD_ANALOG_ID_NOT_UNIQUE"):
            cli._candidate(document, "D-test")
        with self.assertRaisesRegex(ValueError, "CANDIDATE_ID_NOT_UNIQUE"):
            cli._candidate(document, "foreign")
        document["analogs"] = [_analog()]
        document["analogs"][0]["pipeline_status"] = "failed"
        with self.assertRaisesRegex(ValueError, "PIPELINE_NOT_QUALIFIED"):
            cli._candidate(document, "D-test")

    def test_aromatic_all_single_rejected_and_kekule_accepted(self):
        candidate = _candidate()
        candidate["canonical_smiles"] = "c1ccccc1CCO"
        candidate["mapped_smiles"] = (
            "[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1"
            "[CH2:2001][CH2:2002][OH:4001]"
        )
        candidate["atom_roles"] = {"warhead_maps": [1, 2, 3, 4, 5, 6],
                                   "linker_maps": [2001, 2002],
                                   "recruiter_maps": [4001]}
        candidate["attachment_metadata"] = {
            "warhead_linker_bond": {"attachment_atom_map": 6,
                                     "partner_atom_map": 2001,
                                     "bond_type": "SINGLE"},
            "recruiter_linker_bond": {"attachment_atom_map": 4001,
                                       "partner_atom_map": 2002,
                                       "bond_type": "SINGLE"}}
        graph = novel._mapped_graph(candidate)
        names = {int(key): value for key, value in graph["atom_map_to_predicted_atom_name"].items()}
        all_single = [(1, 2, "SING"), (2, 3, "SING"), (3, 4, "SING"),
                      (4, 5, "SING"), (5, 6, "SING"), (6, 1, "SING"),
                      (6, 2001, "SING"), (2001, 2002, "SING"),
                      (2002, 4001, "SING")]
        block = _chem_block(self.root / "single.cif", names, all_single)
        with self.assertRaisesRegex(novel.NovelTernaryError, "BOND_ORDER_MISMATCH"):
            novel._validate_chem_comp_bonds(block, graph, "LIG")
        kekule = [(1, 2, "DOUB"), (2, 3, "SING"), (3, 4, "DOUB"),
                  (4, 5, "SING"), (5, 6, "DOUB"), (6, 1, "SING"),
                  (6, 2001, "SING"), (2001, 2002, "SING"),
                  (2002, 4001, "SING")]
        block = _chem_block(self.root / "kekule.cif", names, kekule)
        self.assertEqual("semantic_match",
                         novel._validate_chem_comp_bonds(block, graph, "LIG"))

    def test_interface_clashes_full_count_bounded_examples(self):
        proteins = {}
        for index in range(300):
            proteins[("target", index, "C")] = {"xyz": [0.0, 0.0, 0.0], "element": "C"}
        for index in range(300):
            proteins[("e3", index, "C")] = {"xyz": [0.0, 0.0, 0.0], "element": "C"}
        result = novel._interface_clashes(proteins)
        self.assertEqual(90000, result["count"])
        self.assertEqual(256, len(result["examples"]))
        self.assertTrue(result["examples_truncated"])

    def test_run_setup_failures_are_all_seed_failures(self):
        plan = self.plan()
        result = novel.run_plan(plan, output_root=self.root / "runs",
                                executable=self.root / "missing-boltz",
                                checkpoint=self.root / "missing-checkpoint",
                                cache=self.root / "cache", timeout=1)
        self.assertEqual(0, result["completed_count"])
        self.assertEqual(2, result["failed_count"])
        self.assertFalse(result["actual_computation"])

    def test_cancellation_is_not_scientific_failure(self):
        plan = self.plan()
        executable, checkpoint = self.root / "boltz", self.root / "checkpoint"
        executable.write_text("stub", encoding="utf-8")
        checkpoint.write_text("stub", encoding="utf-8")
        environment = {"boltz_version": "2.2.1",
                       "rdkit_version": plan["rdkit"]["planning_version"]}
        with mock.patch.object(novel, "_probe_boltz_environment", return_value=environment), \
             mock.patch.object(novel.worker, "inspect_help", return_value={"has_seed": True}):
            result = novel.run_plan(plan, output_root=self.root / "cancelled",
                                    executable=executable, checkpoint=checkpoint,
                                    cache=self.root / "cache", cancel=lambda: True)
        self.assertEqual(len(plan["seeds"]), result["cancelled_count"])
        self.assertEqual(0, result["failed_count"])


if __name__ == "__main__":
    unittest.main()
