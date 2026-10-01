import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rdkit import Chem

from packages.platform import design_panel
from packages.science import acceptance_evidence
from packages.science.dual_e3 import validate


class AcceptanceEvidenceTests(unittest.TestCase):
    def _tree(self, tamper=False):
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name) / "cases" / "acceptance_sources"
        root.mkdir(parents=True)
        files = {}
        for name in acceptance_evidence.EXPECTED_FILES:
            raw = json.dumps({"name": name, "value": 1}).encode()
            (root / name).write_bytes(raw + (b"x" if tamper and name == "route_evidence.json" else b""))
            files[name] = hashlib.sha256(raw).hexdigest()
        (root / "manifest.json").write_text(json.dumps({"version": "test", "files": files}))
        return temporary, root

    def test_missing_evidence_is_explicit(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(acceptance_evidence, "SOURCE", Path(directory) / "missing"), mock.patch.object(acceptance_evidence, "MANIFEST", Path(directory) / "missing" / "manifest.json"):
            result = acceptance_evidence.load_acceptance_evidence()
            self.assertEqual(result["status"], "not_configured")
            self.assertTrue(all(row["status"] == "not_configured" for row in result["sources"].values()))

    def test_tamper_fails(self):
        temporary, root = self._tree(tamper=True)
        try:
            with mock.patch.object(acceptance_evidence, "SOURCE", root), mock.patch.object(acceptance_evidence, "MANIFEST", root / "manifest.json"):
                with self.assertRaisesRegex(ValueError, "HASH_MISMATCH"):
                    acceptance_evidence.load_acceptance_evidence()
        finally:
            temporary.cleanup()

    def test_declared_missing_file_fails_closed(self):
        temporary, root = self._tree()
        try:
            (root / "preparation.json").unlink()
            with mock.patch.object(acceptance_evidence, "SOURCE", root), \
                    mock.patch.object(acceptance_evidence, "MANIFEST", root / "manifest.json"):
                with self.assertRaisesRegex(ValueError, "SOURCE_MISSING"):
                    acceptance_evidence.load_acceptance_evidence()
        finally:
            temporary.cleanup()

    def test_binary_supporting_file_is_hash_bound_without_json_decode(self):
        temporary, root = self._tree()
        try:
            support = b"\x00\xffSDF/raw evidence\n"
            relative = "raw_rcsb/6HAZ.cif.gz"
            path = root / relative
            path.parent.mkdir()
            path.write_bytes(support)
            manifest = json.loads((root / "manifest.json").read_text())
            manifest["files"][relative] = hashlib.sha256(support).hexdigest()
            self.assertNotIn("supporting_files", manifest)
            (root / "manifest.json").write_text(json.dumps(manifest))
            with mock.patch.object(acceptance_evidence, "SOURCE", root), \
                    mock.patch.object(acceptance_evidence, "MANIFEST", root / "manifest.json"):
                loaded = acceptance_evidence.load_acceptance_evidence()
                binding = acceptance_evidence.source_binding()
            self.assertEqual(
                loaded["supporting_files"][relative]["status"],
                "configured_hash_verified",
            )
            self.assertEqual(
                binding["supporting_files"][relative]["sha256"],
                hashlib.sha256(support).hexdigest(),
            )
            self.assertEqual(
                set(binding["files"]),
                set(acceptance_evidence.EXPECTED_FILES) | {relative},
            )
            self.assertEqual(binding["files"][relative]["sha256"], hashlib.sha256(support).hexdigest())
        finally:
            temporary.cleanup()

    def test_unsafe_supporting_paths_are_rejected(self):
        for relative in (
            "../receipt.txt", "C:/receipt.txt", "C:\\receipt.txt",
            "//server/share/file", "receipts\\file.txt", "a:b",
        ):
            temporary, root = self._tree()
            try:
                manifest = json.loads((root / "manifest.json").read_text())
                manifest["files"][relative] = "0" * 64
                (root / "manifest.json").write_text(json.dumps(manifest))
                with mock.patch.object(acceptance_evidence, "SOURCE", root), \
                        mock.patch.object(acceptance_evidence, "MANIFEST", root / "manifest.json"):
                    with self.assertRaisesRegex(ValueError, "SOURCE_PATH"):
                        acceptance_evidence.load_acceptance_evidence()
            finally:
                temporary.cleanup()

    def test_legacy_and_canonical_manifest_entries_must_be_disjoint(self):
        temporary, root = self._tree()
        try:
            relative = "supporting/raw/source.sdf"
            path = root / relative
            path.parent.mkdir(parents=True)
            path.write_bytes(b"source")
            expected = hashlib.sha256(b"source").hexdigest()
            manifest = json.loads((root / "manifest.json").read_text())
            manifest["files"][relative] = expected
            manifest["supporting_files"] = {relative: expected}
            (root / "manifest.json").write_text(json.dumps(manifest))
            with mock.patch.object(acceptance_evidence, "SOURCE", root), \
                    mock.patch.object(acceptance_evidence, "MANIFEST", root / "manifest.json"):
                with self.assertRaisesRegex(ValueError, "ACCEPTANCE_MANIFEST_SCHEMA"):
                    acceptance_evidence.load_acceptance_evidence()
        finally:
            temporary.cleanup()

    def test_all_eighteen_linker_ids_are_accepted(self):
        from packages.science.linker_assessment import expanded_library
        ids = [row["id"] for row in expanded_library()["templates"]]
        self.assertEqual(len(ids), 18)
        self.assertEqual(validate({"linker_ids": ids, "panel_size": 10})["linker_ids"], ids)


class PipelinePolicyTests(unittest.TestCase):
    def _run(self, dock_enabled=True, parent_pass=True):
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)
        source = root / "cases" / "design_sources"
        source.mkdir(parents=True)
        contracts = root / "contracts" / "drafts"
        contracts.mkdir(parents=True)
        (contracts / "design_panel.schema.json").write_text("{}")

        parent = Chem.MolFromSmiles("[CH3:1][OH:2]")
        conformer = Chem.Conformer(parent.GetNumAtoms())
        conformer.Set3D(True)
        for index in range(parent.GetNumAtoms()):
            conformer.SetAtomPosition(index, (float(index), 0.0, 0.0))
        parent.AddConformer(conformer)
        (source / "SMARCA2-neutral-design.sdf").write_text(
            Chem.MolToMolBlock(parent)
        )
        (source / "6HAZ.cif").write_text("mock")

        def row(identifier, smiles):
            return {
                "id": identifier,
                "mapped_smiles": smiles,
                "canonical_smiles": smiles,
                "rule_id": "LH_OH",
                "transformation_class": "ring_expansion",
                "tier": 1,
                "modified_atom_maps": [1],
                "added_atom_maps": [],
                "cheap_filter": {"valid": True},
                "stereochemistry": {"unresolved_stereo_count": 0},
            }

        rows = [
            row("W-pass", "[CH3:1][OH:2]"),
            row("W-fail", "[CH3:1][NH2:2]"),
        ]
        generated = {
            "analogs": copy.deepcopy(rows),
            "rejections": [],
            "rule_catalog": [{
                "rule_id": "LH_OH",
                "transformation_class": "ring_expansion",
            }],
        }
        evidence_sources = {
            name: {
                "status": "not_configured",
                "data": None,
                "relative_path": name,
                "sha256": None,
            }
            for name in acceptance_evidence.EXPECTED_FILES
        }
        catalog = {
            "warhead": {
                "id": "parent",
                "sar": {},
                "protected_maps": [1],
            },
            "acceptance_evidence": {
                "status": "not_configured",
                "manifest": {"version": None},
                "sources": evidence_sources,
            },
            "linkers": {"templates": [{"id": "L1"}]},
            "recruiters": {"recruiters": [{"e3_type": "CRBN"}]},
        }
        calls = []
        dock_results = [
            {
                "status": "completed_with_limits",
                "files": {},
                "results": {"pose_preservation": {
                    "docking_pose_preserved": True,
                    "all_poses_diagnostics": [{
                        "status": "pass",
                        "core_RMSD_A_in_receptor_frame": 0.5,
                    }],
                }},
            },
            {
                "status": "completed_with_limits",
                "files": {},
                "results": {"pose_preservation": {
                    "docking_pose_preserved": False,
                    "all_poses_diagnostics": [{
                        "status": "fail",
                        "core_RMSD_A_in_receptor_frame": 8.0,
                    }],
                }},
            },
            {
                "status": "completed_with_limits",
                "files": {},
                "results": {"pose_preservation": {
                    "docking_pose_preserved": parent_pass,
                    "all_poses_diagnostics": [{
                        "status": "pass" if parent_pass else "fail",
                        "core_RMSD_A_in_receptor_frame": 0.4,
                    }],
                }},
            },
        ]

        def fake_dock(molecule, parent_molecule, protein, protected, folder,
                      check_active=None):
            folder.mkdir(parents=True, exist_ok=True)
            calls.append(("dock", Chem.MolToSmiles(molecule)))
            return copy.deepcopy(dock_results[len(calls) - 1])

        def fake_select(pool, limit):
            calls.append(("select", [record["id"] for record in pool]))
            selected = []
            for record in pool[:limit]:
                selected_row = copy.deepcopy(record)
                selected_row["cluster_id"] = 7
                selected_row["cluster_size"] = 1
                selected.append(selected_row)
            return selected

        assembled = []
        def fake_assemble(analog, recruiter, linker):
            assembled.append(analog["id"])
            return {
                "candidate_id": "C-" + analog["id"],
                "mapped_smiles": analog["mapped_smiles"],
                "e3_type": recruiter["e3_type"],
                "atom_roles": {"linker_maps": []},
            }

        def put(name, data, content_type=None):
            return {"name": name}

        parameters = {
            "dock": dock_enabled,
            "panel_size": 2,
            "exploratory": False,
            "use_api": False,
            "linker_ids": ["L1"],
        }
        validated = dict(parameters)

        patches = [
            mock.patch.object(design_panel, "SOURCE", source),
            mock.patch.object(design_panel, "validate", return_value=validated),
            mock.patch.object(design_panel, "catalog", return_value=catalog),
            mock.patch.object(design_panel, "input_binding", return_value={
                "format": "test", "source_files": {}, "digest": "0" * 64,
            }),
            mock.patch.object(design_panel, "attachment_options", return_value=[{"map": 2}]),
            mock.patch.object(design_panel, "assemble", side_effect=fake_assemble),
            mock.patch.object(design_panel, "_stereoisomers", return_value=copy.deepcopy(rows)),
            mock.patch.object(design_panel, "_source_parent_profile", return_value={
                "interactions": [],
            }),
            mock.patch.object(design_panel, "_pose_directional_profiles"),
            mock.patch.object(design_panel, "_scientific_context", return_value={
                "ligand_hydrogen_preparation": {"status": "mock"},
            }),
            mock.patch.object(design_panel, "_calibration", return_value={
                "CRBN": {"status": "not_configured"},
            }),
            mock.patch.object(design_panel, "_files", return_value={}),
            mock.patch("packages.science.structures.atom_sites", return_value=[]),
            mock.patch("packages.science.warhead_sites.analyze_sites", return_value={
                "atoms": [{"atom_map": 1, "state": "PROTECTED"}, {
                    "atom_map": 2, "state": "MODIFIABLE",
                }],
            }),
            mock.patch("packages.science.warhead_sites.funnel", return_value=[]),
            mock.patch("packages.science.analog_generation.generate", return_value=generated),
            mock.patch("packages.science.analog_generation.select_diverse", side_effect=fake_select),
            mock.patch("packages.science.analog_filters.cheap_filter", return_value={"valid": True}),
            mock.patch("packages.science.analog_filters.count_by_site_and_broad_family", return_value={}),
            mock.patch("packages.science.design_docking.dock", side_effect=fake_dock),
            mock.patch("packages.science.molecules.identity", return_value={"properties": {}}),
            mock.patch("packages.science.e3_benchmark.reference_packet", return_value={}),
            mock.patch("packages.science.linker_assessment.assess_linker", return_value={}),
            mock.patch("packages.science.linker_assessment.route_evidence", return_value={}),
            mock.patch("jsonschema.Draft202012Validator"),
        ]
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)
        self.addCleanup(temporary.cleanup)

        result = design_panel.run_panel(parameters, put)
        return result, calls, assembled

    def test_actual_pipeline_docks_every_eligible_before_selection(self):
        result, calls, assembled = self._run()
        select_index = next(
            index for index, call in enumerate(calls) if call[0] == "select"
        )
        self.assertEqual(
            [call[0] for call in calls[:select_index]],
            ["dock", "dock", "dock"],
        )
        self.assertEqual(calls[select_index], ("select", ["W-pass"]))
        self.assertEqual(result["summary"]["docking_attempted"], 2)
        self.assertEqual(result["summary"]["qualified_analogs"], 1)
        self.assertEqual(assembled, ["W-pass"])

        by_id = {row["id"]: row for row in result["analogs"]}
        self.assertTrue(by_id["W-pass"]["selected"])
        self.assertEqual(by_id["W-pass"]["cluster_id"], 7)
        self.assertFalse(by_id["W-fail"]["selected"])
        self.assertEqual(by_id["W-fail"]["pipeline_status"], "failed_pose_filter")
        self.assertIn("docking", by_id["W-fail"])

    def test_parent_redock_failure_blocks_assembly(self):
        result, calls, assembled = self._run(parent_pass=False)
        select_call = next(call for call in calls if call[0] == "select")
        self.assertEqual(select_call, ("select", []))
        self.assertEqual(result["summary"]["qualified_analogs"], 0)
        self.assertEqual(result["summary"]["protac_count"], 0)
        self.assertEqual(assembled, [])

    def test_no_dock_has_zero_qualified_and_explicit_no_assembly(self):
        result, calls, assembled = self._run(dock_enabled=False)
        self.assertFalse(any(call[0] == "dock" for call in calls))
        self.assertEqual(result["summary"]["qualified_analogs"], 0)
        self.assertEqual(result["summary"]["docking_attempted"], 0)
        self.assertEqual(result["summary"]["protac_count"], 0)
        self.assertEqual(assembled, [])
        self.assertEqual(
            result["summary"]["stage_counts"]["assemble_qualified"]["status"],
            "not_run_docking_disabled_explicit_no_assembly",
        )
        self.assertTrue(
            all("docking" in row for row in result["analogs"])
        )


if __name__ == "__main__":
    unittest.main()
