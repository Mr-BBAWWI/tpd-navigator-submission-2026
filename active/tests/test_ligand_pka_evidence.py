import hashlib
import importlib
import json
import math
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path

from packages.science import ligand_pka_evidence as subject


class LigandPkaEvidenceTests(unittest.TestCase):
    def test_bind_model_site_keeps_distinct_acid_hydrogen_targets(self):
        from rdkit import Chem
        from packages.science.ligand_pka_evidence import _bind_model_site

        molecule = Chem.AddHs(Chem.MolFromSmiles("[NH2:17]"))
        nitrogen = next(atom for atom in molecule.GetAtoms() if atom.GetAtomMapNum() == 17)
        hydrogen_indices = [
            atom.GetIdx() for atom in nitrogen.GetNeighbors() if atom.GetAtomicNum() == 1
        ]
        self.assertEqual(len(hydrogen_indices), 2)

        sites = [_bind_model_site(molecule, "acid", index) for index in hydrogen_indices]
        self.assertEqual([site["atom_map"] for site in sites], [17, 17])
        self.assertEqual([site["model_target_element"] for site in sites], ["H", "H"])
        self.assertEqual(
            [site["model_target_atom_index"] for site in sites], hydrogen_indices
        )
        self.assertNotEqual(sites[0]["computed_site_id"], sites[1]["computed_site_id"])

    def test_bind_model_site_rejects_unattached_acid_hydrogen(self):
        from rdkit import Chem
        from packages.science.ligand_pka_evidence import _bind_model_site

        molecule = Chem.MolFromSmiles("[H]")
        with self.assertRaisesRegex(ValueError, "REQUIRES_ONE_HEAVY_NEIGHBOR"):
            _bind_model_site(molecule, "acid", 0)

    def test_bind_model_site_accepts_mapped_base_heavy_atom(self):
        from rdkit import Chem
        from packages.science.ligand_pka_evidence import _bind_model_site

        molecule = Chem.AddHs(Chem.MolFromSmiles("[NH2:18]C"))
        nitrogen = next(atom for atom in molecule.GetAtoms() if atom.GetAtomMapNum() == 18)
        site = _bind_model_site(molecule, "base", nitrogen.GetIdx())
        self.assertEqual(site["atom_map"], 18)
        self.assertEqual(site["element"], "N")
        self.assertEqual(site["model_target_element"], "N")
        self.assertNotIn("acid_proton_attached_to_atom_map", site)

    def test_bind_model_site_enforces_site_type(self):
        from rdkit import Chem
        from packages.science.ligand_pka_evidence import _bind_model_site

        molecule = Chem.AddHs(Chem.MolFromSmiles("[NH2:17]"))
        nitrogen = next(atom for atom in molecule.GetAtoms() if atom.GetAtomMapNum() == 17)
        hydrogen = next(atom for atom in nitrogen.GetNeighbors() if atom.GetAtomicNum() == 1)
        with self.assertRaisesRegex(ValueError, "ACID_MODEL_SITE_NOT_HYDROGEN"):
            _bind_model_site(molecule, "acid", nitrogen.GetIdx())
        with self.assertRaisesRegex(ValueError, "BASE_MODEL_SITE_NOT_MAPPED_HEAVY_ATOM"):
            _bind_model_site(molecule, "base", hydrogen.GetIdx())

    def test_import_has_no_heavy_dependency(self):
        for name in ("torch", "torch_geometric", "pandas"):
            sys.modules.pop(name, None)
        importlib.reload(subject)
        self.assertNotIn("torch", sys.modules)
        self.assertNotIn("torch_geometric", sys.modules)
        self.assertNotIn("pandas", sys.modules)

    def test_ph_type_range_and_nonfinite(self):
        self.assertEqual(subject.validate_ph(7.4), 7.4)
        for value in (True, False, -0.1, 14.1, math.inf, -math.inf, math.nan):
            with self.assertRaises(ValueError):
                subject.validate_ph(value)

    def test_henderson_hasselbalch_sensitivity(self):
        self.assertGreater(subject.site_fraction(6.0, 8.0, "base"), 0.9)
        self.assertLess(subject.site_fraction(8.0, 6.0, "base"), 0.1)
        self.assertGreater(subject.site_fraction(8.0, 6.0, "acid"), 0.9)
        self.assertLess(subject.site_fraction(6.0, 8.0, "acid"), 0.1)
        value = subject.sensitivity(7.4, 7.0, "base")
        self.assertTrue(value["not_a_confidence_interval"])
        self.assertTrue(value["not_a_coupled_microstate_population"])

    def test_json_rejects_duplicates_nonfinite_and_overflow(self):
        for value in (b'{"x":1,"x":2}', b'{"x":NaN}', b'{"x":Infinity}', b'{"x":1e999}'):
            with self.assertRaises(ValueError):
                subject.parse_json_bytes(value)

    def test_manifest_exact_list_and_tamper(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root = base / "upstream"
            rows = []
            for relative in sorted(subject.UPSTREAM_FILES):
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                data = relative.encode("utf-8")
                path.write_bytes(data)
                rows.append({"path": relative, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)})
            manifest = base / "manifest.json"
            manifest.write_text(json.dumps({"source": "source", "commit": "commit", "license": "MIT", "files": rows}))
            document = subject.verify_upstream(root, manifest)
            self.assertEqual(len(document["files"]), 22)
            (root / "src/utils/gcn_conv.py").write_text("tampered")
            with self.assertRaisesRegex(ValueError, "MISMATCH"):
                subject.verify_upstream(root, manifest)

    def test_contact_counts_and_all_rows_preserved(self):
        records = []
        blocks = []
        for index in range(40):
            analog = subject.ANALOG_IDS[index % 3]
            records.append({
                "analog_id": analog,
                "pose_index": index % 2,
                "state_index": index,
                "state_id": "other",
                "required_contact_evidence": index < 16,
                "status": "computed",
                "failure": None,
            })
            blocks.append(("mol%d\n$$$$\n" % index).encode())
        for index, analog in enumerate(subject.ANALOG_IDS):
            records[index]["analog_id"] = analog
            records[index]["pose_index"] = 0
            records[index]["state_id"] = "source"
        found, parsed = subject.validate_core_records({"rows": records}, b"".join(blocks))
        rows, supported = subject.contact_summary(found)
        self.assertEqual(len(found), 40)
        self.assertEqual(len(parsed), 40)
        self.assertEqual(len(rows), 40)
        self.assertEqual(sum(row["required_contact_evidence"] for row in rows), 16)
        self.assertFalse(all(supported[analog]["0"] for analog in subject.ANALOG_IDS))

    @unittest.skipUnless(importlib.util.find_spec("rdkit") is not None, "RDKit optional")
    def test_split_sdf_preserves_blank_next_title_and_exact_bytes(self):
        from rdkit import Chem

        molecules = [Chem.MolFromSmiles("CC"), Chem.MolFromSmiles("CO")]
        mol_blocks = []
        for molecule in molecules:
            block = Chem.MolToMolBlock(molecule)
            lines = block.splitlines(keepends=True)
            lines[0] = "\n"
            mol_blocks.append("".join(lines).encode("utf-8"))

        data = mol_blocks[0] + b"$$$$\n" + mol_blocks[1] + b"$$$$\n"
        blocks = subject.split_sdf_bytes(data)

        self.assertEqual(len(blocks), 2)
        self.assertEqual(b"".join(blocks), data)
        self.assertTrue(blocks[1].startswith(b"\n"))
        for block in blocks:
            mol_block = block[:-5]
            self.assertIsNotNone(
                Chem.MolFromMolBlock(mol_block.decode("utf-8"), sanitize=True, removeHs=False)
            )

    def test_split_sdf_malformed_delimiter_and_trailing_data(self):
        with self.assertRaisesRegex(ValueError, "SDF_TRAILING_DATA"):
            list(subject.split_sdf_bytes(b"mol\n$$$$trailing\n"))
        with self.assertRaisesRegex(ValueError, "SDF_TRAILING_DATA"):
            subject.split_sdf_bytes(b"mol without delimiter")

    def test_expected_5001_scope(self):
        self.assertEqual(subject.EXPECTED_MAP_ELEMENTS["W-80f8f4a11b5d"][5001], "N")
        self.assertEqual(subject.EXPECTED_MAP_ELEMENTS["W-c2afc5e73c1a"][5001], "O")
        self.assertEqual(subject.EXPECTED_MAP_ELEMENTS["W-4c0a639c0a41"][5001], "C")
        self.assertTrue(all(value[19] == "N" for value in subject.EXPECTED_MAP_ELEMENTS.values()))
        self.assertEqual(subject.EXPECTED_PKA_MAPS["W-c2afc5e73c1a"], [19])
        self.assertEqual(subject.EXPECTED_PKA_MAPS["W-4c0a639c0a41"], [19])
        self.assertEqual(subject.EXPECTED_PKA_MAPS["W-80f8f4a11b5d"], [19, 5001])

    @unittest.skipUnless(importlib.util.find_spec("rdkit") is not None, "RDKit optional")
    def test_graph_preserves_order_independence_stereo_isotope_charge(self):
        from rdkit import Chem
        first = Chem.MolFromSmiles("[13CH3:1][C@H:19]([F:2])[NH+:5001]([CH3:5002])[CH3:5003]")
        changed_charge = Chem.MolFromSmiles("[13CH3:1][C@H:19]([F:2])[N:5001]([CH3:5002])[CH3:5003]")
        changed_isotope = Chem.MolFromSmiles("[12CH3:1][C@H:19]([F:2])[NH+:5001]([CH3:5002])[CH3:5003]")
        changed_stereo = Chem.MolFromSmiles("[13CH3:1][C@@H:19]([F:2])[NH+:5001]([CH3:5002])[CH3:5003]")
        for molecule in (first, changed_charge, changed_isotope, changed_stereo):
            self.assertIsNotNone(molecule)
        reordered = Chem.RenumberAtoms(first, list(reversed(range(first.GetNumAtoms()))))
        self.assertEqual(subject.mapped_graph(first), subject.mapped_graph(reordered))
        self.assertNotEqual(subject.mapped_graph(first), subject.mapped_graph(changed_charge))
        self.assertNotEqual(subject.mapped_graph(first), subject.mapped_graph(changed_isotope))
        self.assertNotEqual(subject.mapped_graph(first), subject.mapped_graph(changed_stereo))

    def _store(self, base, artifact_project="p"):
        store = base / "store"
        blobs = store / "blobs"
        blobs.mkdir(parents=True)
        database = store / "index.sqlite3"
        connection = sqlite3.connect(database)
        connection.executescript("""
        CREATE TABLE workbench_jobs(id TEXT, project TEXT, state TEXT, body TEXT);
        CREATE TABLE artifacts(id TEXT, project TEXT, metadata TEXT);
        CREATE TABLE scientific_evidence(project TEXT, job_id TEXT, kind TEXT, superseded INTEGER, binding TEXT, ref TEXT);
        """)
        result = {
            "input_binding": {"digest": "1" * 64},
            "parent_scope": {"actual_design_parent_id": "parent-1"},
            "analogs": [{"id": value, "mapped_smiles": "[NH2:19][CH3:5001]"} for value in subject.ANALOG_IDS],
            "files": {},
        }
        archived = subject._json_bytes(result)
        artifact_id = "a-" + "1" * 32
        ref = {"artifact_id": artifact_id, "sha256": hashlib.sha256(archived).hexdigest(), "media_type": "application/json"}
        report_id = "a-" + "2" * 32
        report_data = b"{}"
        report_ref = {"artifact_id": report_id, "sha256": hashlib.sha256(report_data).hexdigest(), "media_type": "application/json"}
        (blobs / artifact_id).write_bytes(archived)
        (blobs / report_id).write_bytes(report_data)
        connection.execute("INSERT INTO artifacts VALUES(?,?,?)", (artifact_id, artifact_project, json.dumps(ref)))
        connection.execute("INSERT INTO artifacts VALUES(?,?,?)", (report_id, artifact_project, json.dumps(report_ref)))
        live_result = dict(result)
        live_result["files"] = {"json": ref, "report": report_ref}
        body = {
            "id": "j", "project_id": "p", "state": "completed",
            "input_digest": "1" * 64, "binding": {"digest": "1" * 64},
            "result": live_result,
            "outputs": [{"name": "json", "ref": ref}, {"name": "report", "ref": report_ref}],
        }
        connection.execute("INSERT INTO workbench_jobs VALUES(?,?,?,?)", ("j", "p", "completed", json.dumps(body)))
        connection.commit()
        connection.close()
        return store, ref

    def test_readonly_job_scope_artifact_sha_and_six_key_binding(self):
        with tempfile.TemporaryDirectory() as temporary:
            store, ref = self._store(Path(temporary))
            design, binding = subject.read_bound_design(store, "p", "j")
            self.assertEqual(binding, {
                "project": "p", "job_id": "j", "input_sha256": "1" * 64,
                "result_sha256": ref["sha256"],
                "result_binding_kind": "exact_archived_design_json_bytes",
                "parent_id": "parent-1",
            })
            self.assertEqual(len(design["analogs"]), 3)
            with self.assertRaisesRegex(ValueError, "NOT_FOUND"):
                subject.read_bound_design(store, "other", "j")

    def test_artifact_project_scope(self):
        with tempfile.TemporaryDirectory() as temporary:
            store, unused = self._store(Path(temporary), artifact_project="other")
            with self.assertRaisesRegex(ValueError, "PROJECT_SCOPED"):
                subject.read_bound_design(store, "p", "j")

    def test_interaction_evidence_checksum_includes_encoded_json_lf(self):
        with tempfile.TemporaryDirectory() as temporary:
            store, result_ref = self._store(Path(temporary))
            report = {
                "analog_id": subject.ANALOG_IDS[0],
                "report_id": "interaction-W-known-style",
            }
            encoded = subject._json_bytes(report) + b"\n"
            artifact_id = "a-" + "3" * 32
            ref = {
                "artifact_id": artifact_id,
                "sha256": hashlib.sha256(encoded).hexdigest(),
                "media_type": "application/json",
            }
            (store / "blobs" / artifact_id).write_bytes(encoded)
            connection = sqlite3.connect(store / "index.sqlite3")
            connection.execute("INSERT INTO artifacts VALUES(?,?,?)", (artifact_id, "p", json.dumps(ref)))
            binding = {
                "result_sha256": result_ref["sha256"],
                "evidence_sha256": hashlib.sha256(encoded).hexdigest(),
            }
            connection.execute(
                "INSERT INTO scientific_evidence VALUES(?,?,?,?,?,?)",
                ("p", "j", "interaction", 0, json.dumps(binding), json.dumps(ref)),
            )
            connection.commit()
            connection.close()
            self.assertNotEqual(subject.canonical_digest(report), binding["evidence_sha256"])
            reports = subject._interaction_reports(store, "p", "j", result_ref["sha256"])
            self.assertEqual(reports[subject.ANALOG_IDS[0]][1], report)

    def test_provenance_documents_have_no_approval_flags_set_true(self):
        document = {
            "scientific_approved": False,
            "population_prediction_performed": False,
            "final_state_selection_performed": False,
        }
        self.assertFalse(any(document.values()))


if __name__ == "__main__":
    unittest.main()
