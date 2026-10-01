"""Offline unittest coverage for native deposited-parent docking integrity rules."""

from __future__ import annotations

import hashlib
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem

from packages.science import native_parent_docking as native


def _mapped_halogen(symbol="Cl"):
    mol = Chem.MolFromSmiles(f"[CH3:1][{symbol}:8]")
    if mol is None:
        raise AssertionError("test molecule failed to parse")
    return mol


def _protein():
    return [{
        "label_asym_id": "A", "label_comp_id": "GLY", "label_atom_id": "CA",
        "type_symbol": "C", "auth_seq_id": "1380",
        "pdbx_PDB_ins_code": None, "occupancy": "1.0",
        "B_iso_or_equiv": "12.0", "xyz": [1.0, 2.0, 3.0],
    }]


def _metadata():
    return {
        "pdb": "9D12", "ccd": "A1A1P",
        "source_record_id": "9D12:E:A1A1P", "source_chain": "A",
        "atom_mapping": [
            {"ccd_atom_id": "C1", "atom_map": 1, "rdkit_index_zero_based": 0},
            {"ccd_atom_id": "CL", "atom_map": 8, "rdkit_index_zero_based": 1},
        ],
    }


def _ligand_sites(second_name="CL"):
    return [
        {"id": "1", "label_asym_id": "E", "label_comp_id": "A1A1P",
         "label_atom_id": "C1", "type_symbol": "C", "auth_seq_id": "1",
         "xyz": [4.0, 5.0, 6.0]},
        {"id": "2", "label_asym_id": "E", "label_comp_id": "A1A1P",
         "label_atom_id": second_name, "type_symbol": "Cl", "auth_seq_id": "1",
         "xyz": [7.0, 8.0, 9.0]},
    ]


def _lys(distance, complete=False):
    names = ["N", "CA", "C", "O", "CB", "CG"]
    if complete:
        names += ["CD", "CE", "NZ"]
    return [{
        "label_asym_id": "A", "label_comp_id": "LYS", "label_atom_id": name,
        "type_symbol": "N" if name in {"N", "NZ"} else "O" if name == "O" else "C",
        "auth_seq_id": "1385", "pdbx_PDB_ins_code": None,
        "xyz": [distance, index * 0.01, 0.0],
    } for index, name in enumerate(names)]


class NativeParentDockingTests(unittest.TestCase):
    def test_duplicate_seeds_and_existing_output_are_refused(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            native.validate_seeds([23, 23, 61])
        with tempfile.TemporaryDirectory() as root:
            output = Path(root) / "existing"
            output.mkdir()
            with self.assertRaises(FileExistsError):
                native.run_native_probe("unused.cif", output, seeds=[23, 41, 61])

    def test_exact_mapped_cl_to_br_graph_only(self):
        native.validate_halogen_pair(_mapped_halogen("Cl"), _mapped_halogen("Br"))
        wrong = Chem.MolFromSmiles("[CH2:1]([F:2])[Br:8]")
        with self.assertRaises(ValueError):
            native.validate_halogen_pair(_mapped_halogen("Cl"), wrong)

    def test_native_selection_uses_exact_mapping_and_coordinates(self):
        positioned, receipt = native._validate_native_selection(
            _mapped_halogen(), _metadata(), _ligand_sites()
        )
        conformer = positioned.GetConformer()
        self.assertEqual(list(conformer.GetAtomPosition(0)), [4.0, 5.0, 6.0])
        self.assertEqual(list(conformer.GetAtomPosition(1)), [7.0, 8.0, 9.0])
        self.assertEqual([row["atom_map"] for row in receipt], [1, 8])
        self.assertEqual(native._canonical_mapped(positioned),
                         native._canonical_mapped(_mapped_halogen()))
        with self.assertRaisesRegex(ValueError, "mapping differ"):
            native._validate_native_selection(
                _mapped_halogen(), _metadata(), _ligand_sites("WRONG")
            )

    def test_cif_hash_mismatch_precedes_atom_selection(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            cif = root / "9D12.cif"
            cif.write_text("not used", encoding="ascii")
            metadata = {"source_input_hashes": {"source_pdb_mmcif_sha256": "0" * 64}}
            with patch("packages.science.reference_parents.load_reference_parent",
                       return_value=(_mapped_halogen(), metadata)), patch(
                "packages.science.structures.atom_sites",
                side_effect=AssertionError("atom_sites must not run"),
            ):
                with self.assertRaisesRegex(ValueError, "hash differs"):
                    native.load_native_context(
                        "SMARCA2-9D12-A1A1P", cif, root / "output"
                    )

    def test_source_pair_uses_actual_verified_source_schema(self):
        parent_id = "SMARCA2-9D12-A1A1P"
        source = {"status": "configured_hash_verified", "sha256": "f" * 64,
                  "data": {"pairs": []}}
        receipt = {"result": {
            "selected_parent_measured_evidence": {"selected_parent_id": parent_id},
            "source_catalog": {"sources": {"medchem_evidence.json": source}},
        }}
        with tempfile.TemporaryDirectory() as root:
            export = Path(root) / "export"
            export.mkdir()
            (export / "manifest.json").write_text("{}", encoding="utf-8")
            with patch("packages.science.parent_sar_probe.verify_parent_export",
                       return_value=receipt) as verifier:
                verified, selected, actual, reference = native.verify_source_pair(
                    export, parent_id
                )
            verifier.assert_called_once_with(export / "manifest.json", export.parent)
            self.assertIs(actual, source)
            self.assertIs(verified["result"], receipt["result"])
            self.assertEqual(selected["selected_parent_id"], parent_id)
            self.assertEqual(reference["sha256"], "f" * 64)

    def test_source_pair_refuses_fake_legacy_boolean_flag(self):
        parent_id = "SMARCA2-9D12-A1A1P"
        source = {"configured_hash_verified": True, "sha256": "f" * 64, "data": {}}
        receipt = {"result": {
            "selected_parent_measured_evidence": {"selected_parent_id": parent_id},
            "source_catalog": {"sources": {"medchem_evidence.json": source}},
        }}
        with tempfile.TemporaryDirectory() as root:
            export = Path(root) / "export"
            export.mkdir()
            (export / "manifest.json").write_text("{}", encoding="utf-8")
            with patch("packages.science.parent_sar_probe.verify_parent_export",
                       return_value=receipt):
                with self.assertRaisesRegex(ValueError, "not hash verified"):
                    native.verify_source_pair(export, parent_id)

    def test_actual_site_row_schemas_only(self):
        rows = [
            {"atom_map": 3, "state": "PROTECTED",
             "evidence": {"contacts": {"nearest_heavy_atom_A": 2.75}}},
            {"atom_map": 8, "state": "UNKNOWN", "classification": "PROTECTED",
             "nearest_distance_A": 1.0},
        ]
        self.assertEqual(native._protected_maps_from_sites(rows), [3])
        self.assertEqual(native._map_nearest_distance(rows, 3), 2.75)
        with self.assertRaises(ValueError):
            native._map_nearest_distance(rows, 8)
        rows[0]["evidence"]["contacts"]["nearest_heavy_atom_A"] = float("nan")
        with self.assertRaisesRegex(ValueError, "finite"):
            native._map_nearest_distance(rows, 3)

    def test_canonical_native_pdb_preserves_chain_and_coordinates(self):
        text = native.native_receptor_pdb(_protein())
        atom = text.splitlines()[0]
        self.assertTrue(atom.startswith("ATOM  "))
        self.assertEqual(atom[21], "A")
        self.assertEqual(atom[30:38].strip(), "1.000")
        self.assertEqual(atom[38:46].strip(), "2.000")
        self.assertEqual(atom[46:54].strip(), "3.000")
        self.assertEqual(atom[76:78].strip(), "C")
        rows = _protein()
        rows[0]["label_comp_id"] = "HOH"
        with self.assertRaisesRegex(ValueError, "noncanonical"):
            native.native_receptor_pdb(rows)

    def test_residue_deletions_are_explicit_typed_and_never_automatic(self):
        ligand = [{"xyz_A": [0.0, 0.0, 0.0]}]
        token, records = native.validate_remote_incomplete_residue_exclusions(
            _lys(21.55), ligand, []
        )
        self.assertIsInstance(token, native._ValidatedResidueDeletions)
        self.assertEqual(list(token), [])
        self.assertEqual(records, [])
        token, records = native.validate_remote_incomplete_residue_exclusions(
            _lys(21.55), ligand, ["A:1385"]
        )
        self.assertEqual(list(token), ["A:1385"])
        self.assertEqual(records[0]["missing_heavy_atom_names"], ["CD", "CE", "NZ"])
        with self.assertRaisesRegex(ValueError, "within the 20 A"):
            native.validate_remote_incomplete_residue_exclusions(
                _lys(19.0), ligand, ["A:1385"]
            )
        with self.assertRaisesRegex(ValueError, "Complete residue"):
            native.validate_remote_incomplete_residue_exclusions(
                _lys(25.0, complete=True), ligand, ["A:1385"]
            )
        with tempfile.TemporaryDirectory() as root:
            with self.assertRaisesRegex(ValueError, "explicit validation"):
                native.prepare_native_receptor(_protein(), Path(root), ["A:1385"])

    def test_receptor_preparation_hard_fails_without_output(self):
        class Polymer:
            @staticmethod
            def from_pdb_string(*args, **kwargs):
                if kwargs["allow_bad_res"] is not False or kwargs["residues_to_delete"] is not None:
                    raise AssertionError("unsafe Meeko flags")
                raise RuntimeError("bad residue")

        meeko = types.SimpleNamespace(
            ResidueChemTemplates=types.SimpleNamespace(
                create_from_defaults=lambda: object()),
            MoleculePreparation=type("Preparation", (), {}),
            Polymer=Polymer,
            PDBQTWriterLegacy=types.SimpleNamespace(
                write_from_polymer=lambda polymer: (_ for _ in ()).throw(
                    AssertionError("writer must not run"))),
        )
        with tempfile.TemporaryDirectory() as root, patch.dict(
            sys.modules, {"meeko": meeko}
        ):
            output = Path(root) / "receptor"
            with self.assertRaisesRegex(RuntimeError, "bad residue"):
                native.prepare_native_receptor(_protein(), output)
            self.assertFalse((output / "native-receptor.pdbqt").exists())

    def test_native_pose_coordinates_are_not_aligned(self):
        parent = _mapped_halogen()
        pose = Chem.Mol(parent)
        for mol, offset in ((parent, 0.0), (pose, 5.0)):
            conformer = Chem.Conformer(mol.GetNumAtoms())
            conformer.Set3D(True)
            for index in range(mol.GetNumAtoms()):
                conformer.SetAtomPosition(index, (index + offset, 0.0, 0.0))
            mol.AddConformer(conformer)
        from packages.science.warhead_sites import pose_preservation
        result = pose_preservation(parent, [pose], _protein(), [1, 8])
        diagnostic = result["all_poses_diagnostics"][0]
        self.assertIs(diagnostic["alignment_applied"], False)
        self.assertAlmostEqual(diagnostic["core_RMSD_A_in_receptor_frame"], 5.0)

    def test_manifest_hashes_exact_bytes(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            (root / "a.txt").write_bytes(b"abc\n")
            self.assertEqual(native._manifest(root)["files"], [{
                "path": "a.txt", "bytes": 4,
                "sha256": hashlib.sha256(b"abc\n").hexdigest(),
            }])


if __name__ == "__main__":
    unittest.main()
