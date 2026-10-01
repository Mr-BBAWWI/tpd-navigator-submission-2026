from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rdkit import Chem
from rdkit.Chem import AllChem

from scripts import prepare_contact_comparison as comparison_script


class ContactComparisonTests(unittest.TestCase):
    def _molecule(self, smiles: str, seed: int = 7) -> Chem.Mol:
        molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
        self.assertEqual(AllChem.EmbedMolecule(molecule, randomSeed=seed), 0)
        return molecule

    def _write_one(self, path: Path, molecule: Chem.Mol) -> None:
        writer = Chem.SDWriter(str(path))
        try:
            writer.write(molecule)
        finally:
            writer.close()

    def _make_writable(self, root: Path) -> None:
        if not root.exists():
            return
        root.chmod(0o755)
        for path in root.rglob("*"):
            path.chmod(0o755 if path.is_dir() else 0o644)

    def test_stream_sdf_writer_is_exclusive_and_binary_forward_readable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "poses.sdf"
            first = self._molecule("CO")
            second = self._molecule("CN")
            comparison_script._write_sdf(path, [("first", first), ("second", second)])

            with path.open("rb") as stream:
                records = [
                    molecule
                    for molecule in Chem.ForwardSDMolSupplier(stream, removeHs=False)
                    if molecule is not None
                ]
            self.assertEqual([record.GetProp("_Name") for record in records], ["first", "second"])
            self.assertEqual(path.read_bytes().count(b"$$$$\n"), 2)
            with self.assertRaises(FileExistsError):
                comparison_script._write_sdf(path, [("replacement", first)])

    def test_parent_identity_ignores_charge_and_hydrogen_but_not_heavy_graph(self) -> None:
        protonated = self._molecule("C[NH3+]")
        neutral = self._molecule("CN")
        different = self._molecule("CO")

        same_identity = comparison_script._parent_identity(protonated, neutral)
        self.assertTrue(same_identity["same_parent_heavy_graph"])
        self.assertTrue(same_identity["same_parent_compound"])
        self.assertEqual(
            same_identity["scope"],
            "constitutional heavy graph only; stereochemistry not established",
        )
        self.assertEqual(same_identity["stereo_identity_status"], "not_established")
        self.assertTrue(same_identity["chemical_state_changed"])
        self.assertTrue(same_identity["protonation_or_formal_charge_state_changed"])
        self.assertNotEqual(
            same_identity["before_state"]["chemical_state_sha256"],
            same_identity["after_state"]["chemical_state_sha256"],
        )

        different_identity = comparison_script._parent_identity(protonated, different)
        self.assertFalse(different_identity["same_parent_heavy_graph"])
        self.assertFalse(different_identity["same_parent_compound"])
        self.assertFalse(
            different_identity["protonation_or_formal_charge_state_changed"]
        )
        self.assertNotEqual(
            different_identity["before_heavy_graph_sha256"],
            different_identity["after_heavy_graph_sha256"],
        )

    def test_chemical_state_key_is_atom_map_invariant(self) -> None:
        unmapped = self._molecule("C[C@H](O)F")
        mapped = Chem.Mol(unmapped)
        for atom in mapped.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)

        self.assertEqual(
            comparison_script._chemical_state_key(unmapped),
            comparison_script._chemical_state_key(mapped),
        )
        identity = comparison_script._parent_identity(unmapped, mapped)
        self.assertTrue(identity["same_parent_heavy_graph"])
        self.assertFalse(identity["chemical_state_changed"])
        self.assertFalse(identity["protonation_or_formal_charge_state_changed"])

    def test_main_records_computed_raw_and_optimized_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            before_path = root / "before.sdf"
            after_path = root / "after.sdf"
            protein_path = root / "protein.pdb"
            executable_path = root / "pdb2pqr"
            output = root / "result"
            before_input = self._molecule("CN", seed=11)
            after_input = Chem.Mol(before_input)
            after_conformer = after_input.GetConformer()
            for atom_index in range(after_input.GetNumAtoms()):
                position = after_conformer.GetAtomPosition(atom_index)
                position.x += 1.0
                after_conformer.SetAtomPosition(atom_index, position)
            self._write_one(before_path, before_input)
            self._write_one(after_path, after_input)
            protein_path.write_text("HEADER TEST\n", encoding="utf-8")
            executable_path.write_text("#!/bin/sh\n", encoding="utf-8")

            loaded_before = comparison_script._one_sdf(before_path)
            loaded_after = comparison_script._one_sdf(after_path)
            optimized_before = Chem.Mol(loaded_before)
            optimized_after = Chem.Mol(loaded_after)
            for molecule, displacement in ((optimized_before, 0.25), (optimized_after, 0.5)):
                conformer = molecule.GetConformer()
                hydrogen_index = next(
                    atom.GetIdx()
                    for atom in molecule.GetAtoms()
                    if atom.GetAtomicNum() == 1
                )
                position = conformer.GetAtomPosition(hydrogen_index)
                position.z += displacement
                conformer.SetAtomPosition(hydrogen_index, position)

            optimized_results = [
                ({"hydrogens_optimized": optimized_before}, {"method": "test-before"}),
                ({"hydrogens_optimized": optimized_after}, {"method": "test-after"}),
            ]
            receipt = {
                "heavy_atom_mapping": {
                    "missing_from_output": [],
                    "maximum_displacement_A": 0.0,
                }
            }
            prepared = {
                "protein_atoms": [],
                "protein_hydrogens": [],
                "uncertainties": [],
            }
            profiles = [
                {"interactions": [{"identifier": "raw-before"}]},
                {"interactions": [{"identifier": "raw-after"}]},
                {"interactions": [{"identifier": "retained"}]},
                {"interactions": [{"identifier": "retained"}, {"identifier": "gained"}]},
            ]
            report = {
                "lost_identifiers": [],
                "gained_identifiers": ["gained"],
                "retained_identifiers": ["retained"],
            }
            argv = [
                "prepare_contact_comparison.py",
                "--before-sdf",
                str(before_path),
                "--after-sdf",
                str(after_path),
                "--protein-pdb",
                str(protein_path),
                "--pdb2pqr",
                str(executable_path),
                "--output",
                str(output),
            ]

            try:
                with mock.patch.object(sys, "argv", argv), mock.patch.object(
                    comparison_script, "run_pdb2pqr", return_value=receipt
                ), mock.patch.object(
                    comparison_script, "read_prepared_protein", return_value=prepared
                ), mock.patch.object(
                    comparison_script,
                    "optimize_ligand_hydrogens",
                    side_effect=optimized_results,
                ), mock.patch.object(
                    comparison_script, "interaction_profile", side_effect=profiles
                ), mock.patch.object(
                    comparison_script, "before_after_report", return_value=report
                ):
                    self.assertEqual(comparison_script.main(), 0)

                document = json.loads(
                    (output / "contact_comparison.json").read_text(encoding="utf-8")
                )
                metadata = document["ligand_metadata"]
                self.assertIn("before_crystal", metadata)
                self.assertIn("parent_redock", metadata)
                self.assertTrue(
                    document["parent_compound_identity"]["same_parent_compound"]
                )

                expected_before_source = hashlib.sha256(before_path.read_bytes()).hexdigest()
                expected_after_source = hashlib.sha256(after_path.read_bytes()).hexdigest()
                self.assertEqual(
                    metadata["before_crystal"]["source_sdf_sha256"],
                    expected_before_source,
                )
                self.assertEqual(
                    metadata["parent_redock"]["source_sdf_sha256"],
                    expected_after_source,
                )
                self.assertEqual(
                    metadata["before_crystal"]["raw_pose_sha256"],
                    comparison_script._molecule_sha256(loaded_before),
                )
                self.assertEqual(
                    metadata["parent_redock"]["raw_pose_sha256"],
                    comparison_script._molecule_sha256(loaded_after),
                )

                with (output / "ligand_coordinates.sdf").open("rb") as stream:
                    written = [
                        molecule
                        for molecule in Chem.ForwardSDMolSupplier(stream, removeHs=False)
                        if molecule is not None
                    ]
                self.assertEqual(len(written), 4)
                self.assertEqual(
                    [molecule.GetProp("_Name") for molecule in written],
                    ["before_crystal_raw", "before_crystal_hydrogens_optimized", "parent_redock_raw", "parent_redock_hydrogens_optimized"],
                )
                self.assertEqual(
                    metadata["before_crystal"]["hydrogens_optimized_pose_sha256"],
                    comparison_script._molecule_sha256(written[1]),
                )
                self.assertEqual(
                    metadata["parent_redock"]["hydrogens_optimized_pose_sha256"],
                    comparison_script._molecule_sha256(written[3]),
                )
                effects = document["hydrogen_effects"]
                self.assertEqual(effects["before_crystal"]["raw_profile"], profiles[0])
                self.assertEqual(effects["before_crystal"]["optimized_profile"], profiles[2])
                self.assertEqual(effects["parent_redock"]["raw_profile"], profiles[1])
                self.assertEqual(effects["parent_redock"]["optimized_profile"], profiles[3])
                self.assertEqual(document["before_profile"], profiles[2])
                self.assertEqual(document["after_profile"], profiles[3])
                self.assertNotEqual(
                    metadata["before_crystal"]["raw_pose_sha256"],
                    metadata["before_crystal"]["hydrogens_optimized_pose_sha256"],
                )
                self.assertEqual(
                    (output / "contact_comparison.json").read_bytes()[-1:], b"\n"
                )
            finally:
                self._make_writable(output)


if __name__ == "__main__":
    unittest.main()
