import hashlib
import json
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from packages.science.reference_parents import (
    _checked_catalog_file,
    _find_catalog,
    _safe_relative,
    _source_descriptor,
    apply_transform,
    global_sequence_alignment,
    kabsch_transform,
    load_reference_parent,
)


def _json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def _mapped_3d(smiles="CCO"):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=17) != 0:
        raise AssertionError("RDKit embedding failed")
    mol = Chem.RemoveHs(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    return mol


def _sdf_bytes(mol):
    from io import StringIO

    stream = StringIO()
    writer = Chem.SDWriter(stream)
    writer.write(mol)
    writer.close()
    return stream.getvalue().encode()


def _support_tree(tmp_path, mol=None):
    root = tmp_path / "support"
    parents = root / "parents"
    parents.mkdir(parents=True)
    mol = mol or _mapped_3d()
    parent_id = "SMARCA2-TEST-LIG"
    sdf_path = f"parents/{parent_id}.sdf"
    metadata_path = f"parents/{parent_id}.metadata.json"
    metadata = {
        "id": parent_id,
        "target": "SMARCA2",
        "pdb": "TEST",
        "ccd": "LIG",
        "sar": {},
        "protected_maps": [],
        "canonical_isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
    }
    index = {
        "format_version": "reference-parents-v1",
        "default_parent_id": None,
        "parents": [{"id": parent_id, "sdf": sdf_path, "metadata": metadata_path}],
    }
    data = {
        "index.json": _json_bytes(index),
        sdf_path: _sdf_bytes(mol),
        metadata_path: _json_bytes(metadata),
    }
    manifest = {
        "format_version": "reference-parents-v1",
        "hash_algorithm": "sha256",
        "files": [
            {"path": name, "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
            for name, content in sorted(data.items())
        ],
    }
    for name, content in data.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    (root / "manifest.json").write_bytes(_json_bytes(manifest))
    return root, parent_id


class ReferenceParentTests(unittest.TestCase):
    def test_kabsch_recovers_known_proper_transform_and_preserves_graph(self):
        moving = np.array([
            [0.0, 0.0, 0.0],
            [1.2, 0.1, 0.0],
            [0.2, 1.4, 0.3],
            [-0.1, 0.4, 1.7],
            [1.1, 1.0, 0.8],
        ])
        angle = 0.71
        rotation = np.array([
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ])
        translation = np.array([3.2, -1.7, 4.4])
        reference = moving @ rotation + translation
        result = kabsch_transform(reference, moving)
        fitted = apply_transform(moving, result)
        self.assertLess(result["rmsd_A"], 1e-12)
        np.testing.assert_allclose(fitted, reference, atol=1e-12)
        self.assertAlmostEqual(np.linalg.det(np.asarray(result["rotation_row_vectors"])), 1.0)

        mol = _mapped_3d("C1CCNCC1")
        before = Chem.MolToSmiles(mol, isomericSmiles=True)
        transformed = Chem.Mol(mol)
        xyz = apply_transform(transformed.GetConformer().GetPositions(), result)
        for index, point in enumerate(xyz):
            transformed.GetConformer().SetAtomPosition(index, point)
        self.assertEqual(Chem.MolToSmiles(transformed, isomericSmiles=True), before)
        self.assertEqual(
            [atom.GetAtomMapNum() for atom in transformed.GetAtoms()],
            list(range(1, transformed.GetNumAtoms() + 1)),
        )

    def test_kabsch_rejects_collinear_scientifically_underdetermined_points(self):
        points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        with self.assertRaisesRegex(ValueError, "Collinear"):
            kabsch_transform(points + 1.0, points)

    def test_global_alignment_uses_full_sequences_and_identity_correspondence(self):
        reference = "ACDEFGHIKLMNPQRSTVWY" * 4
        moving = "ACDEFGHIKLMNPQRSTVWY" * 2 + "ACDEYGHIKLMNPQRSTVWY" + "ACDEFGHIKLMNPQRSTVWY"
        result = global_sequence_alignment(reference, moving)
        self.assertEqual(result["reference_length"], 80)
        self.assertEqual(result["moving_length"], 80)
        self.assertEqual(result["aligned_residue_pair_count"], 80)
        self.assertEqual(result["identical_residue_pair_count"], 79)
        self.assertAlmostEqual(result["sequence_identity_over_aligned_residue_pairs"], 79 / 80)
        mismatch = [pair for pair in result["pairs"] if not pair["identity_match"]]
        self.assertEqual(len(mismatch), 1)
        self.assertEqual(
            mismatch[0]["reference_label_seq_id"],
            mismatch[0]["moving_label_seq_id"],
        )

    def test_global_alignment_keeps_stable_core_despite_terminal_tag_ambiguity(self):
        core = "ACDEFGHIKLMNPQRSTVWY" * 4
        result = global_sequence_alignment(core, "AAAA" + core)
        self.assertTrue(result["optimal_alignment_ambiguous"])
        self.assertGreaterEqual(result["identical_residue_pair_count"], 60)
        self.assertTrue(all(pair["stable_across_all_optimal_alignments"] for pair in result["pairs"]))
        self.assertTrue(result["excluded_moving_positions"])

    def test_global_alignment_repeated_core_does_not_choose_arbitrary_tie(self):
        repeat = "ACDEFGHIKLMNPQRSTVWY"
        result = global_sequence_alignment(repeat, repeat + repeat)
        self.assertTrue(result["optimal_alignment_ambiguous"])
        self.assertEqual(result["pairs"], [])
        self.assertEqual(result["aligned_residue_pair_count"], 0)

    def test_global_alignment_rejects_unknown_residues_instead_of_guessing(self):
        with self.assertRaisesRegex(ValueError, "unambiguous"):
            global_sequence_alignment("ACDX", "ACDE")

    def test_loader_verifies_manifest_and_returns_unknown_sar_unchanged(self):
        with TemporaryDirectory() as directory:
            root, parent_id = _support_tree(Path(directory))
            mol, meta = load_reference_parent(parent_id, root)
            self.assertEqual(mol.GetNumAtoms(), 3)
            self.assertEqual(meta["sar"], {})
            self.assertEqual(meta["protected_maps"], [])

    def test_loader_rejects_wrong_parent_hash(self):
        with TemporaryDirectory() as directory:
            root, parent_id = _support_tree(Path(directory))
            with (root / "parents" / f"{parent_id}.sdf").open("ab") as stream:
                stream.write(b"tamper")
            with self.assertRaisesRegex(ValueError, "Manifest verification failed"):
                load_reference_parent(parent_id, root)

    def test_loader_rejects_missing_manifest_map_and_path_injection(self):
        with TemporaryDirectory() as directory:
            root, parent_id = _support_tree(Path(directory))
            manifest = json.loads((root / "manifest.json").read_text())
            manifest["files"] = [
                entry for entry in manifest["files"]
                if not entry["path"].endswith(".metadata.json")
            ]
            (root / "manifest.json").write_bytes(_json_bytes(manifest))
            with self.assertRaisesRegex(ValueError, "does not cover"):
                load_reference_parent(parent_id, root)
            with self.assertRaisesRegex(ValueError, "Invalid parent id"):
                load_reference_parent("../metadata", root)

    def test_different_graph_is_not_autoapproved_even_with_valid_hashes(self):
        with TemporaryDirectory() as directory:
            root, parent_id = _support_tree(Path(directory))
            metadata_path = root / "parents" / f"{parent_id}.metadata.json"
            metadata = json.loads(metadata_path.read_text())
            metadata["canonical_isomeric_smiles"] = "CCN"
            metadata_data = _json_bytes(metadata)
            metadata_path.write_bytes(metadata_data)

            manifest = json.loads((root / "manifest.json").read_text())
            for entry in manifest["files"]:
                if entry["path"].endswith(".metadata.json"):
                    entry["bytes"] = len(metadata_data)
                    entry["sha256"] = hashlib.sha256(metadata_data).hexdigest()
            (root / "manifest.json").write_bytes(_json_bytes(manifest))
            with self.assertRaisesRegex(ValueError, "graph disagrees"):
                load_reference_parent(parent_id, root)

    def test_manifest_relative_path_escape_fails_before_read(self):
        with TemporaryDirectory() as directory:
            root, parent_id = _support_tree(Path(directory))
            manifest = json.loads((root / "manifest.json").read_text())
            manifest["files"].append({
                "path": "../outside",
                "bytes": 0,
                "sha256": "0" * 64,
            })
            (root / "manifest.json").write_bytes(_json_bytes(manifest))
            with self.assertRaisesRegex(ValueError, "Unsafe|escapes"):
                load_reference_parent(parent_id, root)

    def test_safe_relative_rejects_colons_backslashes_unc_and_symlinks(self):
        with TemporaryDirectory() as directory:
            root = Path(directory) / "root"
            root.mkdir()
            (root / "ok.txt").write_text("ok")
            for unsafe in ("C:/outside", "folder\\file", "//server/share", "a/../ok.txt"):
                with self.subTest(unsafe=unsafe):
                    with self.assertRaisesRegex(ValueError, "Unsafe"):
                        _safe_relative(root, unsafe)

            target = root / "ok.txt"
            link = root / "link.txt"
            try:
                os.symlink(target, link)
            except (OSError, NotImplementedError):
                return
            with self.assertRaisesRegex(ValueError, "symbolic link"):
                _safe_relative(root, "link.txt")

    def test_actual_catalog_schema_uses_fixed_catalog_and_source_snapshot_root(self):
        with TemporaryDirectory() as directory:
            catalog_root = Path(directory) / "collector"
            catalog_root.mkdir()
            source_hash = hashlib.sha256(b"source-cif").hexdigest()
            source_relative = f"{source_hash}/6HAZ.cif"
            source_path = catalog_root / "source_snapshots" / source_relative
            source_path.parent.mkdir(parents=True)
            source_path.write_bytes(b"source-cif")

            catalog = {
                "format_version": "collector-v2",
                "source_root": "source_snapshots",
                "cards": [{
                    "group_id": "ligand-example",
                    "domain_partition": "bromodomain",
                    "bound_sdf": {
                        "relative_path": "bound_ligands/example.sdf",
                        "sha256": "0" * 64,
                    },
                    "canonical_isomeric_smiles": "CCO",
                    "representative_record_id": "6HAZ:C:FX5",
                    "experimental_versions": [{
                        "record_id": "6HAZ:C:FX5",
                        "pdb": "6HAZ",
                        "ccd": "FX5",
                        "source_hashes": {"pdb_mmcif_sha256": source_hash},
                        "contacts": {"by_exact_target_chain": [{
                            "target_label_asym_id": "A",
                            "contact_pair_count": 1,
                        }]},
                        "identity": {
                            "canonical_isomeric_smiles": "CCO",
                            "formal_charge": 0,
                        },
                        "atom_mapping": [],
                    }],
                }],
                "sources": [{
                    "bytes": len(b"source-cif"),
                    "relative_path": source_relative,
                    "retrieval": "immutable_cache_reuse",
                    "sha256": source_hash,
                    "url": "https://files.rcsb.org/download/6HAZ.cif",
                }],
            }
            (catalog_root / "catalog.json").write_bytes(_json_bytes(catalog))
            irrelevant = catalog_root / "irrelevant" / "large.json"
            irrelevant.parent.mkdir()
            irrelevant.write_text("not a catalog")

            root, catalog_path, loaded = _find_catalog(catalog_root)
            self.assertEqual(root, catalog_root.resolve())
            self.assertEqual(catalog_path, (catalog_root / "catalog.json").resolve())
            descriptor = _source_descriptor(loaded, source_hash, "6HAZ.cif")
            checked, actual_hash = _checked_catalog_file(
                root / "source_snapshots", descriptor
            )
            self.assertEqual(checked, source_path.resolve())
            self.assertEqual(actual_hash, source_hash)


if __name__ == "__main__":
    unittest.main()
