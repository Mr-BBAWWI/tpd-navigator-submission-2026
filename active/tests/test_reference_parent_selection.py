"""Reference-parent selection, integrity, and default-behavior regressions."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from rdkit import Chem
from rdkit.Geometry import Point3D

from packages.platform.design_panel import input_binding
from packages.science.dual_e3 import SOURCE, validate
from packages.science.reference_parents import (
    load_reference_parent,
    reference_parent_catalog,
    reference_source_binding,
)


def _json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def _sha(value):
    return hashlib.sha256(value).hexdigest()


def _reference_fixture(root):
    root = Path(root)
    parents = root / "parents"
    parents.mkdir(parents=True)
    parent_id = "SMARCA2-TEST-LIG"
    mol = Chem.MolFromSmiles("[CH3:1][OH:2]")
    conformer = Chem.Conformer(mol.GetNumAtoms())
    conformer.SetAtomPosition(0, Point3D(1.0, 2.0, 3.0))
    conformer.SetAtomPosition(1, Point3D(2.0, 2.0, 3.0))
    mol.AddConformer(conformer)
    sdf_path = parents / (parent_id + ".sdf")
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()
    metadata = {
        "format_version": "reference-parents-v1",
        "id": parent_id,
        "target": "SMARCA2",
        "pdb": "TEST",
        "ccd": "LIG",
        "canonical_isomeric_smiles": "CO",
        "sar": {},
        "protected_maps": [],
        "suitability": {"strict_design": False, "exploratory_design": True},
        "limitations": ["Synthetic integrity fixture; not scientific evidence."],
    }
    metadata_path = parents / (parent_id + ".metadata.json")
    metadata_path.write_bytes(_json_bytes(metadata))
    index = {
        "format_version": "reference-parents-v1",
        "default_parent_id": None,
        "reference": {"pdb": "6HAZ", "chain": "A", "sha256": "0" * 64},
        "parents": [{
            "id": parent_id,
            "sdf": "parents/" + parent_id + ".sdf",
            "metadata": "parents/" + parent_id + ".metadata.json",
            "ready": True,
            "suitability": metadata["suitability"],
            "limitations": metadata["limitations"],
        }],
        "not_ready": [{"id": "SMARCA2-NOTREADY-X", "reason": "fixture exclusion"}],
        "limitations": [],
    }
    index_path = root / "index.json"
    index_path.write_bytes(_json_bytes(index))
    entries = []
    for path in (index_path, sdf_path, metadata_path):
        relative = path.relative_to(root).as_posix()
        data = path.read_bytes()
        entries.append({"path": relative, "bytes": len(data), "sha256": _sha(data)})
    manifest = {
        "format_version": "reference-parents-v1",
        "hash_algorithm": "sha256",
        "files": sorted(entries, key=lambda row: row["path"]),
    }
    (root / "manifest.json").write_bytes(_json_bytes(manifest))
    return parent_id, metadata_path


class ReferenceParentSelectionTests(unittest.TestCase):
    def test_default_validation_and_neutral_parent_are_unchanged(self):
        parameters = validate({"dock": False})
        self.assertNotIn("parent_id", parameters)
        parent = Chem.MolFromMolBlock((SOURCE / "SMARCA2-neutral-design.sdf").read_text())
        self.assertEqual(Chem.GetFormalCharge(parent), 0)
        self.assertEqual(validate({"parent_id": "SMARCA2-FX5"})["parent_id"], "SMARCA2-FX5")

    def test_actual_selected_parent_and_invalid_parent(self):
        with tempfile.TemporaryDirectory() as folder:
            parent_id, _ = _reference_fixture(folder)
            catalog = reference_parent_catalog(folder)
            self.assertEqual(catalog["status"], "configured_hash_verified")
            self.assertIn(parent_id, {row["id"] for row in catalog["available"]})
            self.assertIn("SMARCA2-NOTREADY-X", {row.get("id") for row in catalog["not_ready"]})
            mol, metadata = load_reference_parent(parent_id, folder)
            self.assertEqual(metadata["id"], parent_id)
            self.assertEqual(mol.GetNumAtoms(), 2)
            with patch("packages.science.dual_e3.REFERENCE_SOURCE", Path(folder)):
                self.assertEqual(validate({"parent_id": parent_id})["parent_id"], parent_id)
                with self.assertRaisesRegex(ValueError, "Unknown or unavailable parent_id"):
                    validate({"parent_id": "SMARCA2-NOTREADY-X"})

    def test_stale_source_binding_detects_manifest_covered_change(self):
        with tempfile.TemporaryDirectory() as folder:
            _, metadata_path = _reference_fixture(folder)
            before = reference_source_binding(folder)
            metadata_path.write_text("{}\n")
            with self.assertRaisesRegex(ValueError, "Manifest verification failed"):
                reference_source_binding(folder)
            with patch("packages.science.dual_e3.REFERENCE_SOURCE", Path(folder)):
                with self.assertRaises(ValueError):
                    input_binding()
            self.assertEqual(before["status"], "configured_hash_verified")


if __name__ == "__main__":
    unittest.main()
