from __future__ import annotations

import hashlib
import importlib.util
import os
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "inspect_core_state_robustness.py"
SPEC = importlib.util.spec_from_file_location("inspect_core_state_robustness", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(module)


class ClassificationTests(unittest.TestCase):
    def test_hbond_alone_is_insufficient(self):
        self.assertFalse(module.classify_row(True, set(), 0.2))

    def test_hydrophobes_alone_are_insufficient(self):
        self.assertFalse(module.classify_row(False, set(module.HYDROPHOBIC_REQUIRED), 0.2))

    def test_rmsd_does_not_fake_hbond(self):
        self.assertFalse(
            module.classify_row(False, set(module.HYDROPHOBIC_REQUIRED), 0.0)
        )

    def test_contacts_pass_even_when_auxiliary_rmsd_fails(self):
        self.assertTrue(
            module.classify_row(True, set(module.HYDROPHOBIC_REQUIRED), 4.0)
        )

    def test_same_pose_across_all_states_is_required(self):
        rows = []
        for pose in range(5):
            for state in range(2):
                rows.append({
                    "pose_index": pose,
                    "state_index": state,
                    "required_contact_evidence": (pose == state),
                })
        self.assertEqual(module.supported_pose_indices(rows, 2), [])
        for row in rows:
            if row["pose_index"] == 3:
                row["required_contact_evidence"] = True
        self.assertEqual(module.supported_pose_indices(rows, 2), [3])

    def test_missing_state_blocks_support(self):
        rows = [{
            "pose_index": 0,
            "state_index": 0,
            "required_contact_evidence": True,
        }]
        self.assertEqual(module.supported_pose_indices(rows, 2), [])

    def test_duplicate_state_blocks_support(self):
        rows = [
            {"pose_index": 0, "state_index": 0, "required_contact_evidence": True},
            {"pose_index": 0, "state_index": 1, "required_contact_evidence": True},
            {"pose_index": 0, "state_index": 1, "required_contact_evidence": True},
        ]
        self.assertEqual(module.supported_pose_indices(rows, 2), [])


class IntegrityTests(unittest.TestCase):
    @staticmethod
    def digest(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def test_hash_tamper_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact = root / "artifact.bin"
            artifact.write_bytes(b"registered")
            mapping = [{"path": "artifact.bin", "sha256": self.digest(artifact)}]
            module.verify_mappings(root, mapping)
            artifact.write_bytes(b"tampered")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                module.verify_mappings(root, mapping)

    def test_parent_path_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapping = [{"path": "../x", "sha256": "0" * 64}]
            with self.assertRaisesRegex(ValueError, "unsafe artifact path"):
                module.verify_mappings(root, mapping)

    def test_absolute_path_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact = root / "artifact.bin"
            artifact.write_bytes(b"data")
            mapping = [{"path": str(artifact.resolve()), "sha256": self.digest(artifact)}]
            with self.assertRaisesRegex(ValueError, "unsafe artifact path"):
                module.verify_mappings(root, mapping)

    @unittest.skipIf(os.name == "nt", "symlink creation is not consistently permitted on Windows")
    def test_symlink_traversal_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as outside_dir:
            root = Path(temporary)
            outside = Path(outside_dir) / "artifact.bin"
            outside.write_bytes(b"outside")
            (root / "link").symlink_to(Path(outside_dir), target_is_directory=True)
            mapping = [{"path": "link/artifact.bin", "sha256": self.digest(outside)}]
            with self.assertRaisesRegex(ValueError, "symlink traversal"):
                module.verify_mappings(root, mapping)

    def test_recorded_comparison_hash_is_checked(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            comparison = root / "comparison.json"
            comparison.write_text("{}", encoding="utf-8")
            manifest = {
                "source_files": [],
                "artifacts": [{
                    "path": "comparison.json",
                    "sha256": self.digest(comparison),
                }],
            }
            module.verify_comparison_manifest(root, manifest)
            comparison.write_text('{"tampered":true}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                module.verify_comparison_manifest(root, manifest)


if __name__ == "__main__":
    unittest.main()
