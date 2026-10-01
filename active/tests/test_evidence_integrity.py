import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from packages.science.warhead_evidence import (
    SCHEMA_VERSION,
    _matching_curated,
    load_catalog,
)


class EvidenceIntegrityTests(unittest.TestCase):
    def test_curated_record_without_identity_does_not_match_every_ligand(self):
        entries = [{
            "source": "paper", "locator": "table 1", "target_condition": "assay",
            "target_accession": "P51531",
        }]
        self.assertEqual(_matching_curated(entries, "1ABC", "LIG", "CCO"), [])

    def test_curated_match_requires_exact_pair_or_exact_structure(self):
        entries = [
            {"pdb": "1ABC", "ccd": "LIG", "source": "x"},
            {"canonical_isomeric_smiles": "C[C@H](O)F", "source": "y"},
        ]
        self.assertEqual(len(_matching_curated(entries, "1ABC", "LIG", "CCO")), 1)
        self.assertEqual(len(_matching_curated(entries, "9XYZ", "OTHER", "C[C@H](O)F")), 1)
        self.assertEqual(_matching_curated(entries, "1ABC", "OTHER", "CCO"), [])

    def test_catalog_rejects_parent_traversal(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            escaped = root.parent / "escaped-evidence-test.bin"
            escaped.write_bytes(b"outside")
            try:
                catalog = {
                    "schema_version": SCHEMA_VERSION,
                    "sources": [],
                    "artifacts": [{
                        "relative_path": "../escaped-evidence-test.bin",
                        "sha256": hashlib.sha256(b"outside").hexdigest(),
                    }],
                }
                path = root / "catalog.json"
                path.write_text(json.dumps(catalog), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "Unsafe relative path"):
                    load_catalog(path)
            finally:
                escaped.unlink(missing_ok=True)

    def test_catalog_rejects_symlink_artifact(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            real = root / "real.bin"
            real.write_bytes(b"evidence")
            link = root / "linked.bin"
            try:
                link.symlink_to(real)
            except (OSError, NotImplementedError):
                self.skipTest("symlinks are unavailable")
            catalog = {
                "schema_version": SCHEMA_VERSION,
                "sources": [],
                "artifacts": [{
                    "relative_path": "linked.bin",
                    "sha256": hashlib.sha256(b"evidence").hexdigest(),
                }],
            }
            path = root / "catalog.json"
            path.write_text(json.dumps(catalog), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Symlinks are not permitted"):
                load_catalog(path)


if __name__ == "__main__":
    unittest.main()
