from __future__ import annotations

import stat
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from packages.science import parent_funnel as funnel


class ParentFunnelTests(unittest.TestCase):
    def test_exact_derivation_rejects_tampered_count(self) -> None:
        reference = {
            "status": "not_requested",
            "independently_reverified_from_bundle": False,
            "exact_identity_consistent": None,
        }
        derived = {
            "parent_id": "P1",
            "actual_counts": {"cheap_passed": 2, "selected": 0},
            "actual_families": {"cheap_passed": ["substitution"]},
            "source_proof": {"input_binding": {"sha256": "a" * 64}},
            "reference_parent_verification": reference,
        }
        funnel._assert_exact_parent_record(dict(derived), derived)
        tampered = dict(derived)
        tampered["actual_counts"] = {"cheap_passed": 3, "selected": 0}
        with self.assertRaisesRegex(
            funnel.FunnelError,
            "PARENT_SCIENTIFIC_DERIVATION_MISMATCH",
        ):
            funnel._assert_exact_parent_record(tampered, derived)

    def test_bundle_relative_paths_reject_cross_platform_unsafe_forms(self) -> None:
        relative_paths = [
            "",
            "/exports/P1/manifest.json",
            "C:/exports/P1/manifest.json",
            "exports/P1/../manifest.json",
            "exports//P1/manifest.json",
            "exports\\P1\\manifest.json",
            "exports/./P1/manifest.json",
        ]
        for relative in relative_paths:
            with self.subTest(relative=relative):
                with self.assertRaises(funnel.FunnelError):
                    funnel._bundle_parts(relative)

    def test_parent_proof_paths_are_exact(self) -> None:
        copied = {
            "artifact_blob_closure": [
                {
                    "artifact_id": "blob-1",
                    "sha256": "b" * 64,
                    "bytes": 7,
                    "version": 1,
                }
            ]
        }
        proof = funnel._expected_copy_proof(copied, "P1")
        self.assertEqual(proof["manifest_path"], "exports/P1/manifest.json")
        self.assertEqual(proof["receipt_path"], "exports/P1/job-receipt.json")
        self.assertEqual(
            proof["artifact_blobs"][0]["relative_path"],
            "exports/P1/blobs/blob-1",
        )

    def test_onedrive_cloud_reparse_is_not_treated_as_link(self) -> None:
        class FakePath:
            def lstat(self) -> SimpleNamespace:
                return SimpleNamespace(
                    st_mode=stat.S_IFREG,
                    st_file_attributes=0x400,
                    st_reparse_tag=0x9000001A,
                )

            def __str__(self) -> str:
                return "OneDrive-placeholder"

        with mock.patch.object(
            stat,
            "FILE_ATTRIBUTE_REPARSE_POINT",
            0x400,
            create=True,
        ):
            self.assertFalse(funnel._is_link_or_reparse(FakePath()))  # type: ignore[arg-type]

    def test_name_surrogate_and_junction_reparses_are_rejected(self) -> None:
        class FakePath:
            def __init__(self, tag: int) -> None:
                self.tag = tag

            def lstat(self) -> SimpleNamespace:
                return SimpleNamespace(
                    st_mode=stat.S_IFREG,
                    st_file_attributes=0x400,
                    st_reparse_tag=self.tag,
                )

            def __str__(self) -> str:
                return "reparse"

        with mock.patch.object(
            stat,
            "FILE_ATTRIBUTE_REPARSE_POINT",
            0x400,
            create=True,
        ):
            for tag in (0xA000000C, 0xA0000003):
                with self.subTest(tag=hex(tag)):
                    self.assertTrue(
                        funnel._is_link_or_reparse(FakePath(tag))  # type: ignore[arg-type]
                    )

    def test_verified_copy_preserves_source_and_detects_expected_hash(self) -> None:
        with tempfile.TemporaryDirectory(prefix="parent-funnel-copy-") as directory:
            root = Path(directory)
            source_root = root / "source"
            source_root.mkdir()
            source = source_root / "receipt.json"
            source.write_bytes(b"private-source-data\n")
            target_root = root / "target"
            target_root.mkdir()
            target = target_root / "receipt.json"
            digest = funnel.sha256_file(source)

            funnel._copy_verified_file(
                source,
                target,
                source_root,
                "TEST_SOURCE",
                digest,
                source.stat().st_size,
            )

            self.assertTrue(source.exists())
            self.assertEqual(source.read_bytes(), b"private-source-data\n")
            self.assertEqual(target.read_bytes(), source.read_bytes())

    def test_actual_5dkh_reference_ignores_map_only_stereo(self) -> None:
        reference_root = Path("cases/reference_parents")
        mol, metadata = funnel.load_reference_parent("SMARCA2-5DKH-5C0", reference_root)

        check = funnel._reference_check(
            "SMARCA2-5DKH-5C0",
            {"warhead": metadata},
            reference_root,
        )

        self.assertEqual(check["status"], "verified_at_build_time_external_catalog")
        self.assertTrue(check["exact_identity_consistent"])
        self.assertEqual(
            check["reference_map_free_canonical_isomeric_smiles"],
            funnel._map_free_canonical_isomeric_smiles(mol),
        )

    def test_reference_check_still_rejects_opposite_genuine_stereo(self) -> None:
        reference_root = Path("cases/reference_parents")
        _, metadata = funnel.load_reference_parent("SMARCA2-5DKH-5C0", reference_root)
        tampered = dict(metadata)
        tampered["canonical_isomeric_smiles"] = metadata[
            "canonical_isomeric_smiles"
        ].replace("[C@H]2", "[C@@H]2", 1)

        with self.assertRaisesRegex(
            funnel.FunnelError,
            "REFERENCE_PARENT_IDENTITY_OR_MAP_MISMATCH:SMARCA2-5DKH-5C0",
        ):
            funnel._reference_check(
                "SMARCA2-5DKH-5C0",
                {"warhead": tampered},
                reference_root,
            )

    def test_bool_is_not_accepted_as_integer_count(self) -> None:
        with self.assertRaisesRegex(funnel.FunnelError, "COUNT_INVALID"):
            funnel._exact_int(False, 0, "COUNT")


if __name__ == "__main__":
    unittest.main()
