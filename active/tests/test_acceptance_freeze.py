from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "freeze_acceptance_sources.py"
SPEC = importlib.util.spec_from_file_location("freeze_acceptance_sources", SCRIPT)
freeze = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(freeze)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True), encoding="utf-8")


class AcceptanceFreezeTests(unittest.TestCase):
    def fixture(self, root: Path) -> argparse.Namespace:
        catalog = root / "catalog"
        (catalog / "source_snapshots" / digest(b"cif")).mkdir(parents=True)
        (catalog / "source_snapshots" / digest(b"cif") / "x.cif").write_bytes(b"cif")
        (catalog / "bound_ligands").mkdir()
        (catalog / "bound_ligands" / "x.sdf").write_bytes(b"sdf")
        write_json(catalog / "warhead_catalog.json", {
            "sources": [{"relative_path": digest(b"cif") + "/x.cif", "sha256": digest(b"cif")}],
            "artifacts": [{"relative_path": "bound_ligands/x.sdf", "sha256": digest(b"sdf")}],
        })

        medchem = root / "medchem"
        medchem.mkdir()
        write_json(medchem / "medchem_evidence.json", {"sources": {}, "scope": "not approval"})
        write_json(medchem / "route_evidence.json", {"evidence_records": []})
        (medchem / "parent_ligands.sdf").write_bytes(b"parents")

        preparation_dir = root / "preparation"
        preparation_dir.mkdir()
        write_json(preparation_dir / "contact_comparison.json", {
            "report": str((preparation_dir / "before.pdb").absolute()),
            "prepared_structure": str((preparation_dir / "after.pqr").absolute()),
            "executable": "C:\\private\\pdb2pqr.exe",
            "executable_sha256": "a" * 64,
            "command": [
                "C:\\private\\pdb2pqr.exe",
                "--pdb-output=" + str((preparation_dir / "before.pdb").absolute()),
                str((preparation_dir / "after.pqr").absolute()),
            ],
            "stderr_tail": "wrote " + str((preparation_dir / "after.pqr").absolute()),
            "computed_pKa": {"value": 7.1, "measurement_status": "computed_not_measured"},
        })
        (preparation_dir / "before.pdb").write_bytes(b"ATOM\n")
        (preparation_dir / "after.pqr").write_bytes(b"ATOM\n")
        (preparation_dir / "pdb2pqr.stdout.bin").write_bytes(b"PDB2PQR stdout\n")
        (preparation_dir / "pdb2pqr.stderr.bin").write_bytes(b"PDB2PQR stderr\n")
        (preparation_dir / "after.pka").write_bytes(b"pKa output\n")

        reference = root / "6BOY.cif"
        reference.write_bytes(b"reference")
        reference_hash = digest(b"reference")
        seeds = []
        for seed, ligand, e3, contact in ((23, 1.542, 5.2, 0.56), (41, 12.388, 39.2, 0.54), (61, 1.178, 4.8, 0.61)):
            seed_dir = root / f"seed-{seed}"
            output = seed_dir / "boltz_output" / "boltz_results" / "processed" / "msa"
            output.mkdir(parents=True)
            msa = output / "x.npz"
            msa.write_bytes(b"same-msa")
            yaml = seed_dir / "input.yaml"
            yaml.write_bytes(b"name: shared-controlled-input\n")
            log = seed_dir / "inference.log"
            log.write_bytes(f"log {seed}\n".encode())
            output_rel = "boltz_results/processed/msa/x.npz"
            inspection = {
                "status": "success",
                "comparison": {"metrics": {
                    "ligand_heavy_atom_RMSD_after_target_alignment_A": ligand,
                    "e3_CA_RMSD_after_target_alignment_A": e3,
                    "contact_jaccard": contact,
                }},
            }
            receipt = {
                "seed": seed,
                "seed_application": {"actually_applied": True, "actual_seed": seed},
                "status": "failure" if seed == 23 else "success",
                "exit_code": 0,
                "process_state": "completed",
                "inspection": {"status": "inspection_failed"} if seed == 23 else inspection,
                "settings": {"msa_mode": "server"},
                "input": {"msa_mode": "server"},
                "tool_environment": {
                    "boltz_version": "2.2.1",
                    "python_path": "C:\\private\\.venv-boltz\\Scripts\\python.exe",
                },
                "hashes": {
                    "checkpoint_sha256": "0" * 64,
                    "input_yaml_sha256": digest(yaml.read_bytes()),
                    "reference_sha256": reference_hash,
                    "inference_log_sha256": digest(log.read_bytes()),
                    "processed_output_files_sha256": {output_rel: digest(msa.read_bytes())},
                },
            }
            receipt_path = seed_dir / "receipt.json"
            write_json(receipt_path, receipt)
            if seed == 23:
                write_json(seed_dir / "reinspection.json", {
                    "original_receipt_sha256": digest(receipt_path.read_bytes()),
                    "actual_execution_success": True,
                    "structural_assessment_completed": True,
                    "assessment": inspection,
                    "verified_processed_output_files_sha256": receipt["hashes"]["processed_output_files_sha256"],
                })
            seeds.append(seed_dir)
        return argparse.Namespace(
            catalog_dir=catalog,
            medchem_dir=medchem,
            preparation=preparation_dir / "contact_comparison.json",
            seed_dir=seeds,
            extra_source=[],
            reference=reference,
            output=root / "bundle",
        )

    def test_build_is_portable_hash_bound_and_receipt_untouched(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = self.fixture(root)
            receipt = args.seed_dir[0] / "receipt.json"
            original = receipt.read_bytes()
            freeze.build(args)
            self.assertEqual(receipt.read_bytes(), original)
            self.assertEqual(
                (args.output / "supporting/seeds/seed-23/receipt.json").read_bytes(),
                original,
            )
            manifest = json.loads((args.output / "manifest.json").read_text())
            self.assertEqual(manifest["version"], freeze.VERSION)
            for relative, expected in manifest["files"].items():
                self.assertEqual(freeze.sha256_file(args.output / relative), expected)
            self.assertEqual(set(freeze.TOP_FILES).issubset(manifest["files"]), True)
            calibration = json.loads((args.output / "crbn_calibration.json").read_text())
            self.assertEqual(calibration["metrics"]["ligand_rmsd_values_A"], [1.542, 12.388, 1.178])
            self.assertTrue(calibration["new_MSA_baseline"]["same_processed_msa_hashes_across_seeds"])
            self.assertEqual(calibration["seed_receipts"][0]["reinspection_ref"], "supporting/seeds/seed-23/reinspection.json")
            derived = b"".join((args.output / name).read_bytes() for name in freeze.TOP_FILES)
            self.assertNotIn(str(root).encode(), derived)
            self.assertNotIn(b".localdata", derived)
            preparation = json.loads((args.output / "preparation.json").read_text())
            self.assertEqual(preparation["computed_pKa"]["measurement_status"], "computed_not_measured")
            self.assertEqual(preparation["executable"]["status"], "not_packaged_execution_environment")
            self.assertEqual(preparation["executable_sha256"], "a" * 64)
            self.assertEqual(preparation["command"][0]["status"], "not_packaged_execution_environment")
            self.assertEqual(preparation["command"][1], "--pdb-output=supporting/preparation/before.pdb")
            self.assertEqual(preparation["command"][2], "supporting/preparation/after.pqr")
            self.assertTrue((args.output / preparation["command"][1].split("=", 1)[1]).is_file())
            self.assertTrue((args.output / preparation["command"][2]).is_file())
            self.assertEqual(preparation["stderr_tail"], "wrote supporting/preparation/after.pqr")
            self.assertNotIn("supporting/supporting/", preparation["stderr_tail"])
            self.assertEqual(preparation["prepared_structure"], "supporting/preparation/after.pqr")
            self.assertEqual((args.output / "supporting/preparation/after.pka").read_bytes(), b"pKa output\n")
            self.assertEqual((args.output / "supporting/preparation/pdb2pqr.stdout.bin").read_bytes(), b"PDB2PQR stdout\n")
            self.assertEqual((args.output / "supporting/preparation/pdb2pqr.stderr.bin").read_bytes(), b"PDB2PQR stderr\n")
            environment = calibration["seed_receipts"][0]["tool_environment"]
            self.assertEqual(environment["python_path"]["status"], "not_packaged_execution_environment")

    def test_human_filename_is_not_a_scientific_path_but_known_path_field_is(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            builder = freeze.FreezeBuilder(Path(temporary) / "out")
            builder.assert_scientific_refs_packaged({"message": "missing.pdb"})
            with self.assertRaisesRegex(ValueError, "UNRESOLVED_SCIENTIFIC_DATA_FILE"):
                builder.assert_scientific_refs_packaged({"source_path": "missing.pdb"})

    def test_command_unknown_scientific_data_path_fails_not_externalized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            builder = freeze.FreezeBuilder(Path(temporary) / "out")
            command = ["C:\\private\\tool.exe", "C:\\private\\missing.pdb"]
            with self.assertRaisesRegex(ValueError, "UNRESOLVED_PRIVATE_PATH"):
                builder.rewrite({"command": command})

    def test_human_text_rewrite_is_idempotent_and_does_not_replace_path_suffixes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source" / "after.pqr"
            source.parent.mkdir()
            source.write_bytes(b"ATOM\n")
            builder = freeze.FreezeBuilder(root / "out")
            builder.copy(source, "supporting/preparation/after.pqr", (source.name,))
            first = builder.rewrite_human_text("wrote " + str(source.absolute()))
            second = builder.rewrite_human_text(first)
            self.assertEqual(first, "wrote supporting/preparation/after.pqr")
            self.assertEqual(second, first)
            suffix = "unbound/directory/after.pqr"
            self.assertEqual(builder.rewrite_human_text(suffix), suffix)

    def test_cwd_relative_exact_binding_accepts_copied_path_but_unknown_still_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            previous = Path.cwd()
            os.chdir(root)
            try:
                source = Path(".localdata/acceptance-20260930/crbn-server/seed-23/boltz_output/predictions/confidence.json")
                source.parent.mkdir(parents=True)
                source.write_bytes(b"{}")
                builder = freeze.FreezeBuilder(root / "out")
                builder.copy(source, "supporting/seeds/seed-23/boltz_output/predictions/confidence.json")
                supplied = str(source).replace("/", "\\")
                rewritten = builder.rewrite({"prediction": supplied})
                self.assertEqual(
                    rewritten["prediction"],
                    "supporting/seeds/seed-23/boltz_output/predictions/confidence.json",
                )
                with self.assertRaisesRegex(ValueError, "UNRESOLVED_PRIVATE_PATH"):
                    builder.rewrite({"prediction": ".localdata\\unknown\\confidence.json"})
            finally:
                os.chdir(previous)

    def test_preparation_scientific_schema_keys_are_verified(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            builder = freeze.FreezeBuilder(Path(temporary) / "out")
            for document in (
                {"source_paths": ["missing.sdf"]},
                {"before_sdf": "missing.sdf"},
                {"after_sdf": "missing.sdf"},
                {"protein_pdb": "missing.pdb"},
            ):
                with self.subTest(document=document):
                    with self.assertRaisesRegex(ValueError, "UNRESOLVED_SCIENTIFIC_DATA_FILE"):
                        builder.assert_scientific_refs_packaged(document)

    def test_same_hash_ambiguous_alias_uses_deterministic_copy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first" / "shared.sdf"
            second = root / "second" / "shared.sdf"
            first.parent.mkdir()
            second.parent.mkdir()
            first.write_bytes(b"same")
            second.write_bytes(b"same")
            builder = freeze.FreezeBuilder(root / "out")
            builder.copy(first, "supporting/a/shared.sdf", (first.name,))
            builder.copy(second, "supporting/b/shared.sdf", (second.name,))
            self.assertEqual(
                builder.rewrite({"source_filename": "shared.sdf", "sha256": digest(b"same")})["source_filename"],
                "supporting/a/shared.sdf",
            )

    def test_different_hash_ambiguous_alias_requires_disambiguation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first" / "shared.sdf"
            second = root / "second" / "shared.sdf"
            first.parent.mkdir()
            second.parent.mkdir()
            first.write_bytes(b"first")
            second.write_bytes(b"second")
            builder = freeze.FreezeBuilder(root / "out")
            builder.copy(first, "supporting/a/shared.sdf", (first.name,))
            builder.copy(second, "supporting/b/shared.sdf", (second.name,))
            with self.assertRaisesRegex(ValueError, "AMBIGUOUS_SOURCE_BINDING"):
                builder.rewrite({"source_filename": "shared.sdf"})
            resolved = builder.rewrite({
                "source_filename": "shared.sdf",
                "sha256": digest(b"second"),
            })
            self.assertEqual(resolved["source_filename"], "supporting/b/shared.sdf")

    def test_bad_seed_output_hash_rejected_and_failed_scratch_retained(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = self.fixture(root)
            receipt_path = args.seed_dir[1] / "receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["hashes"]["processed_output_files_sha256"]["boltz_results/processed/msa/x.npz"] = "f" * 64
            write_json(receipt_path, receipt)
            with self.assertRaisesRegex(ValueError, "SEED_OUTPUT_HASH_MISMATCH"):
                freeze.build(args)
            self.assertFalse(args.output.exists())
            failed = list(root.glob(".bundle.failed-*"))
            self.assertEqual(len(failed), 1)
            self.assertTrue((failed[0] / "BUILD_FAILED.txt").is_file())

    def test_reinspection_must_match_immutable_original_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = self.fixture(root)
            receipt_path = args.seed_dir[0] / "receipt.json"
            original = receipt_path.read_bytes()
            reinspection_path = args.seed_dir[0] / "reinspection.json"
            reinspection = json.loads(reinspection_path.read_text())
            reinspection["original_receipt_sha256"] = "0" * 64
            write_json(reinspection_path, reinspection)
            with self.assertRaisesRegex(ValueError, "REINSPECTION_ORIGINAL_RECEIPT_HASH_MISMATCH"):
                freeze.build(args)
            self.assertEqual(receipt_path.read_bytes(), original)

    def test_path_traversal_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "TRAVERSAL"):
            freeze.safe_relative("../private/file.json")
        for unsafe in ("a\\b", "a:b", "C:/private/file.json"):
            with self.assertRaisesRegex(ValueError, "INVALID"):
                freeze.safe_relative(unsafe)

    def test_unresolved_private_checkout_path_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            builder = freeze.FreezeBuilder(Path(temporary) / "out")
            with self.assertRaisesRegex(ValueError, "UNRESOLVED_PRIVATE_PATH"):
                builder.rewrite({"prediction": "C:\\private\\checkout\\.localdata\\missing.cif"})

    def test_unresolved_scientific_reference_cannot_fake_existence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            builder = freeze.FreezeBuilder(Path(temporary) / "out")
            fake = builder.destination("supporting/extra/unbound.sdf")
            fake.write_bytes(b"not bound")
            with self.assertRaisesRegex(ValueError, "UNRESOLVED_SCIENTIFIC_DATA_FILE"):
                builder.assert_scientific_refs_packaged({"relative_path": "supporting/extra/unbound.sdf"})

    def test_seed_application_and_metric_bounds_are_enforced(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = self.fixture(root)
            receipt_path = args.seed_dir[1] / "receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["seed_application"]["actually_applied"] = False
            write_json(receipt_path, receipt)
            with self.assertRaisesRegex(ValueError, "SEED_NOT_ACTUALLY_APPLIED"):
                freeze.build(args)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = self.fixture(root)
            receipt_path = args.seed_dir[1] / "receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["inspection"]["comparison"]["metrics"]["contact_jaccard"] = 1.1
            write_json(receipt_path, receipt)
            with self.assertRaisesRegex(ValueError, "SEED_CONTACT_JACCARD_INVALID"):
                freeze.build(args)

    def test_non_log_bin_and_source_ancestor_symlink_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            weights = root / "weights.bin"
            weights.write_bytes(b"weights")
            with self.assertRaisesRegex(ValueError, "FORBIDDEN_SOURCE"):
                freeze.check_source_path(weights)
            real = root / "real"
            real.mkdir()
            (real / "x.pdb").write_bytes(b"ATOM\n")
            link = root / "linked"
            try:
                link.symlink_to(real, target_is_directory=True)
            except OSError:
                self.skipTest("symlinks unavailable")
            with self.assertRaisesRegex(ValueError, "SYMLINK_FORBIDDEN"):
                freeze.regular(link / "x.pdb", "SOURCE")


if __name__ == "__main__":
    unittest.main()
