from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import sys
import tempfile
import types
import unittest
from contextlib import contextmanager
from pathlib import Path
from unittest import mock

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "export_expert_evidence.py"
SPEC = importlib.util.spec_from_file_location("export_expert_evidence", SCRIPT_PATH)
exporter = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(exporter)


class Args:
    store = None
    job_id = None
    design_json = None


class ExportExpertEvidenceTests(unittest.TestCase):
    def test_correct_defaults_are_actual_fixed_sources(self):
        self.assertEqual(exporter.DEFAULT_CHARGED.name, "SMARCA2-parent.sdf")
        self.assertEqual(exporter.DEFAULT_PROTEIN.name, "6HAZ.cif")
        self.assertEqual(exporter.DEFAULT_NEUTRAL.name, "SMARCA2-neutral-design.sdf")
        for path in (
            exporter.DEFAULT_CHARGED,
            exporter.DEFAULT_NEUTRAL,
            exporter.DEFAULT_PROTEIN,
        ):
            self.assertTrue(path.is_file(), path)
            self.assertFalse(path.is_symlink(), path)
        self.assertEqual(exporter._sha(exporter.DEFAULT_CHARGED), exporter.EXPECTED_SOURCE_HASHES["SMARCA2-parent.sdf"])
        self.assertEqual(exporter._sha(exporter.DEFAULT_NEUTRAL), exporter.EXPECTED_SOURCE_HASHES["SMARCA2-neutral-design.sdf"])
        self.assertEqual(exporter._sha(exporter.DEFAULT_PROTEIN), exporter.EXPECTED_SOURCE_HASHES["6HAZ.cif"])

    def test_actual_protein_cif_chain_a_is_an_l_polypeptide(self):
        atoms = exporter._protein(exporter.DEFAULT_PROTEIN, "A")
        self.assertTrue(atoms)
        self.assertEqual({row["label_asym_id"] for row in atoms}, {"A"})
        self.assertEqual({row["auth_asym_id"] for row in atoms}, {"A"})

    def test_unicode_mol_uses_bytes_and_preserves_atom_maps(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "분자 자료"
            directory.mkdir()
            path = directory / "부모 구조.sdf"
            molecule = Chem.MolFromSmiles("[CH3:7][OH:8]")
            stream = io.StringIO()
            writer = Chem.SDWriter(stream)
            writer.write(molecule)
            writer.flush()
            writer.close()
            path.write_text(stream.getvalue(), encoding="utf-8")
            loaded = exporter._mol(path)
            self.assertEqual([atom.GetAtomMapNum() for atom in loaded.GetAtoms()], [7, 8])

    def test_invalid_or_missing_input_rejected_before_output_creation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design = root / "design.json"
            design.write_text("{}", encoding="utf-8")
            output = root / "packet"
            with self.assertRaises(FileNotFoundError):
                exporter.main([
                    "--design-json", str(design),
                    "--output", str(output),
                    "--neutral-sdf", str(root / "missing.sdf"),
                ])
            self.assertFalse(output.exists())

            with self.assertRaisesRegex(ValueError, "MAX_STATES_OUT_OF_RANGE"):
                exporter.main([
                    "--design-json", str(design),
                    "--output", str(output),
                    "--max-states", "17",
                ])
            self.assertFalse(output.exists())

    def test_wrong_scope_candidate_set_and_source_hash_leave_no_folder(self):
        base = {
            "parent_scope": {
                "actual_design_parent_id": exporter.PARENT_ID,
                "receptor_frame": exporter.STRUCTURAL_FRAME,
            },
            "analogs": [{"id": value, "selected": True, "qualified": True} for value in exporter.ANALOG_IDS],
            "protac_candidates": [],
        }
        wrong_scope = dict(base)
        wrong_scope["parent_scope"] = {
            "actual_design_parent_id": "OTHER",
            "receptor_frame": exporter.STRUCTURAL_FRAME,
        }
        with self.assertRaisesRegex(ValueError, "PARENT_SCOPE"):
            exporter._representatives(wrong_scope)

        wrong_ids = dict(base)
        wrong_ids["analogs"] = [{"id": "arbitrary", "selected": True}]
        with self.assertRaisesRegex(ValueError, "SELECTED_ANALOG_SET"):
            exporter._representatives(wrong_ids)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = {}
            for name in exporter.EXPECTED_SOURCE_HASHES:
                path = root / name
                path.write_bytes(b"wrong")
                paths[name] = path
            binding = {
                "design_sources/" + name: {
                    "status": "configured_hash_verified",
                    "sha256": digest,
                }
                for name, digest in exporter.EXPECTED_SOURCE_HASHES.items()
            }
            with self.assertRaisesRegex(ValueError, "ACTUAL_SOURCE_HASH_MISMATCH"):
                exporter._validate_source_bindings(binding, paths)
            self.assertFalse((root / "packet").exists())

    def test_store_mode_calls_real_service_shape_and_preserves_archive_bytes(self):
        archive_document = {
            "protac_candidates": [],
            "x": "exact spacing",
            "input_binding": {
                "source_files": {
                    "design_sources/" + name: {
                        "status": "configured_hash_verified",
                        "sha256": digest,
                    }
                    for name, digest in exporter.EXPECTED_SOURCE_HASHES.items()
                }
            },
        }
        archive = json.dumps(archive_document, separators=(",", ":")).encode("utf-8") + b"\n"
        calls = {}

        class FakeDB:
            def __enter__(self):
                calls["db"] = self
                return self

            def __exit__(self, exc_type, exc, traceback):
                return False

        class FakeStore:
            def __init__(self, path):
                calls["store_path"] = path

            def db(self):
                return FakeDB()

        class FakePort:
            def read(self, path):
                calls["read"] = path
                return archive

        class FakeService:
            def __init__(self, store, project):
                calls["project"] = project
                self.store = store
                self.port = FakePort()

            def _verified_result(self, db, job_id, require_current=False):
                calls["verified"] = (db, job_id, require_current)
                return (
                    {"id": job_id, "state": "stored"},
                    {"files": {"json": "immutable/archive.json"}},
                    {"input_sha256": "i", "result_sha256": "r"},
                    {"runtime": "current-source-fields"},
                )

        platform_module = types.ModuleType("packages.platform")
        acceptance_module = types.ModuleType("packages.platform.scientific_acceptance")
        store_module = types.ModuleType("packages.platform.store")
        acceptance_module.ScientificAcceptanceService = FakeService
        store_module.Store = FakeStore
        args = Args()
        args.store = "/store"
        args.job_id = "job-7"
        with mock.patch.dict(sys.modules, {
            "packages.platform": platform_module,
            "packages.platform.scientific_acceptance": acceptance_module,
            "packages.platform.store": store_module,
        }):
            document, raw, metadata, source = exporter._verified_design(args)
        self.assertEqual(raw, archive)
        self.assertEqual(document["x"], "exact spacing")
        self.assertEqual(calls["project"], "local-research")
        self.assertEqual(calls["read"], "immutable/archive.json")
        self.assertEqual(calls["verified"][1:], ("job-7", True))
        self.assertEqual(metadata["binding_status"], "server_verified_stored_job")
        self.assertEqual(metadata["state"], "stored")
        self.assertNotIn("job", metadata)
        self.assertNotIn("source", metadata)
        self.assertRegex(metadata["source_fingerprint"], r"^[0-9a-f]{64}$")
        self.assertEqual(source, archive_document["input_binding"]["source_files"])

    def test_offline_archive_requires_source_binding_and_is_not_authenticated(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "design.json"
            path.write_text(json.dumps({"protac_candidates": []}), encoding="utf-8")
            args = Args()
            args.design_json = str(path)
            with self.assertRaisesRegex(ValueError, "OFFLINE_INPUT_BINDING_REQUIRED"):
                exporter._verified_design(args)

            document = {
                "protac_candidates": [],
                "input_binding": {
                    "source_files": {
                        "design_sources/" + name: {
                            "status": "configured_hash_verified",
                            "sha256": digest,
                        }
                        for name, digest in exporter.EXPECTED_SOURCE_HASHES.items()
                    }
                },
            }
            raw = json.dumps(document, separators=(",", ":")).encode("utf-8")
            path.write_bytes(raw)
            _, observed, metadata, _ = exporter._verified_design(args)
            self.assertEqual(observed, raw)
            self.assertEqual(metadata["binding_status"], "offline_source_verified")
            self.assertNotEqual(metadata["binding_status"], "server_verified_stored_job")

    def test_source_copies_and_manifest_cover_every_nonmanifest_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            packet = Path(temporary) / "packet"
            sources = packet / "sources"
            sources.mkdir(parents=True)
            payloads = {
                "SMARCA2-parent.sdf": b"charged",
                "SMARCA2-neutral-design.sdf": b"neutral",
                "6HAZ.cif": b"protein",
                "design.json": b'{"portable":true}',
            }
            artifacts = []
            for name, value in payloads.items():
                path = sources / name
                path.write_bytes(value)
                artifacts.append({
                    "path": path.relative_to(packet).as_posix(),
                    "sha256": hashlib.sha256(value).hexdigest(),
                    "bytes": len(value),
                })
            (packet / "index.html").write_text("ok", encoding="utf-8")
            index = packet / "index.html"
            artifacts.append({
                "path": "index.html",
                "sha256": exporter._sha(index),
                "bytes": index.stat().st_size,
            })
            manifest = {"artifacts": artifacts}
            (packet / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

            loaded = json.loads((packet / "manifest.json").read_text(encoding="utf-8"))
            indexed = {row["path"]: row for row in loaded["artifacts"]}
            files = {
                path.relative_to(packet).as_posix(): path
                for path in packet.rglob("*")
                if path.is_file() and path.name != "manifest.json"
            }
            self.assertEqual(set(indexed), set(files))
            for relative, path in files.items():
                self.assertEqual(indexed[relative]["sha256"], exporter._sha(path))
                self.assertEqual(indexed[relative]["bytes"], path.stat().st_size)

    def test_complete_mock_cli_has_safe_links_bounded_names_and_six_schemes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            neutral = root / "SMARCA2-neutral-design.sdf"
            charged = root / "SMARCA2-parent.sdf"
            protein = root / "6HAZ.cif"
            for path in (neutral, charged, protein):
                path.write_bytes(path.name.encode("ascii"))

            candidates = []
            for analog in exporter.ANALOG_IDS:
                for e3 in exporter.E3_ORDER:
                    candidates.append({
                        "candidate_id": f"candidate-{analog}-{e3}<unsafe>",
                        "warhead_analog_id": analog,
                        "e3_type": e3,
                        "linker_id": exporter.LINKER_ID,
                    })
            document = {
                "parent_scope": {
                    "actual_design_parent_id": exporter.PARENT_ID,
                    "receptor_frame": exporter.STRUCTURAL_FRAME,
                },
                "analogs": [
                    {"id": value, "selected": True, "qualified": True, "pose_passed": True}
                    for value in exporter.ANALOG_IDS
                ],
                "protac_candidates": candidates,
                "input_binding": {
                    "source_files": {
                        "design_sources/" + name: {
                            "status": "configured_hash_verified",
                            "sha256": digest,
                        }
                        for name, digest in exporter.EXPECTED_SOURCE_HASHES.items()
                    }
                },
            }
            design = root / "design.json"
            design_raw = json.dumps(document, separators=(",", ":")).encode("utf-8")
            design.write_bytes(design_raw)
            output = root / "expert-packet-v2"

            hashes = dict(exporter.EXPECTED_SOURCE_HASHES)
            real_sha = exporter._sha

            def fake_sha(path):
                path = Path(path)
                if path.name in hashes:
                    return hashes[path.name]
                return real_sha(path)

            def fake_candidate(doc, identifier):
                return next(row for row in doc["protac_candidates"] if row["candidate_id"] == identifier)

            def fake_graph(candidate):
                return {"graph_sha256": hashlib.sha256(candidate["candidate_id"].encode()).hexdigest()}

            state_contact = {
                "protein_atom_id": "A:101:GLY:CA",
                "chain": "A",
                "auth_seq_id": "101",
                "residue": "GLY",
                "protein_atom": "CA",
                "ligand_atom_map": 7,
                "interaction_kind": "heavy_atom_contact",
                "heavy_distance_A": 3.2,
                "directional_hydrogen": None,
            }
            proposal_contact = {
                "protein_atom_id": "A:202:ASP:OD1",
                "chain": "A",
                "auth_seq_id": "202",
                "residue": "ASP",
                "protein_atom": "OD1",
                "ligand_atom_map": 11,
                "interaction_kind": "possible_hbond_contact",
                "heavy_distance_A": 2.9,
                "directional_hydrogen": None,
            }

            def fake_parent(*args, **kwargs):
                parent_dir = Path(kwargs["output_dir"])
                parent_dir.mkdir()
                (parent_dir / "index.html").write_text(
                    "<!doctype html><a href='fingerprint.json'>fingerprint</a>", encoding="utf-8"
                )
                packet = {
                    "format": "parent-fingerprint/20261002.2",
                    "protein_hydrogens_supplied": False,
                    "pH_context": 7.4,
                    "parents": [{
                        "parent_role": "source_charged_parent",
                        "states": [{"state_index": 1, "contacts": [state_contact]}],
                        "common_across_enumerated_states": [proposal_contact],
                        "charge_state_sensitive_exploratory_differences": [],
                        "proposed_protected_interaction_set": {
                            "status": "calculated_proposal_not_expert_approved",
                            "auto_site_role_promotion": False,
                            "auto_state_promotion": False,
                            "contacts": [proposal_contact],
                        },
                    }],
                }
                (parent_dir / "fingerprint.json").write_text(
                    json.dumps(packet), encoding="utf-8"
                )
                return packet

            fake_mol = Chem.MolFromSmiles("CC")
            with mock.patch.object(exporter, "_sha", side_effect=fake_sha), \
                 mock.patch.object(exporter, "_mol", return_value=fake_mol), \
                 mock.patch.object(exporter, "_protein", return_value=[{"label_asym_id": "A"}]), \
                 mock.patch.object(exporter, "_candidate", side_effect=fake_candidate), \
                 mock.patch.object(exporter.novel_ternary, "_mapped_graph", side_effect=fake_graph), \
                 mock.patch.object(exporter, "build_parent_fingerprint", side_effect=fake_parent), \
                 mock.patch.object(exporter, "propose_synthesis", return_value={"status": "proposal_only"}), \
                 mock.patch.object(exporter, "scheme_svg", return_value="<svg xmlns='http://www.w3.org/2000/svg'></svg>"):
                result = exporter.main([
                    "--design-json", str(design),
                    "--output", str(output),
                    "--neutral-sdf", str(neutral),
                    "--charged-sdf", str(charged),
                    "--protein-cif", str(protein),
                ])
            self.assertEqual(result, 0)
            self.assertEqual((output / "sources" / "design.json").read_bytes(), design_raw)
            self.assertEqual(len(list((output / "synthesis").glob("scheme-*.svg"))), 6)
            contacts_csv = (output / "parent-residue-contacts.csv").read_text(encoding="utf-8")
            self.assertIn("GLY", contacts_csv)
            self.assertIn("ASP", contacts_csv)
            self.assertIn("proposed_protected_interaction_set", contacts_csv)
            self.assertNotIn("proposed_common", contacts_csv)

            manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["binding_status"], "offline_source_verified")
            self.assertEqual(manifest["claims"]["readiness_state"], "science_review_pending")
            self.assertFalse(manifest["claims"]["expert_approved"])
            self.assertFalse(manifest["claims"]["scientifically_approved"])
            self.assertEqual(len(manifest["representatives"]), 6)
            for row in manifest["representatives"]:
                self.assertLessEqual(len(Path(row["report"]).stem), 64)
                self.assertRegex(row["report"], r"^[A-Za-z0-9_-]+\.json$")
                self.assertNotIn("<", row["report"])

            page = (output / "index.html").read_text(encoding="utf-8")
            self.assertNotIn("<unsafe>", page)
            self.assertIn("&lt;unsafe&gt;", page)
            exporter._validate_local_links(output)

            indexed = {row["path"]: row for row in manifest["artifacts"]}
            actual = {
                path.relative_to(output).as_posix(): path
                for path in output.rglob("*")
                if path.is_file() and path.name != "manifest.json"
            }
            self.assertEqual(set(indexed), set(actual))
            for relative, path in actual.items():
                self.assertEqual(indexed[relative]["sha256"], fake_sha(path))
                self.assertEqual(indexed[relative]["bytes"], path.stat().st_size)


if __name__ == "__main__":
    unittest.main()
