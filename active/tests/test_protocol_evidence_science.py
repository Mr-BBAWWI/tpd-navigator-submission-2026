from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from packages.science import protocol_evidence as evidence


class ProtocolEvidenceHelpersTest(unittest.TestCase):
    def test_safe_path_rejects_traversal_absolute_drive_and_backslash(self):
        bad = ["../x.json", "/x.json", "C:/x.json", "a\\b.json", "a/./b.json"]
        for value in bad:
            with self.subTest(value=value), self.assertRaises(evidence.ProtocolEvidenceError):
                evidence._safe_relative(value)

    def test_safe_path_accepts_normalized_allowed_path(self):
        self.assertEqual(evidence._safe_relative("evidence/a/plan.json"),
                         "evidence/a/plan.json")

    def test_json_rejects_duplicate_keys(self):
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError, "JSON_DUPLICATE_KEY"):
            evidence._json_bytes(b'{"a":1,"a":2}', "duplicate.json")

    def test_json_rejects_nan_and_infinity(self):
        for token in (b"NaN", b"Infinity", b"-Infinity", b"1e400"):
            with self.subTest(token=token), self.assertRaisesRegex(
                    evidence.ProtocolEvidenceError, "JSON_NONFINITE_NUMBER"):
                evidence._json_bytes(b'{"x":' + token + b"}", "number.json")

    def test_selection_is_confidence_only_with_low_index_tie(self):
        models = [
            {"model_index": index, "confidence": {"confidence_score": score}}
            for index, score in enumerate([0.8, 0.9, 0.9, 0.7, 0.6])
        ]
        # Geometry fields are deliberately irrelevant to the selector.
        models[4]["metrics"] = {"target_CA_RMSD_A": 0.0}
        self.assertEqual(evidence._selected_index(models), 1)

    def test_selection_rejects_bool_model_index(self):
        models = [
            {"model_index": index, "confidence": {"confidence_score": 0.5}}
            for index in range(5)
        ]
        models[1]["model_index"] = True
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError, "MODEL_INDEX_INVALID"):
            evidence._selected_index(models)

    def test_selection_rejects_missing_model(self):
        models = [
            {"model_index": index, "confidence": {"confidence_score": 0.5}}
            for index in range(4)
        ]
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                    "ALL_FIVE_MODELS_REQUIRED"):
            evidence._selected_index(models)

    def test_protocol_id_sources_are_independent(self):
        plan_hash = "a" * 64
        digest = "b" * 64
        self.assertNotEqual("calibration:" + plan_hash, "novel-msa:" + digest)
        self.assertTrue(("calibration:" + plan_hash).startswith("calibration:"))
        self.assertTrue(("novel-msa:" + digest).startswith("novel-msa:"))

    def test_strict_scalar_settings_reject_bool_for_integer(self):
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError, "SETTING"):
            evidence._strict_scalar(True, 1, "SETTING")
        evidence._strict_scalar(False, False, "SETTING")

    def test_calibration_receipt_accepts_actual_offline_msa_shape(self):
        settings = {"msa_mode": "cached_processed_offline"}
        original = dict(settings)
        evidence._verify_calibration_receipt_msa_settings(settings)
        self.assertEqual(settings, original)

        settings["use_msa_server"] = False
        evidence._verify_calibration_receipt_msa_settings(settings)
        settings["use_msa_server"] = True
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                    "CALIBRATION_RECEIPT_SETTINGS_INVALID:use_msa_server"):
            evidence._verify_calibration_receipt_msa_settings(settings)

        receipt = {"effective_command": [
            "boltz", "predict", "input.yaml",
            "--seed", "23", "--recycling_steps", "10",
            "--sampling_steps", "200", "--diffusion_samples", "5",
            "--max_parallel_samples", "1", "--max_msa_seqs", "256",
            "--model", "boltz2",
        ]}
        evidence._verify_effective_command(receipt, 23, "COMMAND")
        self.assertNotIn("--use_msa_server", receipt["effective_command"])

    def test_effective_command_rejects_duplicate_scalar_option(self):
        receipt = {"effective_command": [
            "boltz", "predict", "input.yaml",
            "--seed", "23", "--seed", "23",
            "--recycling_steps", "10", "--sampling_steps", "200",
            "--diffusion_samples", "5", "--max_parallel_samples", "1",
            "--max_msa_seqs", "256", "--model", "boltz2",
        ]}
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError, "COMMAND"):
            evidence._verify_effective_command(receipt, 23, "COMMAND")

    def test_seed_lists_and_receipt_seed_require_exact_ints(self):
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError, "SEEDS"):
            evidence._exact_int_list([True, 41, 61, 79, 97], evidence.SEEDS, "SEEDS")
        evidence._exact_int_list(list(evidence.SEEDS), evidence.SEEDS, "SEEDS")
        command = {"effective_command": [
            "boltz", "predict", "input.yaml",
            "--seed", "23", "--recycling_steps", "10",
            "--sampling_steps", "200", "--diffusion_samples", "5",
            "--max_parallel_samples", "1", "--max_msa_seqs", "256",
            "--model", "boltz2",
        ]}
        evidence._verify_effective_command(command, 23, "COMMAND")


class TrustedPlanSourceBindingTest(unittest.TestCase):
    @staticmethod
    def _source(path: str, role: str, sequence: str) -> dict:
        return {
            "path": path,
            "expected_sha256": "a" * 64,
            "label_asym_id": "A",
            "auth_asym_id": "A",
            "entity_id": "1",
            "role": role,
            "description": "synthetic test metadata; not a scientific fixture",
            "description_terms": ["synthetic"],
            "canonical_sequence": sequence,
            "sequence_length": len(sequence),
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "sequence_source": "_entity_poly.pdbx_seq_one_letter_code_can",
            "expected_deposit": "TEST",
            "observed_entry_id": "TEST",
            "excluded_proteins": [],
        }

    def _plan(self) -> tuple[dict, dict, dict]:
        graph = {
            "candidate_id": "candidate-test",
            "e3_type": "CRBN",
            "graph_sha256": "g" * 64,
            "actualmapped_smiles": "[C:1]",
        }
        binding = {
            "project": "project", "job_id": "job", "input_sha256": "i" * 64,
            "result_sha256": "r" * 64,
        }
        yaml_text = "marker: safe\n"
        plan = {
            "format": evidence.novel_ternary.PLAN_FORMAT,
            "candidate_graph": graph,
            "scientific_scope": {
                "reference_free": True,
                "experimental_ternary_reference": None,
            },
            "sources": {
                "target": self._source("/hostile/absolute/target.cif", "target", "AAAA"),
                "e3": self._source("C:/hostile/e3.cif", "e3", "CCCC"),
            },
            "bindings": {"job": {
                **binding,
                "candidate_graph_sha256": graph["graph_sha256"],
            }},
            "seeds": list(evidence.SEEDS),
            "msa_mode": "server",
            "boltz_input": {
                "yaml": yaml_text,
                "sha256": evidence._sha_bytes(yaml_text.encode("utf-8")),
                "templates": [],
                "constraints": [],
                "potentials": False,
            },
        }
        plan["plan_digest"] = evidence.novel_ternary._digest(plan)
        return plan, binding, graph

    def test_hostile_declared_paths_are_inert_and_fixed_paths_are_validated(self):
        plan, binding, graph = self._plan()
        calls = []

        def validate(source, expected):
            calls.append((source, expected))
            return dict(source)

        with mock.patch.object(evidence, "_trusted_regular_source",
                               side_effect=lambda path, root: path.resolve()), \
                mock.patch.object(evidence.novel_ternary, "_validate_source",
                                  side_effect=validate), \
                mock.patch.object(evidence.novel_ternary, "_mapped_graph",
                                  return_value=graph), \
                mock.patch.object(evidence.novel_ternary, "_document",
                                  return_value={"marker": "safe"}):
            evidence._verify_plan_portable(plan, binding, graph)

        trusted_root = Path(evidence.__file__).resolve().parents[2]
        self.assertEqual([item[1] for item in calls], ["target", "CRBN"])
        self.assertEqual(Path(calls[0][0]["path"]),
                         trusted_root / "cases" / "design_sources" / "6HAZ.cif")
        self.assertEqual(Path(calls[1][0]["path"]),
                         trusted_root / "cases" / "design_sources" / "6BOY.cif")
        self.assertNotIn("hostile", calls[0][0]["path"])
        self.assertNotIn("hostile", calls[1][0]["path"])

    def test_wrong_declared_canonical_sequence_is_rejected_against_trusted_result(self):
        plan, binding, graph = self._plan()

        def validate(source, expected):
            result = dict(source)
            if expected == "target":
                result["canonical_sequence"] = "TRUSTED"
            return result

        with mock.patch.object(evidence, "_trusted_regular_source",
                               side_effect=lambda path, root: path.resolve()), \
                mock.patch.object(evidence.novel_ternary, "_validate_source",
                                  side_effect=validate), \
                mock.patch.object(evidence.novel_ternary, "_mapped_graph",
                                  return_value=graph):
            with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                        "NOVEL_SOURCE_BINDING_MISMATCH:target:canonical_sequence"):
                evidence._verify_plan_portable(plan, binding, graph)

    def test_scientific_scope_must_remain_reference_free(self):
        plan, binding, graph = self._plan()
        plan["scientific_scope"]["reference_free"] = False
        plan["plan_digest"] = evidence.novel_ternary._digest(
            {key: value for key, value in plan.items() if key != "plan_digest"}
        )
        with mock.patch.object(evidence.novel_ternary, "_mapped_graph", return_value=graph):
            with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                        "NOVEL_SCIENTIFIC_SCOPE_INVALID"):
                evidence._verify_plan_portable(plan, binding, graph)

    def test_trusted_source_symlink_is_rejected_where_supported(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            target = root / "target.cif"
            link = root / "link.cif"
            target.write_text("synthetic non-scientific test data", encoding="utf-8")
            try:
                link.symlink_to(target)
            except (OSError, NotImplementedError):
                self.skipTest("symlinks are unavailable on this platform")
            with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                        "TRUSTED_SOURCE_SYMLINK_OR_JUNCTION_FORBIDDEN"):
                evidence._trusted_regular_source(link, root)


class ManifestPreflightTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    @staticmethod
    def _digest(data: bytes) -> str:
        return hashlib.sha256(data).hexdigest()

    def _write_manifest(self, records):
        value = {
            "format": evidence.MANIFEST_FORMAT,
            "files": records,
            "manifest_excludes_itself": True,
        }
        encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
        (self.root / evidence.MANIFEST_NAME).write_bytes(encoded)
        return self._digest(encoded)

    def test_hash_pin_and_complete_file_set(self):
        content = b"{}"
        (self.root / "item.json").write_bytes(content)
        pin = self._write_manifest([{
            "path": "item.json", "sha256": self._digest(content), "bytes": len(content)
        }])
        _manifest, files, normalized = evidence._verify_manifest(self.root, pin)
        self.assertEqual(set(files), {"item.json", evidence.MANIFEST_NAME})
        self.assertEqual({item["path"] for item in normalized}, set(files))
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                    "MANIFEST_HASH_PIN_MISMATCH"):
            evidence._verify_manifest(self.root, "0" * 64)

    def test_missing_manifest_member_rejected(self):
        pin = self._write_manifest([{
            "path": "missing.json", "sha256": "0" * 64, "bytes": 2
        }])
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                    "PACK_FILE_SET_MISMATCH"):
            evidence._verify_manifest(self.root, pin)

    def test_extra_file_rejected(self):
        content = b"{}"
        (self.root / "item.json").write_bytes(content)
        (self.root / "extra.json").write_bytes(content)
        pin = self._write_manifest([{
            "path": "item.json", "sha256": self._digest(content), "bytes": len(content)
        }])
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                    "PACK_FILE_SET_MISMATCH"):
            evidence._verify_manifest(self.root, pin)

    def test_casefold_duplicate_manifest_paths_rejected(self):
        pin = self._write_manifest([
            {"path": "A.json", "sha256": "0" * 64, "bytes": 0},
            {"path": "a.json", "sha256": "0" * 64, "bytes": 0},
        ])
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                    "MANIFEST_PATH_DUPLICATE"):
            evidence._verify_manifest(self.root, pin)

    def test_manifest_traversal_and_absolute_paths_rejected(self):
        for name in ("../bad.json", "/bad.json", "C:/bad.json"):
            with self.subTest(name=name):
                for child in self.root.iterdir():
                    child.unlink()
                pin = self._write_manifest([{
                    "path": name, "sha256": "0" * 64, "bytes": 0
                }])
                with self.assertRaises(evidence.ProtocolEvidenceError):
                    evidence._verify_manifest(self.root, pin)

    def test_wrong_declared_hash_rejected(self):
        content = b"{}"
        (self.root / "item.json").write_bytes(content)
        pin = self._write_manifest([{
            "path": "item.json", "sha256": "0" * 64, "bytes": len(content)
        }])
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                    "FILE_HASH_MISMATCH"):
            evidence._verify_manifest(self.root, pin)


class BindingAndProvenanceTest(unittest.TestCase):
    def test_wrong_job_binding_and_parent_graph_rejected(self):
        expected_binding = {
            "project": "p", "job_id": "j", "input_sha256": "i",
            "result_sha256": "r", "parent_id": "parent",
        }
        graph = {"candidate_id": "c", "warhead_analog_id": "wrong",
                 "graph_sha256": "g"}
        plan = {
            "format": "wrong",
            "candidate_graph": graph,
        }
        with self.assertRaises(evidence.ProtocolEvidenceError):
            evidence._verify_plan_portable(plan, expected_binding, graph)

    def test_expected_binding_requires_exact_fields(self):
        binding = {
            "project": "p", "job_id": "j", "input_sha256": "i" * 64,
            "result_sha256": "r" * 64, "parent_id": "x",
            "result_binding_kind": "exact_archived_design_json_bytes", "extra": "forbidden",
        }
        with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                    "EXPECTED_BINDING_FIELDS_INVALID"):
            evidence._validate_expected_inputs(binding, {"a": {}, "b": {}})

    def test_provenance_hash_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "evidence.json"
            source.write_text("{}", encoding="utf-8")
            provenance = root / evidence.PROVENANCE_NAME
            provenance.write_text(json.dumps({
                "record": {
                    "portable_pack_relative_path": "evidence.json",
                    "source_sha256": "0" * 64,
                }
            }), encoding="utf-8")
            compact = {
                "source_provenance": {
                    "portable_pack_relative_path": evidence.PROVENANCE_NAME,
                    "sha256": evidence._sha_file(provenance),
                },
                "sources": {
                    "x": {
                        "portable_pack_relative_path": "evidence.json",
                        "source_sha256": "0" * 64,
                    }
                },
            }
            with self.assertRaisesRegex(evidence.ProtocolEvidenceError,
                                        "SOURCE_PROVENANCE_HASH_MISMATCH"):
                evidence._verify_provenance(compact, {
                    evidence.PROVENANCE_NAME: provenance,
                    "evidence.json": source,
                })


if __name__ == "__main__":
    unittest.main()
