from __future__ import annotations

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

from scripts import run_novel_msa_protocol as protocol


def completed_receipt(seed, iptm=0.7, endpoint=3.0):
    return {
        "seed": seed,
        "status": "completed",
        "actual_computation": True,
        "settings": copy.deepcopy(protocol.EXPECTED_SETTINGS),
        "selected_model_inspection": {
            "ligand_iptm": 0.01,
            "model_confidence": {"ligand_iptm": 0.02, "iptm": iptm},
            "descriptive_metrics": {
                "other_endpoint": 999.0,
                "linker_endpoint_distance_A": endpoint,
                "contacts_by_protein_role_and_ligand_group": {
                    "target": {"warhead": {"count": 1}},
                    "e3": {"recruiter": {"count": 1}},
                },
            },
        },
    }


class NovelMsaProtocolTests(unittest.TestCase):
    def test_make_plan_original_unchanged_and_hash_recomputed(self):
        old = {
            "sources": {
                "target": {"canonical_sequence": "TARGET"},
                "e3": {"canonical_sequence": "ESEQUENCE"},
            },
            "candidate_graph": {"actualmapped_smiles": "CC"},
            "msa_mode": "single_sequence",
            "boltz_input": {
                "yaml": "old\n",
                "sha256": hashlib.sha256(b"old\n").hexdigest(),
                "templates": ["old"],
                "constraints": ["old"],
                "potentials": True,
            },
            "plan_digest": "old-digest",
        }
        original = copy.deepcopy(old)
        verified = []

        with patch.object(
            protocol.novel,
            "verify_plan",
            side_effect=lambda value: verified.append(copy.deepcopy(value)),
        ), patch.object(
            protocol.novel,
            "_document",
            side_effect=lambda target, e3, smiles, mode: {
                "target": target,
                "e3": e3,
                "smiles": smiles,
                "mode": mode,
            },
        ), patch.object(protocol.novel, "_digest", return_value="new-digest"):
            result = protocol.make_plan(old)

        self.assertEqual(old, original)
        self.assertIsNot(result, old)
        self.assertEqual(result["msa_mode"], "server")
        self.assertEqual(
            result["boltz_input"]["sha256"],
            hashlib.sha256(
                result["boltz_input"]["yaml"].encode("utf-8")
            ).hexdigest(),
        )
        self.assertEqual(
            yaml.safe_load(result["boltz_input"]["yaml"])["mode"], "server"
        )
        self.assertEqual(
            verified[-1]["boltz_input"]["sha256"],
            result["boltz_input"]["sha256"],
        )

    def test_csv_query_exact_sequence_and_gap_semantics(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            csv_path = Path(temporary_directory) / "query.csv"
            csv_path.write_text("sequence,key\nAC-D E,query\n", encoding="utf-8")
            result = protocol.csv_query(csv_path, "ACDE")
            self.assertEqual(result["query_sequence"], "ACDE")

            with self.assertRaisesRegex(
                protocol.ProtocolError, "MSA_QUERY_SEQUENCE_MISMATCH"
            ):
                protocol.csv_query(csv_path, "BRD4")

    def test_expert_metrics_uses_exact_iptm_and_endpoint_paths(self):
        receipts = [
            completed_receipt(seed, 0.70 + index * 0.01, 3.0 + index * 0.1)
            for index, seed in enumerate(protocol.SEEDS)
        ]
        result = protocol.expert_metrics(receipts)
        self.assertEqual(result["finite_selected_ipTM_values"][0], 0.70)
        self.assertEqual(result["finite_selected_endpoint_values_A"][0], 3.0)
        self.assertIs(result["quantitative_criteria_met"], True)
        self.assertIs(result["automatic_approval"], False)

    def test_failure_missing_or_duplicate_seed_cannot_pass(self):
        receipts = [completed_receipt(seed) for seed in protocol.SEEDS]
        receipts[0]["status"] = "failed"
        self.assertIs(
            protocol.expert_metrics(receipts)["quantitative_criteria_met"], False
        )

        missing = [completed_receipt(seed) for seed in protocol.SEEDS[:-1]]
        self.assertIs(
            protocol.expert_metrics(missing)["quantitative_criteria_met"], False
        )

        duplicate = [completed_receipt(seed) for seed in protocol.SEEDS]
        duplicate[-1]["seed"] = protocol.SEEDS[0]
        self.assertIs(
            protocol.expert_metrics(duplicate)["quantitative_criteria_met"], False
        )

    def test_boolean_and_out_of_range_expert_values_are_rejected(self):
        receipts = [completed_receipt(seed) for seed in protocol.SEEDS]
        receipts[0]["selected_model_inspection"]["model_confidence"]["iptm"] = True
        self.assertIs(
            protocol.expert_metrics(receipts)["quantitative_criteria_met"], False
        )
        receipts[0]["selected_model_inspection"]["model_confidence"]["iptm"] = 1.1
        self.assertIs(
            protocol.expert_metrics(receipts)["quantitative_criteria_met"], False
        )

    def test_failed_inspection_is_not_marked_valid(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            tmp_path = Path(temporary_directory)
            cif = tmp_path / "novel_ternary_model_0.cif"
            confidence = tmp_path / "confidence_novel_ternary_model_0.json"
            cif.write_text("data_test\n", encoding="utf-8")
            confidence.write_text(
                json.dumps({"confidence_score": 0.5}), encoding="utf-8"
            )

            with patch.object(
                protocol.novel,
                "assess_novel_outputs",
                return_value={"status": "failed", "actual_computation": False},
            ):
                with self.assertRaisesRegex(
                    protocol.ProtocolError,
                    "INSPECTION_NOT_COMPUTED_HYPOTHESIS",
                ):
                    protocol.inspect_model(
                        tmp_path / "run", 0, cif, confidence, {}
                    )

    def test_inspection_requires_actual_computation(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            tmp_path = Path(temporary_directory)
            cif = tmp_path / "model.cif"
            confidence = tmp_path / "confidence.json"
            cif.write_text("data_test\n", encoding="utf-8")
            confidence.write_text("{}\n", encoding="utf-8")

            with patch.object(
                protocol.novel,
                "assess_novel_outputs",
                return_value={
                    "status": "computed_hypothesis",
                    "actual_computation": False,
                },
            ):
                with self.assertRaisesRegex(
                    protocol.ProtocolError,
                    "INSPECTION_ACTUAL_COMPUTATION_REQUIRED",
                ):
                    protocol.inspect_model(
                        tmp_path / "run", 0, cif, confidence, {}
                    )

    def test_fixed_settings_all_seeds_and_no_automatic_approval(self):
        self.assertEqual(protocol.SEEDS, [23, 41, 61, 79, 97])
        self.assertEqual(len(set(protocol.SEEDS)), 5)
        self.assertEqual(
            protocol.EXPECTED_SETTINGS,
            {
                "accelerator": "gpu",
                "recycling_steps": 10,
                "sampling_steps": 200,
                "diffusion_samples": 5,
                "max_parallel_samples": 1,
                "max_msa_seqs": 256,
                "use_msa_server": False,
                "msa_pairing_strategy": "greedy",
                "preprocessing_threads": 1,
                "potentials": False,
                "templates": [],
                "constraints": [],
            },
        )
        metrics = protocol.expert_metrics(
            [completed_receipt(seed) for seed in protocol.SEEDS]
        )
        self.assertIs(metrics["automatic_approval"], False)

    def test_allowed_manifest_change_is_exact(self):
        self.assertIs(protocol.allowed_manifest_change("processed/manifest.json"), True)
        self.assertIs(
            protocol.allowed_manifest_change("processed/other-manifest.json"), False
        )
        self.assertIs(
            protocol.allowed_manifest_change("processed/subdir/manifest.json"), False
        )
        self.assertIs(protocol.allowed_manifest_change("msa/manifest.json"), False)


if __name__ == "__main__":
    unittest.main()
