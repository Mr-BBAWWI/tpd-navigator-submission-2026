from __future__ import annotations

import importlib.util
import json
import math
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest import mock

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_calibration_protocol.py"
spec = importlib.util.spec_from_file_location("run_calibration_protocol", MODULE_PATH)
protocol = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(protocol)


def models(scores, rmsds=None):
    rmsds = rmsds or [999] * 5
    return [{"model_index": index, "confidence": {"confidence_score": score},
             "reference_rmsd_forbidden_to_selection": rmsds[index]}
            for index, score in enumerate(scores)]


class CalibrationProtocolTests(unittest.TestCase):
    def test_selection_uses_confidence_not_rmsd_and_tie_smaller_index(self):
        value = models([0.1, 0.9, 0.9, 0.2, 0.3], [0.01, 50, 0.001, 2, 3])
        self.assertEqual(protocol.select_model(value), 1)

    def test_selection_rejects_nonfinite_or_nonnumeric_confidence(self):
        for bad in (float("nan"), float("inf"), float("-inf"), None, True, "0.9"):
            with self.subTest(bad=bad):
                value = models([0.1, 0.2, 0.3, 0.4, 0.5])
                value[2]["confidence"]["confidence_score"] = bad
                with self.assertRaisesRegex(Exception, "FINITE_CONFIDENCE_SCORE_REQUIRED"):
                    protocol.select_model(value)

    def test_all_four_boundaries_and_conjunction(self):
        exact = {
            "target_CA_RMSD_A": 1.5,
            "ligand_heavy_atom_RMSD_after_target_alignment_A": 3.0,
            "e3_CA_RMSD_after_target_alignment_A": 10.0,
            "contact_jaccard": 0.4,
        }
        self.assertTrue(protocol.metrics_pass(exact))
        for key in exact:
            with self.subTest(key=key):
                failed = dict(exact)
                failed[key] = exact[key] - 0.001 if key == "contact_jaccard" else exact[key] + 0.001
                self.assertFalse(protocol.metrics_pass(failed))
        self.assertFalse(protocol.metrics_pass({**exact, "contact_jaccard": math.nan}))

    def test_metrics_reject_booleans_negative_values_and_jaccard_over_one(self):
        valid = {
            "target_CA_RMSD_A": 1.0,
            "ligand_heavy_atom_RMSD_after_target_alignment_A": 2.0,
            "e3_CA_RMSD_after_target_alignment_A": 5.0,
            "contact_jaccard": 0.5,
        }
        for key in valid:
            with self.subTest(boolean_key=key):
                changed = dict(valid)
                changed[key] = True
                self.assertFalse(protocol.metrics_pass(changed))
            with self.subTest(negative_key=key):
                changed = dict(valid)
                changed[key] = -0.01
                self.assertFalse(protocol.metrics_pass(changed))
        self.assertFalse(protocol.metrics_pass({**valid, "contact_jaccard": 1.01}))

    def test_fixed_exact_seed_coverage_and_failure_cannot_pass(self):
        self.assertEqual(protocol.SEEDS, [23, 41, 61, 79, 97])
        receipts = [{"seed": seed, "status": "success", "selected_model_pass": index < 3}
                    for index, seed in enumerate(protocol.SEEDS)]
        self.assertTrue(protocol.protocol_pass(receipts))
        receipts[4]["status"] = "failure"
        self.assertFalse(protocol.protocol_pass(receipts))
        self.assertFalse(protocol.protocol_pass(receipts[:4]))

    def test_protocol_rejects_duplicate_seed_even_with_five_receipts(self):
        receipts = [{"seed": seed, "status": "success", "selected_model_pass": True}
                    for seed in [23, 41, 61, 79, 79]]
        self.assertFalse(protocol.protocol_pass(receipts))

    def test_effective_command_is_fixed_offline_and_has_no_forbidden_features(self):
        base = ["boltz", "predict", "input.yaml", "--recycling_steps", "3",
                "--sampling_steps", "200", "--diffusion_samples", "1",
                "--max_parallel_samples", "1", "--max_msa_seqs", "256",
                "--use_msa_server", "--seed", "23"]
        command = protocol.effective_command(base, "options: --model --seed")
        self.assertEqual(command[command.index("--recycling_steps") + 1], "10")
        self.assertEqual(command[command.index("--diffusion_samples") + 1], "5")
        self.assertNotIn("--use_msa_server", command)
        self.assertEqual(command[-2:], ["--model", "boltz2"])
        self.assertFalse(any(word in " ".join(command).lower()
                             for word in ("potential", "template", "constraint")))

    def test_plan_immutability_rejects_existing_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "already-there"
            output.mkdir()
            rc = protocol.main([
                "--output", str(output), "--frozen-seed-directory", str(root / "frozen"),
                "--boltz-executable", str(root / "boltz"),
                "--checkpoint", str(root / "checkpoint"),
                "--cache", str(root / "cache"), "--prepare-only",
            ])
            self.assertEqual(rc, 2)
            self.assertEqual(list(output.iterdir()), [])

    def test_plan_freezes_source_receipt_runner_hash_and_start_time(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            frozen = root / "frozen"
            frozen.mkdir()
            receipt = frozen / "receipt.json"
            receipt.write_text('{"source": true}\n', encoding="utf-8")
            paths = {name: root / name for name in
                     ("reference", "metadata", "checkpoint", "boltz_executable")}
            for name, path in paths.items():
                path.write_text(name, encoding="utf-8")
            paths["frozen_seed_directory"] = frozen
            plan = protocol._plan(
                Namespace(timeout_seconds=1.0), paths,
                {"sha256": "help-hash"}, {"tool": "environment"},
            )
            self.assertIn("started_utc", plan)
            self.assertEqual(plan["sources"]["runner_module_sha256"], protocol.runner_hash())
            self.assertEqual(
                plan["sources"]["frozen_source_receipt"]["sha256"],
                protocol.sha256_file(receipt),
            )

    def test_original_source_settings_are_strict_and_separate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ("reference", "metadata", "boltz", "checkpoint"):
                (root / name).write_text(name, encoding="utf-8")
            captured = {}

            def prepare(reference, metadata, destination, msa_mode):
                destination.write_text("official", encoding="utf-8")
                return {"input_sha256": protocol.worker.FROZEN_CRBN_INPUT_SHA256}

            def install(**kwargs):
                captured["settings"] = kwargs["settings"]
                raise RuntimeError("stop after provenance check")

            run_dir = root / "run"
            paths = {
                "reference": root / "reference", "metadata": root / "metadata",
                "boltz_executable": root / "boltz", "checkpoint": root / "checkpoint",
                "cache": root / "cache", "frozen_seed_directory": root / "frozen",
            }
            args = Namespace(timeout_seconds=1.0)
            with mock.patch.object(protocol.worker, "prepare_input", side_effect=prepare), \
                    mock.patch.object(protocol.worker, "_install_frozen_seed_cache", side_effect=install):
                receipt = protocol._run_seed(
                    23, run_dir, paths, args,
                    {"sha256": "help", "text": "", "has_seed": True}, {},
                )
            self.assertEqual(receipt["status"], "failure")
            self.assertEqual(captured["settings"], protocol.SOURCE_SETTINGS)
            self.assertEqual(captured["settings"]["recycling_steps"], 3)
            self.assertEqual(captured["settings"]["diffusion_samples"], 1)
            self.assertEqual(protocol.PROTOCOL_SETTINGS["recycling_steps"], 10)
            self.assertEqual(protocol.PROTOCOL_SETTINGS["diffusion_samples"], 5)
            stored = json.loads((run_dir / "receipt.json").read_text(encoding="utf-8"))
            self.assertIs(stored["selected_model_pass"], False)


if __name__ == "__main__":
    unittest.main()
