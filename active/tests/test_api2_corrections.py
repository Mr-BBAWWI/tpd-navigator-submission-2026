from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

from packages.science import benchmark_distribution as distribution
from packages.science import boltz_worker
from packages.science.benchmark_distribution import BenchmarkDistributionError
from scripts import run_expert_priority_batch as batch


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def success(
    seed: int,
    e3: float,
    ligand: float,
    contact: float,
    confidence: float = 0.8,
) -> dict:
    return {
        "seed": seed,
        "actual_execution_success": True,
        "hashes": {
            "input_yaml_sha256": "i",
            "checkpoint_sha256": "c",
            "processed_msa_files_sha256": {"x/processed/msa/a.npz": "m"},
        },
        "settings": {"sampling_steps": 200},
        "tool_environment": {
            "boltz_version": "2.2.1",
            "torch_version": "2.7.1",
            "cuda_available": True,
            "cuda_runtime": "12.8",
            "cuda_devices": ["GPU"],
        },
        "inspection": {
            "status": "success",
            "model_confidence": {
                "confidence_score": confidence,
                "nested": {"iptm": confidence - 0.1},
            },
            "comparison": {
                "metrics": {
                    "e3_CA_RMSD_after_target_alignment_A": e3,
                    "ligand_heavy_atom_RMSD_after_target_alignment_A": ligand,
                    "contact_jaccard": contact,
                }
            },
        },
    }


def failed(seed: int) -> dict:
    return {
        "seed": seed,
        "actual_execution_success": False,
        "status": "failure",
        "failure_reason": "NO_OUTPUT",
    }


def candidate_document() -> dict:
    rows = []
    for analog in batch.ANALOG_IDS:
        for e3 in batch.E3_ORDER:
            rows.append(
                {
                    "candidate_id": analog + "-" + e3,
                    "warhead_analog_id": analog,
                    "e3_type": e3,
                    "linker_id": batch.LINKER_ID,
                }
            )
    return {
        "parent_scope": {
            "actual_design_parent_id": "SMARCA2-FX5",
            "receptor_frame": "6HAZ chain A",
        },
        "protac_candidates": rows,
    }


@contextmanager
def frozen_fixture(root: Path, extra_path: str | None = None):
    source = root / "source"
    output = source / "boltz_output"
    files = {
        "boltz_results_6BOY_CRBN/msa/raw.csv": b"raw",
        "boltz_results_6BOY_CRBN/processed/msa/a.npz": b"msa",
        "boltz_results_6BOY_CRBN/processed/records/6BOY_CRBN.json": (
            b'{"case":"6BOY_CRBN","msa":true}'
        ),
        "boltz_results_6BOY_CRBN/processed/manifest.json": b"manifest",
        "boltz_results_6BOY_CRBN/predictions/model.cif": b"prediction",
        "boltz_results_6BOY_CRBN/lightning_logs/hparams.yaml": b"log",
    }
    if extra_path:
        files[extra_path] = b"bad"

    hashes = {}
    for name, data in files.items():
        path = output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        hashes[name] = digest(data)

    input_yaml = root / "input.yaml"
    input_yaml.write_bytes(b"input")
    (source / "6BOY_CRBN.yaml").write_bytes(b"input")

    reference = root / "reference.cif"
    checkpoint = root / "checkpoint"
    executable = root / "boltz"
    reference.write_bytes(b"reference")
    checkpoint.write_bytes(b"checkpoint")
    executable.write_bytes(b"executable")

    receipt = {
        "exit_code": 0,
        "process_state": "completed",
        "status": "failure",
        "inspection": {"status": "inspection_failed"},
        "settings": {"x": 1},
        "tool_environment": {"boltz_version": "2.2.1"},
        "tool": {"help_sha256": "help"},
        "hashes": {
            "input_yaml_sha256": digest(b"input"),
            "checkpoint_sha256": digest(b"checkpoint"),
            "reference_sha256": digest(b"reference"),
            "executable_sha256": digest(b"executable"),
            "processed_output_files_sha256": hashes,
        },
    }
    (source / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")

    with patch.object(
        boltz_worker, "FROZEN_CRBN_INPUT_SHA256", digest(b"input")
    ), patch.object(
        boltz_worker, "FROZEN_CRBN_CHECKPOINT_SHA256", digest(b"checkpoint")
    ):
        yield source, input_yaml, reference, checkpoint, executable


@contextmanager
def installed_fixture(root: Path, extra_path: str | None = None):
    with frozen_fixture(root, extra_path) as fixture:
        source, input_yaml, reference, checkpoint, executable = fixture
        destination = root / "destination"
        result = boltz_worker._install_frozen_seed_cache(
            source_directory=source,
            destination=destination,
            input_yaml=input_yaml,
            reference=reference,
            checkpoint=checkpoint,
            executable=executable,
            environment={"boltz_version": "2.2.1"},
            settings={"x": 1},
            help_sha256="help",
        )
        yield source, destination, result


class TestApi2Corrections(unittest.TestCase):
    def test_distribution_retains_poor_seed_41_and_uses_actual_metrics(self):
        result = distribution.build_benchmark_distribution(
            [],
            [
                success(23, 5, 1, 0.5),
                success(41, 39, 12, 0.1),
                success(61, 3, 1, 0.6),
            ],
        )
        self.assertEqual([row["seed"] for row in result["seeds"]], [23, 41, 61])
        self.assertEqual(
            result["statistics"][
                "comparison.metrics.e3_CA_RMSD_after_target_alignment_A"
            ]["median"],
            5,
        )
        self.assertIs(result["poor_replicates_retained"], True)

    def test_distribution_iqr_is_linear_and_does_not_drop_outlier(self):
        result = distribution.build_benchmark_distribution(
            [],
            [
                success(1, 1, 1, 0.1),
                success(2, 2, 1, 0.2),
                success(3, 3, 1, 0.3),
                success(4, 100, 1, 0.4),
            ],
        )
        stats = result["statistics"][
            "comparison.metrics.e3_CA_RMSD_after_target_alignment_A"
        ]
        self.assertEqual(stats["count"], 4)
        self.assertAlmostEqual(stats["iqr"], 25.5)

    def test_distribution_retains_failed_seed_with_unknown_comparability(self):
        result = distribution.build_benchmark_distribution(
            [], [success(1, 2, 3, 0.4), failed(2)]
        )
        self.assertEqual(result["seed_count"], 2)
        self.assertEqual(result["technical_failure_count"], 1)
        self.assertIs(result["comparability"]["comparability_unknown"], True)

    def test_distribution_rejects_mixed_computed_input(self):
        left = success(1, 2, 3, 0.4)
        right = success(2, 2, 3, 0.4)
        right["hashes"]["input_yaml_sha256"] = "other"
        with self.assertRaisesRegex(BenchmarkDistributionError, "MIXED_INPUT"):
            distribution.build_benchmark_distribution([], [left, right])

    def test_distribution_rejects_mixed_processed_msa(self):
        left = success(1, 2, 3, 0.4)
        right = success(2, 2, 3, 0.4)
        right["hashes"]["processed_msa_files_sha256"] = {
            "x/processed/msa/a.npz": "other"
        }
        with self.assertRaisesRegex(BenchmarkDistributionError, "MIXED_MSA"):
            distribution.build_benchmark_distribution([], [left, right])

    def test_distribution_rejects_numeric_execution_flag(self):
        item = failed(1)
        item["actual_execution_success"] = 0
        with self.assertRaisesRegex(
            BenchmarkDistributionError, "ACTUAL_EXECUTION_SUCCESS_NOT_BOOLEAN"
        ):
            distribution.build_benchmark_distribution([], [item])

    def test_forged_in_memory_frozen_base_receipt_is_rejected(self):
        item = success(23, 5, 1, 0.5)
        item.pop("actual_execution_success")
        item["execution_success"] = True
        item["original_receipt_ref"] = "receipt.json"
        item["reinspection_ref"] = "reinspection.json"
        base = {
            "seed_receipts": [item],
            "new_MSA_baseline": {
                "processed_msa_hashes_by_seed": {"23": {"a.npz": "m"}}
            },
        }
        with self.assertRaisesRegex(
            BenchmarkDistributionError,
            "FROZEN_BASE_EXECUTION_PROVENANCE_INVALID",
        ):
            distribution.build_benchmark_distribution(base, [])

    def test_seed_23_requires_exit0_and_source_bound_reinspection(self):
        with tempfile.TemporaryDirectory() as temporary:
            source_root = Path(temporary) / "cases" / "acceptance_sources"
            seed_root = source_root / "supporting" / "seeds" / "seed-23"
            seed_root.mkdir(parents=True)
            original = {
                "exit_code": 0,
                "process_state": "completed",
                "status": "failure",
                "inspection": {"status": "inspection_failed"},
            }
            original_path = seed_root / "receipt.json"
            original_path.write_text(json.dumps(original), encoding="utf-8")
            original_hash = distribution._sha256(original_path)
            reinspection_path = seed_root / "reinspection.json"
            reinspection_path.write_text(
                json.dumps({
                    "source_receipt_sha256": original_hash,
                    "assessment": {"status": "success"},
                }),
                encoding="utf-8",
            )
            row = success(23, 5, 1, 0.5)
            row.pop("actual_execution_success")
            row.update({
                "execution_success": True,
                "original_receipt_ref": "supporting/seeds/seed-23/receipt.json",
                "original_receipt_sha256": original_hash,
                "reinspection_ref": "supporting/seeds/seed-23/reinspection.json",
            })
            base_document = {
                "seed_receipts": [row],
                "new_MSA_baseline": {
                    "processed_msa_hashes_by_seed": {"23": {"a.npz": "m"}}
                },
            }
            base_path = source_root / "crbn_calibration.json"
            base_path.write_text(json.dumps(base_document), encoding="utf-8")
            with patch.object(
                distribution, "KNOWN_BASE_SHA256", distribution._sha256(base_path)
            ):
                result = distribution.build_benchmark_distribution(base_path, [])
            self.assertEqual(result["technical_success_count"], 1)
            self.assertEqual(
                result["seeds"][0]["execution_status"],
                "frozen_base_hash_verified_reinspection",
            )

            reinspection_path.write_text(
                json.dumps({
                    "note": original_hash,
                    "assessment": {"status": "success"},
                }),
                encoding="utf-8",
            )
            with patch.object(
                distribution, "KNOWN_BASE_SHA256", distribution._sha256(base_path)
            ):
                with self.assertRaisesRegex(
                    BenchmarkDistributionError,
                    "REINSPECTION_SOURCE_BINDING_INVALID",
                ):
                    distribution.build_benchmark_distribution(base_path, [])

    def test_seeds_41_and_61_accept_verified_original_without_reinspection(self):
        with tempfile.TemporaryDirectory() as temporary:
            source_root = Path(temporary) / "cases" / "acceptance_sources"
            rows = []
            baseline = {}
            for seed in (41, 61):
                seed_root = source_root / "supporting" / "seeds" / f"seed-{seed}"
                seed_root.mkdir(parents=True)
                original = success(seed, 39 if seed == 41 else 3, 12 if seed == 41 else 1, 0.1)
                original.update({
                    "exit_code": 0,
                    "process_state": "completed",
                    "status": "success",
                })
                original_path = seed_root / "receipt.json"
                original_path.write_text(json.dumps(original), encoding="utf-8")
                row = success(seed, 39 if seed == 41 else 3, 12 if seed == 41 else 1, 0.1)
                row.pop("actual_execution_success")
                row.update({
                    "execution_success": True,
                    "original_receipt_ref": f"supporting/seeds/seed-{seed}/receipt.json",
                    "original_receipt_sha256": distribution._sha256(original_path),
                    "reinspection_ref": None,
                })
                rows.append(row)
                baseline[str(seed)] = {"a.npz": "m"}
            base_path = source_root / "crbn_calibration.json"
            base_path.write_text(
                json.dumps({
                    "seed_receipts": rows,
                    "new_MSA_baseline": {
                        "processed_msa_hashes_by_seed": baseline
                    },
                }),
                encoding="utf-8",
            )
            with patch.object(
                distribution, "KNOWN_BASE_SHA256", distribution._sha256(base_path)
            ):
                result = distribution.build_benchmark_distribution(base_path, [])
            self.assertEqual(result["technical_success_count"], 2)
            self.assertEqual(
                {row["execution_status"] for row in result["seeds"]},
                {"frozen_base_hash_verified_original_inspection"},
            )

    def test_optional_real_five_seed_distribution_smoke(self):
        root = Path(__file__).resolve().parents[1]
        base_path = root / "cases" / "acceptance_sources" / "crbn_calibration.json"
        receipt_paths = [
            root
            / ".localdata"
            / "expert-closure-20260930"
            / "crbn-extra"
            / f"seed-{seed}"
            / "receipt.json"
            for seed in (79, 97)
        ]
        if not base_path.is_file() or not all(path.is_file() for path in receipt_paths):
            return
        receipts = [
            json.loads(path.read_text(encoding="utf-8")) for path in receipt_paths
        ]
        result = distribution.build_benchmark_distribution(base_path, receipts)
        self.assertEqual(
            [row["seed"] for row in result["seeds"]],
            [23, 41, 61, 79, 97],
        )

    def test_standalone_execution_success_boolean_is_not_trusted(self):
        item = success(1, 2, 3, 0.4)
        item.pop("actual_execution_success")
        item["execution_success"] = True
        with self.assertRaisesRegex(
            BenchmarkDistributionError,
            "FROZEN_BASE_EXECUTION_PROVENANCE_INVALID",
        ):
            distribution.build_benchmark_distribution([], [item])

    def test_frozen_predictions_and_lightning_logs_are_never_copied(self):
        with tempfile.TemporaryDirectory() as temporary:
            with installed_fixture(Path(temporary)) as (_source, destination, result):
                self.assertEqual(list(destination.rglob("*.cif")), [])
                self.assertFalse(
                    any(
                        "lightning_logs" in name
                        for name in result["copied_files_sha256"]
                    )
                )
                self.assertTrue(
                    (
                        destination
                        / "boltz_results_6BOY_CRBN/processed/msa/a.npz"
                    ).is_file()
                )

    def test_frozen_receipt_traversal_is_rejected_before_prefix_filter(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(
                boltz_worker.BoltzWorkerError, "PATH_INVALID|TRAVERSAL"
            ):
                with installed_fixture(
                    Path(temporary),
                    "boltz_results_6BOY_CRBN/../predictions/bad.cif",
                ):
                    pass

    def test_frozen_source_symlink_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with frozen_fixture(root) as fixture:
                source, input_yaml, reference, checkpoint, executable = fixture
                link = source / "boltz_output/boltz_results_6BOY_CRBN/msa/link.csv"
                try:
                    link.symlink_to(source / "6BOY_CRBN.yaml")
                except OSError:
                    self.skipTest("symlinks unavailable")
                with self.assertRaisesRegex(boltz_worker.BoltzWorkerError, "SYMLINK"):
                    boltz_worker._install_frozen_seed_cache(
                        source_directory=source,
                        destination=root / "dest",
                        input_yaml=input_yaml,
                        reference=reference,
                        checkpoint=checkpoint,
                        executable=executable,
                        environment={"boltz_version": "2.2.1"},
                        settings={"x": 1},
                        help_sha256="help",
                    )

    def test_frozen_unexpected_source_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with frozen_fixture(root) as fixture:
                source, input_yaml, reference, checkpoint, executable = fixture
                unexpected = (
                    source
                    / "boltz_output/boltz_results_6BOY_CRBN/msa/hidden.dat"
                )
                unexpected.write_bytes(b"unexpected")
                with self.assertRaisesRegex(
                    boltz_worker.BoltzWorkerError, "FILE_SET_MISMATCH"
                ):
                    boltz_worker._install_frozen_seed_cache(
                        source_directory=source,
                        destination=root / "dest",
                        input_yaml=input_yaml,
                        reference=reference,
                        checkpoint=checkpoint,
                        executable=executable,
                        environment={"boltz_version": "2.2.1"},
                        settings={"x": 1},
                        help_sha256="help",
                    )

    def test_frozen_original_yaml_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with frozen_fixture(root) as fixture:
                source, input_yaml, reference, checkpoint, executable = fixture
                (source / "6BOY_CRBN.yaml").write_bytes(b"tampered")
                with self.assertRaisesRegex(
                    boltz_worker.BoltzWorkerError, "ORIGINAL_INPUT_MISMATCH"
                ):
                    boltz_worker._install_frozen_seed_cache(
                        source_directory=source,
                        destination=root / "dest",
                        input_yaml=input_yaml,
                        reference=reference,
                        checkpoint=checkpoint,
                        executable=executable,
                        environment={"boltz_version": "2.2.1"},
                        settings={"x": 1},
                        help_sha256="help",
                    )

    def test_frozen_source_cache_remains_unchanged(self):
        with tempfile.TemporaryDirectory() as temporary:
            with installed_fixture(Path(temporary)) as (source, destination, result):
                before = result["source_output_files_sha256"].copy()
                boltz_worker._verify_frozen_post_run(result, destination)
                after = {
                    name: boltz_worker.sha256_file(path)
                    for name, path in boltz_worker._tree_files(
                        source / "boltz_output"
                    ).items()
                }
                self.assertEqual(before, after)
                self.assertIs(result["source_cache_unchanged"], True)

    def test_manifest_rewrite_is_reported_but_msa_rewrite_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            with installed_fixture(Path(temporary)) as (_source, destination, result):
                manifest = (
                    destination
                    / "boltz_results_6BOY_CRBN/processed/manifest.json"
                )
                manifest.write_bytes(b"official rewrite")
                boltz_worker._verify_frozen_post_run(result, destination)
                self.assertTrue(
                    result["post_run_manifest_changes"][0]["change_source"].startswith(
                        "official Boltz"
                    )
                )
                msa = (
                    destination
                    / "boltz_results_6BOY_CRBN/processed/msa/a.npz"
                )
                msa.write_bytes(b"changed")
                with self.assertRaisesRegex(
                    boltz_worker.BoltzWorkerError, "CHANGED_AFTER_RUN"
                ):
                    boltz_worker._verify_frozen_post_run(result, destination)

    def test_priority_selection_requires_exact_parent(self):
        document = candidate_document()
        document["parent_scope"]["actual_design_parent_id"] = "OTHER"
        with self.assertRaisesRegex(ValueError, "PARENT_SCOPE"):
            batch.select_priority_candidates(document)

    def test_priority_selection_requires_exact_6haz_frame(self):
        document = candidate_document()
        document["parent_scope"]["receptor_frame"] = "6HAZ chain B"
        with self.assertRaisesRegex(ValueError, "STRUCTURAL_FRAME"):
            batch.select_priority_candidates(document)

    def test_priority_selection_rejects_missing_e3(self):
        document = candidate_document()
        document["protac_candidates"] = [
            row
            for row in document["protac_candidates"]
            if not (
                row["warhead_analog_id"] == batch.ANALOG_IDS[0]
                and row["e3_type"] == "CRBN"
            )
        ]

        def find_candidate(doc, candidate_id):
            return next(
                row
                for row in doc["protac_candidates"]
                if row["candidate_id"] == candidate_id
            )

        def mapped_graph(row):
            return {"graph_sha256": digest(row["candidate_id"].encode("utf-8"))}

        with patch.object(
            batch, "_candidate", side_effect=find_candidate
        ), patch.object(
            batch.novel_ternary, "_mapped_graph", side_effect=mapped_graph
        ):
            with self.assertRaisesRegex(ValueError, "NOT_UNIQUE"):
                batch.select_priority_candidates(document)

    def test_priority_selection_rejects_duplicate_graphs(self):
        document = candidate_document()

        def find_candidate(doc, candidate_id):
            return next(
                row
                for row in doc["protac_candidates"]
                if row["candidate_id"] == candidate_id
            )

        with patch.object(batch, "_candidate", side_effect=find_candidate), patch.object(
            batch.novel_ternary,
            "_mapped_graph",
            return_value={"graph_sha256": "duplicate"},
        ):
            with self.assertRaisesRegex(ValueError, "GRAPH_DUPLICATE"):
                batch.select_priority_candidates(document)

    def test_six_representatives_use_the_same_five_seeds(self):
        document = candidate_document()

        def find_candidate(doc, candidate_id):
            return next(
                row
                for row in doc["protac_candidates"]
                if row["candidate_id"] == candidate_id
            )

        def mapped_graph(row):
            return {"graph_sha256": digest(row["candidate_id"].encode("utf-8"))}

        with patch.object(batch, "_candidate", side_effect=find_candidate), patch.object(
            batch.novel_ternary, "_mapped_graph", side_effect=mapped_graph
        ):
            selected = batch.select_priority_candidates(document)

        seeds = batch._seeds("23,41,61,79,97")
        self.assertEqual(len(selected), 6)
        self.assertTrue(
            all(seeds == [23, 41, 61, 79, 97] for _candidate in selected)
        )

    def test_policy_binding_converts_integer_zero_revision_to_string(self):
        class Service:
            @staticmethod
            def _policy_digest(policy):
                return "policy-digest"

        binding = batch._policy_binding(Service(), {"revision": 0})
        self.assertEqual(binding["revision"], "0")
        self.assertEqual(binding["digest"], "policy-digest")
        self.assertEqual(binding["module"], batch.novel_ternary.__name__)

    def test_policy_binding_rejects_boolean_revision(self):
        class Service:
            @staticmethod
            def _policy_digest(policy):
                return "policy-digest"

        with self.assertRaisesRegex(ValueError, "POLICY_REVISION_REQUIRED"):
            batch._policy_binding(Service(), {"revision": False})

    def test_resume_plan_comparison_ignores_creation_timestamp_only(self):
        first = {
            "created_utc": "one",
            "plan_digest": "a",
            "seeds": [1, 2, 3, 4, 5],
            "x": 1,
        }
        second = {
            "created_utc": "two",
            "plan_digest": "b",
            "seeds": [1, 2, 3, 4, 5],
            "x": 1,
        }
        self.assertEqual(
            batch._immutable_plan_payload(first),
            batch._immutable_plan_payload(second),
        )
        second["seeds"] = [1, 2, 3, 4, 6]
        self.assertNotEqual(
            batch._immutable_plan_payload(first),
            batch._immutable_plan_payload(second),
        )


if __name__ == "__main__":
    unittest.main()
