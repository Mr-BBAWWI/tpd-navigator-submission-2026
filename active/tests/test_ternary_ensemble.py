from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

import scripts.inspect_ternary_ensemble as ensemble
from scripts.inspect_ternary_ensemble import (
    InspectionError,
    _verified_tree,
    complete_link_clusters,
    jaccard,
    protein_contact_signature,
    safe_output_key,
)


def _protein_geometry(transform=None):
    points = {
        ("target", 10, "CA"): np.array([0.0, 0.0, 0.0]),
        ("target", 11, "CA"): np.array([20.0, 0.0, 0.0]),
        ("e3", 30, "CA"): np.array([5.0, 0.0, 0.0]),
        ("e3", 31, "CA"): np.array([30.0, 0.0, 0.0]),
    }
    if transform is not None:
        rotation, translation = transform
        points = {
            key: point @ rotation + translation
            for key, point in points.items()
        }
    return {
        key: {"xyz": point.tolist(), "element": "C"}
        for key, point in points.items()
    }


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _candidate(candidate_id: str) -> dict:
    seeds = []
    for seed in ensemble.SEEDS:
        seeds.append(
            {
                "seed": seed,
                "protein_protein": {
                    "count": 1,
                    "keys": ["10|30"],
                    "full_contact_array": [],
                },
                "target_warhead_site": {"count": 0, "keys": [], "full_contact_array": []},
                "e3_recruiter_site": {"count": 0, "keys": [], "full_contact_array": []},
                "interface_clashes": {
                    "vdw_pair_count": 1,
                    "raw_under_2A_pair_count": 0,
                    "total_protein_heavy_atom_count": 10,
                    "vdw_pairs_per_total_heavy_atom": 0.1,
                    "closest_20_atom_pairs": [],
                },
            }
        )
    return {
        "candidate_id": candidate_id,
        "e3_type": "test-e3",
        "seeds": seeds,
        "pairwise_contact_jaccard": [],
        "contact_frequencies": {
            "protein_protein": {
                "all_contacts": [{"key": "10|30", "seed_count": 5, "frequency": 1.0}],
                "shared_at_least_4_of_5": [{"key": "10|30", "seed_count": 5}],
            }
        },
        "complete_link_contact_jaccard_sensitivity": [
            {"threshold": threshold, "memberships": [list(ensemble.SEEDS)]}
            for threshold in ensemble.THRESHOLDS
        ],
    }


class TernaryEnsembleTests(unittest.TestCase):
    def test_contact_boundary_and_rigid_transform_invariance(self):
        initial = protein_contact_signature(_protein_geometry())
        rotation = np.array(
            [
                [0.0, -1.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        moved = protein_contact_signature(
            _protein_geometry((rotation, np.array([7.0, 4.0, -2.0])))
        )
        self.assertEqual(initial["keys"], ["10|30"])
        self.assertEqual(moved["keys"], initial["keys"])
        self.assertAlmostEqual(
            initial["full_contact_array"][0]["minimum_distance_A"],
            5.0,
        )
        self.assertAlmostEqual(
            moved["full_contact_array"][0]["minimum_distance_A"],
            5.0,
        )

    def test_empty_union_jaccard_is_none(self):
        self.assertIsNone(jaccard([], []))
        self.assertEqual(jaccard([], ["x"]), 0.0)

    def test_complete_link_does_not_chain_merge(self):
        signatures = {
            1: {"a", "b"},
            2: {"a", "b", "c"},
            3: {"b", "c"},
        }
        clusters = complete_link_clusters(signatures, 0.5)
        self.assertEqual(clusters, [[1, 2], [3]])
        self.assertTrue(
            all(not ({1, 2, 3} <= set(cluster)) for cluster in clusters)
        )

    def test_output_key_path_traversal_and_absolute_rejected(self):
        invalid = [
            "../escape.cif",
            "a/../../escape.cif",
            "/absolute.cif",
            "C:\\absolute.cif",
            "a\\escape.cif",
        ]
        for key in invalid:
            with self.subTest(key=key):
                with self.assertRaises(InspectionError):
                    safe_output_key(key)

    def test_verified_tree_accepts_exact_tree_and_detects_modified_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            output.mkdir()
            prediction = output / "prediction.cif"
            prediction.write_text("original\n", encoding="utf-8")
            recorded = {"prediction.cif": _sha(prediction)}

            hashes, selected = _verified_tree(output, recorded, root)
            self.assertEqual(hashes, recorded)
            self.assertEqual(selected, prediction)

            prediction.write_text("modified\n", encoding="utf-8")
            with self.assertRaisesRegex(InspectionError, "OUTPUT_HASH_MISMATCH"):
                _verified_tree(output, recorded, root)

    def test_verified_tree_detects_missing_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            output.mkdir()
            recorded = {
                "prediction.cif": hashlib.sha256(b"missing").hexdigest(),
            }
            with self.assertRaisesRegex(InspectionError, "OUTPUT_FILE_SET_MISMATCH"):
                _verified_tree(output, recorded, root)

    def test_verified_tree_detects_invalid_recorded_path(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            output.mkdir()
            prediction = output / "prediction.cif"
            prediction.write_text("data\n", encoding="utf-8")
            recorded = {"../prediction.cif": _sha(prediction)}
            with self.assertRaises(InspectionError):
                _verified_tree(output, recorded, root)

    def test_verified_tree_rejects_output_outside_root(self):
        with tempfile.TemporaryDirectory() as root_temporary:
            with tempfile.TemporaryDirectory() as output_temporary:
                root = Path(root_temporary)
                output = Path(output_temporary)
                prediction = output / "prediction.cif"
                prediction.write_text("data\n", encoding="utf-8")
                with self.assertRaisesRegex(InspectionError, "PATH_TRAVERSAL_REJECTED"):
                    _verified_tree(output, {"prediction.cif": _sha(prediction)}, root)

    def test_verify_plan_result_is_not_assigned(self):
        plans = []
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for index in range(6):
                candidate_id = f"test-candidate-{index}"
                plan = {
                    "seeds": list(ensemble.SEEDS),
                    "candidate_graph": {
                        "candidate_id": candidate_id,
                        "e3_type": "test-e3",
                    },
                    "plan_digest": f"test-digest-{index}",
                }
                plans.append(plan)
                plan_path = root / f"{candidate_id}.plan.json"
                plan_path.write_text(json.dumps(plan), encoding="utf-8")
                run_dir = root / candidate_id
                run_dir.mkdir()
                for seed in ensemble.SEEDS:
                    (run_dir / f"seed-{seed}").mkdir()

            with mock.patch.object(ensemble.novel, "verify_plan", return_value=None) as verify:
                with self.assertRaises(InspectionError):
                    ensemble._preflight(root)
            verify.assert_called_once()
            self.assertEqual(
                verify.call_args.args[0].get("candidate_graph"),
                plans[0]["candidate_graph"],
            )

    def test_summary_has_explicit_review_status_and_preserves_all_30_results(self):
        candidates = [_candidate(f"D-test-{index}") for index in range(6)]
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "report"
            with mock.patch.object(ensemble, "_preflight", return_value=[{}] * 6), mock.patch.object(
                ensemble,
                "_analyze_candidate",
                side_effect=candidates,
            ):
                summary = ensemble.inspect(Path(temporary), output)

            self.assertIs(summary["scientific_approved"], False)
            self.assertIs(summary["human_review_performed"], False)
            self.assertEqual(summary["raw_seed_result_count"], 30)
            self.assertEqual(len(summary["candidates"]), 6)
            self.assertTrue(
                all(
                    [seed["seed"] for seed in candidate["seeds"]] == ensemble.SEEDS
                    for candidate in summary["candidates"]
                )
            )
            stored = json.loads((output / "summary.json").read_text(encoding="utf-8"))
            self.assertIs(stored["scientific_approved"], False)
            self.assertIs(stored["human_review_performed"], False)
            report = (output / "report.html").read_text(encoding="utf-8")
            self.assertIn("Clash sensitivity", report)
            self.assertIn("Complete-link clustering sensitivity", report)
            self.assertIn("Total protein heavy atoms", report)


if __name__ == "__main__":
    unittest.main()
