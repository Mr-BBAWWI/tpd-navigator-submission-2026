import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from scripts.build_science_closure_pack import (
    assert_close_metrics,
    claimed_output,
    load_json,
    safe_join,
    selected_index,
    validate_pka,
    verify_file,
)


def write_json(path: Path, value):
    path.write_text(json.dumps(value, separators=(",", ":")) + "\n", encoding="utf-8")


class ScienceClosurePackHelpersTest(unittest.TestCase):
    def test_safe_join_accepts_regular_file_and_rejects_cross_platform_escapes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "inside.txt").write_text("evidence", encoding="utf-8")
            self.assertEqual(safe_join(root, "inside.txt"), root.absolute() / "inside.txt")
            for unsafe in ("../outside.txt", "/absolute.txt", r"C:\foreign\file.txt",
                           r"folder\..\file.txt", "folder/name:stream"):
                with self.subTest(unsafe=unsafe), self.assertRaises(ValueError):
                    safe_join(root, unsafe, must_exist=False)

    def test_safe_join_rejects_symlink_when_supported(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "target.txt"
            target.write_text("target", encoding="utf-8")
            link = root / "link.txt"
            try:
                link.symlink_to(target)
            except (OSError, NotImplementedError):
                self.skipTest("symlinks are unavailable")
            with self.assertRaises(ValueError):
                safe_join(root, "link.txt")

    def test_claimed_output_uses_the_actual_path_and_rejects_foreign_prefix(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "seed-101"
            prediction = root / "boltz_output" / "run" / "predictions" / "model.cif"
            prediction.parent.mkdir(parents=True)
            prediction.write_text("model", encoding="utf-8")
            source, relative = claimed_output(root, str(prediction.absolute()))
            self.assertEqual(source, prediction.absolute())
            self.assertEqual(relative, "run/predictions/model.cif")
            foreign = Path(temporary) / "foreign" / "boltz_output" / "run" / "predictions" / "model.cif"
            foreign.parent.mkdir(parents=True)
            foreign.write_text("foreign", encoding="utf-8")
            with self.assertRaises(ValueError):
                claimed_output(root, str(foreign.absolute()))

    def test_load_json_rejects_duplicates_and_all_nonfinite_numbers(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "value.json"
            for document in ('{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}', '{"a":1e999}'):
                path.write_text(document, encoding="utf-8")
                with self.subTest(document=document), self.assertRaises(ValueError):
                    load_json(path)

    def test_metrics_and_indices_are_strict_finite_types(self):
        assert_close_metrics({"x": 1.0}, {"x": 1.00000001}, ["x"])
        for value in (True, "1", float("nan"), float("inf")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                assert_close_metrics({"x": value}, {"x": 1.0}, ["x"])
        self.assertEqual(selected_index(2), 2)
        for value in (True, (2,), {"selected_model_index": 2}):
            with self.subTest(value=value), self.assertRaises(ValueError):
                selected_index(value)

    def test_validate_pka_uses_manifest_and_actual_raw_prediction_shape(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            compact = {
                "scientific_approved": False,
                "population_prediction_performed": False,
                "analog_results": {"A": {}, "B": {}, "W-80f8f4a11b5d": {}},
                "missing_required_maps": {"A": [], "B": [], "W-80f8f4a11b5d": [5001]},
            }
            raw = {
                "predictions": {
                    "A": {"sites": list(range(8))},
                    "B": {"sites": list(range(9))},
                    "W-80f8f4a11b5d": {"sites": list(range(9))},
                }
            }
            write_json(root / "compact-summary.json", compact)
            write_json(root / "raw-predictions.json", raw)
            files = []
            for name in ("compact-summary.json", "raw-predictions.json"):
                path = root / name
                files.append({"path": name, "bytes": path.stat().st_size,
                              "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
            manifest = {
                "scientific_approved": False,
                "population_prediction_performed": False,
                "files": files,
            }
            write_json(root / "manifest.json", manifest)
            loaded_manifest, loaded_compact, manifest_path = validate_pka(root)
            self.assertEqual(loaded_manifest, manifest)
            self.assertEqual(loaded_compact, compact)
            self.assertEqual(manifest_path, root / "manifest.json")

            raw["predictions"]["A"]["sites"].pop()
            write_json(root / "raw-predictions.json", raw)
            manifest["files"][1] = {
                "path": "raw-predictions.json",
                "bytes": (root / "raw-predictions.json").stat().st_size,
                "sha256": hashlib.sha256((root / "raw-predictions.json").read_bytes()).hexdigest(),
            }
            write_json(root / "manifest.json", manifest)
            with self.assertRaisesRegex(ValueError, "expected 26"):
                validate_pka(root)

    def test_verify_file_checks_hash_and_size(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "artifact.bin"
            payload = b"verified evidence"
            path.write_bytes(payload)
            digest = hashlib.sha256(payload).hexdigest()
            verify_file(path, digest, len(payload))
            with self.assertRaises(ValueError):
                verify_file(path, "0" * 64, len(payload))
            with self.assertRaises(ValueError):
                verify_file(path, digest, len(payload) + 1)


if __name__ == "__main__":
    unittest.main()
