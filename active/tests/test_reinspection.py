import hashlib
import os
import tempfile
import unittest
from pathlib import Path

from scripts.reinspect_crbn_calibration import (
    _execution_metadata,
    _verify_processed_outputs,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ReinspectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.output = self.root / "boltz_output"
        self.output.mkdir()
        prediction = self.output / "results" / "prediction.cif"
        confidence = self.output / "results" / "confidence_prediction.json"
        prediction.parent.mkdir()
        prediction.write_text("data_prediction\n", encoding="utf-8")
        confidence.write_text('{"confidence": 0.5}\n', encoding="utf-8")
        self.prediction = prediction
        self.confidence = confidence
        self.receipt = {
            "seed": 0,
            "seed_application": {
                "seed": 0,
                "method": "official_cli_--seed",
                "actually_applied": True,
            },
            "settings": {
                "accelerator": "gpu",
                "recycling_steps": 3,
                "sampling_steps": 200,
                "diffusion_samples": 1,
                "max_parallel_samples": 1,
                "num_workers": 0,
                "preprocessing_threads": 1,
                "no_kernels": True,
                "use_potentials": False,
                "msa_mode": "server",
                "max_msa_seqs": 256,
            },
            "tool_environment": {
                "boltz_version": "2.2.1",
                "torch_version": "2.7.1+cu128",
                "cuda_available": True,
                "cuda_runtime": "12.8",
                "cuda_device_count": 1,
                "cuda_devices": ["Fixture GPU"],
            },
            "hashes": {
                "processed_output_files_sha256": {
                    "results/prediction.cif": _sha256(prediction),
                    "results/confidence_prediction.json": _sha256(confidence),
                }
            },
        }

    def tearDown(self):
        self.temp.cleanup()

    def test_actual_receipt_schema_and_applied_seed_zero(self):
        verified = _verify_processed_outputs(self.receipt, self.output)
        self.assertEqual(
            verified,
            self.receipt["hashes"]["processed_output_files_sha256"],
        )

        metadata = _execution_metadata(self.receipt)
        self.assertEqual(metadata["applied_seed"], 0)
        self.assertEqual(metadata["applied_settings"], self.receipt["settings"])
        self.assertEqual(metadata["model_version"], "2.2.1")
        self.assertEqual(metadata["cuda"]["runtime"], "12.8")
        self.assertEqual(metadata["cuda"]["devices"], ["Fixture GPU"])

    def test_output_tamper_is_rejected(self):
        self.prediction.write_text("data_tampered\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "OUTPUT_HASH_MISMATCH"):
            _verify_processed_outputs(self.receipt, self.output)

    def test_unrecorded_output_file_is_rejected(self):
        (self.output / "unexpected.bin").write_bytes(b"unexpected")
        with self.assertRaisesRegex(ValueError, "OUTPUT_FILE_SET_MISMATCH"):
            _verify_processed_outputs(self.receipt, self.output)

    def test_requested_seed_is_not_treated_as_applied(self):
        self.receipt["seed_application"]["actually_applied"] = False
        with self.assertRaisesRegex(ValueError, "CONFIRMED_SEED_APPLICATION_MISSING"):
            _execution_metadata(self.receipt)

    def test_seed_application_must_match_receipt_seed(self):
        self.receipt["seed_application"]["seed"] = 1
        with self.assertRaisesRegex(ValueError, "CONFIRMED_SEED_APPLICATION_MISMATCH"):
            _execution_metadata(self.receipt)

    def test_output_hash_path_traversal_is_rejected(self):
        digest = self.receipt["hashes"]["processed_output_files_sha256"].pop(
            "results/prediction.cif"
        )
        self.receipt["hashes"]["processed_output_files_sha256"]["../prediction.cif"] = digest
        with self.assertRaisesRegex(ValueError, "OUTPUT_PATH_TRAVERSAL"):
            _verify_processed_outputs(self.receipt, self.output)

    def test_symlink_output_root_is_rejected_before_resolution(self):
        link = self.root / "linked-output"
        try:
            os.symlink(self.output, link, target_is_directory=True)
        except (OSError, NotImplementedError) as error:
            self.skipTest("directory symlinks unavailable: " + str(error))
        with self.assertRaisesRegex(ValueError, "OUTPUT_DIR_SYMLINK_FORBIDDEN"):
            _verify_processed_outputs(self.receipt, link)


if __name__ == "__main__":
    unittest.main()

