import copy
from pathlib import Path
import tempfile
import unittest

from packages.science.boltz_run_record import record
from packages.science.handoff import sha


class RunRecordTests(unittest.TestCase):
    def options(self):
        return dict(run_id="metadata-test", compound_id="C01", requested_samples=1, exit_code=1,
                    settings={"model": "boltz2", "recycling_steps": 3, "sampling_steps": 200, "step_scale": 1.5,
                              "max_parallel_samples": 1, "seed": 23, "use_potentials": False, "no_kernels": True,
                              "msa_mode": "server", "msa_pairing_strategy": "greedy", "msa_server_url": "https://example.invalid"},
                    runtime={"boltz": "2.2.1", "rdkit": "2025.3.6", "torch": "test-only", "cuda": None, "gpu": "test-only"})

    def test_failed_before_weights_or_msa_retains_absence(self):
        result = record(**self.options())
        self.assertIsNone(result["weights_sha256"])
        self.assertIsNone(result["log_sha256"])
        self.assertEqual(result["msa_files_sha256"], {})
        self.assertEqual(result["exit_code"], 1)

    def test_success_requires_actual_weight_file(self):
        options = self.options()
        options["exit_code"] = 0
        with self.assertRaisesRegex(ValueError, "SCHEMA"):
            record(**options)

    def test_hashes_full_weight_log_and_nested_msa_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            weight = root / "weights-test-only"
            weight.write_bytes(b"a" * 1048577)
            log = root / "log-test-only"
            log.write_bytes(b"failure")
            msa = root / "msa" / "nested"
            msa.mkdir(parents=True)
            (msa / "test.a3m").write_bytes(b"test-only")
            options = self.options()
            before = copy.deepcopy(options)
            result = record(**options, checkpoint=weight, log=log, msa_dir=root / "msa")
            self.assertEqual(result["weights_sha256"], sha(weight.read_bytes()))
            self.assertEqual(result["msa_files_sha256"], {"nested/test.a3m": sha(b"test-only")})
            self.assertEqual(result["log_sha256"], sha(b"failure"))
            self.assertEqual(options, before)

    def test_missing_explicit_path_is_error_not_silently_null(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "REGULAR_FILE"):
                record(**self.options(), checkpoint=Path(tmp) / "absent")
            with self.assertRaisesRegex(ValueError, "MSA_DIRECTORY"):
                record(**self.options(), msa_dir=Path(tmp) / "absent")
