"""Core-only archive tests using an explicitly synthetic rejected request."""
import tempfile
import unittest
from pathlib import Path

from packages.science.evidence_common import file_index, seal, write_new_directory
from packages.science.handoff import encoded, sha
from packages.science.ligand_preparation import ARCHIVE, REQUEST, _report, read_preparation, validate_request


def archive(root):
    files = {"raw/model.cif": b"synthetic", "raw/ligand.sdf": b"invalid synthetic SDF", "raw/input.yaml": b"synthetic"}
    refs = {k: {"path": p, "sha256": sha(files[p])} for k, p in zip(("source_structure", "source_molecule", "model_input"), files)}
    request = {"format": REQUEST, "data_mode": "synthetic_test", **refs,
               "binding": {"compound_id": "synthetic", "molecule_id": "synthetic:" + sha(b"N"), "run_id": "synthetic",
                           "model_rank": 0, "input_sha256": refs["model_input"]["sha256"], "structure_sha256": refs["source_structure"]["sha256"]},
               "expected_isomeric_smiles": "N", "chemistry_policy": "preserve_input_no_ph_prediction",
               "selection": {"auth_chain_id": "X", "auth_seq_id": 1, "component_id": "UNK", "insertion_code": ""},
               "atom_mapping": [{"atom_map": 1, "source_atom_name": "N"}]}
    outcome = {"status": "not_prepared", "error": "LIGAND_SINGLE_VALID_SDF_REQUIRED", "chemistry": None, "runtime": {}}
    report = _report(request, outcome, files)
    files.update({"request.json": encoded(request), "outcome.json": encoded(outcome), "report.json": encoded(report)})
    write_new_directory(root, {**files, "manifest.json": encoded(seal({"format": ARCHIVE, "files": file_index(files), "report_digest": report["digest"]}))})
    return request, report, files


def reseal(root, report, files):
    for name, data in files.items():
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).write_bytes(data)
    (root / "manifest.json").write_bytes(encoded(seal({"format": ARCHIVE, "files": file_index(files), "report_digest": report["digest"]})))


class LigandArchiveTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / "archive"
        self.request, self.report, self.files = archive(self.root)

    def test_rejected_request_roundtrips_without_scientific_imports(self):
        import sys
        before = set(sys.modules)
        self.assertEqual(read_preparation(self.root, allow_synthetic=True), self.report)
        self.assertFalse({"rdkit", "gemmi"} & (set(sys.modules) - before))

    def test_synthetic_archive_not_accepted_as_real(self):
        with self.assertRaisesRegex(ValueError, "SYNTHETIC_LIGAND"):
            read_preparation(self.root)

    def test_resealed_report_cannot_override_not_assessed(self):
        self.report["contact_quality_status"] = "passed"
        self.files["report.json"] = encoded(seal({k: v for k, v in self.report.items() if k != "digest"}))
        reseal(self.root, self.report, self.files)
        with self.assertRaisesRegex(ValueError, "LIGAND_REPORT_PROJECTION"):
            read_preparation(self.root, allow_synthetic=True)

    def test_resealed_outcome_cannot_claim_success_without_structure(self):
        self.files["outcome.json"] = encoded({"status": "prepared_for_review", "error": None, "chemistry": {}, "runtime": {}})
        reseal(self.root, self.report, self.files)
        with self.assertRaisesRegex(ValueError, "LIGAND_OUTPUT_STATUS"):
            read_preparation(self.root, allow_synthetic=True)

    def test_extra_unindexed_file_rejected(self):
        (self.root / "extra.json").write_bytes(b"{}")
        with self.assertRaisesRegex(ValueError, "LIGAND_ARCHIVE_UNINDEXED"):
            read_preparation(self.root, allow_synthetic=True)

    def test_extra_indexed_file_rejected(self):
        self.files["extra.json"] = b"{}"; reseal(self.root, self.report, self.files)
        with self.assertRaisesRegex(ValueError, "LIGAND_ARCHIVE_UNEXPECTED"):
            read_preparation(self.root, allow_synthetic=True)

    def test_binding_hash_mismatch_rejected(self):
        self.request["binding"]["structure_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "LIGAND_SOURCE_BINDING"):
            validate_request(self.request, self.files, allow_synthetic=True)

    def test_new_chemical_policy_requires_new_contract(self):
        self.request["chemistry_policy"] = "choose_protonation_at_ph7"
        with self.assertRaisesRegex(ValueError, "SCHEMA_INVALID"):
            validate_request(self.request, self.files, allow_synthetic=True)
