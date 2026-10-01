"""Focused integration tests for computed ligand-pKa pack import."""
from __future__ import annotations

import copy
import hashlib
import io
import json
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

from packages.contracts import ContractError, encoded
from packages.platform.computed_ligand_pka import (
    ALLOWED_PATHS, import_computed_ligand_pka,
)
from scripts import import_computed_ligand_pka as cli

try:
    import test_microstate_population_service as fixtures
except ImportError:
    from tests import test_microstate_population_service as fixtures


class ComputedLigandPkaTests(unittest.TestCase):
    project = fixtures.MicrostatePopulationServiceTests.project
    setUp = fixtures.MicrostatePopulationServiceTests.setUp
    _completed_job = fixtures.MicrostatePopulationServiceTests._completed_job
    _replace_job = fixtures.MicrostatePopulationServiceTests._replace_job
    snapshot = fixtures.MicrostatePopulationServiceTests.snapshot
    report = fixtures.MicrostatePopulationServiceTests.report
    register_interaction = fixtures.MicrostatePopulationServiceTests.register_interaction
    source_document = fixtures.MicrostatePopulationServiceTests.source_document

    def prepare(self, mutate_document=None, extra=None):
        result, _, policy = self.snapshot()
        interaction = self.register_interaction(self.report(), policy)
        document = self.source_document(result, interaction)
        record = document["records"][0]
        record.pop("state_population")
        record["pKa"] = 7.2
        record["source_kind"] = "computed"
        document.update({"scientific_approved": False,
                         "state_selection_performed": False})
        if mutate_document:
            mutate_document(document)
        self.pack_temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.pack_temp.cleanup)
        root = Path(self.pack_temp.name)
        files = {name: (b"synthetic fixture: " + name.encode())
                 for name in ALLOWED_PATHS}
        files["quantitative-source.json"] = encoded(document)
        manifest = {
            "format": "ligand-pka-evidence/20261001.2",
            "files": [{"path": name, "sha256": hashlib.sha256(files[name]).hexdigest(),
                       "bytes": len(files[name])} for name in sorted(files)],
            "final_state_selection_performed": False,
            "population_prediction_performed": False,
            "scientific_approved": False,
        }
        manifest_raw = encoded(manifest)
        for name, raw in {**files, "manifest.json": manifest_raw}.items():
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)
        if extra:
            (root / extra).write_bytes(b"extra")
        return root, hashlib.sha256(manifest_raw).hexdigest(), files, document

    def call(self, root, digest, **changes):
        return import_computed_ligand_pka(
            self.store, changes.get("project", self.project),
            changes.get("job_id", self.job_id), root,
            expected_manifest_sha256=digest)

    def counts(self):
        with self.store.db() as db:
            evidence = db.execute("SELECT COUNT(*) FROM scientific_evidence").fetchone()[0]
            artifacts = db.execute("SELECT COUNT(*) FROM artifacts").fetchone()[0]
            human = tuple(db.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                          for table in ("scientific_assessments", "scientific_decisions",
                                        "scientific_policy_events"))
        blobs = sorted(path.name for path in (self.store.root / "blobs").iterdir())
        return evidence, artifacts, human, blobs

    def test_success_archive_and_repeat_are_stable(self):
        root, digest, files, _ = self.prepare()
        before_human = self.counts()[2]
        made = self.call(root, digest)
        self.assertTrue(made["created"])
        self.assertIsNone(made["binding"]["actor_id"])
        proof = self.service.port.read(made["binding"]["proof_archive_ref"])
        with zipfile.ZipFile(io.BytesIO(proof)) as archive:
            self.assertEqual(archive.namelist(), sorted([*files, "manifest.json"]))
            for name in archive.namelist():
                self.assertEqual(archive.read(name), (root / name).read_bytes())
        snapshot = self.counts()
        repeated = self.call(root, digest)
        self.assertFalse(repeated["created"])
        self.assertEqual(made["id"], repeated["id"])
        self.assertEqual(snapshot, self.counts())
        self.assertEqual(before_human, self.counts()[2])

    def test_stale_runtime_is_preserved_but_stale_inputs_fail(self):
        root, digest, _, _ = self.prepare()
        self._replace_job(lambda body: body.__setitem__("runtime", {"old": True}))
        made = self.call(root, digest)
        self.assertFalse(made["binding"]["source_runtime_current"])

        other = ComputedLigandPkaTests(methodName="runTest")
        other.setUp()
        self.addCleanup(other.doCleanups)
        root2, digest2, _, _ = other.prepare()
        before = other.counts()
        with mock.patch("packages.platform.scientific_acceptance.input_binding",
                        return_value={"digest": "otherdigest"}):
            with self.assertRaisesRegex(
                    ContractError, "SCIENTIFIC_DESIGN_SOURCE_CHANGED"):
                other.call(root2, digest2)
        self.assertEqual(before, other.counts())

    def test_scope_graph_ph_and_nonfinite_rejected(self):
        mutations = [
            lambda d: d["records"][0].__setitem__("job_id", "wrong"),
            lambda d: d["records"][0].__setitem__("project_id", "wrong"),
            lambda d: d["records"][0].__setitem__(
                "state_canonical_isomeric_graph", "CCC"),
            lambda d: d["records"][0].__setitem__("pH", 99),
        ]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                root, digest, _, _ = self.prepare(mutation)
                with self.assertRaises(ContractError):
                    self.call(root, digest)
        root, digest, _, _ = self.prepare()
        path = root / "quantitative-source.json"
        raw = path.read_text().replace("7.2", "NaN")
        path.write_text(raw)
        manifest = json.loads((root / "manifest.json").read_text())
        for item in manifest["files"]:
            if item["path"] == "quantitative-source.json":
                data = path.read_bytes(); item["bytes"] = len(data)
                item["sha256"] = hashlib.sha256(data).hexdigest()
        manifest_raw = encoded(manifest); (root / "manifest.json").write_bytes(manifest_raw)
        with self.assertRaises(ContractError):
            self.call(root, hashlib.sha256(manifest_raw).hexdigest())

    def test_pack_corruption_extra_and_traversal_rejected(self):
        root, digest, _, _ = self.prepare()
        (root / "README.md").write_bytes(b"corrupt")
        with self.assertRaises(ContractError): self.call(root, digest)
        root, digest, _, _ = self.prepare(extra="unexpected.txt")
        with self.assertRaises(ContractError): self.call(root, digest)
        root, digest, _, _ = self.prepare()
        manifest = json.loads((root / "manifest.json").read_text())
        manifest["files"][0]["path"] = "../escape"
        raw = encoded(manifest); (root / "manifest.json").write_bytes(raw)
        with self.assertRaises(ContractError):
            self.call(root, hashlib.sha256(raw).hexdigest())

    def test_duplicate_manifest_keys_and_junction_are_rejected(self):
        root, _, _, _ = self.prepare()
        path = root / "manifest.json"
        text = path.read_text(encoding="utf-8")
        text = text.replace(
            '"format":',
            '"format":"ligand-pka-evidence/20261001.2","format":', 1)
        raw = text.encode("utf-8")
        path.write_bytes(raw)
        with self.assertRaises(ContractError):
            self.call(root, hashlib.sha256(raw).hexdigest())

        root, digest, _, _ = self.prepare()
        junction = root / "README.md"
        with mock.patch.object(Path, "is_junction", create=True,
                               new=lambda path: path == junction):
            with self.assertRaisesRegex(ContractError,
                                        "COMPUTED_LIGAND_PKA_LINK"):
                self.call(root, digest)

    def test_cli_help_is_cwd_independent_and_receipt_preflights(self):
        script = Path(__file__).resolve().parents[1] / "scripts" / \
            "import_computed_ligand_pka.py"
        with tempfile.TemporaryDirectory() as directory:
            completed = subprocess.run(
                [sys.executable, str(script), "--help"], cwd=directory,
                capture_output=True, text=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)

            root = Path(directory)
            receipt = root / "receipt.json"
            receipt.write_text("existing", encoding="utf-8")
            argv = ["--store", str(root / "missing-store"),
                    "--project", "p", "--job-id", "j",
                    "--pack", str(root / "pack"),
                    "--expected-manifest-sha256", "0" * 64,
                    "--receipt", str(receipt)]
            with mock.patch.object(cli, "Store") as store, \
                    mock.patch.object(cli, "import_computed_ligand_pka") as importer, \
                    mock.patch("sys.stderr", new=io.StringIO()):
                self.assertEqual(cli.main(argv), 2)
            store.assert_not_called()
            importer.assert_not_called()

    def test_same_document_changed_pack_and_human_identity_collide(self):
        root, digest, _, document = self.prepare()
        self.call(root, digest)
        (root / "README.md").write_bytes(b"different valid pack")
        manifest = json.loads((root / "manifest.json").read_text())
        for item in manifest["files"]:
            if item["path"] == "README.md":
                raw = (root / "README.md").read_bytes()
                item.update(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
        raw = encoded(manifest); (root / "manifest.json").write_bytes(raw)
        with self.assertRaises(ContractError):
            self.call(root, hashlib.sha256(raw).hexdigest())

        other = ComputedLigandPkaTests(methodName="runTest"); other.setUp()
        self.addCleanup(other.doCleanups)
        root2, digest2, _, doc2 = other.prepare()
        canonical = hashlib.sha256(encoded(doc2)).hexdigest()
        with other.store.db() as db:
            _, result, binding, _ = other.service._verified_result(db, other.job_id)
            made = other.service._register_evidence(
                db, [], other.job_id, "microstate_population_source", canonical,
                {**binding, "actor_id": "human", "source_upload_sha256": canonical,
                 "source_document_sha256": canonical}, doc2, provenance="source")
        with self.assertRaises(ContractError): other.call(root2, digest2)

    def test_failure_after_proof_write_rolls_back_artifact_and_blob(self):
        root, digest, _, _ = self.prepare()
        before = self.counts()
        with mock.patch.object(self.service.__class__, "_register_evidence",
                               side_effect=RuntimeError("injected")):
            with self.assertRaises(RuntimeError): self.call(root, digest)
        self.assertEqual(before, self.counts())


if __name__ == "__main__":
    unittest.main()
