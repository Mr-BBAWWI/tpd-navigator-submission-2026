from __future__ import annotations

import hashlib
import importlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


REVIEWER_DELIVERY = Path(__file__).resolve().parents[1] / "reviewer_delivery"
if str(REVIEWER_DELIVERY) not in sys.path:
    sys.path.insert(0, str(REVIEWER_DELIVERY))
package = importlib.import_module("build_acceptance_package")


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def fingerprint(files: dict[str, str]) -> dict:
    encoded = json.dumps(files, sort_keys=True, separators=(",", ":")).encode()
    return {
        "algorithm": "sha256-file-map-v1",
        "sha256": hashlib.sha256(encoded).hexdigest(),
        "files": files,
    }


class ReleaseIntegrityTests(unittest.TestCase):
    def make_stage(self, root: Path) -> tuple[Path, dict]:
        stage = root / "stage"
        source = stage / "program" / "apps" / "main.py"
        source.parent.mkdir(parents=True)
        source.write_bytes(b"print('staged')\n")
        files = {"apps/main.py": digest(source.read_bytes())}
        source_fingerprint = fingerprint(files)
        manifest = {
            "format": "tpd-acceptance-delivery/2",
            "source_fingerprint": source_fingerprint,
            "files": {"program/apps/main.py": files["apps/main.py"]},
        }
        (stage / "PACKAGE-MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")
        return stage, manifest

    def test_manifest_rejects_unexpected_python_shadow_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            stage, _ = self.make_stage(Path(temporary))
            (stage / "program" / "apps.py").write_text("raise SystemExit\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "PACKAGE_UNEXPECTED_FILE"):
                package.verify_manifest(stage)

    def test_manifest_rejects_explicit_parent_path(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            stage, manifest = self.make_stage(Path(temporary))
            manifest["files"]["../outside"] = "0" * 64
            (stage / "PACKAGE-MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "PACKAGE_MANIFEST_PATH_INVALID"):
                package.verify_manifest(stage)

    def test_staged_source_bytes_must_match_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            stage, manifest = self.make_stage(Path(temporary))
            manifest["files"]["program/apps/main.py"] = digest(b"different\n")
            with self.assertRaisesRegex(RuntimeError, "PACKAGE_STAGED_SOURCE_MISMATCH"):
                package.verify_staged_source_binding(stage, manifest)

    def test_current_source_fingerprint_must_match_stage(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            stage, manifest = self.make_stage(Path(temporary))
            current = fingerprint({"apps/main.py": digest(b"new source\n")})
            with self.assertRaisesRegex(RuntimeError, "PACKAGE_SOURCE_MISMATCH"):
                package.verify_staged_source_binding(stage, manifest, current)

    def test_manifest_rejects_symlink_in_immutable_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage, _ = self.make_stage(root)
            target = root / "outside.py"
            target.write_text("pass\n", encoding="utf-8")
            link = stage / "program" / "shadow.py"
            try:
                link.symlink_to(target)
            except (OSError, NotImplementedError):
                self.skipTest("symlinks unavailable")
            with self.assertRaisesRegex(RuntimeError, "PACKAGE_LINK_NOT_ALLOWED"):
                package.verify_manifest(stage)

    def test_linker_library_is_part_of_source_fingerprint_when_present(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            active = Path(temporary)
            linker = active / "cases" / "linker_library_20260930.json"
            linker.parent.mkdir(parents=True)
            linker.write_text("{}", encoding="utf-8")
            with mock.patch.object(package, "ACTIVE", active):
                result = package.source_fingerprint()
            self.assertEqual(result["files"]["cases/linker_library_20260930.json"], digest(b"{}"))


if __name__ == "__main__":
    unittest.main()
