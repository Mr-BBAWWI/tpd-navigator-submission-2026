from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.package_native_parent_review import build_pack, verify_manifest_tree


def _write_manifest(root: Path, files: dict[str, bytes]) -> None:
    root.mkdir(parents=True)
    rows = []
    for name, data in files.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        rows.append({"path": name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
    (root / "manifest.json").write_text(json.dumps({"files": rows}), encoding="utf-8")


def _j(value: object) -> bytes:
    return json.dumps(value).encode()


def _fixtures(root: Path, smiles: str = "CCO"):
    baseline = root / "baseline"
    _write_manifest(baseline, {
        "native-probe-result.json": _j({"runs": []}),
        "protocol.json": _j({}),
    })
    cases = []
    protocol_cases = []
    panel_files: dict[str, bytes] = {}
    for i in range(12):
        cid = f"N3-LH-{i + 1:02d}"
        pass_count = 2 if i < 9 else 0
        cases.append({
            "id": cid,
            "candidate_id": f"candidate-{i}",
            "mapped_smiles": smiles,
            "candidate_sdf": "candidate.sdf",
            "status": "docked_with_limits" if pass_count else "filtered_not_replaced",
            "docking": {"status": "completed_with_limits"},
            "pose_metrics": {
                "pose_preservation_pass_count": pass_count,
                "best_geometry_core_RMSD_A": 1.2,
                "best_geometry_contact_retention": 0.8,
            },
            "assemblies": [],
        })
        protocol_cases.append({"id": cid})
        panel_files[f"cases/{cid}/candidate.sdf"] = b"sdf\n"
        panel_files[f"cases/{cid}/dock-receipt.json"] = b"{}"
    runs = [
        {"ligand": "native_parent_Cl", "seed": seed, "status": "done", "pose_pass_count": passed}
        for seed, passed in ((23, 1), (41, 0), (61, 1))
    ] + [
        {"ligand": "measured_probe_Br", "seed": seed, "status": "done", "pose_pass_count": passed}
        for seed, passed in ((23, 1), (41, 0), (61, 0))
    ]
    panel_files.update({
        "panel.json": _j({
            "status": "completed_with_limits",
            "requested_count": 12,
            "retained_count": 12,
            "cases": cases,
        }),
        "protocol.json": _j({"requested_candidate_count": 12, "cases": protocol_cases}),
        "compact-summary.json": _j({"native_parent_and_Br_runs": runs}),
    })
    panel = root / "panel"
    _write_manifest(panel, panel_files)
    audit = root / "audit"
    _write_manifest(audit, {"evidence-audit.json": b"{}"})
    export = root / "export"
    export.mkdir()
    (export / "receipt.json").write_text("{}", encoding="utf-8")
    cif = root / "9D12.cif"
    csv = root / "source.csv"
    cif.write_bytes(b"cif")
    csv.write_bytes(b"csv")
    return baseline, panel, audit, export, cif, csv


class PackageNativeParentReviewTests(unittest.TestCase):
    def test_manifest_rejects_tamper_traversal_duplicate_and_windows_drive(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tree = root / "tree"
            _write_manifest(tree, {"a.txt": b"one"})
            (tree / "a.txt").write_bytes(b"two")
            with self.assertRaises(ValueError):
                verify_manifest_tree(tree)

            traversal = root / "traversal"
            _write_manifest(traversal, {"a.txt": b"x"})
            digest = hashlib.sha256(b"x").hexdigest()
            for unsafe in ("../x", "C:/x"):
                manifest = {"files": [{"path": unsafe, "bytes": 1, "sha256": digest}]}
                (traversal / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
                with self.assertRaises(ValueError):
                    verify_manifest_tree(traversal)

            duplicate = root / "duplicate"
            _write_manifest(duplicate, {"a.txt": b"x"})
            data = json.loads((duplicate / "manifest.json").read_text(encoding="utf-8"))
            data["files"].append(dict(data["files"][0]))
            (duplicate / "manifest.json").write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_manifest_tree(duplicate)

    def test_manifest_rejects_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tree = root / "tree"
            tree.mkdir()
            target = root / "target"
            target.write_text("x", encoding="utf-8")
            try:
                (tree / "a.txt").symlink_to(target)
            except (OSError, NotImplementedError):
                self.skipTest("symlinks unavailable")
            row = {"path": "a.txt", "bytes": 1, "sha256": hashlib.sha256(b"x").hexdigest()}
            (tree / "manifest.json").write_text(json.dumps({"files": [row]}), encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_manifest_tree(tree)

    def test_build_pack_html_escapes_and_reports_distinct_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline, panel, audit, export, cif, csv = _fixtures(root, "C[CH2:7]O")
            panel_json = json.loads((panel / "panel.json").read_text(encoding="utf-8"))
            panel_json["cases"][0]["candidate_id"] = '<script>alert("x")</script>'
            files = {
                p.relative_to(panel).as_posix(): p.read_bytes()
                for p in panel.rglob("*")
                if p.is_file() and p.name != "manifest.json"
            }
            files["panel.json"] = _j(panel_json)
            for path in sorted(panel.rglob("*"), reverse=True):
                if path.is_file():
                    path.unlink()
                elif path.is_dir():
                    path.rmdir()
            panel.rmdir()
            _write_manifest(panel, files)

            output = root / "pack"
            native = lambda _: {
                "manifest_sha256": verify_manifest_tree(baseline)["manifest_sha256"],
                "result_sha256": "r",
                "protocol_sha256": "p",
            }
            with patch("scripts.package_native_parent_review._molecule_svg", return_value="<svg></svg>"):
                build_pack(
                    output, panel, baseline, audit, export, cif, csv, None,
                    source_pair_verify=lambda _: ({"ok": True}, {}, {}, {}),
                    native_probe_verify=native,
                )
            page = (output / "index.html").read_text(encoding="utf-8")
            self.assertIn("<title>Native parent N3 12-case review</title>", page)
            self.assertNotIn('<script>alert("x")</script>', page)
            self.assertIn("&lt;script&gt;alert", page)
            self.assertIn("<b>9/12</b>pose-pass cases", page)
            self.assertIn("<b>18</b>total pass poses", page)
            self.assertIn("Parent 2/3, Br 1/3 seed-level pass", page)
            manifest = verify_manifest_tree(output)
            self.assertIn("README.md", manifest["files"])
            self.assertIn("index.html", manifest["files"])
            self.assertNotIn("manifest.json", manifest["files"])
            with self.assertRaises(FileExistsError):
                build_pack(output, panel, baseline, audit, export, cif, csv, None)

    def test_build_rejects_output_containment(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline, panel, audit, export, cif, csv = _fixtures(root)
            with self.assertRaises(ValueError):
                build_pack(panel / "pack", panel, baseline, audit, export, cif, csv, None)
            with self.assertRaises(ValueError):
                build_pack(root, panel, baseline, audit, export, cif, csv, None)

    def test_failure_receipt_preserves_partial_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline, panel, audit, export, cif, csv = _fixtures(root)
            output = root / "pack"
            native = lambda _: {
                "manifest_sha256": verify_manifest_tree(baseline)["manifest_sha256"],
                "result_sha256": "r",
                "protocol_sha256": "p",
            }

            def remove_csv(_: Path):
                csv.unlink()
                return {"ok": True}

            with patch("scripts.package_native_parent_review._molecule_svg", return_value="<svg></svg>"):
                with self.assertRaises(ValueError):
                    build_pack(
                        output, panel, baseline, audit, export, cif, csv, None,
                        source_pair_verify=remove_csv,
                        native_probe_verify=native,
                    )
            self.assertTrue(output.is_dir())
            self.assertTrue((output / "native-baseline").is_dir())
            failure = json.loads((output / "packaging-failure.json").read_text(encoding="utf-8"))
            self.assertEqual(failure["exception_type"], "ValueError")
            self.assertTrue(failure["exception_text"])


if __name__ == "__main__":
    unittest.main()
