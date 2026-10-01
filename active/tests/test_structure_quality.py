"""Synthetic diagnostic edge cases; these fixtures are not molecular predictions."""
import copy
import tempfile
from pathlib import Path
import unittest

from packages.science.handoff import ACTIVE, encoded, sha
from packages.science.evidence_common import file_index, seal
from packages.science.structure_quality import (
    REQUEST, collect_quality, project_request, read_quality, strict_json,
)
from packages.science.quality_evidence import project, markdown, export_quality_evidence, read_quality_evidence


def fixture(chains=("X", "Y"), *, missing=False, clash=True):
    files = {}
    def add(path, value):
        data = value if isinstance(value, bytes) else encoded(value)
        files[path] = data
        return {"path": path, "sha256": sha(data)}
    original = add("raw/original.cif", b"synthetic source; never used for molecular validation")
    hydrogen = add("raw/hydrogen.cif", b"synthetic prepared structure")
    atoms = [{"chainID": chains[0], "resID": 301, "iCode": "", "resName": "ALA", "atomName": "O", "alt": ""},
             {"chainID": chains[1], "resID": 702, "iCode": "", "resName": "LYS", "atomName": "NZ", "alt": ""}]
    inventory_rows = [{"atom": a, "element": e} for a, e in zip(atoms, ("O", "N"))]
    inventory_rows += [{"atom": {**a, "atomName": "H"}, "element": "H"} for a in atoms]
    inventory = add("raw/inventory.json", {"format": "tpd-quality-atom-inventory/0.1.0-draft", "structure_sha256": hydrogen["sha256"], "atoms": inventory_rows})
    def record(path, command, output):
        return add(path, {"exit_code": 0, "command": command, "target_exists": True, "target_sha256": output["sha256"]})
    reduce = record("raw/reduce.json", ["mmtbx.reduce2", "use_neutron_distances=False", "add_flip_movers=False"], hydrogen)
    log = b"Number of hydrogen atoms added to the input model: 2\n"
    if missing:
        log += b"The following H atoms were not placed because they could not be parameterized\n(not enough restraints information)\n H7 UNK Z   5 \n\n"
    directions = {}
    for name, source, target in [("forward", *atoms), ("reverse", *reversed(atoms))]:
        rows = [{"group": "1->2", "src": source, "target": target, "type": t, "gap": -0.55, "dotCount": 3}
                for t in (["bo", "hb"] if clash else ["hb"])]
        output = add(f"raw/{name}.json", {"flat_results": rows})
        command = ["mmtbx.probe2", "/reported/hydrogen.cif", "use_neutron_distances=False", "approach=once",
                   "output.format=json", "output.condensed=True", "source_selection=chain " + source["chainID"],
                   "target_selection=chain " + target["chainID"]]
        directions[name] = {"execution": record(f"raw/{name}-execution.json", command, output), "output": output}
    request = {"format": REQUEST, "data_mode": "synthetic_test",
               "binding": {"compound_id": "synthetic-compound", "molecule_id": "synthetic-molecule", "run_id": "synthetic-run",
                           "model_rank": 0, "input_sha256": "1" * 64, "structure_sha256": original["sha256"]},
               "source_structure": original,
               "method": {"adapter": "cctbx-reduce2-probe2-condensed/1", "tool_version": "2025.11", "settings_note": "Synthetic test only"},
               "preparation_evidence": [add("raw/preparation.json", {"data_mode": "synthetic_test"})],
               "preparation_notes": ["Synthetic fixtures: no scientific result."],
               "analyses": [{"analysis_id": "interface-xray", "hydrogen_convention": "xray",
                             "scope": {"description": "Synthetic chain interface", "source_chains": [chains[0]], "target_chains": [chains[1]],
                                       "atom_typing": "standard_protein", "limitations": []},
                             "hydrogen_model": hydrogen, "hydrogen_inventory": inventory,
                             "hydrogen_model_source_path": "/reported/hydrogen.cif", "reduce_execution": reduce,
                             "reduce_log": add("raw/reduce.log", log), "directions": directions}]}
    return request, files


def replace_ref(files, ref, value):
    data = value if isinstance(value, bytes) else encoded(value)
    files[ref["path"]] = data
    ref["sha256"] = sha(data)


def change_output(request, files, direction, value):
    entry = request["analyses"][0]["directions"][direction]
    replace_ref(files, entry["output"], value)
    receipt = strict_json(files[entry["execution"]["path"]])
    receipt["target_sha256"] = entry["output"]["sha256"]
    replace_ref(files, entry["execution"], receipt)


def synthetic_evidence(request):
    b = request["binding"]
    model = {"rank": 0, "status": "parsed", "summary": {"warnings": []}, "comparison_file": {"path": "comparison.json"}}
    run = {"run_id": b["run_id"], "compound_id": b["compound_id"], "molecule_id": b["molecule_id"],
           "input_sha256": b["input_sha256"], "origin": "synthetic_test", "models": [model]}
    evidence = {"digest": "2" * 64, "case_id": "synthetic-case", "data_mode": "synthetic_test",
                "candidates": [{"compound": {"compound_id": b["compound_id"], "molecule_id": b["molecule_id"]},
                                "prediction": {"status": "available_unreviewed"}},
                               {"compound": {"compound_id": "unrun-compound", "molecule_id": "unrun-molecule"},
                                "prediction": {"status": "not_run"}}]}
    files = {"collection.json": encoded({"runs": [run]}), "manifest.json": b"synthetic base manifest",
             "comparison.json": encoded({"source_manifest": {"files": {"structure": {"sha256": b["structure_sha256"]}}}})}
    return evidence, files


class StructureQualityTests(unittest.TestCase):
    def test_hbond_does_not_cancel_overlap_and_directions_deduplicate(self):
        request, files = fixture()
        report = project_request(request, files)
        analysis = report["analyses"][0]
        self.assertEqual(analysis["bad_overlap_pair_count"], 1)
        self.assertEqual(analysis["bad_overlap_pairs"][0]["types"], ["bo", "hb"])
        self.assertEqual(set(analysis["bad_overlap_pairs"][0]["directions"]), {"forward", "reverse"})
        self.assertFalse(report["summary"]["automatic_acceptance"])

    def test_chain_residue_and_candidate_names_do_not_control_detection(self):
        request, files = fixture(("Protein17", "Other88"))
        report = project_request(request, files)
        self.assertTrue(report["summary"]["has_findings"])
        self.assertEqual(report["analyses"][0]["bad_overlap_pairs"][0]["atoms"][0][0], "Other88")

    def test_no_clashes_means_only_no_findings_within_scope(self):
        request, files = fixture(clash=False)
        summary = project_request(request, files)["summary"]
        self.assertEqual(summary["interpretation_status"], "no_findings_in_declared_scope")
        self.assertFalse(summary["whole_structure_validated"])

    def test_unplaced_h_and_findings_are_visible_together(self):
        request, files = fixture(missing=True)
        report = project_request(request, files)
        self.assertTrue(report["summary"]["has_findings"])
        self.assertTrue(report["summary"]["has_incomplete_checks"])
        self.assertIn("H7 UNK Z   5", report["analyses"][0]["unplaced_hydrogens"])

    def test_no_result_and_no_h_preparation_never_clean(self):
        request, files = fixture(clash=False)
        a = request["analyses"][0]
        for name in ["hydrogen_model", "hydrogen_inventory", "reduce_execution", "reduce_log"]:
            a[name] = None
        a["directions"] = {d: {"execution": None, "output": None} for d in ["forward", "reverse"]}
        report = project_request(request, files)
        self.assertEqual(report["analyses"][0]["finding_status"], "not_assessed")
        self.assertTrue(report["summary"]["has_unrun_checks"])
        self.assertEqual(report["summary"]["interpretation_status"], "review_required")

    def test_partial_direction_preserves_findings_without_full_coverage(self):
        request, files = fixture()
        request["analyses"][0]["directions"]["reverse"] = {"execution": None, "output": None}
        analysis = project_request(request, files)["analyses"][0]
        self.assertEqual(analysis["execution_status"], "partial")
        self.assertEqual(analysis["coverage"], "incomplete")
        self.assertEqual(analysis["bad_overlap_pair_count"], 1)

    def test_failed_command_and_malformed_output_are_recorded_not_clean(self):
        for malformed in [b"not json", b'{"flat_results":[{"gap":NaN}]}']:
            request, files = fixture(clash=False)
            change_output(request, files, "forward", malformed)
            report = project_request(request, files)
            self.assertTrue(report["summary"]["has_failed_checks"])
            self.assertEqual(report["summary"]["interpretation_status"], "review_required")
        request, files = fixture(clash=False)
        ref = request["analyses"][0]["directions"]["forward"]["execution"]
        receipt = strict_json(files[ref["path"]]); receipt["exit_code"] = 1
        replace_ref(files, ref, receipt)
        self.assertTrue(project_request(request, files)["summary"]["has_failed_checks"])

    def test_out_of_scope_or_nonexistent_atoms_cannot_be_counted(self):
        for key, value in [("chainID", "Z"), ("atomName", "FAKE")]:
            request, files = fixture()
            output = strict_json(files[request["analyses"][0]["directions"]["forward"]["output"]["path"]])
            for row in output["flat_results"]:
                row["src"][key] = value
            change_output(request, files, "forward", output)
            self.assertEqual(project_request(request, files)["analyses"][0]["directions"]["forward"]["status"], "failed")

    def test_empty_selection_cannot_make_empty_output_look_clean(self):
        request, files = fixture(clash=False)
        ref = request["analyses"][0]["hydrogen_inventory"]
        inventory = strict_json(files[ref["path"]])
        inventory["atoms"] = [x for x in inventory["atoms"] if x["atom"]["chainID"] == "X"]
        replace_ref(files, ref, inventory)
        for d in ["forward", "reverse"]:
            change_output(request, files, d, {"flat_results": []})
        report = project_request(request, files)
        self.assertTrue(report["summary"]["has_failed_checks"])
        self.assertTrue(report["summary"]["has_incomplete_checks"])

    def test_H_convention_and_source_selection_mismatch_fail(self):
        for setting in ["use_neutron_distances=True", "source_selection=chain Y"]:
            request, files = fixture()
            ref = request["analyses"][0]["directions"]["forward"]["execution"]
            receipt = strict_json(files[ref["path"]])
            prefix = setting.split("=")[0] + "="
            receipt["command"] = [setting if s.startswith(prefix) else s for s in receipt["command"]]
            replace_ref(files, ref, receipt)
            self.assertTrue(project_request(request, files)["summary"]["has_failed_checks"])

    def test_unknown_reduce_log_is_incomplete_even_with_no_clash(self):
        request, files = fixture(clash=False)
        replace_ref(files, request["analyses"][0]["reduce_log"], b"Different tool output\n")
        report = project_request(request, files)
        self.assertTrue(report["summary"]["has_incomplete_checks"])
        self.assertEqual(report["summary"]["interpretation_status"], "review_required")

    def test_reduce_failure_remains_failed_even_if_probe_output_exists(self):
        request, files = fixture(clash=False)
        ref = request["analyses"][0]["reduce_execution"]
        record = strict_json(files[ref["path"]]); record["exit_code"] = 1
        replace_ref(files, ref, record)
        report = project_request(request, files)
        self.assertEqual(report["analyses"][0]["preparation_status"], "failed")
        self.assertTrue(report["summary"]["has_failed_checks"])

    def test_source_and_inventory_hash_mismatch_rejected_or_flagged(self):
        request, files = fixture()
        request["binding"]["structure_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "STRUCTURE_BINDING"):
            project_request(request, files)
        request, files = fixture()
        files[request["source_structure"]["path"]] += b"modified"
        with self.assertRaisesRegex(ValueError, "INPUT_HASH"):
            project_request(request, files)

    def test_duplicates_and_unknown_format_rejected(self):
        with self.assertRaisesRegex(ValueError, "DUPLICATE_JSON_KEY"):
            strict_json(b'{"flat_results":[],"flat_results":[]}')
        request, files = fixture()
        request["analyses"].append(copy.deepcopy(request["analyses"][0]))
        with self.assertRaisesRegex(ValueError, "DUPLICATE_ANALYSIS"):
            project_request(request, files)

    def test_archive_roundtrip_synthetic_guard_and_report_reprojection(self):
        request, files = fixture(missing=True)
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for name, data in files.items():
                path = directory / name; path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(data)
            (directory / "request.json").write_bytes(encoded(request))
            with self.assertRaisesRegex(ValueError, "SYNTHETIC_QUALITY"):
                collect_quality(directory / "request.json", directory / "blocked")
            out = directory / "assessment"
            report = collect_quality(directory / "request.json", out, allow_synthetic=True)
            self.assertEqual(read_quality(out, allow_synthetic=True), report)
            with self.assertRaisesRegex(ValueError, "OUTPUT_EXISTS"):
                collect_quality(directory / "request.json", out, allow_synthetic=True)
            report["summary"]["has_findings"] = False
            (out / "report.json").write_bytes(encoded(report))
            # Even a resealed archive cannot substitute a changed derived view.
            manifest = strict_json((out / "manifest.json").read_bytes())
            manifest["files"]["report.json"] = file_index({"report.json": encoded(report)})["report.json"]
            manifest = seal({k: v for k, v in manifest.items() if k != "digest"})
            (out / "manifest.json").write_bytes(encoded(manifest))
            with self.assertRaisesRegex(ValueError, "REPORT_PROJECTION"):
                read_quality(out, allow_synthetic=True)

    def test_path_traversal_during_collect_rejected(self):
        request, files = fixture()
        request["source_structure"]["path"] = "../outside.cif"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "request.json"; path.write_bytes(encoded(request))
            with self.assertRaisesRegex(ValueError, "RELATIVE_PATH"):
                collect_quality(path, Path(tmp) / "output", allow_synthetic=True)


class QualityEvidenceProjectionTests(unittest.TestCase):
    def test_stored_cpu_evidence_roundtrip_keeps_unrun_candidates_and_source_bytes(self):
        source = ACTIVE / "outputs/b_evidence_20260923"
        before = {name: (source / name).read_bytes() for name in ["manifest.json", "README.md", "evidence.json"]}
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "검증 연결 자료"
            report = export_quality_evidence(source, [], out)
            self.assertEqual(read_quality_evidence(out), report)
            self.assertTrue(all(not c["samples"] for c in report["candidates"]))
            for name, data in before.items():
                self.assertEqual((source / name).read_bytes(), data)
                self.assertEqual((out / "evidence" / name).read_bytes(), data)
            with self.assertRaisesRegex(ValueError, "OUTPUT_EXISTS"):
                export_quality_evidence(source, [], out)

    def test_changed_quality_report_and_unindexed_file_are_rejected(self):
        source = ACTIVE / "outputs/b_evidence_20260923"
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            export_quality_evidence(source, [], out)
            (out / "README.md").write_text("검증 완료", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "FILE_HASH"):
                read_quality_evidence(out)
            manifest = strict_json((out / "manifest.json").read_bytes())
            data = (out / "README.md").read_bytes()
            manifest["files"]["README.md"] = {"sha256": sha(data), "bytes": len(data)}
            manifest = seal({k: v for k, v in manifest.items() if k != "digest"})
            (out / "manifest.json").write_bytes(encoded(manifest))
            with self.assertRaisesRegex(ValueError, "REPORT_PROJECTION"):
                read_quality_evidence(out)
            (out / "extra.txt").write_text("unindexed")
            with self.assertRaisesRegex(ValueError, "UNINDEXED"):
                read_quality_evidence(out)

    def test_all_candidates_and_missing_quality_remain_visible(self):
        request, files = fixture()
        evidence, base = synthetic_evidence(request)
        result = project(evidence, base, [])
        self.assertEqual(result["candidates"][0]["samples"][0]["quality_status"], "not_assessed")
        self.assertIsNone(result["candidates"][0]["samples"][0]["quality_summary"])
        self.assertEqual(result["candidates"][1]["prediction_status"], "not_run")
        self.assertIn("미검증", markdown(result).decode())

    def test_join_is_exact_and_does_not_promote_authority(self):
        request, files = fixture(missing=True)
        evidence, base = synthetic_evidence(request)
        result = project(evidence, base, [("quality/q0000", project_request(request, files))])
        sample = result["candidates"][0]["samples"][0]
        self.assertTrue(sample["quality_summary"]["has_findings"])
        self.assertTrue(sample["quality_summary"]["has_incomplete_checks"])
        self.assertFalse(result["authority"]["approval_record_created"])
        self.assertIn("H7 UNK Z   5", markdown(result).decode())

    def test_wrong_candidate_molecule_input_or_structure_cannot_join(self):
        request, files = fixture()
        evidence, base = synthetic_evidence(request)
        for field in ["compound_id", "molecule_id", "input_sha256", "structure_sha256"]:
            report = project_request(request, files)
            report["binding"][field] = "wrong"
            with self.assertRaisesRegex(ValueError, "SAMPLE_BINDING"):
                project(evidence, base, [("quality/q0000", report)])

    def test_duplicate_or_orphan_reports_cannot_silently_replace_results(self):
        request, files = fixture()
        evidence, base = synthetic_evidence(request)
        report = project_request(request, files)
        with self.assertRaisesRegex(ValueError, "DUPLICATE_SAMPLE"):
            project(evidence, base, [("q1", report), ("q2", report)])
        report["binding"]["model_rank"] = 1
        with self.assertRaisesRegex(ValueError, "ORPHAN_SAMPLE"):
            project(evidence, base, [("q1", report)])

    def test_synthetic_and_real_cannot_be_relabelled(self):
        request, files = fixture()
        evidence, base = synthetic_evidence(request)
        report = project_request(request, files); report["data_mode"] = "real"
        with self.assertRaisesRegex(ValueError, "ORIGIN_MISMATCH"):
            project(evidence, base, [("q1", report)])

    def test_missing_structure_cannot_receive_quality_report(self):
        request, files = fixture()
        evidence, base = synthetic_evidence(request)
        collection = strict_json(base["collection.json"])
        model = collection["runs"][0]["models"][0]
        model.update(status="missing", summary=None, comparison_file=None)
        base["collection.json"] = encoded(collection)
        with self.assertRaisesRegex(ValueError, "UNPARSED_MODEL"):
            project(evidence, base, [("q1", project_request(request, files))])
