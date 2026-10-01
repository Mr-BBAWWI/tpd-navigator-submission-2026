"""Offline, target-independent consumption of Reduce2/Probe2 diagnostics.

This does not execute or authenticate the tools, repair coordinates, or approve
structures. Tool findings and incomplete coverage are independent dimensions.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import re

from .evidence_common import file_index, safe_relative, schema_check, seal, verify_seal, write_new_directory
from .handoff import check, encoded, sha

SCHEMA = "b_structure_quality.schema.json"
REQUEST = "tpd-probe-quality-request/0.1.0-draft"
REPORT = "tpd-structure-quality/0.1.0-draft"
ARCHIVE = "tpd-structure-quality-archive/0.1.0-draft"
PROBE_TYPES = {"wc", "cc", "so", "bo", "hb", "wh"}


def strict_json(data):
    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            check(key not in result, "QUALITY_DUPLICATE_JSON_KEY")
            result[key] = value
        return result
    def invalid_number(value):
        raise ValueError("QUALITY_NONFINITE_JSON_NUMBER: " + value)
    return json.loads(data, object_pairs_hook=unique_pairs, parse_constant=invalid_number)


def refs(request):
    yield request["source_structure"]
    yield from request["preparation_evidence"]
    for analysis in request["analyses"]:
        for name in ("hydrogen_model", "hydrogen_inventory", "reduce_execution", "reduce_log"):
            if analysis[name] is not None:
                yield analysis[name]
        for direction in analysis["directions"].values():
            for name in ("execution", "output"):
                if direction[name] is not None:
                    yield direction[name]


def load_request_files(request, directory):
    schema_check(request, SCHEMA)
    check(request["format"] == REQUEST, "QUALITY_REQUEST_FORMAT")
    files = {}
    for ref in refs(request):
        data = safe_relative(directory, ref["path"]).read_bytes()
        check(sha(data) == ref["sha256"], "QUALITY_INPUT_HASH_MISMATCH")
        check(ref["path"] not in {"request.json", "report.json", "manifest.json"}, "QUALITY_RESERVED_PATH")
        files[ref["path"]] = data
    return files


def atom_key(atom):
    required = {"chainID", "resID", "iCode", "resName", "atomName", "alt"}
    check(isinstance(atom, dict) and set(atom) == required, "QUALITY_PROBE_ATOM_FIELDS")
    check(type(atom["resID"]) is int, "QUALITY_PROBE_RESIDUE")
    check(all(isinstance(atom[k], str) for k in required - {"resID"}), "QUALITY_PROBE_ATOM_TYPES")
    check(bool(atom["chainID"]) and bool(atom["resName"].strip()) and bool(atom["atomName"].strip()), "QUALITY_PROBE_ATOM_EMPTY")
    return (atom["chainID"], atom["resID"], atom["iCode"], atom["resName"].strip(), atom["atomName"].strip(), atom["alt"])


def inventory_from_cif(data):
    """Scientific producer step; Gemmi is optional for portable report readers."""
    import gemmi
    structure = gemmi.make_structure_from_block(gemmi.cif.read_string(data.decode("utf-8")).sole_block())
    check(len(structure) == 1, "QUALITY_SINGLE_MODEL_REQUIRED")
    rows = []
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                if atom.occ <= 0:
                    continue
                rows.append({"atom": {"chainID": chain.name, "resID": residue.seqid.num,
                             "iCode": residue.seqid.icode.strip(), "resName": residue.name,
                             "atomName": atom.name, "alt": atom.altloc.replace("\x00", "").strip()},
                             "element": atom.element.name})
    return {"format": "tpd-quality-atom-inventory/0.1.0-draft", "structure_sha256": sha(data),
            "producer": "gemmi", "version": gemmi.__version__, "atoms": rows}


def inventory_atoms(data, expected_sha):
    inventory = strict_json(data)
    check(inventory.get("format") == "tpd-quality-atom-inventory/0.1.0-draft"
          and inventory.get("structure_sha256") == expected_sha, "QUALITY_INVENTORY_BINDING")
    rows = inventory["atoms"]
    check(isinstance(rows, list) and bool(rows), "QUALITY_EMPTY_INVENTORY")
    result = {}
    for row in rows:
        key = atom_key(row["atom"])
        check(key not in result and isinstance(row.get("element"), str) and bool(row["element"]), "QUALITY_INVENTORY_ATOM")
        result[key] = row["element"]
    return result


def probe_pairs(data, source_chains, target_chains):
    raw = strict_json(data)
    check(isinstance(raw, dict) and set(raw) == {"flat_results"} and isinstance(raw["flat_results"], list), "QUALITY_PROBE_FORMAT")
    pairs = {}
    for row in raw["flat_results"]:
        check(isinstance(row, dict) and row.get("type") in PROBE_TYPES, "QUALITY_PROBE_CONTACT_TYPE")
        check(row.get("group") == "1->2", "QUALITY_PROBE_DIRECTION")
        source, target = atom_key(row["src"]), atom_key(row["target"])
        check(source[0] in source_chains and target[0] in target_chains, "QUALITY_PROBE_OUTSIDE_SCOPE")
        check(type(row.get("dotCount")) is int and row["dotCount"] > 0, "QUALITY_PROBE_DOT_COUNT")
        check(type(row.get("gap")) in (int, float) and math.isfinite(row["gap"]), "QUALITY_PROBE_GAP")
        key = tuple(sorted((source, target)))
        item = pairs.setdefault(key, {"atoms": [list(x) for x in key], "types": set(), "gap_A": row["gap"]})
        check(abs(item["gap_A"] - row["gap"]) < 0.002, "QUALITY_PROBE_INCONSISTENT_PAIR_GAP")
        item["types"].add(row["type"])
    return pairs


def missing_hydrogens(log):
    counts = re.findall(r"Number of hydrogen atoms added to the input model:\s*(\d+)", log)
    check(len(counts) == 1 and int(counts[0]) > 0, "QUALITY_REDUCE_H_COUNT_UNRECOGNIZED")
    marker = "The following H atoms were not placed because they could not be parameterized"
    missing = []
    if marker in log:
        block = log.split(marker, 1)[1]
        check("(not enough restraints information)" in block, "QUALITY_REDUCE_WARNING_UNRECOGNIZED")
        block = block.split("(not enough restraints information)", 1)[1]
        for line in block.splitlines():
            if not line.strip():
                if missing:
                    break
                continue
            check(re.fullmatch(r"\s*\S+\s+\S+\s+\S+\s*-?\d+\s*", line) is not None, "QUALITY_REDUCE_WARNING_UNRECOGNIZED")
            missing.append(line.strip())
        check(bool(missing), "QUALITY_REDUCE_WARNING_EMPTY")
    return int(counts[0]), missing


def command_value(command, name):
    values = [v.split("=", 1)[1] for v in command if v.startswith(name + "=")]
    check(len(values) == 1, "QUALITY_COMMAND_SETTING: " + name)
    return values[0]


def execution(files, ref, output, tool, convention):
    record = strict_json(files[ref["path"]])
    check(type(record.get("exit_code")) is int and record["exit_code"] == 0, "QUALITY_TOOL_FAILED")
    check(output is not None and record.get("target_exists") is True
          and record.get("target_sha256") == output["sha256"], "QUALITY_TOOL_OUTPUT_BINDING")
    command = record.get("command")
    check(isinstance(command, list) and all(isinstance(c, str) for c in command)
          and any(c.replace("\\", "/").rsplit("/", 1)[-1] == tool for c in command), "QUALITY_TOOL_COMMAND")
    check(command_value(command, "use_neutron_distances") == ("True" if convention == "nuclear" else "False"), "QUALITY_H_CONVENTION_MISMATCH")
    return command


def selected_chains(command, field):
    value = command_value(command, field)
    terms = value.split(" or ")
    check(all(re.fullmatch(r"chain [A-Za-z0-9_]+", term) for term in terms), "QUALITY_UNSUPPORTED_SELECTION")
    chains = [term[6:] for term in terms]
    check(len(set(chains)) == len(chains), "QUALITY_DUPLICATE_SELECTION")
    return set(chains)


def project_analysis(analysis, files):
    scope = analysis["scope"]
    source, target = set(scope["source_chains"]), set(scope["target_chains"])
    check(not source & target, "QUALITY_SCOPE_OVERLAP")
    convention = analysis["hydrogen_convention"]
    missing, count, preparation_error = [], None, None
    atoms = None
    try:
        check(analysis["hydrogen_model"] is not None and analysis["hydrogen_inventory"] is not None, "QUALITY_H_MODEL_NOT_PROVIDED")
        atoms = inventory_atoms(files[analysis["hydrogen_inventory"]["path"]], analysis["hydrogen_model"]["sha256"])
        check(source | target <= {k[0] for k in atoms}, "QUALITY_SELECTION_HAS_NO_ATOMS")
        check(analysis["reduce_execution"] is not None and analysis["reduce_log"] is not None, "QUALITY_H_PREPARATION_NOT_PROVIDED")
        command = execution(files, analysis["reduce_execution"], analysis["hydrogen_model"], "mmtbx.reduce2", convention)
        check(command_value(command, "add_flip_movers") == "False", "QUALITY_UNSUPPORTED_HEAVY_ATOM_FLIPS")
        count, missing = missing_hydrogens(files[analysis["reduce_log"]["path"]].decode("utf-8"))
        check(sum(e in {"H", "D"} for e in atoms.values()) == count, "QUALITY_H_COUNT_MISMATCH")
    except (ValueError, KeyError, TypeError, UnicodeError) as error:
        preparation_error = str(error)
    directions, union = {}, {}
    for name, entry in analysis["directions"].items():
        if entry["execution"] is None and entry["output"] is None:
            directions[name] = {"status": "not_run", "reason": "TOOL_RECORD_NOT_PROVIDED"}
            continue
        try:
            check(entry["execution"] is not None, "QUALITY_EXECUTION_RECORD_MISSING")
            command = execution(files, entry["execution"], entry["output"], "mmtbx.probe2", convention)
            check(command_value(command, "approach") == "once" and command_value(command, "output.format") == "json"
                  and command_value(command, "output.condensed") == "True", "QUALITY_UNSUPPORTED_PROBE_MODE")
            check(analysis["hydrogen_model"] is not None
                  and analysis["hydrogen_model_source_path"] in command, "QUALITY_HYDROGEN_INPUT_BINDING")
            a, b = (source, target) if name == "forward" else (target, source)
            check(atoms is not None and source | target <= {k[0] for k in atoms}, "QUALITY_SELECTION_HAS_NO_ATOMS")
            check(selected_chains(command, "source_selection") == a and selected_chains(command, "target_selection") == b, "QUALITY_SELECTION_MISMATCH")
            pairs = probe_pairs(files[entry["output"]["path"]], a, b)
            check(all(atom in atoms for pair in pairs for atom in pair), "QUALITY_PROBE_ATOM_NOT_IN_MODEL")
            for key, item in pairs.items():
                if key in union:
                    check(abs(union[key]["gap_A"] - item["gap_A"]) < 0.002, "QUALITY_DIRECTION_GAP_MISMATCH")
            for key, item in pairs.items():
                stored = union.setdefault(key, {**item, "types": set(), "directions": []})
                stored["types"].update(item["types"])
                if "bo" in item["types"]:
                    stored["directions"].append(name)
            directions[name] = {"status": "completed", "reason": None}
        except (ValueError, KeyError, TypeError, UnicodeError) as error:
            directions[name] = {"status": "failed", "reason": str(error)}
    completed = sum(d["status"] == "completed" for d in directions.values())
    state = "completed" if completed == 2 else "partial" if completed else "not_run" if all(d["status"] == "not_run" for d in directions.values()) else "failed"
    limitations = list(scope["limitations"])
    if preparation_error:
        limitations.append("HYDROGEN_PREPARATION_UNVERIFIED")
    if missing:
        limitations.append("HYDROGENS_UNPLACED")
    if scope["atom_typing"] == "unverified":
        limitations.append("ATOM_TYPING_UNVERIFIED")
    elif atoms is not None:
        standard = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
        if any(k[3] not in standard for k in atoms if k[0] in source | target):
            limitations.append("NONSTANDARD_RESIDUE_TYPING_UNVERIFIED")
    if completed != 2:
        limitations.append("BIDIRECTIONAL_COVERAGE_INCOMPLETE")
    # A partial result may still contain findings; never erase them or call it clean.
    bad = [{**v, "types": sorted(v["types"])} for k, v in sorted(union.items()) if "bo" in v["types"]]
    coverage = "incomplete" if limitations else "declared_scope_only"
    preparation_status = "completed" if preparation_error is None else "not_run" if analysis["reduce_execution"] is None else "failed"
    return {"analysis_id": analysis["analysis_id"], "scope": copy.deepcopy(scope), "hydrogen_convention": convention,
            "execution_status": state, "coverage": coverage, "directions": directions,
            "preparation_status": preparation_status, "hydrogens_added": count, "unplaced_hydrogens": missing, "preparation_error": preparation_error,
            "finding_status": "findings" if bad else "none_observed" if state == "completed" else "not_assessed",
            "bad_overlap_pair_count": len(bad), "bad_overlap_pairs": bad, "limitations": sorted(set(limitations)),
            "atom_namespace": "Probe output auth chain/residue/atom names; hydrogen-prepared structure, not remapped to a reference"}


def project_request(request, files):
    schema_check(request, SCHEMA)
    check(request["format"] == REQUEST, "QUALITY_REQUEST_FORMAT")
    check(request["binding"]["structure_sha256"] == request["source_structure"]["sha256"], "QUALITY_STRUCTURE_BINDING")
    ids = [a["analysis_id"] for a in request["analyses"]]
    check(len(ids) == len(set(ids)), "QUALITY_DUPLICATE_ANALYSIS")
    for ref in refs(request):
        check(ref["path"] in files and sha(files[ref["path"]]) == ref["sha256"], "QUALITY_INPUT_HASH_MISMATCH")
    analyses = [project_analysis(a, files) for a in request["analyses"]]
    findings = any(a["bad_overlap_pair_count"] for a in analyses)
    incomplete = any(a["coverage"] == "incomplete" for a in analyses)
    summary = {"has_findings": findings, "has_incomplete_checks": incomplete,
               "has_failed_checks": any(a["preparation_status"] == "failed" or any(d["status"] == "failed" for d in a["directions"].values()) for a in analyses),
               "has_unrun_checks": any(a["preparation_status"] == "not_run" or any(d["status"] == "not_run" for d in a["directions"].values()) for a in analyses),
               "interpretation_status": "review_required" if findings or incomplete else "no_findings_in_declared_scope",
               "automatic_acceptance": False, "whole_structure_validated": False}
    report = seal({"format": REPORT, "data_mode": request["data_mode"], "binding": copy.deepcopy(request["binding"]),
                 "request_sha256": sha(encoded(request)), "method": copy.deepcopy(request["method"]),
                 "preparation_notes": copy.deepcopy(request["preparation_notes"]), "analyses": analyses, "summary": summary,
                 "provenance": "caller_reported_tool_execution_and_preparation_not_authenticated",
                 "human_review": "pending", "efficacy_established": False})
    schema_check(report, "b_structure_quality_report.schema.json")
    return report


def collect_quality(request_path, destination, *, allow_synthetic=False):
    request_path = Path(request_path)
    request = strict_json(request_path.read_bytes())
    check(request.get("data_mode") != "synthetic_test" or allow_synthetic, "SYNTHETIC_QUALITY_NOT_ALLOWED")
    files = load_request_files(request, request_path.parent)
    report = project_request(request, files)
    files.update({"request.json": encoded(request), "report.json": encoded(report)})
    files["manifest.json"] = encoded(seal({"format": ARCHIVE, "files": file_index(files), "report_digest": report["digest"]}))
    write_new_directory(destination, files)
    return report


def read_quality(directory, *, allow_synthetic=False):
    directory = Path(directory)
    manifest = strict_json((directory / "manifest.json").read_bytes())
    verify_seal(manifest)
    check(manifest["format"] == ARCHIVE, "QUALITY_ARCHIVE_FORMAT")
    files = {}
    for name, item in manifest["files"].items():
        data = safe_relative(directory, name).read_bytes()
        check(sha(data) == item["sha256"] and len(data) == item["bytes"], "QUALITY_ARCHIVE_HASH")
        files[name] = data
    actual = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file()}
    check(actual == set(files) | {"manifest.json"}, "QUALITY_ARCHIVE_UNINDEXED_FILES")
    request, report = [strict_json(files[n + ".json"]) for n in ("request", "report")]
    check(request.get("data_mode") != "synthetic_test" or allow_synthetic, "SYNTHETIC_QUALITY_NOT_ALLOWED")
    check(encoded(project_request(request, files)) == encoded(report), "QUALITY_REPORT_PROJECTION")
    check(manifest["report_digest"] == report["digest"], "QUALITY_REPORT_DIGEST")
    return report
