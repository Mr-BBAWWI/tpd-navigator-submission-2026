#!/usr/bin/env python3
"""Reinspect fixed-frame analog microstates against registered protein evidence."""
from __future__ import annotations

import argparse
import hashlib
import html
import io
import json
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import Chem
from packages.science import chemical_states
from packages.science import analog_state_comparison

ANALOG_STATES = {
    "W-c2afc5e73c1a": 2,
    "W-4c0a639c0a41": 2,
    "W-80f8f4a11b5d": 4,
}
HYDROPHOBIC_REQUIRED = {"VAL1408", "PHE1409", "ILE1470"}
FORMAT = "core-state-robustness/20261001.1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return value


def safe_file(root: Path, relative: str) -> Path:
    if not isinstance(relative, str) or not relative or "\x00" in relative:
        raise ValueError("invalid artifact path")
    item = Path(relative)
    if item.is_absolute() or any(part in ("", ".", "..") for part in item.parts):
        raise ValueError(f"unsafe artifact path: {relative}")
    root = root.resolve(strict=True)
    current = root
    for part in item.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError(f"symlink traversal rejected: {relative}")
    resolved = current.resolve(strict=True)
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"artifact escapes root: {relative}") from exc
    if not resolved.is_file():
        raise ValueError(f"artifact is not a regular file: {relative}")
    return resolved


def verify_mappings(root: Path, mappings: list[dict[str, Any]]) -> dict[str, str]:
    verified: dict[str, str] = {}
    for row in mappings:
        if not isinstance(row, dict) or "path" not in row or "sha256" not in row:
            raise ValueError("artifact mapping requires path and sha256")
        relative, expected = row["path"], row["sha256"]
        if not isinstance(expected, str) or len(expected) != 64:
            raise ValueError(f"invalid artifact hash for {relative!r}")
        path = safe_file(root, relative)
        actual = sha256(path)
        if actual != expected.lower():
            raise ValueError(f"artifact hash mismatch: {relative}")
        if relative in verified and verified[relative] != actual:
            raise ValueError(f"conflicting artifact mapping: {relative}")
        verified[relative] = actual
    return verified


def verify_comparison_manifest(root: Path, manifest: dict[str, Any]) -> dict[str, str]:
    mappings = []
    for key in ("source_files", "artifacts"):
        rows = manifest.get(key, [])
        if not isinstance(rows, list):
            raise ValueError(f"manifest {key} must be a list")
        mappings.extend(rows)
    verified = verify_mappings(root, mappings)
    if "comparison.json" in verified:
        comparison = safe_file(root, "comparison.json")
        if verified["comparison.json"] != sha256(comparison):
            raise ValueError("recorded comparison.json hash mismatch")
    return verified


def _source_hash(comparison: dict[str, Any]) -> str:
    source = comparison.get("protein_source", {})
    value = source.get("sha256")
    if not isinstance(value, str):
        raise ValueError("comparison protein_source hash is absent")
    return value


def verify_source_binding(comparison: dict[str, Any], evidence: dict[str, Any], root: Path) -> Path:
    expected = _source_hash(comparison)
    evidence_hash = evidence.get("sha256", {}).get("source_cif_file")
    if expected != evidence_hash:
        raise ValueError("comparison/protein-evidence source CIF mismatch")
    configured = comparison.get("verified_source", {}).get("current_input_binding", {})
    configured = configured.get("source_files", {}).get("design_sources/6HAZ.cif", {})
    if configured and configured.get("sha256") != expected:
        raise ValueError("graph binding/source CIF mismatch")
    source = safe_file(root, "sources/6HAZ.cif")
    if sha256(source) != expected:
        raise ValueError("comparison source CIF actual hash mismatch")
    if not evidence.get("state_flags", {}).get("source_heavy_coordinates_accepted"):
        raise ValueError("protein source-frame coordinates were not accepted")
    if evidence.get("failures") != []:
        raise ValueError("protein evidence contains failures")
    frame = evidence.get("source_frame", {})
    if frame.get("coordinate_tolerance_A") != 0.0:
        raise ValueError("protein source-frame tolerance is not exact")
    if frame.get("pdb_cif_maximum_displacement_A") != 0.0:
        raise ValueError("protein PDB/CIF coordinates changed")
    if frame.get("prepared_source_maximum_displacement_A") != 0.0:
        raise ValueError("prepared protein heavy coordinates changed")
    return source


def find_analogs(comparison: dict[str, Any]) -> dict[str, dict[str, Any]]:
    found: dict[str, dict[str, Any]] = {}

    def visit(value: Any) -> None:
        if isinstance(value, dict):
            analog = value.get("analog_id")
            if analog in ANALOG_STATES and isinstance(value.get("poses"), list):
                previous = found.get(analog)
                if previous is not None and previous is not value:
                    raise ValueError(f"duplicate analog record: {analog}")
                found[analog] = value
                return
            for child in value.values():
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)

    visit(comparison)
    if set(found) != set(ANALOG_STATES):
        raise ValueError("comparison does not contain exactly the three required analogs")
    return found


def load_sdf(path: Path) -> Chem.Mol:
    data = path.read_bytes()
    if b"V3000" not in data:
        raise ValueError(f"state artifact is not V3000: {path}")
    molecules = list(Chem.ForwardSDMolSupplier(
        io.BytesIO(data), removeHs=False, sanitize=True, strictParsing=True,
    ))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError(f"state artifact must contain one molecule: {path}")
    return molecules[0]


def heavy_map_coordinates(mol: Chem.Mol) -> dict[int, tuple[float, float, float]]:
    if mol.GetNumConformers() != 1:
        raise ValueError("state molecule requires one conformer")
    conformer = mol.GetConformer()
    result = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        atom_map = atom.GetAtomMapNum()
        if atom_map <= 0 or atom_map in result:
            raise ValueError("heavy atom maps must be unique and positive")
        point = conformer.GetAtomPosition(atom.GetIdx())
        result[atom_map] = (point.x, point.y, point.z)
    return result


def validate_state(before: Chem.Mol, after: Chem.Mol, metadata: dict[str, Any], expected_maps: set[int]) -> None:
    before_xyz = heavy_map_coordinates(before)
    after_xyz = heavy_map_coordinates(after)
    if set(before_xyz) != expected_maps or set(after_xyz) != expected_maps:
        raise ValueError("state SDF heavy atom-map set mismatch")
    if before_xyz != after_xyz:
        raise ValueError("state heavy coordinates were not exactly preserved")
    if Chem.GetFormalCharge(after) != metadata.get("formal_charge"):
        raise ValueError("state SDF formal charge mismatch")
    nitrogen_meta = metadata.get("mapped_nitrogens", {})
    for atom in after.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetAtomMapNum() <= 0:
            continue
        row = nitrogen_meta.get(str(atom.GetAtomMapNum()))
        attached = sum(n.GetAtomicNum() == 1 for n in atom.GetNeighbors())
        hydrogens = attached + int(atom.GetTotalNumHs(includeNeighbors=False))
        if not row or row.get("formal_charge") != atom.GetFormalCharge() or row.get("hydrogen_count") != hydrogens:
            raise ValueError(f"state nitrogen metadata mismatch at map {atom.GetAtomMapNum()}")


def load_protein(evidence: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    atoms = evidence.get("protein_atoms")
    hydrogens = evidence.get("protein_hydrogens")
    if not isinstance(atoms, list) or not atoms or not isinstance(hydrogens, list) or not hydrogens:
        raise ValueError("registered protein atoms and hydrogens are required")
    if any(
        not isinstance(row, dict)
        or row.get("parent_assignment_status") != "validated_source_parent"
        for row in hydrogens
    ):
        raise ValueError("every protein hydrogen requires a validated source parent")
    return atoms, hydrogens


def hydrophobic_residues(summary: list[dict[str, Any]]) -> set[str]:
    result = set()
    for row in summary:
        for atom_id in row.get("protein_atom_ids", []):
            fields = atom_id.split(":")
            if len(fields) == 4:
                result.add(fields[2].upper() + fields[1])
    return result


def classify_row(hbond: bool, hydrophobes: set[str], rmsd: float | None) -> bool:
    """Required contact evidence is independent of the auxiliary RMSD diagnostic."""
    return bool(hbond and HYDROPHOBIC_REQUIRED <= set(hydrophobes))


def supported_pose_indices(rows: list[dict[str, Any]], requested_states: int) -> list[int]:
    supported = []
    for pose in range(5):
        selected = [row for row in rows if row["pose_index"] == pose]
        indices = {row["state_index"] for row in selected}
        if (len(selected) == requested_states and indices == set(range(requested_states))
                and all(row["required_contact_evidence"] for row in selected)):
            supported.append(pose)
    return supported


def write_outputs(output: Path, summary: dict[str, Any], sdf_records: list[bytes], sources: dict[str, str]) -> None:
    output.mkdir(parents=False)
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    sdf_path = output / "all-states.sdf"
    with sdf_path.open("wb") as handle:
        for record in sdf_records:
            handle.write(record)
            if not record.endswith(b"\n"):
                handle.write(b"\n")
    rows = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(key)))}</td>" for key in
        ("analog_id", "pose_index", "state_index", "status", "required_contact_evidence",
         "auxiliary_core_RMSD_within_1A")) + "</tr>" for row in summary["rows"]
    )
    page = ("<!doctype html><meta charset='utf-8'><title>Core-state robustness</title>"
            "<h1>Core-state robustness diagnostic</h1>"
            "<p>Pose 0 is prior Vina order, not a scientific-best selection. pH 7.4 is a preparation scenario, not assay-state proof.</p>"
            "<table border='1'><tr><th>analog</th><th>pose</th><th>state</th><th>status</th>"
            "<th>required contacts</th><th>RMSD ≤1 Å</th></tr>" + rows + "</table>")
    html_path = output / "index.html"
    html_path.write_text(page, encoding="utf-8")
    artifacts = []
    for name in ("summary.json", "index.html", "all-states.sdf"):
        path = output / name
        artifacts.append({"path": name, "sha256": sha256(path), "bytes": path.stat().st_size})
    manifest = {
        "format": "core-state-robustness-manifest/20261001.1",
        "source_inputs": sources,
        "artifacts": artifacts,
        "scientific_approved": False,
        "human_review_performed": False,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def run(args: argparse.Namespace) -> dict[str, Any]:
    root = Path(args.comparison_root).resolve(strict=True)
    if not root.is_dir():
        raise ValueError("comparison root is not a directory")
    output = Path(args.output)
    if output.exists():
        raise FileExistsError("output must be a new directory")
    manifest_path = safe_file(root, "manifest.json")
    comparison_path = safe_file(root, "comparison.json")
    manifest = read_json(manifest_path)
    comparison_artifacts = verify_comparison_manifest(root, manifest)
    comparison = read_json(comparison_path)

    evidence_path = Path(args.protein_evidence).resolve(strict=True)
    if sha256(evidence_path) != args.protein_evidence_sha256.lower():
        raise ValueError("protein evidence SHA-256 mismatch")
    evidence = read_json(evidence_path)
    if evidence.get("format") != "protein-hydrogen-evidence/1.1":
        raise ValueError("unsupported protein evidence format")
    verify_source_binding(comparison, evidence, root)
    protein_atoms, protein_hydrogens = load_protein(evidence)

    rows, sdf_records = [], []
    for analog_id, analog in sorted(find_analogs(comparison).items()):
        poses = analog.get("poses", [])
        if [pose.get("pose_index") for pose in poses] != list(range(5)):
            raise ValueError(f"{analog_id} must retain all five Vina-order poses")
        for pose in poses:
            report = pose.get("report", pose)
            states = report.get("states", [])
            requested = ANALOG_STATES[analog_id]
            if len(states) != requested or {s.get("state_index") for s in states} != set(range(requested)):
                raise ValueError(f"{analog_id} pose {pose['pose_index']} has missing states")
            expected_maps = set(report.get("expected_heavy_atom_maps", []))
            rmsd = report.get("shared_frame_core_rmsd_A")
            if (isinstance(rmsd, bool) or not isinstance(rmsd, (int, float))
                    or not math.isfinite(rmsd)):
                raise ValueError("finite numeric shared-frame core RMSD is required")
            for state in sorted(states, key=lambda row: row["state_index"]):
                artifacts = state.get("artifacts", {})
                before_relative = artifacts.get("before_sdf", "")
                after_relative = artifacts.get("after_sdf", "")
                for relative in (before_relative, after_relative):
                    if not isinstance(relative, str) or relative not in comparison_artifacts:
                        raise ValueError(f"state SDF is not registered in comparison manifest: {relative!r}")
                before_path = safe_file(root, before_relative)
                after_path = safe_file(root, after_relative)
                if sha256(before_path) != comparison_artifacts[before_relative]:
                    raise ValueError(f"state SDF hash mismatch on use: {before_relative}")
                if sha256(after_path) != comparison_artifacts[after_relative]:
                    raise ValueError(f"state SDF hash mismatch on use: {after_relative}")
                before, after = load_sdf(before_path), load_sdf(after_path)
                validate_state(before, after, state, expected_maps)
                sdf_records.append(after_path.read_bytes())
                status, error = "computed_diagnostic", None
                try:
                    profile = chemical_states.interaction_profile(
                        after, protein_atoms, protein_hydrogens=protein_hydrogens
                    )
                    measurement = analog_state_comparison._measurement(
                        after, protein_atoms, 17, "ASN", "1464", "OD1", profile
                    )
                    hydrophobic = analog_state_comparison._hydrophobic_summary(profile)
                except Exception as exc:
                    profile, measurement, hydrophobic = {}, None, []
                    status, error = "failed_review_required", f"{type(exc).__name__}: {exc}"
                hbond = bool(measurement and measurement.get("directional_hbond_observed"))
                residues = hydrophobic_residues(hydrophobic)
                rows.append({
                    "analog_id": analog_id,
                    "pose_index": pose["pose_index"],
                    "pose_order_note": "pose 0 is prior Vina order; not scientific best",
                    "state_index": state["state_index"],
                    "state_id": state.get("state_id"),
                    "status": status,
                    "failure": error,
                    "required_contact_evidence": classify_row(hbond, residues, float(rmsd)),
                    "map17_ASN1464_OD1": measurement,
                    "hydrophobic_contacts_VAL1408_PHE1409_ILE1470": hydrophobic,
                    "observed_required_hydrophobic_residues": sorted(residues),
                    "interaction_profile": profile,
                    "auxiliary_core_RMSD_A": float(rmsd),
                    "auxiliary_core_RMSD_within_1A": float(rmsd) <= 1.0,
                    "heavy_coordinates_exactly_unchanged": True,
                    "protein_hydrogen_orientation_requires_review": True,
                    "terminal_OXT": "pending_review",
                })
    if len(rows) != 40:
        raise ValueError(f"expected 40 retained rows, found {len(rows)}")

    analog_summaries = []
    for analog_id, count in ANALOG_STATES.items():
        analog_rows = [row for row in rows if row["analog_id"] == analog_id]
        supported = supported_pose_indices(analog_rows, count)
        analog_summaries.append({
            "analog_id": analog_id,
            "pose_count": 5,
            "requested_state_count_per_pose": count,
            "retained_row_count": len(analog_rows),
            "pose0_policy": "prior Vina order; explicitly not a scientific-best pose",
            "supported_pose_indices": supported,
            "all_states_robust_at_same_pose": bool(supported),
        })
    pka_rows = evidence.get("actual_pka_entries", [])
    summary = {
        "format": FORMAT,
        "rows": rows,
        "row_count": len(rows),
        "analogs": analog_summaries,
        "required_contact_definition": "directional map17-ASN1464-OD1 H-bond AND contacts to VAL1408, PHE1409, ILE1470",
        "auxiliary_RMSD_policy": "numeric diagnostic only; it cannot create or override contact evidence",
        "protein_hydrogens": "PDB2PQR registered; donor tests available; orientation_requires_review",
        "parent_map17_note": "Ligand-donor map17 to protein acceptor can be observed independently of protein hydrogen orientation.",
        "terminal_OXT": "pending_review",
        "microstate_gate": "pending",
        "state_selection_performed": False,
        "tautomer_selection_performed": False,
        "population_prediction_performed": False,
        "ligand_pKa_prediction_performed": False,
        "protein_actual_pKa_entry_count": len(pka_rows),
        "protein_actual_pKa_entries": pka_rows,
        "protein_registered_artifact_manifest": evidence.get("artifact_manifest", []),
        "pH_7_4_note": "Preparation scenario only; not proof of an assay-state population.",
        "diagnostic_only": True,
        "scientific_approved": False,
        "human_review_performed": False,
    }
    sources = {
        "comparison.json": sha256(comparison_path),
        "comparison_manifest.json": sha256(manifest_path),
        "protein_evidence.json": sha256(evidence_path),
        "protein_source_cif": _source_hash(comparison),
    }
    write_outputs(output, summary, sdf_records, sources)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison-root", required=True)
    parser.add_argument("--protein-evidence", required=True)
    parser.add_argument("--protein-evidence-sha256", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
