"""Source-supported exploratory native-9D12 N3 attachment panel.

This experiment is deliberately separate from strict scientific gates. It preserves
UNKNOWN state at parent map 3, makes no potency, approval, route, or binding claim,
and retains every requested candidate and docking outcome without replacement.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem

from packages.science.analog_filters import cheap_filter
from packages.science.mapped_stereo import same_mapped_tetrahedral_stereo
from packages.science.native_parent_docking import (
    dock_native,
    load_native_context,
    verify_source_pair,
)

PARENT_ID = "SMARCA2-9D12-A1A1P"
SOURCE_DOI = "10.1021/acs.jmedchem.4c01903"
SEED = 23
EXHAUSTIVENESS = 16
FRAGMENTS = (
    ("N3-LH-01", "CCO"),
    ("N3-LH-02", "CCCO"),
    ("N3-LH-03", "CCCCO"),
    ("N3-LH-04", "CCCCCO"),
    ("N3-LH-05", "CCCCCCO"),
    ("N3-LH-06", "CCN"),
    ("N3-LH-07", "CCCN"),
    ("N3-LH-08", "CCCCN"),
    ("N3-LH-09", "CCCCCN"),
    ("N3-LH-10", "COCCO"),
    ("N3-LH-11", "CCOCCO"),
    ("N3-LH-12", "CCOCCN"),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode("utf-8")


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(_json_bytes(value))


def _mapped_atom(mol: Chem.Mol, atom_map: int) -> Chem.Atom:
    atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomMapNum() == atom_map]
    if len(atoms) != 1:
        raise ValueError(f"Expected exactly one atom at map {atom_map}")
    return atoms[0]


def _mapped_smiles(mol: Chem.Mol) -> str:
    copy_mol = Chem.Mol(mol)
    Chem.AssignStereochemistry(copy_mol, cleanIt=True, force=True)
    return Chem.MolToSmiles(copy_mol, canonical=True, isomericSmiles=True)


def _graph_hash(mol: Chem.Mol) -> str:
    copy_mol = Chem.Mol(mol)
    for atom in copy_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    value = Chem.MolToSmiles(copy_mol, canonical=True, isomericSmiles=True)
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _source_hash(value: Any) -> str:
    return hashlib.sha256(_json_bytes(value)).hexdigest()


def _walk(value: Any, path: str = "data") -> Iterable[tuple[str, Any]]:
    yield path, value
    if isinstance(value, dict):
        for key, item in value.items():
            yield from _walk(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from _walk(item, f"{path}[{index}]")


def _valid_pair_locators(value: Any) -> bool:
    return (
        isinstance(value, list)
        and len(value) == 2
        and all(
            isinstance(row, dict)
            and row.get("doi") == SOURCE_DOI
            and isinstance(row.get("filename"), str)
            and bool(row["filename"])
            and isinstance(row.get("compound_id_column"), str)
            and bool(row["compound_id_column"])
            for row in value
        )
    )


def find_source_linkage(selected_evidence: Any, source_data: Any) -> dict[str, Any]:
    """Join the one verified selected-parent SAR row to its exact raw-source pair."""
    selected_snapshot = _source_hash(selected_evidence)
    source_snapshot = _source_hash(source_data)
    if not isinstance(selected_evidence, dict):
        raise ValueError("Selected-parent evidence is malformed")
    matched = selected_evidence.get("matched_record")
    if (
        selected_evidence.get("status") != "exact_measured_parent_join"
        or not isinstance(matched, dict)
        or matched.get("compound_id") != "SMI-6080"
    ):
        raise ValueError("Selected parent is not the exact measured SMI-6080 record")
    selected_rows = selected_evidence.get("source_sar")
    selected_matches = [
        (index, row) for index, row in enumerate(selected_rows if isinstance(selected_rows, list) else [])
        if isinstance(row, dict)
        and row.get("pair_id") == "SMI-6080__SMD-6087"
        and row.get("requested_change_label") == "terminal_N_derivatization_to_PROTAC"
        and row.get("selected_parent_pair_side") == "left"
        and row.get("other_compound_id") == "SMD-6087"
        and row.get("selected_parent_affected_atom_maps") == [3]
        and row.get("source_site_consensus") is True
    ]
    if len(selected_matches) != 1:
        raise ValueError("Selected evidence must contain one exact SMI-6080/SMD-6087 N3 SAR row")
    selected_index, selected_row = selected_matches[0]
    if not _valid_pair_locators(selected_row.get("source_row_locators")):
        raise ValueError("Selected SAR pair has invalid primary-source row locators")

    raw_rows = source_data.get("sar_pair_evidence") if isinstance(source_data, dict) else None
    raw_matches = [
        row for row in raw_rows if isinstance(row, dict)
        and row.get("pair_id") == "SMI-6080__SMD-6087"
        and row.get("left_compound_id") == "SMI-6080"
        and row.get("right_compound_id") == "SMD-6087"
        and row.get("requested_change_label") == "terminal_N_derivatization_to_PROTAC"
    ] if isinstance(raw_rows, list) else []
    if len(raw_matches) != 1:
        raise ValueError("Raw source must contain one exact SMI-6080/SMD-6087 pair")
    raw_row = raw_matches[0]
    if not _valid_pair_locators(raw_row.get("source_row_locators")):
        raise ValueError("Raw SAR pair has invalid primary-source row locators")
    possibilities = raw_row.get("mapping_possibilities")
    raw_sites = [
        row.get("left_affected_parent_atom_maps")
        for row in possibilities if isinstance(row, dict)
    ] if isinstance(possibilities, list) else []
    if not raw_sites or any(value != [19] for value in raw_sites):
        raise ValueError("Raw source SMI-6080 site must have exact map-19 consensus")
    if (
        selected_row.get("source_affected_atom_maps") != [19]
        or selected_row.get("source_row_locators") != raw_row.get("source_row_locators")
    ):
        raise ValueError("Selected SAR row does not match the verified raw-source pair")
    if (
        _source_hash(selected_evidence) != selected_snapshot
        or _source_hash(source_data) != source_snapshot
    ):
        raise ValueError("Source evidence changed during read-only inspection")
    return {
        "doi": SOURCE_DOI,
        "source_full_pointer": f"result.selected_parent_measured_evidence.source_sar[{selected_index}]",
        "source_record": copy.deepcopy(selected_row),
        "raw_source_record": copy.deepcopy(raw_row),
        "locators": copy.deepcopy(selected_row["source_row_locators"]),
        "source_data_sha256": source_snapshot,
        "selected_evidence_sha256": selected_snapshot,
        "source_parent_atom_maps": [19],
        "selected_parent_atom_maps": [3],
        "interpretation": "The original source map 19 corresponds to selected-parent map 3. This supports exploratory terminal N derivatization only, not potency, approval, binding, or a synthetic route.",
    }


def _original_maps(parent: Chem.Mol) -> list[int]:
    maps = [atom.GetAtomMapNum() for atom in parent.GetAtoms() if atom.GetAtomicNum() > 1]
    if any(value <= 0 for value in maps) or len(maps) != len(set(maps)):
        raise ValueError("Parent heavy atoms must have unique positive maps")
    return sorted(maps)


def _validate_parent_n3(parent: Chem.Mol) -> None:
    atom = _mapped_atom(parent, 3)
    degree = sum(neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors())
    if (
        atom.GetAtomicNum() != 7
        or atom.GetFormalCharge() != 0
        or atom.GetIsAromatic()
        or degree != 2
        or atom.GetTotalNumHs() != 1
    ):
        raise ValueError("Parent map 3 must be one neutral, nonaromatic, degree-2 N-H")


def _verify_original_graph(parent: Chem.Mol, candidate: Chem.Mol) -> None:
    parent_maps = _original_maps(parent)
    candidate_by_map = {
        atom.GetAtomMapNum(): atom for atom in candidate.GetAtoms() if atom.GetAtomicNum() > 1
    }
    if not set(parent_maps) <= set(candidate_by_map):
        raise ValueError("Candidate removed an original heavy atom")
    for atom_map in parent_maps:
        old = _mapped_atom(parent, atom_map)
        new = candidate_by_map[atom_map]
        if (
            old.GetAtomicNum(), old.GetIsotope(), old.GetFormalCharge(), old.GetIsAromatic()
        ) != (
            new.GetAtomicNum(), new.GetIsotope(), new.GetFormalCharge(), new.GetIsAromatic()
        ):
            raise ValueError(f"Original atom chemistry changed at map {atom_map}")
        if not same_mapped_tetrahedral_stereo(old, new):
            raise ValueError(f"Original mapped stereochemistry changed at map {atom_map}")
    for bond in parent.GetBonds():
        left = bond.GetBeginAtom().GetAtomMapNum()
        right = bond.GetEndAtom().GetAtomMapNum()
        other = candidate.GetBondBetweenAtoms(
            candidate_by_map[left].GetIdx(), candidate_by_map[right].GetIdx()
        )
        if (
            other is None
            or other.GetBondType() != bond.GetBondType()
            or other.GetIsAromatic() != bond.GetIsAromatic()
            or int(other.GetStereo()) != int(bond.GetStereo())
        ):
            raise ValueError(f"Original bond changed between maps {left} and {right}")
    original_set = set(parent_maps)
    external = []
    for bond in candidate.GetBonds():
        left = bond.GetBeginAtom().GetAtomMapNum()
        right = bond.GetEndAtom().GetAtomMapNum()
        if (left in original_set) != (right in original_set):
            external.append((left, right, bond.GetBondType()))
    if len(external) != 1 or 3 not in external[0][:2] or external[0][2] != Chem.BondType.SINGLE:
        raise ValueError("Candidate must add exactly one N3-C pendant bond")


def build_n3_candidate(parent: Chem.Mol, fragment_id: str, fragment_smiles: str) -> dict[str, Any]:
    """Consume the parent N3-H and attach the fragment through its first carbon."""
    _validate_parent_n3(parent)
    fragment = Chem.MolFromSmiles(fragment_smiles)
    if fragment is None or fragment.GetNumAtoms() < 2:
        raise ValueError("Invalid terminal-handle fragment")
    first = fragment.GetAtomWithIdx(0)
    if first.GetAtomicNum() != 6:
        raise ValueError("Fragment must attach through its first carbon")

    editable = Chem.RWMol(Chem.Mol(parent))
    n3 = _mapped_atom(editable, 3)
    n3.SetNumExplicitHs(0)
    n3.SetNoImplicit(True)
    added_maps = []
    fragment_indices = []
    for offset, atom in enumerate(fragment.GetAtoms()):
        copied = Chem.Atom(atom)
        copied.SetAtomMapNum(5000 + offset)
        copied.SetChiralTag(atom.GetChiralTag())
        fragment_indices.append(editable.AddAtom(copied))
        added_maps.append(5000 + offset)
    for bond in fragment.GetBonds():
        editable.AddBond(
            fragment_indices[bond.GetBeginAtomIdx()],
            fragment_indices[bond.GetEndAtomIdx()],
            bond.GetBondType(),
        )
    editable.AddBond(n3.GetIdx(), fragment_indices[0], Chem.BondType.SINGLE)
    candidate = editable.GetMol()
    Chem.SanitizeMol(candidate)
    Chem.AssignStereochemistry(candidate, cleanIt=True, force=True)
    _verify_original_graph(parent, candidate)

    resulting_n3 = _mapped_atom(candidate, 3)
    if resulting_n3.GetFormalCharge() != 0 or resulting_n3.GetTotalNumHs() != 0:
        raise ValueError("Resulting map 3 must be a neutral tertiary nitrogen")
    endpoint = candidate.GetAtomWithIdx(fragment_indices[-1])
    if endpoint.GetAtomicNum() not in (7, 8) or endpoint.GetFormalCharge() != 0 or endpoint.GetTotalNumHs() < 1:
        raise ValueError("Fragment must introduce a terminal neutral N-H or O-H endpoint")

    original_maps = _original_maps(parent)
    mapped_smiles = _mapped_smiles(candidate)
    record = {
        "id": fragment_id,
        "candidate_id": "N3-" + _graph_hash(candidate)[:16],
        "parent_id": PARENT_ID,
        "parent_warhead": PARENT_ID,
        "rule_id": fragment_id,
        "transformation_class": "linker_handle_introduction",
        "hypothesis_class": "source-extrapolation",
        "fragment_smiles": fragment_smiles,
        "mapped_smiles": mapped_smiles,
        "removed_atom_maps": [],
        "added_atom_maps": added_maps,
        "attachment_site_atom_maps": [3],
        "modified_atom_maps": [3],
        "protected_atom_maps": [value for value in original_maps if value != 3],
        "unchanged_original_graph_maps": original_maps,
        "protected_graph_preserved": True,
        "site_status": [{"atom_map": 3, "state": "UNKNOWN", "scientific_pending": True}],
        "source_state_promotion": False,
        "strict_qualified_for_counts": False,
        "evidence_confidence": "low",
        "reactivity_risk": "unknown",
        "binding_claim": False,
        "potency_claim": False,
        "approval_claim": False,
        "synthetic_route_status": "not_assessed",
    }
    record["cheap_filter"] = cheap_filter(parent, record)
    return record


def enumerate_n3_candidates(parent: Chem.Mol) -> list[dict[str, Any]]:
    records = [build_n3_candidate(parent, fragment_id, smiles) for fragment_id, smiles in FRAGMENTS]
    if len(records) != 12 or len({row["candidate_id"] for row in records}) != 12:
        raise ValueError("The predeclared panel must contain exactly twelve distinct graphs")
    return records


def _manifest_files(root: Path) -> list[dict[str, Any]]:
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError("Symlinks are forbidden in experiment artifacts")
        if path.is_file() and path.name != "manifest.json":
            files.append({
                "path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            })
    return files


def verify_native_probe(native_probe: Path) -> dict[str, Any]:
    """Verify every native-probe artifact and reject traversal or symlink aliases."""
    root = Path(native_probe).resolve()
    manifest_path = root / "manifest.json"
    if not root.is_dir() or not manifest_path.is_file() or manifest_path.is_symlink():
        raise ValueError("Native probe manifest is missing or unsafe")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows = manifest.get("files")
    if not isinstance(rows, list):
        raise ValueError("Native probe manifest has no files list")
    declared = {}
    for row in rows:
        if not isinstance(row, dict) or set(row) != {"path", "bytes", "sha256"}:
            raise ValueError("Malformed native probe manifest row")
        relative = Path(row["path"])
        if relative.is_absolute() or ".." in relative.parts or relative.as_posix() == "manifest.json":
            raise ValueError("Native probe manifest path escapes its root")
        path = root / relative
        cursor = root
        for part in relative.parts:
            cursor = cursor / part
            if cursor.is_symlink():
                raise ValueError("Native probe artifact uses a symlink")
        resolved = path.resolve()
        if not resolved.is_relative_to(root) or not path.is_file():
            raise ValueError("Native probe artifact is outside or absent")
        if path.stat().st_size != row["bytes"] or _sha256(path) != row["sha256"]:
            raise ValueError("Native probe artifact hash or size mismatch")
        normalized = relative.as_posix()
        if normalized in declared:
            raise ValueError("Native probe manifest contains a duplicate path")
        declared[normalized] = row
    required = {"native-probe-result.json", "protocol.json"}
    if not required <= set(declared):
        raise ValueError("Native probe manifest must declare protocol and result files")
    actual = {row["path"] for row in _manifest_files(root)}
    if actual != set(declared):
        raise ValueError("Native probe manifest does not cover all and only output files")
    result_path = root / "native-probe-result.json"
    protocol_path = root / "protocol.json"
    result = json.loads(result_path.read_text(encoding="utf-8"))
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    return {
        "root": str(root),
        "manifest_sha256": _sha256(manifest_path),
        "result_sha256": _sha256(result_path),
        "protocol_sha256": _sha256(protocol_path),
        "result": result,
        "protocol": protocol,
    }


def _preservation_pass_count(run: dict[str, Any]) -> int:
    rows = run.get("results", {}).get("pose_preservation", {}).get("all_poses_diagnostics", [])
    return sum(
        1 for row in rows if isinstance(row, dict)
        and row.get("docking_pose_preserved") is True
        and row.get("status") == "pass"
    ) if isinstance(rows, list) else 0


def validate_parent_baseline(native_receipt: dict[str, Any]) -> dict[str, Any]:
    result = native_receipt["result"]
    runs = [row for row in result.get("runs", []) if row.get("predeclared_ligand") == "native_parent_Cl"]
    by_seed = {row.get("seed"): row for row in runs}
    if set(by_seed) != {23, 41, 61} or len(runs) != 3:
        raise ValueError("Native probe must retain exactly the three declared parent seeds")
    if any(row.get("status") != "completed_with_limits" for row in runs):
        raise ValueError("All three native parent docking runs must be complete")
    pass_counts = {seed: _preservation_pass_count(run) for seed, run in by_seed.items()}
    if pass_counts[23] < 1:
        raise ValueError("Predeclared parent seed 23 has no pose-preservation pass")
    failed_seeds = sorted(seed for seed, count in pass_counts.items() if count == 0)
    return {
        "parent_run_count": 3,
        "parent_complete_count": 3,
        "pass_counts_by_seed": {str(key): value for key, value in sorted(pass_counts.items())},
        "retained_failed_seeds": failed_seeds,
        "selected_baseline_seed": 23,
        "selection_policy": "seed 23 was predeclared; no post-hoc best-seed selection",
        "three_of_three_pass_claim": all(value > 0 for value in pass_counts.values()),
        "experimental_demo_only": True,
    }


def _remote_exclusions(protocol: dict[str, Any]) -> list[str]:
    context = protocol.get("native_context", {})
    rows = context.get("protein_selection", {}).get("explicit_remote_incomplete_residue_omissions", [])
    if not isinstance(rows, list):
        raise ValueError("Old native context has malformed residue omissions")
    result = []
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get("residue_id"), str):
            raise ValueError("Old native context has malformed residue omission record")
        result.append(row["residue_id"])
    return result


def _site_state(site_rows: Any, atom_map: int) -> str:
    rows = [row for row in site_rows if isinstance(row, dict) and row.get("atom_map") == atom_map]
    if len(rows) != 1:
        raise ValueError(f"Native site analysis must contain exactly one map {atom_map} row")
    return str(rows[0].get("state", "UNKNOWN")).upper()


def _write_sdf(path: Path, mol: Chem.Mol) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    stream = path.open("w", encoding="utf-8", newline="\n")
    writer = Chem.SDWriter(stream)
    writer.SetForceV3000(True)
    try:
        writer.write(mol)
    finally:
        writer.close()
        stream.close()


def _find_catalog_item(values: Any, item_id: str | None = None, e3_type: str | None = None) -> dict[str, Any]:
    candidates = []
    for _, value in _walk(values, "catalog"):
        if isinstance(value, dict):
            if item_id is not None and value.get("id") == item_id:
                candidates.append(value)
            elif e3_type is not None and value.get("e3_type") == e3_type and value.get("mapped_smiles"):
                candidates.append(value)
    if not candidates:
        raise ValueError(f"Catalog item unavailable: {item_id or e3_type}")
    return copy.deepcopy(sorted(candidates, key=lambda row: str(row.get("id")))[0])


def _assemble_branches(record: dict[str, Any], case_dir: Path) -> list[dict[str, Any]]:
    from packages.science.dual_e3 import assemble, attachment_options, catalog

    options = attachment_options(record)
    preferred = sorted(
        options,
        key=lambda row: (row.get("element") not in ("N", "O"), int(row.get("atom_map", 10**9))),
    )
    if not preferred:
        raise ValueError("Candidate has no newly introduced terminal N-H/O-H attachment option")
    chosen = preferred[0]
    if chosen.get("source") != "added_atom_maps":
        raise ValueError("Assembly attachment must use a newly introduced endpoint")
    data = catalog()
    linker = _find_catalog_item(data, item_id="alkyl_c6")
    if linker.get("id") != "alkyl_c6":
        raise ValueError("Catalog linker identity mismatch")
    branches = []
    for e3_type in ("CRBN", "VHL"):
        recruiter = _find_catalog_item(data, e3_type=e3_type)
        if recruiter.get("e3_type") != e3_type:
            raise ValueError("Catalog recruiter identity mismatch")
        assembled = assemble(
            record,
            recruiter,
            linker,
            orientation="forward",
            attachment_map=int(chosen["atom_map"]),
        )
        assembled["scientific_pending"] = True
        assembled["strict_qualified_for_counts"] = False
        assembled["binding_claim"] = False
        assembled["synthetic_route_status"] = "not_assessed"
        mol = Chem.MolFromSmiles(assembled["mapped_smiles"])
        if mol is None:
            raise ValueError("Assembled mapped SMILES cannot be parsed")
        sdf_path = case_dir / "assemblies" / f"{e3_type.lower()}-alkyl_c6.sdf"
        _write_sdf(sdf_path, mol)
        assembled["sdf"] = sdf_path.relative_to(case_dir).as_posix()
        branches.append(assembled)
    return branches


def _pose_metrics(result: dict[str, Any]) -> dict[str, Any]:
    results = result.get("results", {})
    poses = results.get("poses", [])
    diagnostics = results.get("pose_preservation", {}).get("all_poses_diagnostics", [])
    pose_scores = {}
    for index, row in enumerate(poses if isinstance(poses, list) else []):
        if isinstance(row, dict):
            score = row.get("score", {}).get("affinity")
            if isinstance(score, (int, float)) and not isinstance(score, bool) and math.isfinite(float(score)):
                pose_scores[index] = float(score)
    rows = []
    for index, row in enumerate(diagnostics if isinstance(diagnostics, list) else []):
        if not isinstance(row, dict):
            continue
        rmsd = row.get("core_RMSD_A_in_receptor_frame")
        retention = row.get("protected_atom_near_residue_preservation", {}).get("retention_fraction")
        score = row.get("caller_supplied_score")
        if not isinstance(score, (int, float)) or isinstance(score, bool) or not math.isfinite(float(score)):
            score = pose_scores.get(row.get("pose_index_zero_based", index))
        rows.append({
            "pose_index_zero_based": row.get("pose_index_zero_based", index),
            "score": float(score) if isinstance(score, (int, float)) and not isinstance(score, bool) and math.isfinite(float(score)) else None,
            "core_RMSD_A": float(rmsd) if isinstance(rmsd, (int, float)) and not isinstance(rmsd, bool) and math.isfinite(float(rmsd)) else None,
            "contact_retention": float(retention) if isinstance(retention, (int, float)) and not isinstance(retention, bool) and math.isfinite(float(retention)) else None,
            "docking_pose_preserved": row.get("docking_pose_preserved") is True and row.get("status") == "pass",
        })
    geometry_rows = [row for row in rows if row["core_RMSD_A"] is not None]
    best_geometry = min(geometry_rows, key=lambda row: row["core_RMSD_A"]) if geometry_rows else None
    scores = [row["score"] for row in rows if row["score"] is not None]
    if not scores:
        scores = list(pose_scores.values())
    return {
        "pose_rows": rows,
        "global_best_raw_score": min(scores) if scores else None,
        "best_geometry_pose_index_zero_based": best_geometry["pose_index_zero_based"] if best_geometry else None,
        "best_geometry_core_RMSD_A": best_geometry["core_RMSD_A"] if best_geometry else None,
        "best_geometry_contact_retention": best_geometry["contact_retention"] if best_geometry else None,
        "best_geometry_score": best_geometry["score"] if best_geometry else None,
        "pose_preservation_pass_count": _preservation_pass_count(result),
    }


def _report(cases: list[dict[str, Any]], baseline: dict[str, Any]) -> str:
    lines = [
        "# Native 9D12 N3 attachment panel",
        "",
        "Exploratory source-extrapolation only. No potency, binding, approval, route, or strict-gate claim is made.",
        "Map 3 remains UNKNOWN and scientific_pending; pose diagnostics do not promote it.",
        "",
        f"Parent seed 23 preservation passes: {baseline['pass_counts_by_seed']['23']}",
        f"Retained parent seeds with no preservation pass: {baseline['retained_failed_seeds']}",
        "",
        "| candidate | status | raw score | core RMSD A | contact retention | assemblies |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for case in cases:
        metrics = case.get("pose_metrics", {})
        lines.append(
            "| {id} | {status} | {score} | {rmsd} | {retention} | {assemblies} |".format(
                id=case["candidate_id"],
                status=case["status"],
                score=metrics.get("global_best_raw_score"),
                rmsd=metrics.get("best_geometry_core_RMSD_A"),
                retention=metrics.get("best_geometry_contact_retention"),
                assemblies=len(case.get("assemblies", [])),
            )
        )
    lines.extend([
        "",
        "All twelve predeclared cases are retained. Failed or filtered cases are not replaced.",
        "Docking scores are not affinity measurements, and assembled branches are graph hypotheses only.",
        "",
    ])
    return "\n".join(lines)


def run_native_attachment_panel(
    native_probe: Path,
    native_cif: Path,
    source_export: Path,
    output: Path,
    seed: int = SEED,
    exhaustiveness: int = EXHAUSTIVENESS,
) -> dict[str, Any]:
    if seed != 23:
        raise ValueError("This predeclared experiment requires exact seed 23")
    if exhaustiveness != 16:
        raise ValueError("This predeclared experiment requires exhaustiveness 16")
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")

    native_receipt = verify_native_probe(native_probe)
    baseline = validate_parent_baseline(native_receipt)
    verified, selected, source_object, medchem_ref = verify_source_pair(source_export, PARENT_ID)
    source_before = _source_hash(source_object)
    source_data = source_object.get("data")
    linkage = find_source_linkage(selected, source_data)
    if _source_hash(source_object) != source_before:
        raise ValueError("Verified source object was mutated")

    exclusions = _remote_exclusions(native_receipt["protocol"])
    parent, protein, context = load_native_context(
        PARENT_ID,
        native_cif,
        output / "receptorfiles",
        exclude_remote_incomplete_residues=exclusions,
    )
    old_hash = native_receipt["protocol"].get("native_context", {}).get("native_cif_sha256")
    if old_hash and context.get("native_cif_sha256") != old_hash:
        raise ValueError("Fresh native context differs from the verified baseline CIF")

    from packages.science.medchem_sar import link_selected_parent
    recomputed_selected = link_selected_parent(parent, PARENT_ID, source_object)
    recomputed_linkage = find_source_linkage(recomputed_selected, source_data)
    linkage_fields = (
        "pair_id", "requested_change_label", "selected_parent_pair_side",
        "other_compound_id", "source_affected_atom_maps",
        "selected_parent_affected_atom_maps", "source_site_consensus",
        "source_row_locators",
    )
    if any(
        recomputed_linkage["source_record"].get(key) != linkage["source_record"].get(key)
        for key in linkage_fields
    ):
        raise ValueError("Fresh native parent does not reproduce the verified selected SAR site pair")

    from packages.science.warhead_sites import analyze_sites
    native_sites = analyze_sites(parent, protein, {}, [])
    if _site_state(native_sites.get("atoms", []), 3) != "UNKNOWN":
        raise ValueError("Native map 3 must remain UNKNOWN")
    protected_maps = sorted(
        row["atom_map"] for row in native_sites.get("atoms", [])
        if isinstance(row, dict) and row.get("state") == "PROTECTED"
    )
    if 3 in protected_maps:
        raise ValueError("Docking PROTECTED mask must exclude UNKNOWN map 3")

    candidates = enumerate_n3_candidates(parent)
    protocol_cases = [{
        "id": row["id"],
        "candidate_id": row["candidate_id"],
        "fragment_smiles": row["fragment_smiles"],
        "source_pointer": linkage["source_full_pointer"],
        "family": "linker_handle_introduction",
        "seed": 23,
    } for row in candidates]
    protocol = {
        "format_version": "native-9D12-N3-attachment-panel-v1",
        "status": "predeclared_before_any_candidate_vina_run",
        "parent_id": PARENT_ID,
        "strict_qualified_for_counts": False,
        "requested_candidate_count": 12,
        "failure_replacement_policy": "none",
        "source_export_receipt": verified,
        "selected_parent_evidence": selected,
        "source_medchem_evidence": medchem_ref,
        "source_linkage": linkage,
        "source_object_sha256": source_before,
        "native_probe": {
            key: native_receipt[key] for key in ("manifest_sha256", "result_sha256", "protocol_sha256")
        },
        "parent_baseline": baseline,
        "native_context": context,
        "explicit_remote_exclusions_reused_exactly": exclusions,
        "site_analysis": native_sites,
        "site_state": {"atom_map": 3, "original": "UNKNOWN", "experiment": "UNKNOWN", "promotion": False},
        "masks": {
            "native_PROTECTED_docking_maps_excluding_N3": protected_maps,
            "unchanged_all_original_graph_maps_identity_only_not_docking_RMSD": _original_maps(parent),
        },
        "threshold_scope": "existing pose_preservation diagnostics; exploratory only",
        "docking": {"seed": 23, "exhaustiveness": 16, "num_modes": 5, "timeout_seconds": 180},
        "assembly": {"linker_id": "alkyl_c6", "e3_types": ["CRBN", "VHL"]},
        "cases": protocol_cases,
        "claims": {"binding": False, "potency": False, "approval": False, "synthetic_route": False},
    }
    _write_json(output / "protocol.json", protocol)

    receptor_pdbqt = Path(context["receptor_preparation"]["pdbqt"])
    outcomes = []
    for record in candidates:
        case_dir = output / "cases" / record["id"]
        case_dir.mkdir(parents=True, exist_ok=False)
        mol = Chem.MolFromSmiles(record["mapped_smiles"])
        if mol is None:
            raise ValueError("Generated candidate cannot be parsed")
        _write_sdf(case_dir / "candidate.sdf", mol)
        _write_json(case_dir / "candidate.json", record)
        outcome = {
            "id": record["id"],
            "candidate_id": record["candidate_id"],
            "mapped_smiles": record["mapped_smiles"],
            "candidate_sdf": "candidate.sdf",
            "cheap_filter": record["cheap_filter"],
            "strict_qualified_for_counts": False,
            "scientific_pending": True,
            "assemblies": [],
        }
        if record["cheap_filter"].get("valid") is not True:
            outcome["status"] = "filtered_not_replaced"
            outcome["docking"] = None
        else:
            ligand = Chem.Mol(mol)
            ligand.RemoveAllConformers()
            docking = dock_native(
                ligand,
                parent,
                protein,
                protected_maps,
                receptor_pdbqt,
                case_dir / "dock",
                seed=23,
                exhaustiveness=16,
                timeout=180,
            )
            _write_json(case_dir / "dock-receipt.json", docking)
            outcome["docking"] = docking
            outcome["pose_metrics"] = _pose_metrics(docking)
            if docking.get("status") == "completed_with_limits":
                outcome["status"] = "docked_with_limits"
                if outcome["pose_metrics"]["pose_preservation_pass_count"] >= 1:
                    try:
                        outcome["assemblies"] = _assemble_branches(record, case_dir)
                        outcome["status"] = "docked_and_assembled_hypotheses"
                    except Exception as exc:
                        outcome["status"] = "assembly_failed_not_replaced"
                        outcome["assembly_error"] = str(exc)
            else:
                outcome["status"] = "docking_failed_not_replaced"
        _write_json(case_dir / "outcome.json", outcome)
        outcomes.append(outcome)

    final = {
        "status": "completed_with_limits",
        "parent_id": PARENT_ID,
        "requested_count": 12,
        "retained_count": len(outcomes),
        "replacement_count": 0,
        "strict_qualified_for_counts": False,
        "source_state_map_3": "UNKNOWN",
        "scientific_pending": True,
        "parent_baseline": baseline,
        "cases": outcomes,
        "claims": {"binding": False, "potency": False, "approval": False, "synthetic_route": False},
    }
    _write_json(output / "panel.json", final)
    (output / "report.md").write_text(_report(outcomes, baseline), encoding="utf-8", newline="\n")

    graphs = []
    for row in outcomes:
        mol = Chem.MolFromSmiles(row["mapped_smiles"])
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        graphs.append((
            Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
            Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        ))
    native_runs = []
    for run in native_receipt["result"].get("runs", []):
        if not isinstance(run, dict):
            continue
        diagnostics = run.get("results", {}).get("pose_preservation", {}).get("all_poses_diagnostics", [])
        core_values = [
            float(row["core_RMSD_A_in_receptor_frame"])
            for row in diagnostics if isinstance(row, dict)
            and isinstance(row.get("core_RMSD_A_in_receptor_frame"), (int, float))
            and not isinstance(row.get("core_RMSD_A_in_receptor_frame"), bool)
        ] if isinstance(diagnostics, list) else []
        preservation = run.get("results", {}).get("pose_preservation", {})
        native_runs.append({
            "ligand": run.get("predeclared_ligand"),
            "seed": run.get("seed"),
            "status": run.get("status"),
            "pose_pass_count": _preservation_pass_count(run),
            "minimum_core_RMSD_A": min(core_values) if core_values else None,
            "core_atom_maps": diagnostics[0].get("core_atom_maps", []) if isinstance(diagnostics, list) and diagnostics and isinstance(diagnostics[0], dict) else [],
            "policy": preservation.get("policy"),
        })
    attempted = [row for row in outcomes if row.get("docking") is not None]
    completed = [row for row in attempted if row["docking"].get("status") == "completed_with_limits"]
    compact = {
        "format_version": "native-9D12-N3-attachment-compact-summary-v1",
        "family_1": "exploratory",
        "requested_count": 12,
        "retained_count": len(outcomes),
        "docking_attempt_count": len(attempted),
        "docking_completed_count": len(completed),
        "docking_failed_count": len(attempted) - len(completed),
        "filtered_count": sum(row.get("status") == "filtered_not_replaced" for row in outcomes),
        "unique_constitutional_map_free_nonstereo_graph_count": len({row[0] for row in graphs}),
        "unique_isomeric_graph_count": len({row[1] for row in graphs}),
        "pose_pass_case_count": sum(row.get("pose_metrics", {}).get("pose_preservation_pass_count", 0) > 0 for row in outcomes),
        "actual_e3_assembly_success": {
            "case_count": sum(bool(row.get("assemblies")) for row in outcomes),
            "branch_count": sum(len(row.get("assemblies", [])) for row in outcomes),
            "by_e3_type": {
                e3: sum(branch.get("e3_type") == e3 for row in outcomes for branch in row.get("assemblies", []))
                for e3 in ("CRBN", "VHL")
            },
        },
        "candidate_rows": [{
            "candidate_id": row["candidate_id"],
            "status": row["status"],
            "best_geometry_pose_index_zero_based": row.get("pose_metrics", {}).get("best_geometry_pose_index_zero_based"),
            "best_geometry_core_RMSD_A": row.get("pose_metrics", {}).get("best_geometry_core_RMSD_A"),
            "best_geometry_contact_retention": row.get("pose_metrics", {}).get("best_geometry_contact_retention"),
            "global_best_raw_score": row.get("pose_metrics", {}).get("global_best_raw_score"),
        } for row in outcomes],
        "native_parent_and_Br_runs": native_runs,
        "native_run_count": len(native_runs),
        "source_site": {"selected_parent_atom_map": 3, "state": "UNKNOWN", "raw_source_atom_map": 19},
        "comparison_policy": "No fixed 6HAZ numeric improvement is claimed because protocols and atom masks differ.",
        "strict_counts": {"family_1": 0, "total": 0},
    }
    _write_json(output / "compact-summary.json", compact)
    _write_json(output / "manifest.json", {
        "format_version": "native-9D12-N3-attachment-panel-v1",
        "hash_algorithm": "sha256",
        "files": _manifest_files(output),
        "manifest_self_hash_policy": "manifest.json is excluded",
    })
    return final


__all__ = [
    "FRAGMENTS",
    "build_n3_candidate",
    "enumerate_n3_candidates",
    "find_source_linkage",
    "run_native_attachment_panel",
    "validate_parent_baseline",
    "verify_native_probe",
]
