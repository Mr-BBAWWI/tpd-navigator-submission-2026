"""Exact-graph synthesis evidence and conservative disconnection review."""
from __future__ import annotations

import json
from typing import Any

from rdkit import Chem

VERSION = "synthesis-review/20261001.1"
_ROLE_KEYS = ("warhead_maps", "linker_maps", "recruiter_maps")
_FLAG_SMARTS = {"amide": "[NX3][CX3](=[OX1])", "imide": "[NX3]([CX3](=[OX1]))[CX3](=[OX1])", "glutarimide": "O=C1CCC(=O)N1", "hydroxyproline": "[C@H]1C[C@H](O)CN1", "labile_ester": "[CX3](=[OX1])[OX2][#6]", "labile_acetal": "[CX4]([OX2])([OX2])", "labile_disulfide": "[SX2][SX2]"}


def _mol(smiles: Any, name: str) -> Chem.Mol:
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError(f"{name} must be a nonempty SMILES string")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError(f"{name} is not a valid single molecular graph")
    return mol


def _canonical(mol: Chem.Mol, isomeric: bool = True) -> str:
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=isomeric)


def _maps(mol: Chem.Mol) -> dict[int, int]:
    result = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        atom_map = atom.GetAtomMapNum()
        if atom_map <= 0 or atom_map in result:
            raise ValueError("all candidate heavy atoms require unique positive atom maps")
        result[atom_map] = atom.GetIdx()
    return result


def _roles(candidate: dict[str, Any], graph_maps: set[int]) -> dict[str, list[int]]:
    raw = candidate.get("atom_roles")
    if not isinstance(raw, dict) or set(raw) != set(_ROLE_KEYS):
        raise ValueError("atom_roles must contain exactly warhead_maps, linker_maps, and recruiter_maps")
    result, seen = {}, set()
    for key in _ROLE_KEYS:
        values = raw[key]
        if not isinstance(values, list) or any(isinstance(x, bool) or not isinstance(x, int) or x <= 0 for x in values):
            raise ValueError(f"{key} must be a positive integer list")
        if len(values) != len(set(values)) or seen & set(values):
            raise ValueError("atom role lists must be disjoint and internally unique")
        seen.update(values)
        result[key] = sorted(values)
    if seen != graph_maps:
        raise ValueError("atom roles must form a complete partition of the mapped heavy graph")
    return result


def _attachment_bonds(candidate: dict[str, Any], mol: Chem.Mol, mapping: dict[int, int], roles: dict[str, list[int]]):
    metadata = candidate.get("attachment_metadata")
    if not isinstance(metadata, dict):
        raise ValueError("attachment_metadata is required")
    role_sets = {key: set(value) for key, value in roles.items()}
    bonds = []
    for name in ("warhead_linker_bond", "recruiter_linker_bond"):
        record = metadata.get(name)
        if not isinstance(record, dict) or record.get("bond_type") != "SINGLE":
            raise ValueError(f"{name} must state a SINGLE bond")
        a, b = record.get("attachment_atom_map"), record.get("partner_atom_map")
        if a not in mapping or b not in mapping or a == b:
            raise ValueError(f"{name} references absent or identical atom maps")
        bond = mol.GetBondBetweenAtoms(mapping[a], mapping[b])
        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
            raise ValueError(f"{name} does not match an actual SINGLE bond")
        outer = role_sets["warhead_maps"] if name.startswith("warhead") else role_sets["recruiter_maps"]
        if not ((a in outer and b in role_sets["linker_maps"]) or (b in outer and a in role_sets["linker_maps"])):
            raise ValueError(f"{name} does not connect the declared role partitions")
        bonds.append((bond.GetIdx(), min(a, b), max(a, b), name))
    if bonds[0][0] == bonds[1][0]:
        raise ValueError("attachment records must identify distinct bonds")
    return bonds


def _cut_report(mol: Chem.Mol, bonds: list[tuple[int, int, int, str]]) -> dict[str, Any]:
    fragments = Chem.GetMolFrags(Chem.FragmentOnBonds(mol, [x[0] for x in bonds], addDummies=False), asMols=True, sanitizeFrags=True)
    if len(fragments) != 3:
        raise ValueError("cutting the stated bonds must produce exactly three fragments")
    endpoint_maps = {value for bond in bonds for value in bond[1:3]}
    records = []
    for fragment in fragments:
        maps = sorted(atom.GetAtomMapNum() for atom in fragment.GetAtoms() if atom.GetAtomicNum() > 1)
        records.append({"mapped_smiles": Chem.MolToSmiles(fragment, canonical=True, isomericSmiles=True), "canonical_smiles": _canonical(fragment), "atom_maps": maps, "end_chemistry": [{"atom_map": atom.GetAtomMapNum(), "element": atom.GetSymbol(), "formal_charge": atom.GetFormalCharge(), "heavy_degree_after_cut": atom.GetDegree(), "total_hydrogens_after_cut": int(atom.GetTotalNumHs(includeNeighbors=True))} for atom in fragment.GetAtoms() if atom.GetAtomicNum() > 1 and atom.GetAtomMapNum() in endpoint_maps]})
    records.sort(key=lambda item: item["atom_maps"])
    return {"status": "computed_graph_disconnection", "cut_bonds": [{"name": name, "atom_maps": [a, b], "bond_type": "SINGLE"} for _, a, b, name in bonds], "fragments": records, "interpretation": "Graph cut and reaction hypothesis only; no recipe is generated."}


def _flags(mol: Chem.Mol) -> list[dict[str, Any]]:
    result = []
    for name, smarts in _FLAG_SMARTS.items():
        query = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(query) if query is not None else ()
        if matches:
            result.append({"flag": name, "match_count": len(matches), "atom_maps": [sorted(mol.GetAtomWithIdx(i).GetAtomMapNum() for i in match) for match in matches], "status": "structural_review_flag"})
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            result.append({"flag": "specified_tetrahedral_stereochemistry", "atom_maps": [[atom.GetAtomMapNum()]], "status": "preserve_and_review"})
    return result


def _constraint_review(mol: Chem.Mol, constraints: Any) -> dict[str, Any]:
    if constraints is None:
        return {"status": "unknown_no_expert_constraints", "blocked": False, "matches": []}
    if not isinstance(constraints, dict):
        raise TypeError("constraints must be a dictionary")
    blocked = constraints.get("blocked_smarts", [])
    if not isinstance(blocked, list):
        raise ValueError("constraints.blocked_smarts must be a list")
    matches = []
    for position, entry in enumerate(blocked):
        if isinstance(entry, str):
            identifier, smarts, reason = f"constraint-{position}", entry, None
        elif isinstance(entry, dict):
            identifier, smarts, reason = entry.get("id", f"constraint-{position}"), entry.get("smarts"), entry.get("reason")
        else:
            raise TypeError("blocked_smarts entries must be strings or dictionaries")
        query = Chem.MolFromSmarts(smarts) if isinstance(smarts, str) else None
        if query is None:
            raise ValueError(f"invalid blocked SMARTS: {smarts}")
        found = mol.GetSubstructMatches(query)
        if found:
            matches.append({"id": identifier, "smarts": smarts, "reason": reason, "match_count": len(found)})
    return {"status": "blocked_by_explicit_expert_smarts" if matches else "unknown_no_blocking_match", "blocked": bool(matches), "matches": matches}


def _meaningful_string(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip()) and value.strip().lower() not in {"unknown", "n/a", "none", "not reported", "paper", "source"}


def _source_valid(source: Any, locator: Any) -> bool:
    if _meaningful_string(source) and _meaningful_string(locator):
        return True
    if not isinstance(source, dict):
        return False
    doi = source.get("doi")
    title = source.get("title")
    return _meaningful_string(doi) and doi.lower().startswith("10.") and _meaningful_string(title) and _meaningful_string(locator)


def _step_valid(step: Any) -> bool:
    if not isinstance(step, dict):
        return False
    reagents = step.get("reagents")
    return isinstance(reagents, list) and bool(reagents) and all(_meaningful_string(x) for x in reagents) and all(_meaningful_string(step.get(field)) for field in ("conditions", "purification", "characterization"))


def _route_records(records: Any, candidate_id: str, canonical: str) -> dict[str, Any]:
    if records is None:
        records = []
    if not isinstance(records, list):
        raise TypeError("exact_route_records must be a list")
    reports, complete = [], 0
    for position, record in enumerate(records):
        if not isinstance(record, dict):
            reports.append({"record_index": position, "status": "invalid_record", "failure_type": "TypeError", "failure": "route record must be a dictionary", "failure_preserved": True, "documentation_complete": False})
            continue
        report = {"record_index": position, "candidate_id": record.get("candidate_id"), "source_record": record, "documentation_complete": False}
        try:
            exact_id = record.get("candidate_id") == candidate_id
            route_graph = _mol(record.get("canonical_smiles"), "route canonical_smiles")
            graph_matches = _canonical(route_graph) == canonical
            steps = record.get("steps")
            step_complete = isinstance(steps, list) and bool(steps) and all(_step_valid(step) for step in steps)
            source_complete = _source_valid(record.get("source"), record.get("locator"))
            fields = {"candidate_id": exact_id, "canonical_smiles": graph_matches, "source_and_locator": source_complete, "steps": step_complete}
            is_complete = all(fields.values())
            complete += int(is_complete)
            report.update({"status": "complete_exact_evidence" if is_complete else "invalid_or_foreign_route", "field_presence": fields, "exact_graph": graph_matches, "documentation_complete": is_complete, "failure_preserved": not is_complete})
        except Exception as exc:
            report.update({"status": "invalid_record", "failure_type": type(exc).__name__, "failure": str(exc), "failure_preserved": True})
        reports.append(report)
    return {"records": reports, "complete_exact_evidence_count": complete, "documented_route_for_graph": complete > 0, "synthesized": False, "approved": False, "C01_C02_yield_policy": "Yields are never transferred to another graph.", "artifact_resolution": "Evidence artifacts must be resolved by the server-side document service, never by a browser."}


def assess_synthesis(candidate, *, exact_route_records=None, constraints=None):
    """Validate one exact candidate graph and preserve route-record failures."""
    if not isinstance(candidate, dict):
        raise TypeError("candidate must be a dictionary")
    candidate_id = candidate.get("candidate_id")
    if not isinstance(candidate_id, str) or not candidate_id:
        raise ValueError("candidate_id is required")
    mapped = _mol(candidate.get("mapped_smiles"), "candidate mapped_smiles")
    mapping = _maps(mapped)
    actual_canonical = _canonical(mapped)
    stated = _mol(candidate.get("canonical_smiles"), "candidate canonical_smiles")
    if _canonical(stated) != actual_canonical:
        raise ValueError("candidate canonical_smiles does not match the full mapped graph")
    roles = _roles(candidate, set(mapping))
    bonds = _attachment_bonds(candidate, mapped, mapping, roles)
    expert = _constraint_review(mapped, constraints)
    routes = _route_records(exact_route_records, candidate_id, actual_canonical)
    result = {"format": VERSION, "candidate_id": candidate_id, "graph_validation": {"valid": True, "canonical_smiles": actual_canonical, "mapped_smiles": Chem.MolToSmiles(mapped, canonical=True, isomericSmiles=True), "heavy_atom_count": len(mapping), "atom_roles": roles, "attachment_metadata_matches_graph": True}, "retrosynthetic_cuts": _cut_report(mapped, bonds), "structural_flags": _flags(mapped), "chemical_compatibility": {"status": "reaction_hypotheses_only", "computed": True, "demonstrated_synthesis": False, "reason": "Structural alerts and graph cuts do not establish chemoselectivity, conditions, route success, or yield."}, "expert_constraints": expert, "route_evidence": routes, "documented_route_for_graph": routes["documented_route_for_graph"], "synthesized": False, "scientific_approved": False, "overall_status": "blocked" if expert["blocked"] else ("documented_exact_route_pending_expert_review" if routes["documented_route_for_graph"] else "route_unknown_pending_expert_review"), "limitations": ["No recipe, condition, purification, characterization, or yield is invented.", "Complete documentation does not prove independent synthesis or grant approval."]}
    json.dumps(result, ensure_ascii=False, allow_nan=False)
    return result


__all__ = ["VERSION", "assess_synthesis"]
