"""Mapped, coordinate-bound interaction diagnostics without approval authority."""
from __future__ import annotations

import hashlib
import json
import math
import re
from typing import Any

import numpy as np
from rdkit import Chem

from .chemical_states import before_after_report, enumerate_microstates, interaction_profile, optimize_ligand_hydrogens

VERSION = "interaction-review/20260930.2"
_ALLOWED_KINDS = {"directional_hbond", "salt_bridge", "aromatic", "hydrophobic", "proximity"}
_AROMATIC_RING_ATOMS = {
    "PHE": frozenset({"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
    "TYR": frozenset({"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
    "TRP": frozenset({"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}),
    "HIS": frozenset({"CG", "ND1", "CD2", "CE1", "NE2"}),
}
_PROFILE_KINDS = {
    "directional_hbond": {"directional_hbond"},
    "salt_bridge": {"ionic_contact"},
    "aromatic": {"aromatic_contact_geometry"},
    "hydrophobic": {"hydrophobic_contact"},
    "proximity": {"directional_hbond", "ionic_contact", "aromatic_contact_geometry", "hydrophobic_contact", "possible_hbond_contact"},
}
_LIGAND_TOKEN = re.compile(r"^L:([1-9][0-9]*)$")
_LIGAND_RING_TOKEN = re.compile(r"^L:ring:(L:[1-9][0-9]*(?:,L:[1-9][0-9]*)*)$")


def _finite(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _hash(value: Any) -> str:
    text = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str, allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _validate_mol(mol: Chem.Mol, name: str) -> tuple[dict[int, int], np.ndarray]:
    if not isinstance(mol, Chem.Mol):
        raise TypeError(f"{name} must be an RDKit Mol")
    if mol.GetNumConformers() != 1:
        raise ValueError(f"{name} requires exactly one conformer")
    conformer = mol.GetConformer()
    if not conformer.Is3D():
        raise ValueError(f"{name} requires an actual 3D conformer")
    xyz = np.asarray(conformer.GetPositions(), dtype=float)
    if xyz.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(xyz).all():
        raise ValueError(f"{name} has missing or nonfinite coordinates")
    mapping: dict[int, int] = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        atom_map = atom.GetAtomMapNum()
        if atom_map <= 0 or atom_map in mapping:
            raise ValueError(f"{name} has missing or duplicate heavy-atom maps")
        mapping[atom_map] = atom.GetIdx()
    if not mapping:
        raise ValueError(f"{name} has no mapped heavy atoms")
    return mapping, xyz


def _protein_id(record: dict[str, Any]) -> str:
    required = ("label_asym_id", "auth_seq_id", "label_comp_id", "label_atom_id")
    if any(key not in record for key in required):
        raise ValueError("protein atom record lacks an exact chain/residue/component/atom identifier")
    return ":".join((str(record["label_asym_id"]), str(record["auth_seq_id"]), str(record["label_comp_id"]).upper(), str(record["label_atom_id"]).upper()))


def _validate_requirements(requirements: Any, protein_ids: set[str]) -> list[dict[str, Any]]:
    if requirements is None:
        return []
    if not isinstance(requirements, list):
        raise TypeError("requirements must be a list")
    result, identifiers = [], set()
    for position, raw in enumerate(requirements):
        if not isinstance(raw, dict):
            raise TypeError(f"requirements[{position}] must be a dictionary")
        missing = {"id", "kind", "ligand_maps", "protein_atom_ids", "required"} - set(raw)
        if missing:
            raise ValueError(f"requirements[{position}] missing fields: {sorted(missing)}")
        identifier, kind = raw["id"], raw["kind"]
        if not isinstance(identifier, str) or not identifier or identifier in identifiers:
            raise ValueError("requirement ids must be unique nonempty strings")
        if kind not in _ALLOWED_KINDS or not isinstance(raw["required"], bool):
            raise ValueError("invalid requirement kind or required flag")
        maps = raw["ligand_maps"]
        if not isinstance(maps, list) or not maps or any(isinstance(x, bool) or not isinstance(x, int) or x <= 0 for x in maps):
            raise ValueError("ligand_maps must be a nonempty positive integer list")
        protein = raw["protein_atom_ids"]
        if not isinstance(protein, list) or not protein or any(not isinstance(x, str) or not x for x in protein):
            raise ValueError("protein_atom_ids must be a nonempty string list")
        if len(maps) != len(set(maps)) or len(protein) != len(set(protein)):
            raise ValueError("requirement participants must be unique")
        if not set(protein) <= protein_ids:
            raise ValueError("protein_atom_ids contain absent exact protein atoms")
        identifiers.add(identifier)
        result.append({"id": identifier, "kind": kind, "ligand_maps": sorted(maps), "protein_atom_ids": sorted(protein), "required": raw["required"]})
    return result


def _parse_participant(value: Any) -> dict[str, Any] | None:
    if not isinstance(value, str):
        return None
    ligand = _LIGAND_TOKEN.fullmatch(value)
    if ligand:
        return {"type": "ligand_atom", "atom_map": int(ligand.group(1))}
    ligand_ring = _LIGAND_RING_TOKEN.fullmatch(value)
    if ligand_ring:
        maps = frozenset(int(token[2:]) for token in ligand_ring.group(1).split(","))
        return {"type": "ligand_ring", "atom_maps": maps} if maps else None
    fields = value.split(":")
    if len(fields) == 4 and all(fields):
        return {"type": "protein_atom", "id": f"{fields[0]}:{fields[1]}:{fields[2].upper()}:{fields[3].upper()}"}
    if len(fields) == 3 and all(fields):
        return {"type": "protein_residue", "chain": fields[0], "residue": fields[1], "component": fields[2].upper()}
    return None


def _protein_participant_matches(parsed: dict[str, Any], requested: str) -> bool:
    chain, residue, component, atom = requested.split(":", 3)
    if parsed["type"] == "protein_atom":
        return parsed["id"] == requested
    return (
        parsed["type"] == "protein_residue"
        and parsed["chain"] == chain
        and parsed["residue"] == residue
        and parsed["component"] == component
        and atom in _AROMATIC_RING_ATOMS.get(component, frozenset())
    )


def _participant_match(item: dict[str, Any], requirement: dict[str, Any]) -> bool:
    if item.get("kind") not in _PROFILE_KINDS[requirement["kind"]]:
        return False
    raw = item.get("participants", [])
    if not isinstance(raw, list):
        return False
    parsed = [token for token in (_parse_participant(value) for value in raw) if token is not None]
    required_maps = frozenset(requirement["ligand_maps"])
    ligand_hit = any(
        token["type"] == "ligand_atom" and token["atom_map"] in required_maps
        or token["type"] == "ligand_ring" and token["atom_maps"] == required_maps
        for token in parsed
    )
    protein_hit = any(
        _protein_participant_matches(token, requested)
        for token in parsed if token["type"] in {"protein_atom", "protein_residue"}
        for requested in requirement["protein_atom_ids"]
    )
    return ligand_hit and protein_hit


def _matches(profile: dict[str, Any], requirement: dict[str, Any]) -> list[dict[str, Any]]:
    interactions = profile.get("interactions", [])
    return [item for item in interactions if isinstance(item, dict) and _participant_match(item, requirement)] if isinstance(interactions, list) else []


def _pending_direction(profile: dict[str, Any], requirement: dict[str, Any], protein_hydrogens: Any) -> bool:
    if requirement["kind"] != "directional_hbond" or protein_hydrogens:
        return False
    possible = {**requirement, "kind": "proximity"}
    return any(isinstance(item, dict) and item.get("kind") == "possible_hbond_contact" and _participant_match(item, possible) for item in profile.get("interactions", []))


def _transfer_heavy_coordinates(state: Chem.Mol, source: Chem.Mol) -> Chem.Mol:
    source_maps, source_xyz = _validate_mol(source, "state coordinate source")
    state_maps: dict[int, int] = {}
    for atom in state.GetAtoms():
        if atom.GetAtomicNum() > 1:
            atom_map = atom.GetAtomMapNum()
            if atom_map <= 0 or atom_map in state_maps:
                raise ValueError("state has missing or duplicate heavy-atom maps")
            state_maps[atom_map] = atom.GetIdx()
    if set(state_maps) != set(source_maps):
        raise ValueError("state heavy maps differ from coordinate source")
    state.RemoveAllConformers()
    conformer = Chem.Conformer(state.GetNumAtoms())
    conformer.Set3D(True)
    center = np.mean(source_xyz[list(source_maps.values())], axis=0)
    for atom in state.GetAtoms():
        point = source_xyz[source_maps[atom.GetAtomMapNum()]] if atom.GetAtomicNum() > 1 else center
        conformer.SetAtomPosition(atom.GetIdx(), tuple(float(x) for x in point))
    state.AddConformer(conformer, assignId=True)
    return state


def _profile_for(mol: Chem.Mol, protein_atoms: Any, protein_hydrogens: Any, seed: int):
    structures, receipt = optimize_ligand_hydrogens(mol, seed=seed)
    profile = interaction_profile(structures["hydrogens_optimized"], protein_atoms, protein_hydrogens=protein_hydrogens)
    return profile, receipt


def _receipt_review(receipt: Any) -> bool:
    return isinstance(receipt, dict) and receipt.get("requires_review") is True


def assess_interactions(parent_mol, pose_mol, protein_atoms, *, protein_hydrogens=None, requirements=None, pH=7.4, seed=23, max_states=8):
    """Assess exact mapped interactions in the supplied fixed receptor frame."""
    parent_maps, parent_xyz = _validate_mol(parent_mol, "parent_mol")
    pose_maps, pose_xyz = _validate_mol(pose_mol, "pose_mol")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    if isinstance(max_states, bool) or not isinstance(max_states, int) or not 1 <= max_states <= 16:
        raise ValueError("max_states must be an integer from 1 to 16")
    context_ph = _finite(pH, "pH")
    if not isinstance(protein_atoms, (list, tuple)) or not protein_atoms:
        raise ValueError("protein_atoms must be nonempty")
    protein_ids = set()
    for record in protein_atoms:
        if not isinstance(record, dict):
            raise TypeError("protein atom records must be dictionaries")
        identifier = _protein_id(record)
        if identifier in protein_ids:
            raise ValueError(f"duplicate protein atom id: {identifier}")
        protein_ids.add(identifier)
    normalized = _validate_requirements(requirements, protein_ids)
    parent_profile, parent_receipt = _profile_for(parent_mol, protein_atoms, protein_hydrogens, seed)
    pose_profile, pose_receipt = _profile_for(pose_mol, protein_atoms, protein_hydrogens, seed)
    requirement_reports, failures = [], []
    unresolved_required_protein_donor = False
    for requirement in normalized:
        required_maps = set(requirement["ligand_maps"])
        in_parent, in_pose = required_maps <= set(parent_maps), required_maps <= set(pose_maps)
        before = _matches(parent_profile, requirement) if in_parent else []
        after = _matches(pose_profile, requirement) if in_pose else []
        pending_h = in_pose and not after and _pending_direction(pose_profile, requirement, protein_hydrogens)
        review = _receipt_review(parent_receipt) or _receipt_review(pose_receipt) or any(item.get("requires_review") is True for item in after)
        if not in_parent:
            status, reason = "pending", "baseline_definition_required"
        elif not before:
            status, reason = "pending", "baseline_definition_required"
        elif not in_pose:
            status, reason = "lost", "required ligand map is absent from the pose"
        elif after and not review:
            status, reason = "preserved", "observed in parent and pose with exact participants"
        elif after:
            status, reason = "pending", "optimization or interaction receipt requires review"
        elif pending_h:
            status, reason = "pending", "required protein-donor hydrogen direction is unresolved"
        else:
            status, reason = "lost", "parent interaction is not observed in the pose"
        passes = status == "preserved"
        unresolved_required_protein_donor |= bool(requirement["required"] and pending_h)
        if requirement["required"] and not passes:
            failures.append({"requirement_id": requirement["id"], "status": status, "reason": reason, "preserved_in_failure_ledger": True})
        requirement_reports.append({**requirement, "status": status, "computed_pass": passes, "parent_maps_present": in_parent, "pose_maps_present": in_pose, "before_matches": before, "after_matches": after, "needs_expert": status == "pending" or review, "reason": reason})
    alternatives = []
    fixed_hash = _hash({"pose_maps": sorted(pose_maps), "pose_xyz": pose_xyz.tolist(), "protein_atoms": protein_atoms, "protein_hydrogens": protein_hydrogens})
    source_hash = _hash(Chem.MolToMolBlock(pose_mol))
    for index, state in enumerate(enumerate_microstates(Chem.Mol(pose_mol), max_states=max_states, pH=context_ph)):
        entry = {"index": index, "origin": state.GetProp("chemical_state_origin") if state.HasProp("chemical_state_origin") else "unknown", "mapped_smiles": Chem.MolToSmiles(state, canonical=True, isomericSmiles=True), "formal_charge": sum(atom.GetFormalCharge() for atom in state.GetAtoms()), "pKa": "unknown", "population": "unknown", "pH_context": context_ph, "population_prediction_performed": False, "source_structure_sha256": source_hash, "fixed_coordinates_sha256": fixed_hash}
        try:
            transferred = _transfer_heavy_coordinates(Chem.Mol(state), pose_mol)
            profile, receipt = _profile_for(transferred, protein_atoms, protein_hydrogens, seed)
            reports = []
            for requirement in normalized:
                observed = bool(_matches(profile, requirement))
                pending_h = not observed and _pending_direction(profile, requirement, protein_hydrogens)
                reports.append({"id": requirement["id"], "observed": observed, "unambiguous": observed and not _receipt_review(receipt), "pending_missing_protein_hydrogen": pending_h})
            entry.update({"status": "computed", "hydrogen_receipt": receipt, "hydrogen_receipt_sha256": _hash(receipt), "interaction_profile": profile, "requirements": reports})
        except Exception as exc:
            entry.update({"status": "failed", "failure_type": type(exc).__name__, "failure": str(exc), "failure_preserved": True})
        alternatives.append(entry)
    if not normalized:
        policy_status = "computed_diagnostic_pending_policy"
    elif failures:
        policy_status = "requirements_not_satisfied"
    else:
        policy_status = "computed_requirements_satisfied_pending_scientific_review"
    result = {"format": VERSION, "coordinate_frame": "caller-bound fixed receptor frame; no alignment applied", "mapped_heavy_atoms": {"parent": sorted(parent_maps), "pose": sorted(pose_maps), "common": sorted(set(parent_maps) & set(pose_maps))}, "profiles": {"parent": parent_profile, "pose": pose_profile}, "before_after": before_after_report(parent_profile, pose_profile), "requirements": requirement_reports, "state_alternatives": alternatives, "hydrogen_receipts": {"parent": parent_receipt, "pose": pose_receipt}, "failures": failures, "policy_status": policy_status, "protein_hydrogen_status": "required_protein_donor_pending" if unresolved_required_protein_donor else ("provided" if protein_hydrogens else "not_required_or_not_implicated"), "pH_context": context_ph, "pKa_prediction_performed": False, "population_prediction_performed": False, "scientific_approved": False, "limitations": ["Heavy-atom proximity is not a hydrogen bond without donor typing and explicit direction.", "Computed contacts do not establish affinity, efficacy, selectivity, or experimental preservation."]}
    json.dumps(result, ensure_ascii=False, allow_nan=False)
    return result


__all__ = ["VERSION", "assess_interactions"]
