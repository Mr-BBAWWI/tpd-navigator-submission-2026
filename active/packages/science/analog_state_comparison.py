"""Analog-specific protonation-state diagnostics in an unchanged docking frame."""
from __future__ import annotations

import hashlib
import io
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem

from .chemical_states import _change_proton, interaction_profile, optimize_ligand_hydrogens

VERSION = "analog-state-comparison/20261001.3"
CORE_MAPS = (1, 2, 3, 4, 5, 6, 7, 9, 10, 13, 14, 15, 16, 17, 18, 20)
_PARENT_MAPS = frozenset(range(1, 21))
_EXPECTED = {
    "W-c2afc5e73c1a": (
        "[cH:1]1[c:2](-[c:4]2[cH:5][c:15]([N:18]3[CH2:9][CH2:12]"
        "[N:19]([CH2:5000][OH:5001])[CH2:8][CH2:11]3)[c:6]([NH2:17])"
        "[n:16][n:7]2)[c:3]([OH:20])[cH:10][cH:13][cH:14]1"
    ),
    "W-4c0a639c0a41": (
        "[cH:1]1[c:2](-[c:4]2[cH:5][c:15]([N:18]3[CH2:9][CH2:12]"
        "[N:19]([CH2:5000][CH2:5001][CH2:5002][OH:5003])[CH2:8][CH2:11]3)"
        "[c:6]([NH2:17])[n:16][n:7]2)[c:3]([OH:20])[cH:10][cH:13][cH:14]1"
    ),
    "W-80f8f4a11b5d": (
        "[cH:1]1[c:2](-[c:4]2[cH:5][c:15]([N:18]3[CH2:9][CH2:12]"
        "[N:19]([CH2:5000][NH2:5001])[CH2:8][CH2:11]3)[c:6]([NH2:17])"
        "[n:16][n:7]2)[c:3]([OH:20])[cH:10][cH:13][cH:14]1"
    ),
}
_STATE_SITES = {
    "W-c2afc5e73c1a": (19,),
    "W-4c0a639c0a41": (19,),
    "W-80f8f4a11b5d": (19, 5001),
}


def _sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha_text(value: str) -> str:
    return _sha_bytes(value.encode("utf-8"))


def _heavy_maps(mol: Chem.Mol, name: str) -> dict[int, int]:
    if not isinstance(mol, Chem.Mol):
        raise TypeError(f"{name} must be an RDKit Mol")
    result: dict[int, int] = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        atom_map = atom.GetAtomMapNum()
        if atom_map <= 0 or atom_map in result:
            raise ValueError(f"{name} requires unique positive heavy-atom maps")
        result[atom_map] = atom.GetIdx()
    if not result:
        raise ValueError(f"{name} has no mapped heavy atoms")
    return result


def _xyz(mol: Chem.Mol, name: str) -> np.ndarray:
    if mol.GetNumConformers() != 1:
        raise ValueError(f"{name} requires exactly one coordinate conformer")
    result = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    if result.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(result).all():
        raise ValueError(f"{name} coordinates are missing or nonfinite")
    return result


def _strip_hydrogens(mol: Chem.Mol, name: str) -> Chem.Mol:
    _xyz(mol, name)
    try:
        result = Chem.RemoveHs(Chem.Mol(mol), sanitize=True)
        Chem.SanitizeMol(result)
        Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    except Exception as exc:
        raise ValueError(f"{name} cannot be safely normalized by removing hydrogens") from exc
    _heavy_maps(result, name)
    _xyz(result, name)
    return result


def _atom_label(atom: Chem.Atom) -> tuple[int, int, int, int]:
    return (
        atom.GetAtomicNum(), atom.GetIsotope(), int(atom.GetChiralTag()),
        atom.GetFormalCharge(),
    )


def _bond_label(bond: Chem.Bond) -> tuple[str, bool, int, int]:
    return (
        str(bond.GetBondType()), bool(bond.GetIsAromatic()),
        int(bond.GetStereo()), int(bond.GetBondDir()),
    )


def _graph(mol: Chem.Mol, maps: set[int] | frozenset[int] | None = None) -> tuple[dict, dict]:
    indices = _heavy_maps(mol, "molecule")
    selected = set(indices) if maps is None else set(maps)
    if not selected <= set(indices):
        raise ValueError("requested graph maps are absent")
    atoms = {m: _atom_label(mol.GetAtomWithIdx(indices[m])) for m in selected}
    bonds = {}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if a in selected and b in selected:
            bonds[tuple(sorted((a, b)))] = _bond_label(bond)
    return atoms, bonds


def _mapped_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _mapless_nonstereogenic_nitrogen(mol: Chem.Mol, atom_map: int) -> bool:
    """Confirm by mapless canonical symmetry that a mapped N is not stereogenic."""
    mapping = _heavy_maps(mol, "stereo symmetry check")
    if atom_map not in mapping:
        return False
    copy = Chem.Mol(mol)
    index = mapping[atom_map]
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
        ranks = list(Chem.CanonicalRankAtoms(
            copy, breakTies=False, includeChirality=False, includeIsotopes=True
        ))
    except Exception as exc:
        raise ValueError("mapless chemical-symmetry stereo check failed") from exc
    atom = copy.GetAtomWithIdx(index)
    if atom.GetAtomicNum() != 7 or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
        return False
    ligand_classes = [ranks[neighbor.GetIdx()] for neighbor in atom.GetNeighbors()]
    ligand_classes.extend(
        [-(1 + int(atom.GetIsotope()))]
        * int(atom.GetTotalNumHs(includeNeighbors=False))
    )
    return len(ligand_classes) != len(set(ligand_classes))


def _verify_heavy_graph_and_stereo(
        source: Chem.Mol, restored: Chem.Mol, context: str) -> None:
    """Verify heavy chemistry while allowing only map-induced pseudo-N stereo."""
    source_maps = _heavy_maps(source, f"{context} source")
    restored_maps = _heavy_maps(restored, f"{context} restored")
    if set(source_maps) != set(restored_maps):
        raise ValueError(f"{context} changed heavy-atom map identity")

    for atom_map, source_index in source_maps.items():
        left = source.GetAtomWithIdx(source_index)
        right = restored.GetAtomWithIdx(restored_maps[atom_map])
        if (
            left.GetAtomicNum(), left.GetIsotope(), left.GetFormalCharge()
        ) != (
            right.GetAtomicNum(), right.GetIsotope(), right.GetFormalCharge()
        ):
            raise ValueError(f"{context} changed heavy-atom labels at map {atom_map}")

        left_tag = left.GetChiralTag()
        right_tag = right.GetChiralTag()
        if left_tag == right_tag:
            continue
        if left_tag != Chem.ChiralType.CHI_UNSPECIFIED:
            raise ValueError(f"{context} changed explicit source stereo at map {atom_map}")
        if (
            atom_map not in {19, 5001}
            or left.GetAtomicNum() != 7
            or not _mapless_nonstereogenic_nitrogen(source, atom_map)
            or not _mapless_nonstereogenic_nitrogen(restored, atom_map)
        ):
            raise ValueError(f"{context} introduced unsupported stereo at map {atom_map}")

    source_bonds = _graph(source)[1]
    restored_bonds = _graph(restored)[1]
    if source_bonds != restored_bonds:
        raise ValueError(f"{context} changed heavy-atom bonds or explicit bond stereo")


def _atom_h_count(atom: Chem.Atom) -> int:
    attached = sum(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors())
    return attached + int(atom.GetTotalNumHs(includeNeighbors=False))


def _nitrogen_metrics(mol: Chem.Mol) -> dict[str, dict[str, int]]:
    mol.UpdatePropertyCache(strict=False)
    result = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetAtomMapNum() > 0:
            result[str(atom.GetAtomMapNum())] = {
                "formal_charge": int(atom.GetFormalCharge()),
                "hydrogen_count": _atom_h_count(atom),
            }
    return dict(sorted(result.items(), key=lambda item: int(item[0])))


def _state(source: Chem.Mol, protonated: tuple[int, ...]) -> Chem.Mol:
    result = Chem.Mol(source)
    for atom_map in protonated:
        mapping = _heavy_maps(result, "state under construction")
        changed = _change_proton(result, mapping[atom_map], 1)
        if changed is None:
            raise ValueError(f"protonation failed at atom map {atom_map}")
        result = changed
    try:
        Chem.SanitizeMol(result)
        Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    except Exception as exc:
        raise ValueError("constructed protonation state is not sanitizable") from exc
    return result


def _sdf_text(mol: Chem.Mol) -> str:
    """Serialize as V3000 and verify maps, graph, and coordinates on load-back."""
    original_maps = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
    original_xyz = _xyz(mol, "SDF export source")
    stream = io.StringIO()
    writer = Chem.SDWriter(stream)
    writer.SetForceV3000(True)
    writer.write(mol)
    writer.flush()
    writer.close()
    text = stream.getvalue()
    loaded = list(Chem.ForwardSDMolSupplier(
        io.BytesIO(text.encode("utf-8")), removeHs=False,
        sanitize=True, strictParsing=True,
    ))
    if len(loaded) != 1 or loaded[0] is None:
        raise ValueError("V3000 SDF export cannot be loaded back")
    restored = loaded[0]
    if [atom.GetAtomMapNum() for atom in restored.GetAtoms()] != original_maps:
        raise ValueError("V3000 SDF export lost or changed atom maps")
    _verify_heavy_graph_and_stereo(mol, restored, "V3000 SDF export")
    restored_xyz = _xyz(restored, "SDF export load-back")
    if restored_xyz.shape != original_xyz.shape or not np.allclose(
            restored_xyz, original_xyz, rtol=0.0, atol=5.1e-5):
        raise ValueError("V3000 SDF export changed source coordinates")
    return text


def _safe_artifact(directory: Path, name: str, text: str) -> str:
    path = directory / name
    if path.exists():
        raise FileExistsError(path)
    path.write_text(text, encoding="utf-8", newline="\n")
    return path.relative_to(directory).as_posix()


def _protein_id(record: dict[str, Any]) -> str:
    return ":".join((
        str(record.get("label_asym_id", "")), str(record.get("auth_seq_id", "")),
        str(record.get("label_comp_id", "")).upper(),
        str(record.get("label_atom_id", "")).upper(),
    ))


def _measurement(mol: Chem.Mol, protein_atoms: Any, ligand_map: int,
                 residue: str, sequence: str, atom_name: str,
                 profile: dict) -> dict | None:
    mapping = _heavy_maps(mol, "optimized state")
    if ligand_map not in mapping:
        return None
    xyz = _xyz(mol, "optimized state")[mapping[ligand_map]]
    matches = [
        row for row in protein_atoms
        if str(row.get("label_comp_id", "")).upper() == residue
        and str(row.get("auth_seq_id", "")) == sequence
        and str(row.get("label_atom_id", "")).upper() == atom_name
    ]
    if not matches:
        return None
    row = matches[0]
    distance = float(np.linalg.norm(xyz - np.asarray(row["xyz"], dtype=float)))
    pid = _protein_id(row)
    interactions = [
        item for item in profile.get("interactions", [])
        if f"L:{ligand_map}" in item.get("participants", [])
        and pid in item.get("participants", [])
    ]
    return {
        "ligand_atom_map": ligand_map,
        "protein_atom_id": pid,
        "heavy_distance_A": round(distance, 4),
        "directional_hbond_observed": any(
            item.get("kind") == "directional_hbond" for item in interactions
        ),
        "matching_profile_interactions": interactions,
    }


def _hydrophobic_summary(profile: dict) -> list[dict]:
    wanted = {("VAL", "1408"), ("PHE", "1409"), ("ILE", "1470")}
    result = []
    for item in profile.get("interactions", []):
        if item.get("kind") != "hydrophobic_contact":
            continue
        proteins = []
        for participant in item.get("participants", []):
            fields = participant.split(":") if isinstance(participant, str) else []
            if len(fields) == 4 and (fields[2].upper(), fields[1]) in wanted:
                proteins.append(participant)
        if proteins:
            result.append({
                "protein_atom_ids": proteins,
                "ligand_participants": [
                    value for value in item.get("participants", [])
                    if isinstance(value, str) and value.startswith("L:")
                ],
                "distance_A": item.get("distance_A"),
                "profile_interaction_id": item.get("id"),
            })
    return result


def _json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("nonfinite report value")
        return value
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    raise TypeError(f"non-JSON report value: {type(value).__name__}")


def build_state_set(analog_id, mapped_pose, parent_mol, protein_atoms, output_dir=None):
    """Build the fixed, analog-specific 2/2/4 state set and real diagnostics.

    No tautomer, pKa, population, alignment, state selection, or approval operation
    is performed. Optimization failures remain represented as review rows.
    """
    if analog_id not in _EXPECTED:
        raise ValueError(f"unknown analog id: {analog_id}")
    if not isinstance(protein_atoms, (list, tuple)) or not protein_atoms:
        raise ValueError("protein_atoms must be a nonempty sequence")

    pose = _strip_hydrogens(mapped_pose, "mapped_pose")
    parent = _strip_hydrogens(parent_mol, "parent_mol")
    expected = Chem.MolFromSmiles(_EXPECTED[analog_id])
    if expected is None:
        raise RuntimeError("internal expected analog graph is invalid")
    Chem.SanitizeMol(expected)

    pose_maps = _heavy_maps(pose, "mapped_pose")
    parent_maps = _heavy_maps(parent, "parent_mol")
    expected_maps = _heavy_maps(expected, "expected analog")
    if set(pose_maps) != set(expected_maps):
        raise ValueError("mapped_pose atom-map identity does not match selected analog")
    if set(parent_maps) != _PARENT_MAPS:
        raise ValueError("parent_mol must contain exactly mapped parent heavy atoms 1..20")
    if _graph(pose) != _graph(expected):
        raise ValueError("mapped_pose graph does not match selected analog")
    if _graph(pose, _PARENT_MAPS) != _graph(expected, _PARENT_MAPS):
        raise ValueError("mapped_pose parent-core topology mismatch")
    expected_parent_atoms, expected_parent_bonds = _graph(expected, _PARENT_MAPS)
    parent_atoms, parent_bonds = _graph(parent, _PARENT_MAPS)
    # Parent map 19 differs only in its source ionization charge/H state; topology
    # and atom identity are compared without formal charge.
    for atom_map in _PARENT_MAPS:
        if parent_atoms[atom_map][:3] != expected_parent_atoms[atom_map][:3]:
            raise ValueError(f"parent atom identity mismatch at map {atom_map}")
    if parent_bonds != expected_parent_bonds:
        raise ValueError("parent-core bond topology or stereo labels mismatch")

    pose_xyz = _xyz(pose, "mapped_pose")
    parent_xyz = _xyz(parent, "parent_mol")
    differences = np.asarray([
        pose_xyz[pose_maps[m]] - parent_xyz[parent_maps[m]] for m in CORE_MAPS
    ])
    shared_rmsd = float(np.sqrt(np.mean(np.sum(differences * differences, axis=1))))
    source_smiles = _mapped_smiles(pose)
    source_metrics = _nitrogen_metrics(pose)

    destination = None
    if output_dir is not None:
        destination = Path(output_dir)
        if destination.exists():
            raise FileExistsError(destination)
        destination.mkdir(parents=True)

    sites = _STATE_SITES[analog_id]
    combinations = [()]
    for site in sites:
        combinations += [tuple((*existing, site)) for existing in list(combinations)]
    rows = []
    source_graph = _graph(pose)
    for index, protonated in enumerate(combinations):
        identifier = "source" if not protonated else "protonated-" + "-".join(map(str, protonated))
        entry: dict[str, Any] = {
            "state_index": index,
            "state_id": identifier,
            "protonated_atom_maps": list(protonated),
            "status": "pending_computation",
            "selected": False,
            "formal_charge": None,
            "mapped_nitrogens": {},
            "mapped_smiles": None,
            "mapped_smiles_sha256": None,
            "artifacts": {},
        }
        state = _state(pose, protonated)
        if _graph(state)[1] != source_graph[1]:
            raise ValueError("protonation changed heavy-atom bond labels")
        state_atoms = _graph(state)[0]
        source_atoms = source_graph[0]
        for atom_map in source_atoms:
            if atom_map not in protonated and state_atoms[atom_map] != source_atoms[atom_map]:
                raise ValueError(f"unrequested atom change at map {atom_map}")
        metrics = _nitrogen_metrics(state)
        smiles = _mapped_smiles(state)
        entry.update({
            "formal_charge": int(Chem.GetFormalCharge(state)),
            "mapped_nitrogens": metrics,
            "mapped_smiles": smiles,
            "mapped_smiles_sha256": _sha_text(smiles),
            "other_mapped_nitrogens_unchanged": all(
                metrics[key] == source_metrics[key]
                for key in source_metrics if int(key) not in protonated
            ),
        })
        if destination is not None:
            entry["artifacts"]["before_sdf"] = _safe_artifact(
                destination, f"state-{index:02d}-before.sdf", _sdf_text(state)
            )
        try:
            structures, receipt = optimize_ligand_hydrogens(state, seed=23 + index)
        except Exception as exc:
            entry.update({
                "status": "failed_review_required",
                "failure_stage": "hydrogen_optimization",
                "failure_type": type(exc).__name__,
                "failure": str(exc),
                "failure_retained": True,
            })
            rows.append(entry)
            continue
        if not isinstance(structures, dict) or "hydrogens_optimized" not in structures:
            raise ValueError("hydrogen optimization returned no optimized structure")
        if not isinstance(receipt, dict):
            raise ValueError("hydrogen optimization receipt must be an object")
        optimized = structures["hydrogens_optimized"]
        _verify_heavy_graph_and_stereo(
            state, optimized, "hydrogen optimization"
        )
        optimized_maps = _heavy_maps(optimized, "optimized state")
        if set(optimized_maps) != set(pose_maps):
            raise ValueError("hydrogen optimization changed heavy-atom map identity")
        optimized_xyz = _xyz(optimized, "optimized state")
        drift = max(float(np.linalg.norm(
            optimized_xyz[optimized_maps[m]] - pose_xyz[pose_maps[m]]
        )) for m in pose_maps)
        try:
            profile = interaction_profile(
                optimized, protein_atoms, protein_hydrogens=None
            )
        except Exception as exc:
            entry.update({
                "status": "failed_review_required",
                "failure_stage": "interaction_profile",
                "failure_type": type(exc).__name__,
                "failure": str(exc),
                "failure_retained": True,
                "hydrogen_optimization": receipt,
                "heavy_coordinate_max_displacement_A": drift,
                "heavy_coordinates_exactly_preserved": drift == 0.0,
            })
            rows.append(entry)
            continue
        if not isinstance(profile, dict):
            raise ValueError("interaction profile must be an object")
        if destination is not None:
            entry["artifacts"]["after_sdf"] = _safe_artifact(
                destination, f"state-{index:02d}-after.sdf", _sdf_text(optimized)
            )
        entry.update({
            "status": "optimization_review" if receipt.get("requires_review") else "computed",
            "hydrogen_optimization": receipt,
            "heavy_coordinate_max_displacement_A": drift,
            "heavy_coordinates_exactly_preserved": drift == 0.0,
            "interaction_profile": profile,
            "map17_ASN1464_OD1": _measurement(
                optimized, protein_atoms, 17, "ASN", "1464", "OD1", profile
            ),
            "map20_TYR1421_OH": _measurement(
                optimized, protein_atoms, 20, "TYR", "1421", "OH", profile
            ),
            "hydrophobic_contacts_VAL1408_PHE1409_ILE1470": _hydrophobic_summary(profile),
            "protein_hydrogen_status": "absent",
            "directional_limitation": (
                "Protein hydrogens were not supplied; protein-donor direction cannot be "
                "confirmed and is never auto-approved. Ligand-donor direction uses only "
                "the explicitly optimized ligand hydrogens."
            ),
        })
        rows.append(entry)

    report = {
        "format": VERSION,
        "analog_id": analog_id,
        "source_mapped_smiles": source_smiles,
        "source_mapped_smiles_sha256": _sha_text(source_smiles),
        "source_mapped_nitrogens": source_metrics,
        "stereo_provenance": {
            "authoritative_source": "input molecular graph",
            "explicit_source_stereo_conserved": True,
            "coordinates_are_measured_configuration": False,
            "allowed_roundtrip_inference": (
                "Only source-unspecified nitrogen maps 19/5001 may acquire a "
                "3D-derived tag, and only when mapless canonical chemical symmetry "
                "confirms equivalent ligands and RDKit stereo cleaning removes it."
            ),
        },
        "expected_heavy_atom_maps": sorted(expected_maps),
        "source_like_aromatic_tautomer_retained": True,
        "tautomer_enumeration_performed": False,
        "state_count": len(rows),
        "states": rows,
        "coordinate_frame": "shared input frame; no alignment or realignment applied",
        "core_maps": list(CORE_MAPS),
        "shared_frame_core_rmsd_A": shared_rmsd,
        "shared_frame_core_rmsd_limit_A": 1.0,
        "shared_frame_core_rmsd_within_limit": shared_rmsd <= 1.0,
        "selected_state": None,
        "microstate_gate": "pending",
        "core_final": "pending",
        "scientific_approved": False,
        "population_prediction_performed": False,
        "pKa_prediction_performed": False,
        "limitations": [
            "States are retained diagnostic structures, not predicted populations.",
            "No state is selected and no interaction establishes affinity, efficacy, or degradation.",
            "Missing protein hydrogens limit directional interpretation and never imply approval.",
        ],
    }
    safe = _json_safe(report)
    json.dumps(safe, ensure_ascii=False, allow_nan=False)
    return safe


__all__ = ["VERSION", "CORE_MAPS", "build_state_set"]
