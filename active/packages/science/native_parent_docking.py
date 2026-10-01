"""Exploratory docking in a deposited parent/receptor coordinate frame.

This module is intentionally separate from the fixed 6HAZ docking adapter. It does
not mutate strict-gate behavior, align native poses into another receptor frame, or
turn exploratory docking diagnostics into scientific approval.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem

from packages.science.design_docking import (
    _VINA_RESULT,
    _box,
    _kill,
    _restore_maps,
    _short_windows_path,
    _validate_inputs,
)

_METHOD = "native-deposited-parent-vina-1.2.7-experiment-v1"
_CANONICAL_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_CANONICAL_HEAVY_ATOMS = {
    "ALA": {"N", "CA", "C", "O", "CB"},
    "ARG": {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "ASN": {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"},
    "ASP": {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"},
    "CYS": {"N", "CA", "C", "O", "CB", "SG"},
    "GLN": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"},
    "GLU": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"},
    "GLY": {"N", "CA", "C", "O"},
    "HIS": {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "ILE": {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"},
    "LEU": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"},
    "LYS": {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"},
    "MET": {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"},
    "PHE": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": {"N", "CA", "C", "O", "CB", "CG", "CD"},
    "SER": {"N", "CA", "C", "O", "CB", "OG"},
    "THR": {"N", "CA", "C", "O", "CB", "OG1", "CG2"},
    "TRP": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "VAL": {"N", "CA", "C", "O", "CB", "CG1", "CG2"},
}
_REMOTE_INCOMPLETE_RADIUS_A = 20.0


class _ValidatedResidueDeletions(tuple):
    """Opaque result of explicit remote-incomplete-residue validation."""


def _residue_id(row: dict[str, Any]) -> str:
    chain = str(row.get("label_asym_id", ""))
    seq = str(row.get("auth_seq_id", ""))
    insertion = row.get("pdbx_PDB_ins_code")
    suffix = "" if insertion in (None, False, "", ".", "?") else str(insertion)
    return f"{chain}:{seq}{suffix}"


def validate_remote_incomplete_residue_exclusions(
    protein: list[dict[str, Any]],
    ligand_atoms: list[dict[str, Any]],
    requested_ids: Any,
) -> tuple[_ValidatedResidueDeletions, list[dict[str, Any]]]:
    """Validate only explicitly named, incomplete canonical residues beyond 20 A."""
    if requested_ids is None:
        requested_ids = []
    if not isinstance(requested_ids, (list, tuple)):
        raise ValueError("Remote incomplete residue exclusions must be a list")
    if any(not isinstance(value, str) or not value for value in requested_ids):
        raise ValueError("Remote incomplete residue IDs must be nonempty strings")
    if len(requested_ids) != len(set(requested_ids)):
        raise ValueError("Duplicate remote incomplete residue IDs are forbidden")

    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in protein:
        grouped.setdefault(_residue_id(row), []).append(row)
    ligand_xyz = []
    for row in ligand_atoms:
        xyz = row.get("xyz_A", row.get("xyz"))
        if not isinstance(xyz, list) or len(xyz) != 3:
            raise ValueError("Native ligand heavy-atom coordinates are malformed")
        ligand_xyz.append(tuple(float(value) for value in xyz))
    if not ligand_xyz:
        raise ValueError("Native ligand heavy-atom selection is empty")

    records = []
    for residue_id in requested_ids:
        rows = grouped.get(residue_id)
        if not rows:
            raise ValueError(f"Requested residue exclusion does not exist: {residue_id}")
        names = {str(row.get("label_atom_id", "")) for row in rows}
        residues = {str(row.get("label_comp_id", "")).upper() for row in rows}
        if len(residues) != 1 or next(iter(residues)) not in _CANONICAL_HEAVY_ATOMS:
            raise ValueError(f"Requested residue is not one canonical residue: {residue_id}")
        residue = next(iter(residues))
        expected = _CANONICAL_HEAVY_ATOMS[residue]
        unexpected = sorted(names - expected)
        missing = sorted(expected - names)
        if unexpected:
            raise ValueError(f"Requested residue has unexpected heavy atoms: {residue_id}")
        if not missing:
            raise ValueError(f"Complete residue cannot be excluded: {residue_id}")
        minimum = min(
            math.dist(tuple(float(value) for value in row["xyz"]), ligand_point)
            for row in rows
            for ligand_point in ligand_xyz
        )
        if minimum <= _REMOTE_INCOMPLETE_RADIUS_A:
            raise ValueError(
                f"Incomplete residue {residue_id} is within the 20 A experimental radius"
            )
        records.append(
            {
                "residue_id": residue_id,
                "residue_name": residue,
                "actual_heavy_atom_names": sorted(names),
                "canonical_expected_heavy_atom_names": sorted(expected),
                "missing_heavy_atom_names": missing,
                "minimum_native_ligand_heavy_atom_distance_A": minimum,
                "developer_experimental_radius_A": _REMOTE_INCOMPLETE_RADIUS_A,
                "scientific_gate": False,
                "explicitly_requested": True,
            }
        )
    return _ValidatedResidueDeletions(requested_ids), records


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")


def _write_json(path: Path, value: Any) -> None:
    path.write_bytes(_json_bytes(value))


def _version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _safe_number(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): _safe_number(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe_number(item) for item in value]
    if hasattr(value, "item"):
        return _safe_number(value.item())
    return value


def _canonical_mapped(mol: Chem.Mol) -> str:
    copy = Chem.Mol(mol)
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)


def _mapping_rows(metadata: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows = metadata.get("atom_mapping")
    if not isinstance(rows, list) or not rows:
        raise ValueError("Native parent metadata has no atom_mapping")
    result: dict[str, dict[str, Any]] = {}
    used_maps: set[int] = set()
    used_indices: set[int] = set()
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("Native parent atom_mapping rows must be dictionaries")
        name = row.get("ccd_atom_id")
        atom_map = row.get("atom_map")
        index = row.get("rdkit_index_zero_based")
        if (
            not isinstance(name, str)
            or not name
            or isinstance(atom_map, bool)
            or not isinstance(atom_map, int)
            or atom_map <= 0
            or isinstance(index, bool)
            or not isinstance(index, int)
            or index < 0
        ):
            raise ValueError("Native parent atom_mapping is malformed")
        if name in result or atom_map in used_maps or index in used_indices:
            raise ValueError("Native parent atom_mapping is not one-to-one")
        result[name] = dict(row)
        used_maps.add(atom_map)
        used_indices.add(index)
    return result


def _validate_native_selection(
    parent: Chem.Mol,
    metadata: dict[str, Any],
    ligand_sites: list[dict[str, Any]],
) -> tuple[Chem.Mol, list[dict[str, Any]]]:
    if metadata.get("source_record_id") != "9D12:E:A1A1P":
        raise ValueError("Native ligand source_record_id must be 9D12:E:A1A1P")
    if metadata.get("source_chain") != "A":
        raise ValueError("Native source protein chain must be label chain A")
    if metadata.get("pdb") != "9D12" or metadata.get("ccd") != "A1A1P":
        raise ValueError("Native parent metadata must identify 9D12/A1A1P")
    if {str(row.get("label_asym_id")) for row in ligand_sites} != {"E"}:
        raise ValueError("Native ligand selection must contain only label chain E")
    if {str(row.get("label_comp_id")) for row in ligand_sites} != {"A1A1P"}:
        raise ValueError("Native ligand selection has an unexpected residue")

    mapping = _mapping_rows(metadata)
    by_name: dict[str, dict[str, Any]] = {}
    for row in ligand_sites:
        name = str(row.get("label_atom_id"))
        if name in by_name:
            raise ValueError("Native ligand selection contains a duplicate CCD atom name")
        by_name[name] = row
    if set(by_name) != set(mapping):
        raise ValueError("Native ligand atom sites and registered CCD mapping differ")
    if len(mapping) != parent.GetNumHeavyAtoms():
        raise ValueError("Native ligand mapping does not cover every parent heavy atom")

    original = Chem.Mol(parent)
    positioned = Chem.Mol(parent)
    positioned.RemoveAllConformers()
    conformer = Chem.Conformer(positioned.GetNumAtoms())
    conformer.Set3D(True)
    receipt = []
    seen_indices: set[int] = set()
    for name, mapped in mapping.items():
        index = mapped["rdkit_index_zero_based"]
        if index >= positioned.GetNumAtoms() or index in seen_indices:
            raise ValueError("Native ligand mapping references an invalid RDKit index")
        seen_indices.add(index)
        atom = positioned.GetAtomWithIdx(index)
        if atom.GetAtomicNum() <= 1:
            raise ValueError("Native ligand mapping references a hydrogen")
        if atom.GetAtomMapNum() != mapped["atom_map"]:
            raise ValueError("Native ligand map number differs from registered metadata")
        site = by_name[name]
        if str(site.get("type_symbol", "")).title() != atom.GetSymbol():
            raise ValueError("Native ligand CCD/atom-site element mismatch")
        xyz = site.get("xyz")
        if not isinstance(xyz, list) or len(xyz) != 3 or not all(
            math.isfinite(float(value)) for value in xyz
        ):
            raise ValueError("Native ligand atom site has malformed coordinates")
        conformer.SetAtomPosition(index, tuple(float(value) for value in xyz))
        receipt.append(
            {
                "atom_map": atom.GetAtomMapNum(),
                "rdkit_index_zero_based": index,
                "ccd_atom_id": name,
                "atom_site_id": str(site.get("id")),
                "element": atom.GetSymbol(),
                "label_asym_id": "E",
                "label_comp_id": "A1A1P",
                "auth_seq_id": str(site.get("auth_seq_id")),
                "xyz_A": [float(value) for value in xyz],
            }
        )
    positioned.AddConformer(conformer, assignId=True)
    if _canonical_mapped(positioned) != _canonical_mapped(original):
        raise ValueError("Assigning native coordinates changed parent graph identity")
    if Chem.GetFormalCharge(positioned) != Chem.GetFormalCharge(original):
        raise ValueError("Assigning native coordinates changed parent charge")
    if [a.GetAtomMapNum() for a in positioned.GetAtoms()] != [
        a.GetAtomMapNum() for a in original.GetAtoms()
    ]:
        raise ValueError("Assigning native coordinates changed parent atom maps")
    return positioned, sorted(receipt, key=lambda row: row["atom_map"])


def _pdb_atom_name(name: str, element: str) -> str:
    value = str(name)[:4]
    if len(value) == 4:
        return value
    if len(element.strip()) == 1 and not value[:1].isdigit():
        return (" " + value).ljust(4)
    return value.ljust(4)


def native_receptor_pdb(protein: list[dict[str, Any]]) -> str:
    """Serialize the exact selected native protein heavy atoms as chain-A PDB."""
    if not isinstance(protein, list) or not protein:
        raise ValueError("Native protein selection is empty")
    lines = []
    for serial, row in enumerate(protein, start=1):
        if str(row.get("label_asym_id")) != "A":
            raise ValueError("Native receptor may contain only label chain A")
        residue = str(row.get("label_comp_id", "")).upper()
        if residue not in _CANONICAL_RESIDUES:
            raise ValueError(
                "Native receptor contains a noncanonical residue, water, metal, or ligand"
            )
        element = str(row.get("type_symbol", "")).strip().title()
        atom_name = str(row.get("label_atom_id", ""))
        if not atom_name or element in {"", "H", "D"}:
            raise ValueError("Native receptor contains an invalid heavy atom")
        try:
            auth_seq = int(str(row.get("auth_seq_id")))
        except ValueError as exc:
            raise ValueError("Native receptor auth_seq_id is not an integer") from exc
        insertion = row.get("pdbx_PDB_ins_code")
        insertion = " " if insertion in (None, False, "", ".", "?") else str(insertion)[0]
        xyz = [float(value) for value in row.get("xyz", [])]
        if len(xyz) != 3 or not all(math.isfinite(value) for value in xyz):
            raise ValueError("Native receptor has malformed coordinates")
        occupancy = float(row.get("occupancy", 1.0))
        bfactor = float(row.get("B_iso_or_equiv", 0.0))
        if not math.isfinite(occupancy) or not math.isfinite(bfactor):
            raise ValueError("Native receptor occupancy or B factor is invalid")
        lines.append(
            f"ATOM  {serial:5d} {_pdb_atom_name(atom_name, element)} "
            f"{residue:>3s} A{auth_seq:4d}{insertion}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"{occupancy:6.2f}{bfactor:6.2f}          {element:>2s}  "
        )
    lines.extend(["TER", "END"])
    return "\n".join(lines) + "\n"


def prepare_native_receptor(
    protein: list[dict[str, Any]],
    output: Path,
    validated_residues_to_delete: _ValidatedResidueDeletions | None = None,
) -> dict[str, Any]:
    """Prepare a rigid native receptor, failing closed on every Meeko error."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy, Polymer, ResidueChemTemplates

    if validated_residues_to_delete is None:
        validated_residues_to_delete = _ValidatedResidueDeletions()
    if not isinstance(validated_residues_to_delete, _ValidatedResidueDeletions):
        raise ValueError("Receptor residue deletions must come from explicit validation")
    delete_ids = list(validated_residues_to_delete)

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    pdb_path = output / "native-receptor.pdb"
    pdbqt_path = output / "native-receptor.pdbqt"
    pdb_text = native_receptor_pdb(protein)
    pdb_path.write_text(pdb_text, encoding="ascii", newline="\n")

    templates = ResidueChemTemplates.create_from_defaults()
    preparation = MoleculePreparation()
    polymer = Polymer.from_pdb_string(
        pdb_text,
        chem_templates=templates,
        mk_prep=preparation,
        allow_bad_res=False,
        residues_to_delete=delete_ids or None,
        bonds_to_delete=None,
        delete_bad_res_from_box_radius=None,
    )
    written = PDBQTWriterLegacy.write_from_polymer(polymer)
    if not isinstance(written, tuple) or len(written) != 2:
        raise RuntimeError("Unexpected Meeko polymer writer result")
    rigid, flex = written
    if not isinstance(rigid, str) or not rigid.strip():
        raise RuntimeError("Meeko produced no rigid receptor PDBQT")
    if flex:
        raise RuntimeError("Native receptor preparation unexpectedly produced flexible residues")
    pdbqt_path.write_text(rigid, encoding="ascii", newline="\n")
    return {
        "pdb": str(pdb_path),
        "pdbqt": str(pdbqt_path),
        "pdb_sha256": _sha256(pdb_path),
        "pdbqt_sha256": _sha256(pdbqt_path),
        "parameters": {
            "Polymer.from_pdb_string": {
                "chem_templates": "ResidueChemTemplates.create_from_defaults",
                "mk_prep": "MoleculePreparation defaults",
                "allow_bad_res": False,
                "residues_to_delete": delete_ids or None,
                "bonds_to_delete": None,
                "delete_bad_res_from_box_radius": None,
            },
            "network_template_use": False,
            "canonical_amino_acids_only": True,
            "bad_residue_deletion": bool(delete_ids),
            "explicit_remote_incomplete_residue_omission": bool(delete_ids),
            "fixed_receptor_fallback": False,
            "hydrogens": "added by receptor preparation; not experimentally measured",
            "pH": "unknown; no pH-based state selection",
        },
        "versions": {
            "meeko": _version("meeko"),
            "rdkit": _version("rdkit"),
        },
    }


def load_native_context(
    parent_id: str,
    cif_path: Path,
    output: Path,
    exclude_remote_incomplete_residues: Any = (),
):
    """Load a registered parent and replace only its coordinates with native 9D12 sites."""
    from packages.science.dual_e3 import REFERENCE_SOURCE
    from packages.science.reference_parents import load_reference_parent
    from packages.science.structures import atom_sites

    cif_path = Path(cif_path).resolve()
    if not cif_path.is_file():
        raise FileNotFoundError(cif_path)
    if parent_id != "SMARCA2-9D12-A1A1P":
        raise ValueError("This native experiment is registered only for SMARCA2-9D12-A1A1P")
    parent, metadata = load_reference_parent(parent_id, REFERENCE_SOURCE)
    actual_hash = _sha256(cif_path)
    expected_hash = metadata.get("source_input_hashes", {}).get(
        "source_pdb_mmcif_sha256"
    )
    if not isinstance(expected_hash, str) or not _SHA256.fullmatch(expected_hash):
        raise ValueError("Registered native CIF hash is missing or malformed")
    if actual_hash != expected_hash:
        raise ValueError("Native CIF hash differs from registered parent provenance")

    ligand_sites = atom_sites(cif_path, ["E"])
    native_parent, mapping_receipt = _validate_native_selection(
        parent, metadata, ligand_sites
    )
    protein = atom_sites(cif_path, ["A"])
    if {str(row.get("label_asym_id")) for row in protein} != {"A"}:
        raise ValueError("Native receptor selection must contain exactly label chain A")
    if any(str(row.get("label_comp_id", "")).upper() not in _CANONICAL_RESIDUES for row in protein):
        raise ValueError("Native receptor selection contains waters, metals, or nonprotein residues")

    validated_deletions, omission_records = validate_remote_incomplete_residue_exclusions(
        protein, mapping_receipt, exclude_remote_incomplete_residues
    )
    omitted = set(validated_deletions)
    evaluation_protein = [row for row in protein if _residue_id(row) not in omitted]

    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    try:
        receptor = prepare_native_receptor(protein, output, validated_deletions)
    except Exception as exc:
        failure = {
            "status": "failed_receptor_preflight",
            "error": str(exc),
            "parent_id": parent_id,
            "native_cif": str(cif_path),
            "native_cif_sha256": actual_hash,
            "original_native_protein_heavy_atom_count": len(protein),
            "explicit_residue_omissions": omission_records,
            "allow_bad_res": False,
            "fixed_receptor_fallback": False,
            "all_runs_complete": False,
        }
        _write_json(output / "native-context-failure.json", failure)
        raise
    metadata_hash = hashlib.sha256(_json_bytes(metadata)).hexdigest()
    native_parent.SetProp("native_parent_id", parent_id)
    native_parent.SetProp("native_cif_sha256", actual_hash)
    native_parent.SetProp("native_metadata_sha256", metadata_hash)
    receipt = {
        "parent_id": parent_id,
        "source_record_id": "9D12:E:A1A1P",
        "native_cif": str(cif_path),
        "native_cif_sha256": actual_hash,
        "registered_source_pdb_mmcif_sha256": expected_hash,
        "registered_parent_metadata_sha256": metadata_hash,
        "ligand_selection": {
            "label_asym_id": "E",
            "label_comp_id": "A1A1P",
            "model": 1,
            "hydrogens": "excluded",
            "altloc": "blank or A",
            "mapping_policy": "registered atom map to exact CCD atom name; no inferred cross-parent mapping",
            "atoms": mapping_receipt,
        },
        "protein_selection": {
            "label_asym_id": "A",
            "model": 1,
            "hydrogens": "excluded",
            "waters": "excluded",
            "metals": "excluded",
            "other_chains": "excluded",
            "noncanonical_residues": "refused",
            "original_heavy_atom_count": len(protein),
            "docking_receptor_heavy_atom_count": len(evaluation_protein),
            "original_native_pdb_retains_all_coordinates": True,
            "explicit_remote_incomplete_residue_omissions": omission_records,
            "automatic_residue_deletion": False,
            "remote_omission_scientific_limitation": bool(omission_records),
        },
        "receptor_preparation": receptor,
        "coordinate_frame": "deposited 9D12 coordinates; no alignment",
        "scientific_scope": "exploratory protocol comparison only; no approval",
    }
    _write_json(output / "native-context.json", receipt)
    return native_parent, evaluation_protein, receipt


def _dock_failure(result: dict[str, Any], status: str, error: Any) -> dict[str, Any]:
    result["status"] = status
    result["results"] = {"error": str(error)}
    return _safe_number(result)


def dock_native(
    mol,
    parent,
    protein,
    protected_maps,
    receptor_pdbqt,
    workdir,
    seed=23,
    exhaustiveness=16,
    timeout=180,
):
    """Run five actual Vina poses in the unchanged native receptor frame."""
    module_root = Path(__file__).resolve().parents[2]
    vina = module_root / "vendor" / "vina" / "vina_1.2.7_win.exe"
    run_dir = Path(workdir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {
        "status": "failed",
        "method": _METHOD,
        "diagnostic_role": "completed_with_limits_not_official_M2_gate_or_strict_acceptance",
        "real_docking": False,
        "native_receptor_frame": True,
        "alignment_applied": False,
        "seed": int(seed),
        "exhaustiveness": int(exhaustiveness),
        "timeout_seconds": int(timeout),
        "cpu": 2,
        "num_modes_requested": 5,
        "files": {
            "ligand_pdbqt": "ligand.pdbqt",
            "receptor_pdbqt": "receptor.pdbqt",
            "poses_pdbqt": "poses.pdbqt",
            "poses_sdf": "poses.sdf",
            "stdout": "vina.stdout.txt",
            "stderr": "vina.stderr.txt",
        },
        "sourcehashes": {},
        "parameters": {},
        "results": {},
        "limitations": [
            "Exploratory deposited-frame docking is not scientific approval.",
            "Hydrogens are computationally added and pH is unknown.",
            "No explicit waters are included.",
            "Scores are not affinity measurements.",
        ],
    }
    try:
        if (
            isinstance(seed, bool)
            or not isinstance(seed, int)
            or seed < 0
            or seed > 0xFFFFFFFF
        ):
            raise ValueError("seed must be a nonnegative 32-bit integer")
        if (
            isinstance(exhaustiveness, bool)
            or not isinstance(exhaustiveness, int)
            or not 1 <= exhaustiveness <= 64
        ):
            raise ValueError("exhaustiveness must be an integer from 1 through 64")
        if int(timeout) != 180:
            raise ValueError("Native experiment timeout is fixed at 180 seconds")
        if not vina.is_file():
            raise FileNotFoundError("Pinned Vina 1.2.7 binary is missing")
        receptor_pdbqt = Path(receptor_pdbqt).resolve()
        if not receptor_pdbqt.is_file() or receptor_pdbqt.stat().st_size == 0:
            raise FileNotFoundError("Prepared native receptor PDBQT is missing")
        protected = _validate_inputs(mol, parent, protected_maps)
        if Chem.GetFormalCharge(mol) != 0:
            raise ValueError("Native docking ligand must be neutral")
        if not parent.HasProp("native_cif_sha256"):
            raise ValueError("Parent lacks verified native-context provenance")
        if not isinstance(protein, list) or not protein:
            raise ValueError("Native protein atom selection is missing")

        result["sourcehashes"] = {
            "vina_1.2.7_win.exe": _sha256(vina),
            "native_receptor_pdbqt": _sha256(receptor_pdbqt),
            "native_cif_sha256": parent.GetProp("native_cif_sha256"),
            "native_parent_metadata_sha256": parent.GetProp(
                "native_metadata_sha256"
            ),
        }
        receptor = run_dir / "receptor.pdbqt"
        shutil.copyfile(receptor_pdbqt, receptor)
        ligand_path = run_dir / "ligand.pdbqt"
        poses_path = run_dir / "poses.pdbqt"
        sdf_path = run_dir / "poses.sdf"

        from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate

        ligand = Chem.AddHs(Chem.Mol(mol))
        params = AllChem.ETKDGv3()
        params.randomSeed = int(seed)
        params.clearConfs = True
        if AllChem.EmbedMolecule(ligand, params) != 0:
            raise RuntimeError("ETKDGv3 ligand embedding failed")
        if not AllChem.MMFFHasAllMoleculeParams(ligand):
            raise RuntimeError("MMFF parameters are unavailable for the ligand")
        optimization = AllChem.MMFFOptimizeMolecule(ligand, maxIters=100)
        setups = MoleculePreparation().prepare(ligand)
        if len(setups) != 1:
            raise RuntimeError("Meeko did not produce exactly one ligand setup")
        pdbqt, success, message = PDBQTWriterLegacy.write_string(setups[0])
        if not success or not pdbqt.strip():
            raise RuntimeError(f"Meeko ligand writing failed: {message}")
        ligand_path.write_text(pdbqt, encoding="ascii", newline="\n")

        center, size = _box(parent)
        command = [
            _short_windows_path(vina),
            "--receptor", "receptor.pdbqt",
            "--ligand", "ligand.pdbqt",
            "--center_x", f"{center[0]:.6f}",
            "--center_y", f"{center[1]:.6f}",
            "--center_z", f"{center[2]:.6f}",
            "--size_x", f"{size[0]:.6f}",
            "--size_y", f"{size[1]:.6f}",
            "--size_z", f"{size[2]:.6f}",
            "--exhaustiveness", str(int(exhaustiveness)),
            "--seed", str(int(seed)),
            "--cpu", "2",
            "--num_modes", "5",
            "--out", "poses.pdbqt",
        ]
        result["parameters"] = {
            "command": command,
            "embedding": "RDKit ETKDGv3",
            "embedding_seed": int(seed),
            "optimization": "MMFF maxIters=100",
            "optimization_return_code": int(optimization),
            "ligand_preparation": "Meeko MoleculePreparation defaults",
            "receptor_frame": "native 9D12; no ligand or receptor alignment",
            "versions": {
                "rdkit": _version("rdkit"),
                "meeko": _version("meeko"),
            },
        }
        creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        with (run_dir / "vina.stdout.txt").open("wb") as stdout, (
            run_dir / "vina.stderr.txt"
        ).open("wb") as stderr:
            process = subprocess.Popen(
                command,
                cwd=str(run_dir),
                stdout=stdout,
                stderr=stderr,
                shell=False,
                creationflags=creationflags,
            )
            result["real_docking"] = True
            started = time.monotonic()
            try:
                while process.poll() is None:
                    if time.monotonic() - started > 180.0:
                        _kill(process)
                        raise TimeoutError("Vina exceeded the 180 second timeout")
                    time.sleep(0.2)
            except BaseException:
                _kill(process)
                raise
            if process.returncode != 0:
                raise RuntimeError(f"Vina failed with exit code {process.returncode}")

        if not poses_path.is_file() or poses_path.stat().st_size == 0:
            raise RuntimeError("Vina produced no pose output")
        text = poses_path.read_text(encoding="utf-8", errors="replace")
        scores = [
            {
                "affinity": float(match.group(1)),
                "rmsd_lb": float(match.group(2)) if match.group(2) else None,
                "rmsd_ub": float(match.group(3)) if match.group(3) else None,
            }
            for match in _VINA_RESULT.finditer(text)
        ]
        blocks = re.findall(r"MODEL\s+\d+.*?ENDMDL", text, flags=re.S)
        if not 1 <= len(blocks) <= 5 or len(scores) != len(blocks):
            raise ValueError("Native docking must produce one through five paired poses and scores")
        actual_pose_count = len(blocks)
        if actual_pose_count < 5:
            result.setdefault("warnings", []).append(
                f"Vina returned {actual_pose_count} mode(s) although five were requested"
            )

        restored = []
        for block in blocks:
            converted = RDKitMolCreate.from_pdbqt_mol(
                PDBQTMolecule(block + "\n", skip_typing=True)
            )
            if len(converted) != 1 or converted[0] is None:
                raise ValueError("Meeko could not restore a native docking pose")
            pose, mapping_count, direct_rmsd = _restore_maps(
                converted[0], mol, parent, protected
            )
            restored.append((pose, mapping_count, direct_rmsd))

        sdf_stream = sdf_path.open("w", encoding="utf-8", newline="\n")
        writer = Chem.SDWriter(sdf_stream)
        writer.SetForceV3000(True)
        pose_records = []
        summaries = []
        try:
            for index, ((pose, mapping_count, direct_rmsd), score) in enumerate(
                zip(restored, scores), start=1
            ):
                pose.SetIntProp("vina_mode", index)
                pose.SetDoubleProp("vina_affinity", score["affinity"])
                pose.SetProp("coordinate_frame", "native_9D12_no_alignment")
                pose.SetBoolProp("alignment_applied", False)
                if direct_rmsd is not None:
                    pose.SetDoubleProp("parent_protected_core_RMSD", direct_rmsd)
                writer.write(pose)
                pose_records.append({"mol": pose, "score": score["affinity"]})
                summaries.append(
                    {
                        "mode": index,
                        "score": score,
                        "mapping_count": mapping_count,
                        "core_RMSD_A_in_native_receptor_frame": direct_rmsd,
                        "alignment_applied": False,
                    }
                )
        finally:
            writer.close()
            sdf_stream.close()

        from packages.science.warhead_sites import pose_preservation

        preservation = pose_preservation(
            parent, pose_records, protein, sorted(protected)
        )
        if any(
            row.get("alignment_applied") is not False
            for row in preservation.get("all_poses_diagnostics", [])
        ):
            raise ValueError("Native pose preservation unexpectedly applied alignment")
        result["status"] = "completed_with_limits"
        result["results"] = {
            "pose_count": actual_pose_count,
            "num_modes_requested": 5,
            "num_modes_actual": actual_pose_count,
            "poses": summaries,
            "pose_preservation": preservation,
            "box_center": center,
            "box_size": size,
            "SDF_format": "V3000",
        }
        return _safe_number(result)
    except TimeoutError as exc:
        return _dock_failure(result, "failed_timeout", exc)
    except Exception as exc:
        return _dock_failure(result, "failed_docking", exc)


def validate_seeds(seeds: Any) -> list[int]:
    if not isinstance(seeds, (list, tuple)) or not seeds:
        raise ValueError("seeds must be a nonempty list")
    result = []
    for value in seeds:
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError("seeds must be integers")
        number = value
        if number < 0 or number > 0xFFFFFFFF:
            raise ValueError("seeds must be nonnegative 32-bit integers")
        result.append(number)
    if len(result) != len(set(result)):
        raise ValueError("Duplicate docking seeds are forbidden")
    return result


def validate_halogen_pair(parent: Chem.Mol, candidate: Chem.Mol) -> None:
    """Require the exact mapped Cl-to-Br graph change at map 8 and nothing else."""
    if not isinstance(parent, Chem.Mol) or not isinstance(candidate, Chem.Mol):
        raise TypeError("Halogen pair must contain RDKit molecules")
    parent_atoms = {a.GetAtomMapNum(): a for a in parent.GetAtoms() if a.GetAtomicNum() > 1}
    candidate_atoms = {
        a.GetAtomMapNum(): a for a in candidate.GetAtoms() if a.GetAtomicNum() > 1
    }
    if set(parent_atoms) != set(candidate_atoms) or 8 not in parent_atoms:
        raise ValueError("Halogen pair mapped atom sets differ")
    for atom_map in sorted(parent_atoms):
        left, right = parent_atoms[atom_map], candidate_atoms[atom_map]
        expected = 35 if atom_map == 8 else left.GetAtomicNum()
        if right.GetAtomicNum() != expected:
            raise ValueError("Halogen pair differs outside the exact map-8 Cl-to-Br change")
        if atom_map == 8 and left.GetAtomicNum() != 17:
            raise ValueError("Native parent map 8 is not chlorine")
        if (
            left.GetFormalCharge() != right.GetFormalCharge()
            or left.GetIsAromatic() != right.GetIsAromatic()
            or left.GetIsotope() != right.GetIsotope()
        ):
            raise ValueError("Halogen pair changes atom chemistry beyond element at map 8")
    for bond in parent.GetBonds():
        first = bond.GetBeginAtom().GetAtomMapNum()
        second = bond.GetEndAtom().GetAtomMapNum()
        other = candidate.GetBondBetweenAtoms(
            candidate_atoms[first].GetIdx(), candidate_atoms[second].GetIdx()
        )
        if other is None or other.GetBondType() != bond.GetBondType():
            raise ValueError("Halogen pair bond graph differs")
    if parent.GetNumBonds() != candidate.GetNumBonds():
        raise ValueError("Halogen pair has extra bonds")


def verify_source_pair(export_dir: Path, expected_parent_id: str):
    """Use only the verifier result and its hash-verified embedded source object."""
    from packages.science.parent_sar_probe import verify_parent_export

    export_dir = Path(export_dir).resolve()
    manifest = export_dir / "manifest.json"
    receipt = verify_parent_export(manifest, export_dir.parent)
    if not isinstance(receipt, dict) or not isinstance(receipt.get("result"), dict):
        raise ValueError("Parent export verification returned no result")
    result = receipt["result"]
    selected = result.get("selected_parent_measured_evidence")
    if not isinstance(selected, dict):
        raise ValueError("Verified export has no selected_parent_measured_evidence")
    if selected.get("selected_parent_id") != expected_parent_id:
        raise ValueError("Measured-evidence export identifies a different parent")
    sources = result.get("source_catalog", {}).get("sources", {})
    source_key = "medchem_evidence.json"
    source_object = sources.get(source_key)
    if not isinstance(source_object, dict):
        raise ValueError("Verified result has no embedded medchem evidence source object")
    source_sha = source_object.get("sha256")
    if (
        source_object.get("status") != "configured_hash_verified"
        or not isinstance(source_sha, str)
        or not _SHA256.fullmatch(source_sha.lower())
        or "data" not in source_object
    ):
        raise ValueError("Embedded medchem evidence source is not hash verified")
    medchem_ref = {
        "path": "result.source_catalog.sources.medchem_evidence.json",
        "sha256": source_sha.lower(),
    }
    return {"result": result}, selected, source_object, medchem_ref


def _candidate_from_builder(value: Any) -> Chem.Mol:
    if isinstance(value, Chem.Mol):
        return value
    if isinstance(value, tuple):
        matches = [item for item in value if isinstance(item, Chem.Mol)]
    elif isinstance(value, dict):
        matches = [value[key] for key in ("mol", "candidate", "br_mol") if isinstance(value.get(key), Chem.Mol)]
    else:
        matches = []
    if len(matches) != 1:
        raise ValueError("Halogen probe builder did not return exactly one RDKit molecule")
    return matches[0]


def _ligand_input(parent: Chem.Mol) -> Chem.Mol:
    result = Chem.Mol(parent)
    result.RemoveAllConformers()
    return result


def _graph_sha256(mol: Chem.Mol) -> str:
    return hashlib.sha256(_canonical_mapped(mol).encode("utf-8")).hexdigest()


def _protected_maps_from_sites(site_atoms: Any) -> list[int]:
    if not isinstance(site_atoms, list):
        raise ValueError("Native site analysis has no atom rows")
    result = []
    for row in site_atoms:
        if not isinstance(row, dict):
            continue
        atom_map = row.get("atom_map")
        if row.get("state") == "PROTECTED":
            if isinstance(atom_map, bool) or not isinstance(atom_map, int) or atom_map <= 0:
                raise ValueError("PROTECTED native site row has no valid atom map")
            result.append(atom_map)
    if not result:
        raise ValueError("Native site analysis selected no PROTECTED atom maps")
    if len(result) != len(set(result)):
        raise ValueError("Native site analysis has duplicate PROTECTED atom maps")
    return sorted(result)


def _map_nearest_distance(site_atoms: Any, atom_map: int) -> float:
    if not isinstance(site_atoms, list):
        raise ValueError("Site atom context is missing")
    rows = [row for row in site_atoms if isinstance(row, dict) and row.get("atom_map") == atom_map]
    if len(rows) != 1:
        raise ValueError(f"Site context must contain exactly one map {atom_map} row")
    row = rows[0]
    evidence = row.get("evidence")
    contacts = evidence.get("contacts") if isinstance(evidence, dict) else None
    value = contacts.get("nearest_heavy_atom_A") if isinstance(contacts, dict) else None
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
    ):
        raise ValueError(f"Site map {atom_map} row has no finite nearest-distance measurement")
    return float(value)


def _aggregate_runs(runs: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [run for run in runs if run.get("status") == "completed_with_limits"]
    affinities = []
    rmsds = []
    contacts = []
    clashes = []
    for run in completed:
        for pose in run.get("results", {}).get("poses", []):
            score = pose.get("score", {}).get("affinity")
            if isinstance(score, (int, float)):
                affinities.append(float(score))
        preservation = run.get("results", {}).get("pose_preservation", {})
        for pose in preservation.get("all_poses_diagnostics", []):
            rmsd = pose.get("core_RMSD_A_in_receptor_frame")
            if isinstance(rmsd, (int, float)):
                rmsds.append(float(rmsd))
            retained = pose.get("protected_atom_near_residue_preservation", {}).get(
                "retention_fraction"
            )
            if isinstance(retained, (int, float)):
                contacts.append(float(retained))
            overlap = pose.get("clashes", {}).get("maximum_heavy_atom_vdw_overlap_A")
            if isinstance(overlap, (int, float)):
                clashes.append(float(overlap))

    def distribution(values: list[float]) -> dict[str, Any]:
        return {
            "count": len(values),
            "minimum": min(values) if values else None,
            "maximum": max(values) if values else None,
            "mean": sum(values) / len(values) if values else None,
            "values": values,
        }

    return {
        "predeclared_run_count": len(runs),
        "completed_run_count": len(completed),
        "failed_run_count": len(runs) - len(completed),
        "affinity_distribution": distribution(affinities),
        "exact_core_RMSD_A_distribution": distribution(rmsds),
        "protected_contact_retention_distribution": distribution(contacts),
        "maximum_clash_overlap_A_distribution": distribution(clashes),
        "all_results_retained": True,
        "all_pass_claimed": False,
    }


def _manifest(output: Path) -> dict[str, Any]:
    files = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            relative = path.relative_to(output).as_posix()
            files.append(
                {"path": relative, "bytes": path.stat().st_size, "sha256": _sha256(path)}
            )
    return {
        "format_version": "native-parent-probe-v1",
        "hash_algorithm": "sha256",
        "files": files,
        "manifest_self_hash_policy": "manifest.json is excluded",
    }


def run_native_probe(
    native_cif: Path,
    output: Path,
    parent_id: str = "SMARCA2-9D12-A1A1P",
    seeds: Any = (23, 41, 61),
    exhaustiveness: int = 16,
    source_export: Path | None = None,
    exclude_remote_incomplete_residues: Any = (),
) -> dict[str, Any]:
    """Run the predeclared parent/Br experiment and retain every run result."""
    seeds = validate_seeds(seeds)
    if len(seeds) != 3:
        raise ValueError("Native protocol requires exactly three distinct seeds")
    if (
        isinstance(exhaustiveness, bool)
        or not isinstance(exhaustiveness, int)
        or not 1 <= exhaustiveness <= 64
    ):
        raise ValueError("exhaustiveness must be an integer from 1 through 64")
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    if source_export is None:
        raise ValueError("A verified parent source export is required for the measured Cl/Br pair")

    context_dir = output / "receptorfiles"
    try:
        parent, protein, context = load_native_context(
            parent_id,
            native_cif,
            context_dir,
            exclude_remote_incomplete_residues=exclude_remote_incomplete_residues,
        )
    except Exception as exc:
        output.mkdir(parents=True, exist_ok=True)
        failure = {
            "status": "failed_receptor_preflight",
            "diagnostic": "no_seed_runs_started",
            "error": str(exc),
            "parent_id": parent_id,
            "all_runs_complete": False,
            "scientific_scope": "remote receptor omission is scientifically limited",
        }
        _write_json(output / "native-probe-result.json", failure)
        _write_json(output / "manifest.json", _manifest(output))
        return failure
    from packages.science.warhead_sites import analyze_sites

    native_sites = analyze_sites(parent, protein, {}, [])
    verified, selected_evidence, medchem_source, medchem_ref = verify_source_pair(
        source_export, parent_id
    )
    from packages.science.parent_sar_probe import build_exact_halogen_probe

    linkage = dict(selected_evidence, sites=native_sites["atoms"])
    proposal = build_exact_halogen_probe(parent, linkage, medchem_source)
    if not isinstance(proposal, dict) or not isinstance(proposal.get("mapped_smiles"), str):
        raise ValueError("Halogen probe builder returned no mapped_smiles proposal")
    candidate = Chem.MolFromSmiles(proposal["mapped_smiles"])
    if candidate is None:
        raise ValueError("Halogen probe mapped_smiles could not be parsed")
    validate_halogen_pair(parent, candidate)
    if Chem.GetFormalCharge(candidate) != 0:
        raise ValueError("Measured Br probe must be neutral")
    clean = Chem.Mol(candidate)
    clean.RemoveAllConformers()
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    if any(str(info.specified) == "Unspecified" for info in Chem.FindPotentialStereo(clean)):
        raise ValueError("Measured Br probe has unspecified stereochemistry")

    parent_protected_maps = _protected_maps_from_sites(native_sites["atoms"])
    candidate_protected_maps = [value for value in parent_protected_maps if value != 8]
    if parent_protected_maps == candidate_protected_maps:
        raise ValueError("Parent and candidate masks must remain explicitly distinct")
    predeclared = [
        {"ligand": ligand, "seed": seed}
        for ligand in ("native_parent_Cl", "measured_probe_Br")
        for seed in seeds
    ]
    protocol = {
        "format_version": "native-parent-probe-v1",
        "status": "predeclared_before_runs",
        "parent_id": parent_id,
        "native_context": context,
        "source_export_receipt": verified,
        "source_medchem_evidence": medchem_ref,
        "source_proof": proposal,
        "site_analysis": native_sites,
        "distinct_context_comparison": {
            "native_deposited_9D12": {
                "map_3_nearest_A": _map_nearest_distance(native_sites["atoms"], 3),
                "sites_source": "native deposited parent conformer",
            },
            "verified_fixed_context": {
                "map_3_nearest_A": _map_nearest_distance(verified["result"]["sites"]["atoms"], 3),
                "sites_source": "verified parent export result",
            },
            "frames_differ": True,
            "scope": "descriptive protocol comparison only; no policy mutation",
        },
        "candidate_rule": "exact measured Cl-to-Br graph; only mapped atom 8 changes element",
        "candidate_excluded_changed_map": 8,
        "mask_record": {
            "parent_native_PROTECTED_maps": parent_protected_maps,
            "candidate_parent_mask_minus_changed_map_8": candidate_protected_maps,
            "strictly_incomparable": True,
        },
        "modifiable_promotion": False,
        "incomparable_masks": "strict",
        "input_graph_sha256": {
            "parent": _graph_sha256(parent),
            "candidate": _graph_sha256(candidate),
        },
        "runs": predeclared,
        "docking": {
            "seeds": seeds,
            "exhaustiveness": int(exhaustiveness),
            "cpu": 2,
            "timeout_seconds": 180,
            "num_modes": 5,
            "coordinate_frame": "native 9D12 without post-hoc alignment",
        },
        "interpretation": "experiment only; not official M2, gate, strict acceptance, or approval",
    }
    _write_json(output / "protocol.json", protocol)
    ligands_dir = output / "ligands"
    ligands_dir.mkdir()
    for name, molecule in (("native-parent-Cl.sdf", parent), ("measured-probe-Br.sdf", candidate)):
        sdf_stream = (ligands_dir / name).open("w", encoding="utf-8", newline="\n")
        writer = Chem.SDWriter(sdf_stream)
        writer.SetForceV3000(True)
        try:
            writer.write(molecule)
        finally:
            writer.close()
            sdf_stream.close()

    receptor_pdbqt = Path(context["receptor_preparation"]["pdbqt"])
    all_runs = []
    for ligand_name, molecule, protected_maps in (
        ("native_parent_Cl", _ligand_input(parent), parent_protected_maps),
        ("measured_probe_Br", _ligand_input(candidate), candidate_protected_maps),
    ):
        for seed in seeds:
            run_dir = output / "runs" / ligand_name / f"seed-{seed}"
            result = dock_native(
                molecule,
                parent,
                protein,
                protected_maps,
                receptor_pdbqt,
                run_dir,
                seed=seed,
                exhaustiveness=exhaustiveness,
                timeout=180,
            )
            result["predeclared_ligand"] = ligand_name
            _write_json(run_dir / "per-seed-result.json", result)
            all_runs.append(result)

    aggregate = _aggregate_runs(all_runs)
    final = {
        "status": "completed_with_limits" if aggregate["completed_run_count"] else "failed",
        "diagnostic": "not_official_M2_gate_or_strict_acceptance",
        "parent_id": parent_id,
        "runs": all_runs,
        "analysis": aggregate,
        "coordinate_frame": "native 9D12",
        "alignment_applied": False,
        "scientific_scope": "protocol comparison only; no approval",
    }
    _write_json(output / "native-probe-result.json", final)
    summary = [
        "# Native 9D12 parent/probe docking experiment",
        "",
        "This is an exploratory deposited-receptor-frame protocol comparison only.",
        "It is not official M2, a gate result, strict acceptance, or scientific approval.",
        "",
        f"- Parent: `{parent_id}`",
        f"- Predeclared runs: {aggregate['predeclared_run_count']}",
        f"- Completed runs: {aggregate['completed_run_count']}",
        f"- Failed runs: {aggregate['failed_run_count']}",
        "- Poses are compared directly in the native 9D12 frame without alignment.",
        "- Every run result, including failures, is retained.",
        "- No UNKNOWN site was promoted to MODIFIABLE.",
        "",
    ]
    (output / "summary.md").write_text("\n".join(summary), encoding="utf-8", newline="\n")
    _write_json(output / "manifest.json", _manifest(output))
    return final


__all__ = [
    "dock_native",
    "load_native_context",
    "native_receptor_pdb",
    "prepare_native_receptor",
    "run_native_probe",
    "validate_halogen_pair",
    "validate_remote_incomplete_residue_exclusions",
    "validate_seeds",
    "verify_source_pair",
]
