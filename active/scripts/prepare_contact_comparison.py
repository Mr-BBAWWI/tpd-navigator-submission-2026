#!/usr/bin/env python3
"""Prepare hash-bound before/after ligand contact profiles against one protein.

No network retrieval, docking, efficacy scoring, or Probe/cctbx equivalence is
performed. The supplied ligand poses must already contain exactly one conformer.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.chemical_states import (
    before_after_report,
    interaction_profile,
    optimize_ligand_hydrogens,
    read_prepared_protein,
    run_pdb2pqr,
)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _one_sdf(path: Path) -> Chem.Mol:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"regular SDF required: {path}")
    with path.open("rb") as stream:
        molecules = [
            molecule
            for molecule in Chem.ForwardSDMolSupplier(stream, removeHs=False)
            if molecule is not None
        ]
    if len(molecules) != 1 or molecules[0].GetNumConformers() != 1:
        raise ValueError(f"exactly one valid single-conformer molecule required: {path}")
    return molecules[0]


def _coordinate_text(value: float) -> str:
    rounded = round(float(value), 4)
    if rounded == 0.0:
        rounded = 0.0
    return f"{rounded:.4f}"


def _molecule_hash_payload(molecule: Chem.Mol) -> bytes:
    if molecule.GetNumConformers() != 1:
        raise ValueError("exactly one conformer is required for a pose hash")
    conformer = molecule.GetConformer()
    atoms: list[dict[str, Any]] = []
    for atom in molecule.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        atoms.append({
            "atomic_number": atom.GetAtomicNum(),
            "formal_charge": atom.GetFormalCharge(),
            "isotope": atom.GetIsotope(),
            "aromatic": atom.GetIsAromatic(),
            "coordinates_A": [
                _coordinate_text(position.x),
                _coordinate_text(position.y),
                _coordinate_text(position.z),
            ],
        })
    bonds = [
        {
            "begin": bond.GetBeginAtomIdx(),
            "end": bond.GetEndAtomIdx(),
            "type": str(bond.GetBondType()),
            "aromatic": bond.GetIsAromatic(),
        }
        for bond in molecule.GetBonds()
    ]
    payload = {"atoms": atoms, "bonds": bonds}
    return json.dumps(
        payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _molecule_sha256(molecule: Chem.Mol) -> str:
    return _sha256_bytes(_molecule_hash_payload(molecule))


def _heavy_graph_key(molecule: Chem.Mol) -> str:
    graph = Chem.RWMol()
    atom_mapping: dict[int, int] = {}
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        graph_atom = Chem.Atom(atom.GetAtomicNum())
        graph_atom.SetFormalCharge(0)
        graph_atom.SetNoImplicit(True)
        graph_atom.SetNumExplicitHs(0)
        graph_atom.SetIsAromatic(atom.GetIsAromatic())
        atom_mapping[atom.GetIdx()] = graph.AddAtom(graph_atom)
    for bond in molecule.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin not in atom_mapping or end not in atom_mapping:
            continue
        graph.AddBond(atom_mapping[begin], atom_mapping[end], bond.GetBondType())
        added_bond = graph.GetBondBetweenAtoms(atom_mapping[begin], atom_mapping[end])
        added_bond.SetIsAromatic(bond.GetIsAromatic())
    return Chem.MolToSmiles(
        graph.GetMol(), canonical=True, isomericSmiles=False, allHsExplicit=True
    )


def _chemical_state_key(molecule: Chem.Mol) -> str:
    normalized = Chem.RemoveHs(Chem.Mol(molecule))
    for atom in normalized.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(
        normalized, canonical=True, isomericSmiles=True, allHsExplicit=True
    )


def _state_summary(molecule: Chem.Mol) -> dict[str, Any]:
    state_key = _chemical_state_key(molecule)
    return {
        "explicit_hydrogen_atom_count": sum(
            1 for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 1
        ),
        "total_formal_charge": sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()),
        "chemical_state_sha256": _sha256_bytes(state_key.encode("utf-8")),
    }


def _parent_identity(before: Chem.Mol, after: Chem.Mol) -> dict[str, Any]:
    before_key = _heavy_graph_key(before)
    after_key = _heavy_graph_key(after)
    same_heavy_graph = before_key == after_key
    before_state = _state_summary(before)
    after_state = _state_summary(after)
    chemical_state_changed = (
        same_heavy_graph
        and before_state["chemical_state_sha256"]
        != after_state["chemical_state_sha256"]
    )
    formal_charge_changed = (
        same_heavy_graph
        and before_state["total_formal_charge"]
        != after_state["total_formal_charge"]
    )
    return {
        "same_parent_heavy_graph": same_heavy_graph,
        "same_parent_compound": same_heavy_graph,
        "scope": "constitutional heavy graph only; stereochemistry not established",
        "stereo_identity_status": "not_established",
        "criterion": "Constitutional heavy-atom elemental connectivity and bond identity only.",
        "before_heavy_graph_sha256": _sha256_bytes(before_key.encode("utf-8")),
        "after_heavy_graph_sha256": _sha256_bytes(after_key.encode("utf-8")),
        "chemical_state_changed": chemical_state_changed,
        "protonation_or_formal_charge_state_changed": formal_charge_changed,
        "before_state": before_state,
        "after_state": after_state,
    }


def _write_sdf(path: Path, molecules: list[tuple[str, Chem.Mol]]) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        writer = Chem.SDWriter(stream)
        try:
            for name, molecule in molecules:
                copy = Chem.Mol(molecule)
                copy.SetProp("_Name", name)
                writer.write(copy)
            writer.flush()
        finally:
            writer.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before-sdf", required=True, type=Path)
    parser.add_argument("--after-sdf", required=True, type=Path)
    parser.add_argument("--protein-pdb", required=True, type=Path)
    parser.add_argument("--pdb2pqr", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ph", type=float, default=7.4)
    args = parser.parse_args()

    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    output.mkdir(parents=True)

    before_source_sha256 = _file_sha256(args.before_sdf)
    after_source_sha256 = _file_sha256(args.after_sdf)
    protein_source_sha256 = _file_sha256(args.protein_pdb)
    before_raw = _one_sdf(args.before_sdf)
    after_raw = _one_sdf(args.after_sdf)
    parent_identity = _parent_identity(before_raw, after_raw)
    after_label = (
        "parent_redock"
        if parent_identity["same_parent_heavy_graph"]
        else "different_compound_pose"
    )

    receipt = run_pdb2pqr(
        args.protein_pdb,
        output / "protein_preparation",
        args.pdb2pqr,
        pH=args.ph,
        preserve_heavy=True,
    )
    mapping = receipt["heavy_atom_mapping"]
    if mapping["missing_from_output"] or mapping["maximum_displacement_A"] != 0.0:
        raise ValueError("PROTEIN_HEAVY_ATOM_DRIFT_PROHIBITS_CONTACT_COMPARISON")
    prepared = read_prepared_protein(receipt)

    before_structures, before_h_receipt = optimize_ligand_hydrogens(before_raw)
    after_structures, after_h_receipt = optimize_ligand_hydrogens(after_raw)
    before_ligand = before_structures["hydrogens_optimized"]
    after_ligand = after_structures["hydrogens_optimized"]
    if before_ligand.GetNumConformers() != 1 or after_ligand.GetNumConformers() != 1:
        raise ValueError("hydrogen optimization must return exactly one conformer")

    before_raw_profile = interaction_profile(
        before_raw, prepared["protein_atoms"], prepared["protein_hydrogens"]
    )
    after_raw_profile = interaction_profile(
        after_raw, prepared["protein_atoms"], prepared["protein_hydrogens"]
    )
    before = interaction_profile(
        before_ligand, prepared["protein_atoms"], prepared["protein_hydrogens"]
    )
    after = interaction_profile(
        after_ligand, prepared["protein_atoms"], prepared["protein_hydrogens"]
    )
    comparison = before_after_report(before, after)
    raw_comparison = before_after_report(before_raw_profile, after_raw_profile)
    raw_hydrogen_uncertainty = {
        "before_crystal": sum(1 for atom in before_raw.GetAtoms() if atom.GetAtomicNum() == 1) == 0,
        after_label: sum(1 for atom in after_raw.GetAtoms() if atom.GetAtomicNum() == 1) == 0,
        "meaning": "True means raw-hydrogen-dependent contacts are unknown, not absent or negative.",
    }

    _write_sdf(
        output / "ligand_coordinates.sdf",
        [
            ("before_crystal_raw", before_raw),
            ("before_crystal_hydrogens_optimized", before_ligand),
            (f"{after_label}_raw", after_raw),
            (f"{after_label}_hydrogens_optimized", after_ligand),
        ],
    )
    ligand_metadata = {
        "source_paths": {
            "before_sdf": args.before_sdf.as_posix(),
            "after_sdf": args.after_sdf.as_posix(),
            "protein_pdb": args.protein_pdb.as_posix(),
        },
        "source_sha256": {
            "before_sdf": before_source_sha256,
            "after_sdf": after_source_sha256,
            "protein_pdb": protein_source_sha256,
        },
        "before_crystal": {
            "source_sdf_sha256": before_source_sha256,
            "raw_pose_sha256": _molecule_sha256(before_raw),
            "hydrogens_optimized_pose_sha256": _molecule_sha256(before_ligand),
        },
        after_label: {
            "source_sdf_sha256": after_source_sha256,
            "raw_pose_sha256": _molecule_sha256(after_raw),
            "hydrogens_optimized_pose_sha256": _molecule_sha256(after_ligand),
        },
    }
    document = {
        "format": "fixed-heavy-contact-comparison/1.0",
        "protein_preparation": receipt,
        "prepared_protein_uncertainties": prepared["uncertainties"],
        "parent_compound_identity": parent_identity,
        "ligand_metadata": ligand_metadata,
        "before_ligand_hydrogen_receipt": before_h_receipt,
        "after_ligand_hydrogen_receipt": after_h_receipt,
        "hydrogen_effects": {
            "before_crystal": {
                "raw_profile": before_raw_profile,
                "optimized_profile": before,
                "raw_vs_optimized_report": before_after_report(before_raw_profile, before),
            },
            after_label: {
                "raw_profile": after_raw_profile,
                "optimized_profile": after,
                "raw_vs_optimized_report": before_after_report(after_raw_profile, after),
            },
            "before_after_report": {
                "raw": raw_comparison,
                "optimized": comparison,
            },
            "raw_missing_hydrogen_uncertainty": raw_hydrogen_uncertainty,
            "interpretation": "Raw missing hydrogens make hydrogen-dependent contacts unknown; reports are separate and are not averaged.",
        },
        "before_profile": before,
        "after_profile": after,
        "comparison": comparison,
        "actual_interaction_counts": {
            "before": len(before["interactions"]),
            "after": len(after["interactions"]),
            "lost": len(comparison["lost_identifiers"]),
            "gained": len(comparison["gained_identifiers"]),
            "retained": len(comparison["retained_identifiers"]),
        },
        "scope": "Generic supplied pose pair; no claim of docking, efficacy, Probe, or cctbx equivalence.",
    }
    with (output / "contact_comparison.json").open(
        "x", encoding="utf-8", newline="\n"
    ) as stream:
        json.dump(document, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")

    for path in output.rglob("*"):
        if path.is_file():
            path.chmod(0o444)
    for path in sorted(
        (item for item in output.rglob("*") if item.is_dir()), reverse=True
    ):
        path.chmod(0o555)
    output.chmod(0o555)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
