"""Atom-map-invariant tetrahedral stereochemistry utilities.

RDKit tetrahedral CW/CCW tags are relative to an atom's neighbor ordering.  They
must therefore be normalized before molecules produced by canonical SMILES,
atom renumbering, or other graph round trips are compared.
"""
from __future__ import annotations

from typing import Any

from rdkit import Chem


_VIRTUAL_HYDROGEN = (1, 0)


class MappedStereoError(ValueError):
    """The mapped atom cannot be compared by the supported stereo policy."""


def _permutation_parity(values: list[tuple[int, int]]) -> int:
    ordered = sorted(values)
    positions = {value: index for index, value in enumerate(ordered)}
    permutation = [positions[value] for value in values]
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return inversions & 1


def _mapped_ligands(atom: Chem.Atom) -> list[tuple[int, int]]:
    ligands: list[tuple[int, int]] = []
    explicit_hydrogen_neighbors = 0
    seen_maps: set[int] = set()
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 1:
            explicit_hydrogen_neighbors += 1
            ligands.append(_VIRTUAL_HYDROGEN)
            continue
        atom_map = neighbor.GetAtomMapNum()
        if atom_map <= 0:
            raise MappedStereoError("TETRAHEDRAL_NEIGHBOR_MAP_MISSING")
        if atom_map in seen_maps:
            raise MappedStereoError("TETRAHEDRAL_NEIGHBOR_MAP_DUPLICATE")
        seen_maps.add(atom_map)
        ligands.append((0, atom_map))

    attached_hydrogens = atom.GetNumImplicitHs() + atom.GetNumExplicitHs()
    for _ in range(attached_hydrogens):
        ligands.append(_VIRTUAL_HYDROGEN)
    if explicit_hydrogen_neighbors + attached_hydrogens > 1:
        raise MappedStereoError("TETRAHEDRAL_MULTIPLE_HYDROGEN_LIGANDS")
    if len(ligands) != 4 or len(set(ligands)) != 4:
        raise MappedStereoError("TETRAHEDRAL_LIGAND_COUNT_UNSUPPORTED")
    return ligands


def mapped_tetrahedral_parity(atom: Chem.Atom) -> tuple[Any, ...]:
    """Return an atom-map-ordered tetrahedral signature.

    Unspecified stereochemistry remains explicitly distinct. Defined CW/CCW
    tags are combined with the parity of the current neighbor order relative to
    atom-map order. One attached hydrogen is represented by a stable virtual
    marker whether RDKit stores it implicitly, as an atom property, or as an
    explicit hydrogen neighbor.

    Non-tetrahedral chiral types and malformed mapped tetrahedral centers raise
    ``MappedStereoError`` so callers can fail closed.
    """
    tag = atom.GetChiralTag()
    if tag == Chem.ChiralType.CHI_UNSPECIFIED:
        return ("UNSPECIFIED",)
    if tag not in {
        Chem.ChiralType.CHI_TETRAHEDRAL_CW,
        Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
    }:
        raise MappedStereoError("NON_TETRAHEDRAL_STEREO_UNSUPPORTED")

    ligands = _mapped_ligands(atom)
    tag_parity = 0 if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW else 1
    normalized_parity = tag_parity ^ _permutation_parity(ligands)
    return ("TETRAHEDRAL", tuple(sorted(ligands)), normalized_parity)


def same_mapped_tetrahedral_stereo(before: Chem.Atom,
                                    after: Chem.Atom) -> bool:
    """Compare mapped tetrahedral stereo, returning false on unsupported data.

    Defined centers include their mapped ligand set in the signature. A changed
    stereo-neighbor set is therefore rejected rather than inferred to preserve
    geometry. Unspecified atoms compare only as unspecified, allowing ordinary
    edits at atoms where no stereochemistry was asserted.
    """
    try:
        return mapped_tetrahedral_parity(before) == mapped_tetrahedral_parity(after)
    except MappedStereoError:
        return False
