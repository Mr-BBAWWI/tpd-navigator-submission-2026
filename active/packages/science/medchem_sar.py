"""Exact selected-parent linkage to hash-verified medicinal-chemistry evidence.

This module exposes measured records and informational source-pair SAR only. It does
not activate SAR sites, protected atoms, generation rules, or claim measurements for
new analogues.
"""
from __future__ import annotations

import copy

from rdkit import Chem


FORMAT = "selected-parent-medchem-linkage/1.0"
MAX_FULL_ISOMORPHISMS = 100000


def _heavy_graph(mol):
    """Return an H-free/map-free graph while preserving every formal charge."""
    if mol is None:
        raise ValueError("PARENT_MOLECULE_REQUIRED")
    graph = Chem.RemoveHs(Chem.Mol(mol), sanitize=True)
    for atom in graph.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(graph)
    return graph


def _connectivity_key(mol):
    graph = _heavy_graph(mol)
    return Chem.MolToSmiles(graph, canonical=True, isomericSmiles=False)


def _formal_charges(mol):
    return sorted(atom.GetFormalCharge() for atom in _heavy_graph(mol).GetAtoms())


def _mapped_source_molecule(smiles, name):
    """Mirror prepare-medchem-evidence's deterministic mapped_molecule procedure."""
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("INVALID_MEDCHEM_SOURCE_SMILES")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=True, includeChirality=True))
    order = sorted(range(mol.GetNumAtoms()), key=lambda index: (ranks[index], index))
    old_to_map = {old: new + 1 for new, old in enumerate(order)}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() > 0:
            atom.SetAtomMapNum(old_to_map[atom.GetIdx()])
    mol.SetProp("_Name", name)
    mol.SetProp("source_smiles", smiles)
    mol.SetProp("atom_map_definition", "canonical-rank heavy-atom map; 1-based")
    return mol


def _records(data):
    catalog = data.get("parent_catalog", {}) if isinstance(data, dict) else {}
    records = catalog.get("records", []) if isinstance(catalog, dict) else []
    return records if isinstance(records, list) else []


def _full_isomorphisms(source, parent):
    source_graph = _heavy_graph(source)
    parent_graph = _heavy_graph(parent)
    if source_graph.GetNumAtoms() != parent_graph.GetNumAtoms():
        return [], None
    # Full graph matching uses elements, aromaticity, bond orders and formal charges.
    matches = parent_graph.GetSubstructMatches(
        source_graph,
        uniquify=False,
        useChirality=False,
        maxMatches=MAX_FULL_ISOMORPHISMS,
    )
    # RDKit truncates at maxMatches. Reaching the cap cannot establish complete
    # all-mapping consensus, even when the true count happens to equal the cap.
    if len(matches) >= MAX_FULL_ISOMORPHISMS:
        return None, "FULL_GRAPH_ISOMORPHISM_LIMIT_REACHED"
    return [
        match for match in matches if len(match) == parent_graph.GetNumAtoms()
    ], None


def _source_site_consensus(pair, source_side):
    possibilities = pair.get("mapping_possibilities")
    if not isinstance(possibilities, list) or not possibilities:
        return None
    key = (
        "right_affected_atom_maps"
        if source_side == "right"
        else "left_affected_parent_atom_maps"
    )
    values = []
    for possibility in possibilities:
        sites = possibility.get(key)
        if not isinstance(sites, list) or not sites or any(type(v) is not int for v in sites):
            return None
        values.append(tuple(sorted(set(sites))))
    return list(values[0]) if all(value == values[0] for value in values) else None


def _parent_site_consensus(source, parent, source_maps):
    # Build the indexed correspondence arrays independently of the map-free
    # identity graphs. Removing H from both molecules first keeps mapping tuple
    # indices aligned even when explicit hydrogens were interspersed originally.
    source_heavy = Chem.RemoveHs(Chem.Mol(source), sanitize=True)
    parent_heavy = Chem.RemoveHs(Chem.Mol(parent), sanitize=True)
    Chem.SanitizeMol(source_heavy)
    Chem.SanitizeMol(parent_heavy)

    mappings, reason = _full_isomorphisms(source_heavy, parent_heavy)
    if reason is not None:
        return None, None, reason
    if not mappings:
        return None, 0, None
    source_index_by_map = {
        atom.GetAtomMapNum(): atom.GetIdx()
        for atom in source_heavy.GetAtoms()
        if atom.GetAtomicNum() > 1
    }
    if any(value not in source_index_by_map for value in source_maps):
        return None, len(mappings), None
    mapped_sets = []
    for mapping in mappings:
        parent_maps = []
        for source_map in source_maps:
            parent_atom = parent_heavy.GetAtomWithIdx(
                mapping[source_index_by_map[source_map]]
            )
            parent_map = parent_atom.GetAtomMapNum()
            if parent_map <= 0:
                return None, len(mappings), None
            parent_maps.append(parent_map)
        mapped_sets.append(tuple(sorted(set(parent_maps))))
    if not all(value == mapped_sets[0] for value in mapped_sets):
        return None, len(mappings), None
    return list(mapped_sets[0]), len(mappings), None


def link_selected_parent(parent, parent_id, medchem_source):
    """Link an actual selected parent to exactly one measured source graph."""
    unavailable = {
        "format": FORMAT,
        "status": "not_joined",
        "selected_parent_id": parent_id,
        "identity_policy": "exact neutral heavy-atom connectivity after removeHs/atom-map removal only; formal charges must match",
        "matched_record": None,
        "measurements": [],
        "source_sar": [],
        "new_analog_measurement_status": "not_proved_for_any_new_analog",
        "changes_design_policy": False,
    }
    if not isinstance(medchem_source, dict) or medchem_source.get("status") != "configured_hash_verified":
        unavailable["reason"] = "MEDCHEM_SOURCE_NOT_HASH_VERIFIED"
        return unavailable
    data = medchem_source.get("data")
    if not isinstance(data, dict):
        unavailable["reason"] = "MEDCHEM_SOURCE_SCHEMA_UNSUPPORTED"
        return unavailable

    parent_graph = _heavy_graph(parent)
    # Identity intentionally uses the map-free graph, while correspondence uses
    # an H-free copy that retains the selected parent's atom maps.
    parent_for_correspondence = Chem.RemoveHs(Chem.Mol(parent), sanitize=True)
    Chem.SanitizeMol(parent_for_correspondence)
    parent_key = _connectivity_key(parent_graph)
    parent_charge = int(Chem.GetFormalCharge(parent_graph))
    matches = []
    source_molecules = {}
    for record in _records(data):
        identity = record.get("identity", {})
        smiles = identity.get("source_smiles") if isinstance(identity, dict) else None
        if not isinstance(smiles, str):
            continue
        try:
            source_mol = _mapped_source_molecule(smiles, str(record.get("compound_id", "source")))
        except ValueError:
            continue
        if (
            _connectivity_key(source_mol) == parent_key
            and int(Chem.GetFormalCharge(_heavy_graph(source_mol))) == parent_charge
            and _formal_charges(source_mol) == _formal_charges(parent_graph)
        ):
            isomorphisms, mapping_reason = _full_isomorphisms(
                source_mol, parent_graph
            )
            if mapping_reason is not None:
                unavailable["reason"] = mapping_reason
                return unavailable
            if isomorphisms:
                matches.append(record)
                source_molecules[record.get("compound_id")] = source_mol
    if len(matches) != 1:
        unavailable["reason"] = "NO_EXACT_GRAPH_AND_CHARGE_MATCH" if not matches else "MULTIPLE_EXACT_SOURCE_MATCHES"
        unavailable["exact_match_count"] = len(matches)
        return unavailable

    matched = matches[0]
    compound_id = matched.get("compound_id")
    source_mol = source_molecules[compound_id]
    measured = copy.deepcopy(matched.get("measurements", []))
    source_sar = []
    for pair in data.get("sar_pair_evidence", []):
        if not isinstance(pair, dict):
            continue
        if pair.get("left_compound_id") == compound_id:
            side = "left"
        elif pair.get("right_compound_id") == compound_id:
            side = "right"
        else:
            continue
        source_sites = _source_site_consensus(pair, side)
        if source_sites is None:
            continue
        parent_sites, isomorphism_count, mapping_reason = _parent_site_consensus(
            source_mol, parent_for_correspondence, source_sites
        )
        if mapping_reason is not None or parent_sites is None:
            continue
        source_sar.append({
            "pair_id": pair.get("pair_id"),
            "requested_change_label": pair.get("requested_change_label"),
            "selected_parent_pair_side": side,
            "other_compound_id": pair.get("right_compound_id") if side == "left" else pair.get("left_compound_id"),
            "selected_parent_measurements": copy.deepcopy(pair.get(side + "_measurements", [])),
            "other_measurements": copy.deepcopy(pair.get(("right" if side == "left" else "left") + "_measurements", [])),
            "source_row_locators": copy.deepcopy(pair.get("source_row_locators", [])),
            "source_affected_atom_maps": source_sites,
            "selected_parent_affected_atom_maps": parent_sites,
            "source_site_consensus": True,
            "full_graph_mapping_unambiguous": isomorphism_count == 1,
            "full_graph_isomorphism_count": isomorphism_count,
            "status": "informational_exact_observed_pair_only",
            "interpretation": "Observed source pair only; it does not allow arbitrary transformations or activate a design site.",
            "protected_atom_activation": False,
        })

    return {
        **unavailable,
        "status": "exact_measured_parent_join",
        "reason": None,
        "exact_match_count": 1,
        "matched_record": {
            "compound_id": compound_id,
            "record_class": matched.get("record_class"),
            "identity": copy.deepcopy(matched.get("identity")),
            "source": copy.deepcopy(matched.get("source")),
            "identity_status": "connectivity_matched_stereochemistry_not_established",
            "formal_charge": parent_charge,
        },
        "measurements": measured,
        "source_sar": source_sar,
    }


__all__ = ["FORMAT", "link_selected_parent"]
