"""Structure-only inexpensive filters and enumeration coverage diagnostics.

These checks are developer assumptions for triage before docking. They are not
medicinal-chemistry verdicts and do not predict activity, safety, or synthesis.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any, Iterable

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

from packages.science.mapped_stereo import (
    MappedStereoError,
    mapped_tetrahedral_parity,
    same_mapped_tetrahedral_stereo,
)


_DEVELOPER_THRESHOLDS = {
    "absolute_MW_delta_soft_A": 180.0,
    "absolute_TPSA_delta_soft_A2": 70.0,
    "absolute_clogP_delta_soft": 2.0,
    "product_MW_soft": 650.0,
}

_ALERTS = {
    "peroxide": Chem.MolFromSmarts("[O;X2]-[O;X2]"),
    "N_N_bond": Chem.MolFromSmarts("[N,n]-[N,n]"),
    "nitro": Chem.MolFromSmarts("[N+](=O)[O-]"),
    "aldehyde": Chem.MolFromSmarts("[CX3H1](=O)[#6,#7,#8,#16]"),
    "Michael_acceptor": Chem.MolFromSmarts("[C,c]=[C,c]-[C,S](=O)"),
    "azide": Chem.MolFromSmarts("[$([N-]=[N+]=N),$([N]=[N+]=[N-])]"),
}


def _mapped_index(mol: Chem.Mol) -> tuple[dict[int, int], list[int]]:
    result: dict[int, int] = {}
    duplicates: list[int] = []
    for atom in mol.GetAtoms():
        atom_map = atom.GetAtomMapNum()
        if atom_map <= 0:
            continue
        if atom_map in result:
            duplicates.append(atom_map)
        result[atom_map] = atom.GetIdx()
    return result, sorted(set(duplicates))


def _protected_atom_signature(atom: Chem.Atom) -> tuple[Any, ...]:
    return (
        atom.GetAtomicNum(), atom.GetIsotope(), atom.GetFormalCharge(),
        atom.GetIsAromatic(), atom.GetTotalNumHs(includeNeighbors=True),
        mapped_tetrahedral_parity(atom),
        sorted((neighbor.GetAtomMapNum(), neighbor.GetAtomicNum()) for neighbor in atom.GetNeighbors()),
    )


def _protected_graph_reasons(parent: Chem.Mol, product: Chem.Mol,
                             protected_maps: set[int]) -> list[str]:
    reasons: list[str] = []
    pidx, _ = _mapped_index(parent)
    qidx, _ = _mapped_index(product)
    for atom_map in sorted(protected_maps):
        if atom_map not in pidx:
            reasons.append(f"PROTECTED_MAP_ABSENT_FROM_PARENT:{atom_map}")
        elif atom_map not in qidx:
            reasons.append(f"PROTECTED_MAP_REMOVED:{atom_map}")
        else:
            try:
                changed = _protected_atom_signature(
                    parent.GetAtomWithIdx(pidx[atom_map])
                ) != _protected_atom_signature(product.GetAtomWithIdx(qidx[atom_map]))
            except MappedStereoError:
                changed = True
            if changed:
                reasons.append(f"PROTECTED_ATOM_GRAPH_CHANGED:{atom_map}")
    for bond in parent.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if (a not in protected_maps and b not in protected_maps) or a not in qidx or b not in qidx:
            continue
        other = product.GetBondBetweenAtoms(qidx[a], qidx[b])
        if other is None or other.GetBondType() != bond.GetBondType() or int(other.GetStereo()) != int(bond.GetStereo()):
            reasons.append(f"PROTECTED_BOND_CHANGED:{min(a,b)}-{max(a,b)}")
    return reasons


def _carbon_stereo_reasons(parent: Chem.Mol, product: Chem.Mol) -> list[str]:
    reasons: list[str] = []
    pidx, _ = _mapped_index(parent)
    qidx, _ = _mapped_index(product)
    for atom_map in sorted(set(pidx) & set(qidx)):
        before = parent.GetAtomWithIdx(pidx[atom_map])
        after = product.GetAtomWithIdx(qidx[atom_map])
        if before.GetAtomicNum() == 6 and not same_mapped_tetrahedral_stereo(before, after):
            reasons.append(f"PARENT_CARBON_STEREO_CHANGED:{atom_map}")
    for bond in parent.GetBonds():
        if int(bond.GetStereo()) == int(Chem.BondStereo.STEREONONE):
            continue
        a, b = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
        if a not in qidx or b not in qidx:
            reasons.append(f"PARENT_STEREOBOND_REMOVED:{min(a,b)}-{max(a,b)}")
            continue
        other = product.GetBondBetweenAtoms(qidx[a], qidx[b])
        if other is None or int(other.GetStereo()) != int(bond.GetStereo()):
            reasons.append(f"PARENT_STEREOBOND_CHANGED:{min(a,b)}-{max(a,b)}")
    return reasons


def _properties(mol: Chem.Mol) -> dict[str, float]:
    return {
        "MW": float(Descriptors.MolWt(mol)),
        "TPSA_A2": float(rdMolDescriptors.CalcTPSA(mol)),
        "clogP": float(Crippen.MolLogP(mol)),
    }


def _sa_score(mol: Chem.Mol) -> tuple[float | None, str]:
    try:
        from rdkit.Contrib.SA_Score import sascorer
        return float(sascorer.calculateScore(mol)), "RDKit Contrib SA_Score"
    except Exception:
        return None, "unavailable; no substitute score inferred"


def cheap_filter(parent: Chem.Mol, record: dict[str, Any]) -> dict[str, Any]:
    """Validate one generated graph and report transparent structural risk flags."""
    hard: list[str] = []
    soft: list[dict[str, Any]] = []
    if not isinstance(parent, Chem.Mol):
        raise TypeError("parent must be an RDKit Mol")
    if not isinstance(record, dict):
        raise TypeError("record must be a dictionary")
    try:
        parent_copy = Chem.Mol(parent)
        Chem.SanitizeMol(parent_copy)
    except Exception as exc:
        raise ValueError("parent does not sanitize") from exc

    smiles = record.get("mapped_smiles")
    product = Chem.MolFromSmiles(str(smiles), sanitize=False) if smiles else None
    if product is None:
        return {
            "valid": False,
            "hard_reasons": ["PRODUCT_PARSE_FAILED"],
            "soft_flags": [],
            "developer_assumptions": dict(_DEVELOPER_THRESHOLDS),
        }
    try:
        Chem.SanitizeMol(product)
    except Exception as exc:
        hard.append("PRODUCT_SANITIZATION_OR_VALENCE_FAILED:" + (str(exc) or exc.__class__.__name__))
        return {
            "valid": False,
            "hard_reasons": hard,
            "soft_flags": [],
            "developer_assumptions": dict(_DEVELOPER_THRESHOLDS),
        }

    if any(atom.GetAtomicNum() == 0 for atom in product.GetAtoms()):
        hard.append("DUMMY_ATOM_PRESENT")
    if len(Chem.GetMolFrags(product)) != 1:
        hard.append("DISCONNECTED_PRODUCT")

    pidx, parent_duplicates = _mapped_index(parent_copy)
    qidx, product_duplicates = _mapped_index(product)
    if any(atom.GetAtomicNum() > 0 and atom.GetAtomMapNum() <= 0 for atom in parent_copy.GetAtoms()):
        hard.append("PARENT_ATOM_MAP_MISSING_OR_ZERO")
    if any(atom.GetAtomicNum() > 0 and atom.GetAtomMapNum() <= 0 for atom in product.GetAtoms()):
        hard.append("PRODUCT_ATOM_MAP_MISSING_OR_ZERO")
    if parent_duplicates:
        hard.append("DUPLICATE_PARENT_MAPS")
    if product_duplicates:
        hard.append("DUPLICATE_PRODUCT_MAPS")
    parent_maps, product_maps = set(pidx), set(qidx)
    declared_removed = {int(value) for value in record.get("removed_atom_maps", [])}
    declared_added = {int(value) for value in record.get("added_atom_maps", [])}
    if parent_maps - product_maps != declared_removed:
        hard.append("REMOVED_MAP_IDENTITY_MISMATCH")
    if product_maps - parent_maps != declared_added:
        hard.append("ADDED_MAP_IDENTITY_MISMATCH")
    if parent_maps & declared_added:
        hard.append("ADDED_MAP_COLLIDES_WITH_PARENT")

    protected = {int(value) for value in record.get("protected_atom_maps", [])}
    hard.extend(_protected_graph_reasons(parent_copy, product, protected))
    hard.extend(_carbon_stereo_reasons(parent_copy, product))
    if hard:
        return {
            "valid": False,
            "hard_reasons": hard,
            "soft_flags": [],
            "developer_assumptions": dict(_DEVELOPER_THRESHOLDS),
        }

    for name, query in _ALERTS.items():
        if query is None or not product.HasSubstructMatch(query):
            continue
        if name == "azide":
            role = "introduced_click_handle" if record.get("rule_id") == "LH_AZIDE" else "unassigned_azide"
            soft.append({"flag": "REACTIVE_ALERT_AZIDE", "role": role})
        elif name == "N_N_bond" and product.HasSubstructMatch(_ALERTS["azide"]):
            continue
        else:
            soft.append({"flag": "REACTIVE_ALERT_" + name.upper(), "role": "structural_risk_for_review"})

    parent_properties = _properties(parent_copy)
    product_properties = _properties(product)
    deltas = {key: product_properties[key] - parent_properties[key] for key in parent_properties}
    if abs(deltas["MW"]) > _DEVELOPER_THRESHOLDS["absolute_MW_delta_soft_A"]:
        soft.append({"flag": "LARGE_MW_DELTA", "value": deltas["MW"]})
    if abs(deltas["TPSA_A2"]) > _DEVELOPER_THRESHOLDS["absolute_TPSA_delta_soft_A2"]:
        soft.append({"flag": "LARGE_TPSA_DELTA", "value": deltas["TPSA_A2"]})
    if abs(deltas["clogP"]) > _DEVELOPER_THRESHOLDS["absolute_clogP_delta_soft"]:
        soft.append({"flag": "LARGE_CLOGP_DELTA", "value": deltas["clogP"]})
    if product_properties["MW"] > _DEVELOPER_THRESHOLDS["product_MW_soft"]:
        soft.append({"flag": "HIGH_PRODUCT_MW", "value": product_properties["MW"]})
    sa_score, sa_method = _sa_score(product)

    return {
        "valid": not hard,
        "hard_reasons": hard,
        "soft_flags": soft,
        "properties": {
            "parent": parent_properties,
            "product": product_properties,
            "delta_product_minus_parent": deltas,
            "synthetic_accessibility_score": sa_score,
            "synthetic_accessibility_method": sa_method,
        },
        "developer_assumptions": {
            **_DEVELOPER_THRESHOLDS,
            "interpretation": "Soft thresholds and alerts are triage assumptions, not medicinal-chemistry, safety, or activity verdicts.",
        },
    }


def _graph_without_maps_or_stereo(record: dict[str, Any]) -> str | None:
    smiles = record.get("mapped_smiles") or record.get("canonical_smiles")
    mol = Chem.MolFromSmiles(str(smiles)) if smiles else None
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in mol.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
        bond.SetBondDir(Chem.BondDir.NONE)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def count_by_site_and_broad_family(records: Iterable[dict[str, Any]],
                                   sites: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Count unique constitutional graphs and expose the 20-per-modifiable-site gap."""
    site_rows = {int(site["atom_map"]): site for site in sites}
    graphs_by_site_family: dict[int, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    accepted_graphs_by_site: dict[int, set[str]] = defaultdict(set)
    all_graphs: set[str] = set()
    invalid_graph_records = 0
    catalog_rule_ids: set[str] = set()
    applied_rule_ids: set[str] = set()

    for record in records:
        graph = _graph_without_maps_or_stereo(record)
        if graph is None:
            invalid_graph_records += 1
            continue
        all_graphs.add(graph)
        rule_id = str(record.get("rule_id") or "unknown")
        catalog_rule_ids.add(rule_id)
        filter_result = record.get("cheap_filter") or record.get("filter_result") or {}
        passed_filter = isinstance(filter_result, dict) and filter_result.get("valid") is True
        if passed_filter:
            applied_rule_ids.add(rule_id)
        family = str(record.get("transformation_class") or "unknown")
        if family in {"ring_expansion", "ring_contraction"}:
            family = "ring_modification"
        statuses = record.get("site_status") or []
        status_by_map = {
            int(row["atom_map"]): str(row.get("state", "UNKNOWN")).upper()
            for row in statuses if isinstance(row, dict) and "atom_map" in row
        }
        touched = record.get("attachment_site_atom_maps") or record.get("modified_atom_maps") or []
        for raw_map in touched:
            atom_map = int(raw_map)
            graphs_by_site_family[atom_map][family].add(graph)
            state = status_by_map.get(atom_map, str(site_rows.get(atom_map, {}).get("state", "UNKNOWN")).upper())
            if state == "MODIFIABLE" and passed_filter:
                accepted_graphs_by_site[atom_map].add(graph)

    counts = {
        str(atom_map): {family: len(graphs) for family, graphs in sorted(families.items())}
        for atom_map, families in sorted(graphs_by_site_family.items())
    }
    gaps = []
    modifiable_site_count = sum(
        str(site.get("state", "UNKNOWN")).upper() == "MODIFIABLE"
        for site in site_rows.values()
    )
    for atom_map, site in sorted(site_rows.items()):
        if str(site.get("state", "UNKNOWN")).upper() != "MODIFIABLE":
            continue
        count = len(accepted_graphs_by_site.get(atom_map, set()))
        gaps.append({
            "atom_map": atom_map,
            "required_unique_graphs": 20,
            "accepted_unique_graphs": count,
            "gap": max(0, 20 - count),
            "meets_requirement": count >= 20,
        })
    return {
        "counts_by_site_and_broad_family": counts,
        "unique_graph_count_without_maps_or_stereo": len(all_graphs),
        "acceptance_gaps_for_modifiable_sites": gaps,
        "unknown_sites_count_toward_acceptance": False,
        "invalid_graph_record_count": invalid_graph_records,
        "rule_coverage": {
            "catalog_rule_count": len(catalog_rule_ids),
            "catalog_rule_ids": sorted(catalog_rule_ids),
            "actually_applied_passing_rule_count": len(applied_rule_ids),
            "actually_applied_passing_rule_ids": sorted(applied_rule_ids),
        },
        "modifiable_site_count": modifiable_site_count,
        "requires_at_least_two_modifiable_sites": True,
        "meets_all_acceptance_requirements": (
            modifiable_site_count >= 2 and bool(gaps)
            and all(row["meets_requirement"] for row in gaps)
        ),
        "counting_policy": "Constitutional graphs are canonicalized without atom maps or stereochemistry; ring expansion and contraction are merged as ring_modification; acceptance includes only records whose cheap-filter result explicitly has valid=true.",
    }
