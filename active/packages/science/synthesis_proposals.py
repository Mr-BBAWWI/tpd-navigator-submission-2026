"""Candidate-specific graph disconnections and synthesis proposals without route approval."""
from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from typing import Any

from rdkit import Chem
from rdkit.Chem import Draw

from .synthesis_review import assess_synthesis

VERSION = "synthesis-proposals/20261002.1"


def _unmapped_canonical(mol: Chem.Mol) -> str:
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)


def _candidate_mol(
    candidate: dict[str, Any], *, expected_canonical: str | None = None
) -> tuple[Chem.Mol, dict[int, int]]:
    if not isinstance(candidate, dict):
        raise TypeError("candidate must be a dictionary")
    mol = Chem.MolFromSmiles(candidate.get("mapped_smiles", ""))
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("candidate mapped_smiles must be one valid graph")

    stated = Chem.MolFromSmiles(candidate.get("canonical_smiles", ""))
    if stated is None or len(Chem.GetMolFrags(stated)) != 1:
        raise ValueError("candidate canonical_smiles must be one valid graph")
    actual_canonical = _unmapped_canonical(mol)
    if _unmapped_canonical(stated) != actual_canonical:
        raise ValueError("candidate canonical_smiles does not match mapped_smiles")
    if expected_canonical is not None and actual_canonical != expected_canonical:
        raise ValueError("candidate graph disagrees with synthesis graph validation")

    mapping = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        number = atom.GetAtomMapNum()
        if number <= 0 or number in mapping:
            raise ValueError("all candidate heavy atoms require unique positive maps")
        mapping[number] = atom.GetIdx()
    return mol, mapping


def _cuts(candidate: dict[str, Any], mol: Chem.Mol, mapping: dict[int, int]) -> list[dict[str, Any]]:
    roles = candidate.get("atom_roles")
    if not isinstance(roles, dict) or set(roles) != {"warhead_maps", "linker_maps", "recruiter_maps"}:
        raise ValueError("exact atom role partition is required")
    role_sets = {}
    seen = set()
    for key in ("warhead_maps", "linker_maps", "recruiter_maps"):
        values = roles[key]
        if not isinstance(values, list) or any(type(x) is not int or x <= 0 for x in values):
            raise ValueError(f"{key} must be a positive integer list")
        current = set(values)
        if len(current) != len(values) or current & seen:
            raise ValueError("atom role maps must be unique and disjoint")
        seen |= current
        role_sets[key] = current
    if seen != set(mapping):
        raise ValueError("atom roles must exactly partition the mapped graph")
    metadata = candidate.get("attachment_metadata")
    if not isinstance(metadata, dict):
        raise ValueError("attachment_metadata is required")
    result = []
    for name, outer_key in (("warhead_linker_bond", "warhead_maps"), ("recruiter_linker_bond", "recruiter_maps")):
        record = metadata.get(name)
        if not isinstance(record, dict):
            raise ValueError(f"{name} is required")
        a, b = record.get("attachment_atom_map"), record.get("partner_atom_map")
        if type(a) is not int or type(b) is not int or a not in mapping or b not in mapping or a == b:
            raise ValueError(f"{name} has absent or ambiguous atom maps")
        if not ((a in role_sets[outer_key] and b in role_sets["linker_maps"]) or (b in role_sets[outer_key] and a in role_sets["linker_maps"])):
            raise ValueError(f"{name} does not cross the declared role boundary")
        bond = mol.GetBondBetweenAtoms(mapping[a], mapping[b])
        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE or record.get("bond_type") != "SINGLE":
            raise ValueError(f"{name} does not identify the exact actual SINGLE bond")
        result.append({"name": name, "bond_index": bond.GetIdx(), "atom_maps": [a, b], "outer_role": outer_key})
    if result[0]["bond_index"] == result[1]["bond_index"]:
        raise ValueError("attachment cuts must be distinct")
    return result


def _fragments(mol: Chem.Mol, cuts: list[dict[str, Any]]) -> tuple[Chem.Mol, list[dict[str, Any]]]:
    labels = [(cut["atom_maps"][0], cut["atom_maps"][1]) for cut in cuts]
    fragmented = Chem.FragmentOnBonds(mol, [cut["bond_index"] for cut in cuts], addDummies=True, dummyLabels=labels)
    Chem.AssignStereochemistry(fragmented, cleanIt=True, force=True)
    records = []
    for fragment in Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=False):
        try:
            Chem.SanitizeMol(fragment)
        except Exception:
            fragment.UpdatePropertyCache(strict=False)
        Chem.AssignStereochemistry(fragment, cleanIt=True, force=True)
        maps = sorted(atom.GetAtomMapNum() for atom in fragment.GetAtoms() if atom.GetAtomicNum() > 1)
        records.append({
            "mapped_smiles_with_attachment_dummies": Chem.MolToSmiles(fragment, canonical=True, isomericSmiles=True),
            "atom_maps": maps,
            "stereochemistry": [{"atom_map": atom.GetAtomMapNum(), "CIP": atom.GetProp("_CIPCode")} for atom in fragment.GetAtoms() if atom.HasProp("_CIPCode")],
            "attachment_dummies": [{"isotope_endpoint_map": atom.GetIsotope(), "neighbor_map": next(iter(atom.GetNeighbors())).GetAtomMapNum() if atom.GetDegree() == 1 else None} for atom in fragment.GetAtoms() if atom.GetAtomicNum() == 0],
        })
    records.sort(key=lambda row: row["atom_maps"])
    return fragmented, records


def _endpoint(mol: Chem.Mol, mapping: dict[int, int], cut: dict[str, Any], roles: dict[str, Any]) -> dict[str, Any]:
    a, b = cut["atom_maps"]
    linker = set(roles["linker_maps"])
    outer_map, linker_map = (b, a) if a in linker else (a, b)
    outer, link = mol.GetAtomWithIdx(mapping[outer_map]), mol.GetAtomWithIdx(mapping[linker_map])
    elements = f"{outer.GetSymbol()}-{link.GetSymbol()}"
    if outer.GetSymbol() == "N" and link.GetSymbol() == "C":
        chemical_class = "N-C_single_bond"
        proposal = "N-C bond-forming proposal, potentially alkylation-type; not an amide recipe"
    elif outer.GetSymbol() == "O" and link.GetSymbol() == "C":
        chemical_class = "O-C_single_bond"
        proposal = "O-C bond-forming proposal, potentially substitution/alkylation-type; graph edge is not a click reaction"
    elif {outer.GetSymbol(), link.GetSymbol()} == {"C", "N"}:
        chemical_class = "C-N_single_bond"
        proposal = "C-N bond-forming proposal; transformation class requires precursor-level review"
    else:
        chemical_class = f"{elements}_single_bond"
        proposal = f"{elements} graph-bond proposal; reaction class unresolved"
    carbonyl_neighbor = any(n.GetAtomicNum() == 6 and any(x.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(n.GetIdx(), x.GetIdx()).GetBondType() == Chem.BondType.DOUBLE for x in n.GetNeighbors()) for n in outer.GetNeighbors())
    return {"cut": cut["name"], "outer_atom_map": outer_map, "linker_atom_map": linker_map, "elements": elements, "chemical_class": chemical_class, "proposal": proposal, "outer_atom_is_internal_amide_or_carbamate_context": bool(outer.GetSymbol() == "N" and carbonyl_neighbor), "status": "graph_space_edge_not_completed_chemical_reaction"}


def _functional_flags(mol: Chem.Mol) -> list[dict[str, Any]]:
    patterns = {
        "ester": "[CX3](=[OX1])[OX2][#6]", "phenol": "[c][OH1]", "alcohol": "[C][OH1]",
        "free_amine": "[N;H1,H2;!$(N-C=O)]", "amide": "[NX3][CX3](=[OX1])",
        "potential_electrophile": "[C;X3](=[O,S])[Cl,Br,I]", "disulfide": "[SX2][SX2]",
    }
    result = []
    for name, smarts in patterns.items():
        query = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(query) if query is not None else ()
        if matches:
            result.append({"flag": name, "match_count": len(matches), "mapped_atom_ids": [sorted(mol.GetAtomWithIdx(i).GetAtomMapNum() for i in match if mol.GetAtomWithIdx(i).GetAtomicNum() > 1) for match in matches], "status": "manual_review"})
    active = [row for row in result if row["flag"] in {"phenol", "alcohol", "free_amine"}]
    if sum(row["match_count"] for row in active) > 1:
        result.append({"flag": "multiple_competing_phenol_amine_or_alcohol_sites", "match_count": sum(row["match_count"] for row in active), "mapped_atom_ids": [maps for row in active for maps in row["mapped_atom_ids"]], "status": "regioselectivity_and_protection_review"})
    if any(row["flag"] == "ester" for row in result):
        result.append({"flag": "ester_handle_policy", "status": "precursor_or_protection_only", "interpretation": "An ester may inform precursor/protection planning but cannot be counted directly as the final attachment reaction."})
    return result


def _sa_score(mol: Chem.Mol) -> dict[str, Any]:
    try:
        from rdkit.Contrib.SA_Score import sascorer
        value = float(sascorer.calculateScore(mol))
        return {"status": "available", "value": value, "method": "RDKit Contrib SA_Score sascorer", "interpretation": "screening descriptor only; not synthetic proof"}
    except (ImportError, ModuleNotFoundError, AttributeError):
        return {"status": "unavailable", "value": None, "method": "RDKit Contrib SA_Score not installed", "interpretation": "No fallback or fabricated score was used."}


def _source_suggestions(candidate: dict[str, Any], endpoints: list[dict[str, Any]]) -> dict[str, Any]:
    evidence = candidate.get("route_evidence")
    records = evidence.get("reusable_precedent_records", []) if isinstance(evidence, dict) else []
    endpoint_classes = {row.get("chemical_class") for row in endpoints}
    candidate_id = candidate.get("candidate_id")
    candidate_canonical = candidate.get("canonical_smiles")
    suggestions = []

    for record in records if isinstance(records, list) else []:
        if not isinstance(record, dict):
            continue
        source = record.get("source_record", record)
        if not isinstance(source, dict):
            continue

        template_class = record.get(
            "template_chemical_class", source.get("template_chemical_class")
        )
        reference = record.get(
            "graph_confirmed_reference", source.get("graph_confirmed_reference")
        )
        if template_class not in endpoint_classes or not isinstance(reference, dict):
            continue
        if (
            reference.get("confirmed") is not True
            or reference.get("scope") != "exact_candidate_graph"
            or reference.get("candidate_id") != candidate_id
            or reference.get("canonical_smiles") != candidate_canonical
        ):
            continue

        reference_mol = Chem.MolFromSmiles(reference.get("canonical_smiles", ""))
        candidate_mol = Chem.MolFromSmiles(candidate_canonical or "")
        if (
            reference_mol is None
            or candidate_mol is None
            or _unmapped_canonical(reference_mol)
            != _unmapped_canonical(candidate_mol)
        ):
            continue

        doi = source.get("doi")
        locator = source.get("locator", record.get("locator"))
        if (
            isinstance(doi, str)
            and doi.startswith("10.")
            and isinstance(locator, str)
            and locator.strip()
        ):
            suggestions.append(
                {
                    "doi": doi,
                    "locator": locator,
                    "template_chemical_class": template_class,
                    "status": "graph_scoped_source_template_suggestion_only",
                    "conditions_transferred": False,
                }
            )

    return {
        "status": (
            "compatible_precise_sources_found"
            if suggestions
            else "unknown_no_graph_confirmed_chemically_compatible_template"
        ),
        "suggestions": suggestions,
        "invented_conditions": False,
        "invented_reagents": False,
        "invented_yield": False,
    }


def scheme_svg(candidate: dict[str, Any]) -> str:
    base = assess_synthesis(candidate)
    mol, mapping = _candidate_mol(
        candidate,
        expected_canonical=base["graph_validation"]["canonical_smiles"],
    )
    cuts = _cuts(candidate, mol, mapping)
    fragmented, _ = _fragments(mol, cuts)
    fragments = list(Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=False))
    for fragment in fragments:
        for atom in fragment.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetProp("atomLabel", f"cut:{atom.GetIsotope()}")
            elif atom.GetAtomMapNum() in {value for cut in cuts for value in cut["atom_maps"]}:
                atom.SetProp("atomLabel", f"{atom.GetSymbol()}:{atom.GetAtomMapNum()}")
    legend = ["actual mapped fragment; attachment endpoints labelled"] * len(fragments)
    svg = Draw.MolsToGridImage(
        fragments,
        molsPerRow=3,
        subImgSize=(430, 330),
        legends=legend,
        useSVG=True,
    )
    if isinstance(svg, bytes):
        svg = svg.decode("utf-8")
    root = ET.fromstring(svg)
    namespace = root.tag.partition("}")[0].lstrip("{")
    metadata_tag = f"{{{namespace}}}metadata" if namespace else "metadata"
    metadata = ET.SubElement(
        root,
        metadata_tag,
        {
            "id": "attachment-map-labels",
            "component-count": str(len(fragments)),
        },
    )
    metadata.text = json.dumps(
        [
            {
                "cut": cut["name"],
                "atom_maps": cut["atom_maps"],
            }
            for cut in cuts
        ],
        sort_keys=True,
    )
    return ET.tostring(root, encoding="unicode")


def propose_synthesis(candidate: dict[str, Any], *, exact_route_records=None, constraints=None) -> dict[str, Any]:
    """Produce a candidate-specific proposal; never a recipe or synthesis proof."""
    base = assess_synthesis(candidate, exact_route_records=exact_route_records, constraints=constraints)
    mol, mapping = _candidate_mol(
        candidate,
        expected_canonical=base["graph_validation"]["canonical_smiles"],
    )
    cuts = _cuts(candidate, mol, mapping)
    _, fragments = _fragments(mol, cuts)
    endpoints = [_endpoint(mol, mapping, cut, candidate["atom_roles"]) for cut in cuts]
    flags = _functional_flags(mol)
    result = {
        "format": VERSION, "candidate_id": candidate["candidate_id"],
        "graph_review": base["graph_validation"],
        "exact_attachment_cuts": [{key: value for key, value in cut.items() if key != "bond_index"} for cut in cuts],
        "actual_fragment_products": fragments,
        "reaction_proposal": {
            "status": "PROPOSAL_requires_expert_review", "endpoints": endpoints,
            "linker_role": "bifunctional graph connector between declared warhead and recruiter partitions",
            "completed_chemical_reaction": False, "reagents_claimed": False, "conditions_claimed": False,
            "interpretation": "The displayed cut is a graph-space edge. Dummy/radical-like fragment endpoints are bookkeeping, not isolated products or a completed reaction.",
        },
        "functional_group_and_competing_site_checks": flags,
        "SA_Score": _sa_score(mol),
        "source_template_suggestions": _source_suggestions(candidate, endpoints),
        "independent_synthesis_disposition": "unverified_proposal_requires_review",
        "synthesis_unverified": True, "unsynthesizable": False,
        "scientific_approved": False, "expert_approved": False,
        "limitations": ["SA_Score is screening, not synthetic proof.", "Graph compatibility does not establish chemoselectivity, route completion, yield, purification, characterization, or independent synthesis.", "An unavailable exact template does not mean unsynthesizable."],
    }
    json.dumps(result, ensure_ascii=False, allow_nan=False)
    return result


__all__ = ["VERSION", "propose_synthesis", "scheme_svg"]
