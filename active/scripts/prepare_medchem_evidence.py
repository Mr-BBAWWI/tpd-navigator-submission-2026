#!/usr/bin/env python3
"""Build source-bound medicinal-chemistry evidence from supplied local files.

The builder performs no network retrieval and does not fill unreported experimental
conditions. The output directory must not already exist.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Any

import gemmi
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFMCS

VERSION = "prepare-medchem-evidence/1.0"
DOI = "10.1021/acs.jmedchem.4c01903"
PARENT_IDS = (
    "SMI-6077", "SMI-6078", "SMI-6079", "SMI-6080", "SMI-6084", "SMI-6409",
    "SMI-6085", "SMI-1069", "SMI-1073", "SMI-1074", "SMI-1089", "SMI-3204",
    "SMI-3108", "SMI-3100", "SMI-3105", "SMI-3106",
)
SAR_PAIRS = (
    ("SMI-6079", "SMI-6080", "central_ring_size_5_to_6"),
    ("SMI-6080", "SMI-6084", "central_ring_size_6_to_7"),
    ("SMI-6080", "SMI-6085", "chlorine_to_bromine"),
    ("SMI-6080", "SMI-1074", "piperidine_substitution_position"),
    ("SMI-6080", "SMD-6087", "terminal_N_derivatization_to_PROTAC"),
)
MEASUREMENT_COLUMNS = ("SMARCA2 IC50 (nM)", "SMARCA4 IC50 (nM)")
ROUTE_CONFIG = {
    "C01": {
        "paper_compound": "2",
        "name": "PROTAC 1",
        "pdb": "6HAY",
        "precursor": "16",
        "heading": "Synthetic methods for the preparation of PROTAC 1 (2)",
    },
    "C02": {
        "paper_compound": "3",
        "name": "PROTAC 2",
        "pdb": "6HAX",
        "precursor": "19",
        "heading": "Synthetic methods for the preparation of PROTAC 2 (3)",
    },
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def source_info(path: Path) -> dict[str, Any]:
    data = path.read_bytes()
    return {
        "path_as_supplied": str(path),
        "filename": path.name,
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
    }


def canonical_connectivity(mol: Chem.Mol) -> str:
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in copy.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
    Chem.SanitizeMol(copy)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=False)


def mapped_molecule(smiles: str, name: str) -> tuple[Chem.Mol, dict[str, Any]]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError(f"invalid connected SMILES for {name}")
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
    identity_copy = Chem.Mol(mol)
    for atom in identity_copy.GetAtoms():
        atom.SetAtomMapNum(0)
    identity = {
        "source_smiles": smiles,
        "canonical_isomeric_smiles": Chem.MolToSmiles(identity_copy, canonical=True, isomericSmiles=True),
        "canonical_connectivity_smiles": canonical_connectivity(mol),
        "inchikey": Chem.MolToInchiKey(Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))),
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "stereocenters": [
            {"atom_index": index, "assignment": assignment}
            for index, assignment in Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        ],
    }
    return mol, identity


def parse_measurement(raw: str, target: str, label: str) -> dict[str, Any]:
    value = raw.strip()
    compact = value.replace(",", "").replace(" ", "")
    match = re.fullmatch(r"([<>]=?)?([0-9]+(?:\.[0-9]+)?)", compact)
    if not match:
        raise ValueError(f"unsupported {label} value for {target}: {raw!r}")
    relation = match.group(1) or "="
    return {
        "target": target,
        "endpoint": "IC50",
        "unit": "nM",
        "relation": relation,
        "value": float(match.group(2)),
        "raw_string": raw,
        "source_column_label": label,
        "assay": "TR-FRET",
        "assay_source_locator": "Supporting Information Figures S1-S3",
        "conditions": {"buffer": "unknown", "temperature": "unknown"},
        "not_reported_as": ["Kd"],
    }


def read_publisher_csv(path: Path) -> tuple[list[dict[str, Any]], dict[str, Chem.Mol]]:
    raw = path.read_bytes()
    text = raw.decode("cp1252")
    normalized = text.replace("\r\r\n", "\r\n").replace("\r\n", "\n").replace("\r", "\n")
    reader = csv.DictReader(normalized.splitlines())
    if reader.fieldnames is None:
        raise ValueError("CSV header missing")
    required = {"Ligands", "Smiles", *MEASUREMENT_COLUMNS}
    if not required.issubset(reader.fieldnames):
        raise ValueError(f"CSV missing columns: {sorted(required - set(reader.fieldnames))}")

    by_id: dict[str, tuple[int, dict[str, str]]] = {}
    for record_number, row in enumerate(reader, start=2):
        compound_id = (row.get("Ligands") or "").strip()
        if compound_id:
            if compound_id in by_id:
                raise ValueError(f"duplicate compound row: {compound_id}")
            by_id[compound_id] = (record_number, row)
    missing = [compound_id for compound_id in PARENT_IDS if compound_id not in by_id]
    if missing:
        raise ValueError(f"missing required parent rows: {', '.join(missing)}")

    records: list[dict[str, Any]] = []
    molecules: dict[str, Chem.Mol] = {}
    for compound_id in PARENT_IDS:
        record_number, row = by_id[compound_id]
        smiles = row["Smiles"]
        if not smiles.strip():
            raise ValueError(f"missing parent SMILES: {compound_id}")
        values = [row.get(column, "") for column in MEASUREMENT_COLUMNS]
        if any(not value.strip() for value in values):
            raise ValueError(f"incomplete parent IC50 row: {compound_id}")
        mol, identity = mapped_molecule(smiles, compound_id)
        molecules[compound_id] = mol
        measurements = [
            parse_measurement(row[column], column.split()[0], column)
            for column in MEASUREMENT_COLUMNS
        ]
        records.append({
            "compound_id": compound_id,
            "record_class": "measured_parent_ligand",
            "identity": identity,
            "measurements": measurements,
            "source": {
                "doi": DOI,
                "filename": path.name,
                "record_number_1_based_including_header": record_number,
                "compound_id_column": "Ligands",
            },
        })
    return records, molecules


def select_diverse_catalog(records: list[dict[str, Any]], molecules: dict[str, Chem.Mol], count: int = 8) -> list[str]:
    ids = [record["compound_id"] for record in records]
    fingerprints = {
        compound_id: AllChem.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(molecules[compound_id])
        for compound_id in ids
    }
    selected = [max(ids, key=lambda item: (molecules[item].GetNumHeavyAtoms(), item))]
    while len(selected) < min(count, len(ids)):
        remaining = [item for item in ids if item not in selected]
        choice = max(
            remaining,
            key=lambda item: (
                min(1.0 - DataStructs.TanimotoSimilarity(fingerprints[item], fingerprints[used]) for used in selected),
                item,
            ),
        )
        selected.append(choice)
    return selected


def atom_description(mol: Chem.Mol, index: int) -> dict[str, Any]:
    atom = mol.GetAtomWithIdx(index)
    return {
        "atom_index": index,
        "atom_map": atom.GetAtomMapNum(),
        "element": atom.GetSymbol(),
        "aromatic": atom.GetIsAromatic(),
        "formal_charge": atom.GetFormalCharge(),
    }


def graph_difference(left: Chem.Mol, right: Chem.Mol, left_match: tuple[int, ...], right_match: tuple[int, ...]) -> dict[str, Any]:
    left_to_right = dict(zip(left_match, right_match))
    right_to_left = {value: key for key, value in left_to_right.items()}
    left_unmatched = sorted(set(range(left.GetNumAtoms())) - set(left_match))
    right_unmatched = sorted(set(range(right.GetNumAtoms())) - set(right_match))
    affected_left: set[int] = set()
    affected_right: set[int] = set()
    bond_changes: list[dict[str, Any]] = []

    for index in left_unmatched:
        for neighbor in left.GetAtomWithIdx(index).GetNeighbors():
            common_left = neighbor.GetIdx()
            if common_left in left_to_right:
                affected_left.add(common_left)
                affected_right.add(left_to_right[common_left])
    for index in right_unmatched:
        for neighbor in right.GetAtomWithIdx(index).GetNeighbors():
            common_right = neighbor.GetIdx()
            if common_right in right_to_left:
                affected_right.add(common_right)
                affected_left.add(right_to_left[common_right])
    for left_bond in left.GetBonds():
        a, b = left_bond.GetBeginAtomIdx(), left_bond.GetEndAtomIdx()
        if a not in left_to_right or b not in left_to_right:
            continue
        right_bond = right.GetBondBetweenAtoms(left_to_right[a], left_to_right[b])
        left_type = str(left_bond.GetBondType())
        right_type = str(right_bond.GetBondType()) if right_bond else None
        if right_type != left_type:
            affected_left.update((a, b))
            affected_right.update((left_to_right[a], left_to_right[b]))
            bond_changes.append({"left_maps": [left.GetAtomWithIdx(a).GetAtomMapNum(), left.GetAtomWithIdx(b).GetAtomMapNum()], "left_type": left_type, "right_type": right_type})
    common_correspondences = sorted(
        (
            left.GetAtomWithIdx(left_index).GetAtomMapNum(),
            right.GetAtomWithIdx(right_index).GetAtomMapNum(),
        )
        for left_index, right_index in left_to_right.items()
    )
    return {
        "left_unmatched_atoms": [atom_description(left, index) for index in left_unmatched],
        "right_unmatched_atoms": [atom_description(right, index) for index in right_unmatched],
        "common_atom_map_correspondences": [
            {"left_atom_map": left_map, "right_atom_map": right_map}
            for left_map, right_map in common_correspondences
        ],
        "left_affected_parent_atom_maps": sorted(left.GetAtomWithIdx(index).GetAtomMapNum() for index in affected_left),
        "right_affected_atom_maps": sorted(right.GetAtomWithIdx(index).GetAtomMapNum() for index in affected_right),
        "changed_common_bonds": bond_changes,
    }


def mcs_pair(left_id: str, right_id: str, label: str, molecules: dict[str, Chem.Mol], records_by_id: dict[str, dict[str, Any]]) -> dict[str, Any]:
    left, right = molecules[left_id], molecules[right_id]
    result = rdFMCS.FindMCS(
        [left, right],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        matchValences=True,
        timeout=2,
    )
    if result.canceled or not result.smartsString:
        raise ValueError(f"MCS failed for {left_id}/{right_id}")
    query = Chem.MolFromSmarts(result.smartsString)
    left_matches = left.GetSubstructMatches(query, uniquify=False, maxMatches=10000)
    right_matches = right.GetSubstructMatches(query, uniquify=False, maxMatches=10000)
    if len(left_matches) >= 10000 or len(right_matches) >= 10000:
        raise ValueError(f"MCS match enumeration limit reached for {left_id}/{right_id}")
    possibilities = []
    seen: set[str] = set()
    for left_match in left_matches:
        for right_match in right_matches:
            difference = graph_difference(left, right, left_match, right_match)
            key = json.dumps(difference, sort_keys=True, separators=(",", ":"))
            if key not in seen:
                seen.add(key)
                possibilities.append(difference)
    definitive = len(possibilities) == 1
    return {
        "pair_id": f"{left_id}__{right_id}",
        "requested_change_label": label,
        "left_compound_id": left_id,
        "right_compound_id": right_id,
        "method": "RDKit MCS with element, exact bond-order, valence, and complete-ring constraints",
        "mcs_smarts": result.smartsString,
        "mcs_atom_count": result.numAtoms,
        "mapping_possibility_count": len(possibilities),
        "mapping_possibilities": possibilities,
        "atom_specific_sar_status": "definitive_for_this_pair" if definitive else "ambiguous_MCS_mapping_review_required",
        "review_required": not definitive,
        "interpretation": "Measured analogue differences support only this observed pair and do not establish an arbitrary new transformation.",
        "left_measurements": records_by_id[left_id]["measurements"],
        "right_measurements": records_by_id[right_id].get("measurements", []),
        "source_row_locators": [records_by_id[left_id]["source"], records_by_id[right_id]["source"]],
        "protected_atom_activation": False,
    }


def safe_mcs_pair(left_id: str, right_id: str, label: str, molecules: dict[str, Chem.Mol], records_by_id: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Return explicit unknown SAR evidence when constrained MCS cannot be resolved."""
    try:
        return mcs_pair(left_id, right_id, label, molecules, records_by_id)
    except Exception as exc:
        return {
            "pair_id": f"{left_id}__{right_id}",
            "requested_change_label": label,
            "left_compound_id": left_id,
            "right_compound_id": right_id,
            "method": "RDKit constrained MCS attempted",
            "mapping_possibility_count": None,
            "mapping_possibilities": [],
            "atom_specific_sar_status": "unknown_due_to_MCS_failure_or_limit",
            "review_required": True,
            "failure_type": exc.__class__.__name__,
            "failure_message": str(exc),
            "interpretation": "No atom-specific SAR assignment is made. Ring-expansion and other symmetric mappings remain under review.",
            "left_measurements": records_by_id[left_id].get("measurements", []),
            "right_measurements": records_by_id[right_id].get("measurements", []),
            "source_row_locators": [records_by_id[left_id]["source"], records_by_id[right_id]["source"]],
            "protected_atom_activation": False,
        }


def _cif_rows(block: gemmi.cif.Block, category: str) -> list[dict[str, str]]:
    values = block.get_mmcif_category(category)
    return [dict(zip(values, row)) for row in zip(*values.values())] if values else []


def read_ccd(path: Path) -> tuple[str, Chem.Mol]:
    block = gemmi.cif.read_file(str(path)).sole_block()
    atoms = _cif_rows(block, "_chem_comp_atom.")
    bonds = _cif_rows(block, "_chem_comp_bond.")
    if not atoms or not bonds:
        raise ValueError("CCD atom and bond tables required")
    rw = Chem.RWMol()
    indices: dict[str, int] = {}
    for row in atoms:
        atom = Chem.Atom(row["type_symbol"].title())
        charge = row.get("charge", "0")
        atom.SetFormalCharge(int(charge) if charge not in {"?", ".", ""} else 0)
        indices[row["atom_id"]] = rw.AddAtom(atom)
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for row in bonds:
        rw.AddBond(indices[row["atom_id_1"]], indices[row["atom_id_2"]], orders[row["value_order"]])
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return block.name, Chem.RemoveHs(mol)


def compare_ccd(path: Path, parent: Chem.Mol) -> dict[str, Any]:
    component_id, ccd = read_ccd(path)
    parent_connectivity = canonical_connectivity(parent)
    ccd_connectivity = canonical_connectivity(ccd)
    connectivity_match = parent_connectivity == ccd_connectivity
    parent_unassigned = any(value == "?" for _, value in Chem.FindMolChiralCenters(parent, includeUnassigned=True))
    ccd_isomeric = Chem.MolToSmiles(ccd, canonical=True, isomericSmiles=True)
    parent_copy = Chem.Mol(parent)
    for atom in parent_copy.GetAtoms():
        atom.SetAtomMapNum(0)
    parent_isomeric = Chem.MolToSmiles(parent_copy, canonical=True, isomericSmiles=True)
    if not connectivity_match:
        stereo_status = "not_comparable_connectivity_mismatch"
    elif parent_unassigned or parent_isomeric != ccd_isomeric:
        stereo_status = "connectivity_match_stereochemistry_not_established_by_csv"
    else:
        stereo_status = "represented_stereochemistry_match"
    return {
        "pdb_id": "9D12",
        "ccd_component_expected": "A1A1P",
        "ccd_component_observed": component_id,
        "compared_parent": "SMI-6080",
        "connectivity_identity_status": "exact_connectivity_match" if connectivity_match else "connectivity_mismatch",
        "stereochemistry_identity_status": stereo_status,
        "parent_canonical_connectivity_smiles": parent_connectivity,
        "ccd_canonical_connectivity_smiles": ccd_connectivity,
        "parent_canonical_isomeric_smiles": parent_isomeric,
        "ccd_canonical_isomeric_smiles": ccd_isomeric,
        "atropisomer_scope": "Axial/atropisomeric identity is not asserted unless explicitly represented by the source graph and RDKit stereochemical model.",
        "no_co_crystal_claim_for_other_parents": True,
    }


def read_text_source(path: Path) -> str:
    if path.suffix.lower() == ".pdf":
        from pypdf import PdfReader
        return "\n\n".join(page.extract_text() or "" for page in PdfReader(str(path)).pages)
    raw = path.read_bytes()
    for encoding in ("utf-8", "cp1252"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            pass
    raise ValueError(f"cannot decode text source: {path}")


def paragraphs(text: str) -> list[dict[str, Any]]:
    chunks = [re.sub(r"\s+", " ", item).strip() for item in re.split(r"(?:\r?\n\s*){2,}", text) if item.strip()]
    return [{"paragraph_index_1_based": index, "text": value, "sha256": sha256_bytes(value.encode("utf-8"))} for index, value in enumerate(chunks, 1)]


def find_paragraph(items: list[dict[str, Any]], predicate: Any, start: int = 0, stop: int | None = None) -> dict[str, Any]:
    for item in items[start:stop]:
        if predicate(item["text"]):
            return item
    raise ValueError("required route paragraph not found")


def optional_paragraph(
    items: list[dict[str, Any]], predicate: Any, start: int = 0, stop: int | None = None
) -> dict[str, Any] | None:
    try:
        return find_paragraph(items, predicate, start, stop)
    except ValueError:
        return None


def paragraph_ref(item: dict[str, Any], summary: str) -> dict[str, Any]:
    return {"summary": summary, "paragraph_index_1_based": item["paragraph_index_1_based"], "paragraph_sha256": item["sha256"]}


def route_record(compound_id: str, config: dict[str, str], items: list[dict[str, Any]], source_url: str, source_hash: str) -> dict[str, Any]:
    heading = find_paragraph(items, lambda text: config["heading"] in text)
    section_start = heading["paragraph_index_1_based"]
    next_section = optional_paragraph(
        items,
        lambda text: text.startswith("Synthetic method") and config["heading"] not in text,
        section_start,
    )
    section_stop = next_section["paragraph_index_1_based"] - 1 if next_section else len(items)
    final_heading = find_paragraph(
        items,
        lambda text: text.endswith(f"({config['paper_compound']})") and "carboxamide" in text,
        section_start,
        section_stop,
    )
    final_start = final_heading["paragraph_index_1_based"]
    procedure = find_paragraph(
        items,
        lambda text: f"{config['precursor']} (" in text and "NaBH(OAc)3" in text,
        final_start,
        section_stop,
    )
    characterization_items = items[procedure["paragraph_index_1_based"]:section_stop]
    nmr = next((item for item in characterization_items if re.match(r"(?:1H|13C) NMR\b", item["text"])), None)
    hrms = next((item for item in characterization_items if re.search(r"\bHRMS\b", item["text"], re.IGNORECASE)), None)
    proc = procedure["text"]
    required_tokens = [config["precursor"], "HCl/THF", "(1).2HCl", "NEt3", "NaBH(OAc)3", "MgSO4"]
    missing = [token for token in required_tokens if token not in proc]
    if missing:
        raise ValueError(f"route {compound_id} procedure lacks: {missing}")
    yield_match = re.search(r"\((\d+(?:\.\d+)?\s*mg),\s*(\d+(?:\.\d+)?\s*% yield)\)", proc)
    if not yield_match:
        raise ValueError(f"route {compound_id} yield not found")
    purification_match = re.search(
        r"purified by (.+?)(?=\s+(?:to (?:give|obtain|afford)|yielding|affording)\b|$)",
        proc,
        flags=re.IGNORECASE,
    )
    if not purification_match:
        raise ValueError(f"route {compound_id} purification not found")
    materials = [config["precursor"], "0.5M HCl/THF", "(1).2HCl", "NEt3", "NaBH(OAc)3", "MgSO4"]
    solvents = "DMF" if compound_id == "C01" else "DCE/DMSO"
    if solvents not in proc:
        raise ValueError(f"route {compound_id} solvent not found")
    materials.append(solvents)
    return {
        "record_id": f"original-route-{compound_id}",
        "compound_id": compound_id,
        "paper_compound_id": config["paper_compound"],
        "paper_name": config["name"],
        "pdb_id": config["pdb"],
        "source_url": source_url,
        "locator": f"{config['heading']}; final compound ({config['paper_compound']})",
        "materials": materials,
        "steps": [
            paragraph_ref(procedure, f"Acidic acetal cleavage of precursor {config['precursor']} to the aldehyde; LCMS reported full deprotection."),
            paragraph_ref(procedure, "Reductive amination with SMARCABD ligand (1).2HCl and sodium triacetoxyborohydride."),
        ],
        "yield": {"isolated_mass": yield_match.group(1), "reported_yield": yield_match.group(2), "source": paragraph_ref(procedure, "Final isolated yield.")},
        "purification": {"summary": purification_match.group(1), "source": paragraph_ref(procedure, "Final purification method.")},
        "characterization": [
            (
                {"type": "NMR", **paragraph_ref(nmr, re.sub(r"\s+", " ", nmr["text"][:180]).rstrip() + "…")}
                if nmr
                else {"type": "NMR", "status": "unknown_not_found_in_bounded_final_compound_section"}
            ),
            (
                {"type": "HRMS", **paragraph_ref(hrms, hrms["text"])}
                if hrms
                else {"type": "HRMS", "status": "unknown_not_found_in_bounded_final_compound_section"}
            ),
        ],
        "reusable_precedent": True,
        "scope": "Original reported route for this exact source compound only; generated CRBN/SAR/linker candidates require a new route assessment.",
        "source_document_sha256": source_hash,
    }


def item_after(items: list[dict[str, Any]], reference: dict[str, Any], text: str) -> bool:
    return any(item["text"] == text and item["paragraph_index_1_based"] > reference["paragraph_index_1_based"] for item in items)


def extract_routes(path: Path, source_url: str) -> list[dict[str, Any]]:
    text = read_text_source(path)
    items = paragraphs(text)
    digest = sha256_bytes(path.read_bytes())
    return [route_record(compound_id, config, items, source_url, digest) for compound_id, config in ROUTE_CONFIG.items()]


def write_sdf(path: Path, records: list[dict[str, Any]], molecules: dict[str, Chem.Mol]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        writer = Chem.SDWriter(stream)
        try:
            for record in records:
                mol = Chem.Mol(molecules[record["compound_id"]])
                mol.SetProp("record_class", "measured_parent_ligand")
                for measurement in record["measurements"]:
                    key = f"{measurement['target']}_IC50_nM"
                    mol.SetProp(key, measurement["raw_string"])
                writer.write(mol)
        finally:
            writer.close()


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")


def build(args: argparse.Namespace) -> None:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"output directory already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        csv_path = args.csv.resolve()
        binding_path = args.binding_si.resolve()
        route_path = args.route_si.resolve()
        records, molecules = read_publisher_csv(csv_path)
        records_by_id = {record["compound_id"]: record for record in records}

        # Add only the source PROTAC required for the terminal-N comparison.
        raw_text = csv_path.read_bytes().decode("cp1252").replace("\r\r\n", "\n").replace("\r", "\n")
        rows = {row["Ligands"].strip(): row for row in csv.DictReader(raw_text.splitlines()) if row.get("Ligands")}
        if "SMD-6087" not in rows or not rows["SMD-6087"].get("Smiles", "").strip():
            raise ValueError("SMD-6087 source row and SMILES required for requested SAR comparison")
        protac_mol, protac_identity = mapped_molecule(rows["SMD-6087"]["Smiles"], "SMD-6087")
        molecules["SMD-6087"] = protac_mol
        records_by_id["SMD-6087"] = {
            "compound_id": "SMD-6087",
            "record_class": "complete_PROTAC_not_parent",
            "identity": protac_identity,
            "measurements": [],
            "source": {"doi": DOI, "filename": csv_path.name, "compound_id_column": "Ligands"},
        }

        sar = [safe_mcs_pair(left, right, label, molecules, records_by_id) for left, right, label in SAR_PAIRS]
        ccd_result = compare_ccd(args.ccd.resolve(), molecules["SMI-6080"]) if args.ccd else {
            "status": "not_supplied",
            "required_match": "PDB 9D12 CCD A1A1P against SMI-6080",
            "no_identity_claim": True,
        }
        routes = extract_routes(route_path, args.route_url)
        binding_text = read_text_source(binding_path)
        if "TR-FRET Assay" not in binding_text or "SMARCA2" not in binding_text or "SMARCA4" not in binding_text:
            raise ValueError("binding SI does not contain required SMARCA2/4 TR-FRET evidence")

        evidence = {
            "schema_version": 1,
            "builder_version": VERSION,
            "scope": "Source-bound measured parent catalogue and pairwise SAR evidence; no efficacy, Kd, protected-atom activation, or unmatched co-crystal claim.",
            "sources": {
                "publisher_csv": source_info(csv_path),
                "binding_supporting_information": source_info(binding_path),
                "route_supporting_note": source_info(route_path),
                "ccd": source_info(args.ccd.resolve()) if args.ccd else None,
            },
            "parent_catalog": {
                "record_count": len(records),
                "all_record_ids": [record["compound_id"] for record in records],
                "diverse_display_subset_ids": select_diverse_catalog(records, molecules),
                "selection_method": "deterministic greedy MaxMin Morgan-radius-2 fingerprint diversity",
                "records": records,
                "excluded_complete_protac_ids": [compound_id for compound_id in rows if compound_id.startswith("SMD-")],
            },
            "sar_pair_evidence": sar,
            "crystal_identity": ccd_result,
        }
        write_json(temporary / "medchem_evidence.json", evidence)
        write_json(temporary / "route_evidence.json", {"schema_version": 1, "evidence_records": routes})
        write_sdf(temporary / "parent_ligands.sdf", records, molecules)
        outputs = {}
        for path in sorted(temporary.iterdir()):
            outputs[path.name] = {"sha256": sha256_bytes(path.read_bytes()), "size_bytes": path.stat().st_size}
        write_json(temporary / "manifest.json", {"builder_version": VERSION, "outputs": outputs})
        os.replace(temporary, output)
        for path in output.iterdir():
            path.chmod(0o444)
        output.chmod(0o555)
    except Exception as exc:
        # Preserve the computed scratch directory for audit and forensic review.
        # Inputs are never modified and the requested output path remains absent.
        if temporary.exists():
            marker = temporary / "BUILD_FAILED.txt"
            try:
                marker.write_text(f"{exc.__class__.__name__}: {exc}\n", encoding="utf-8")
            except OSError:
                pass
        raise


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--csv", required=True, type=Path, help="publisher jm4c01903_si_006.csv")
    result.add_argument("--binding-si", required=True, type=Path, help="binding SI PDF or extracted text")
    result.add_argument("--route-si", required=True, type=Path, help="original 2019 supplementary note PDF or text")
    result.add_argument("--route-url", required=True, help="source URL supplied by the caller")
    result.add_argument("--ccd", type=Path, help="optional local A1A1P CCD CIF for PDB 9D12 identity comparison")
    result.add_argument("--output", required=True, type=Path, help="new output directory")
    return result


def main() -> None:
    build(parser().parse_args())


if __name__ == "__main__":
    main()
