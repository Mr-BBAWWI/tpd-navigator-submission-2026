"""Evidence-bound atom-level exit-vector and pose-preservation analysis."""

import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFreeSASA


_METHOD_VERSION = "exit-vector-developer-heuristics-v1"
_ALLOWED_SAR_STATUS = {
    "known_attachment_precedent",
    "substitution_tolerated",
    "essential",
    "unknown",
}
_EVIDENCE_TIERS = {
    "cocrytal": 0,
    "co_crystal": 0,
    "cocrystal": 0,
    "experimentalbinding": 1,
    "experimental_binding": 1,
    "knownpro_tac": 2,
    "knownprotac": 2,
    "known_protac": 2,
    "dockingsupported": 3,
    "docking_supported": 3,
    "unknown": 4,
}

_HBD_QUERY = Chem.MolFromSmarts(
    "[$([N;!H0;v3,v4&+1]),$([n;H1;+0]),$([O,S;H1;+0])]"
)
_HBA_QUERY = Chem.MolFromSmarts(
    "[$([O,S;H0;v2]),$([O,S;-]),"
    "$([N;v3;!$(N-*=[O,N,P,S]);!$(N-[!#6;!#1]);!$(N~[O,N,S])]),"
    "$([nH0,o,s;+0])]"
)


def _finite_float(value, name):
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _require_one_conformer(mol, name):
    if not isinstance(mol, Chem.Mol):
        raise TypeError(f"{name} must be an RDKit Mol")
    if mol.GetNumConformers() != 1:
        raise ValueError(f"{name} must contain exactly one conformer")
    positions = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    if positions.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(positions).all():
        raise ValueError(f"{name} has nonfinite or malformed coordinates")
    return positions


def _atom_maps(mol, require_all=False):
    mapping = {}
    for atom in mol.GetAtoms():
        atom_map = atom.GetAtomMapNum()
        if atom_map < 0:
            raise ValueError("Negative atom-map numbers are invalid")
        if atom_map == 0:
            if require_all and atom.GetAtomicNum() > 1:
                raise ValueError("Every heavy ligand atom must have a positive atom-map number")
            continue
        if atom_map in mapping:
            raise ValueError(f"Duplicate atom-map number: {atom_map}")
        mapping[atom_map] = atom.GetIdx()
    return mapping


def _normalize_map(value, name="atom map"):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive integer")
    try:
        number = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a positive integer") from exc
    if str(value).strip() not in {str(number), f"+{number}"} and not isinstance(value, int):
        raise ValueError(f"{name} must be a positive integer")
    if number <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return number


def _normalize_map_set(values, name):
    if values is None:
        return set()
    if isinstance(values, (str, bytes, dict)):
        raise TypeError(f"{name} must be an iterable of atom-map numbers")
    return {_normalize_map(value, name) for value in values}


def _protein_arrays(protein):
    if not isinstance(protein, (list, tuple)) or not protein:
        raise ValueError("protein must be a nonempty list of atom records")
    records = []
    coordinates = []
    elements = []
    periodic = Chem.GetPeriodicTable()
    required = {
        "xyz",
        "type_symbol",
        "label_comp_id",
        "auth_seq_id",
        "label_atom_id",
        "label_asym_id",
    }
    for position, record in enumerate(protein):
        if not isinstance(record, dict):
            raise TypeError(f"protein record {position} must be a dictionary")
        missing = sorted(required - set(record))
        if missing:
            raise ValueError(f"protein record {position} is missing: {', '.join(missing)}")
        symbol = str(record["type_symbol"]).strip().title()
        if symbol in {"H", "D"}:
            continue
        try:
            atomic_number = periodic.GetAtomicNumber(symbol)
        except Exception as exc:
            raise ValueError(f"Unsupported protein element at record {position}: {symbol}") from exc
        if atomic_number <= 1:
            raise ValueError(f"Unsupported protein element at record {position}: {symbol}")
        xyz = np.asarray(record["xyz"], dtype=float)
        if xyz.shape != (3,) or not np.isfinite(xyz).all():
            raise ValueError(f"Protein record {position} has nonfinite or malformed coordinates")
        clean = dict(record)
        clean["type_symbol"] = symbol
        clean["xyz"] = [float(value) for value in xyz]
        records.append(clean)
        coordinates.append(xyz)
        elements.append(symbol)
    if not records:
        raise ValueError("protein contains no eligible heavy atoms")
    return records, np.asarray(coordinates, dtype=float), elements


def _vdw_radius(symbol):
    try:
        radius = float(Chem.GetPeriodicTable().GetRvdw(symbol))
    except Exception as exc:
        raise ValueError(f"No RDKit van der Waals radius for element {symbol}") from exc
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError(f"Invalid RDKit van der Waals radius for element {symbol}")
    return radius


def _point_molecule(elements, coordinates):
    if len(elements) != len(coordinates):
        raise ValueError("Element and coordinate counts differ")
    editable = Chem.RWMol()
    conformer = Chem.Conformer(len(elements))
    radii = []
    for index, (symbol, xyz) in enumerate(zip(elements, coordinates)):
        atom = Chem.Atom(symbol)
        if atom.GetAtomicNum() <= 1:
            raise ValueError("SASA point molecule must contain heavy atoms only")
        editable.AddAtom(atom)
        point = np.asarray(xyz, dtype=float)
        if point.shape != (3,) or not np.isfinite(point).all():
            raise ValueError("Nonfinite SASA coordinate")
        conformer.SetAtomPosition(index, tuple(float(value) for value in point))
        radii.append(_vdw_radius(symbol))
    result = editable.GetMol()
    result.AddConformer(conformer)
    return result, radii


def _sasa_values(elements, coordinates, protein_elements=None, protein_coordinates=None):
    isolated, ligand_radii = _point_molecule(elements, coordinates)
    options = rdFreeSASA.SASAOpts()
    options.probeRadius = 1.4
    rdFreeSASA.CalcSASA(isolated, ligand_radii, opts=options)
    isolated_values = [
        float(isolated.GetAtomWithIdx(index).GetDoubleProp("SASA"))
        for index in range(len(elements))
    ]

    all_elements = list(elements)
    all_coordinates = list(coordinates)
    if protein_elements is not None:
        all_elements.extend(protein_elements)
        all_coordinates.extend(protein_coordinates)
    complex_mol, all_radii = _point_molecule(all_elements, all_coordinates)
    rdFreeSASA.CalcSASA(complex_mol, all_radii, opts=options)
    complex_values = [
        float(complex_mol.GetAtomWithIdx(index).GetDoubleProp("SASA"))
        for index in range(len(elements))
    ]
    if not np.isfinite(isolated_values).all() or not np.isfinite(complex_values).all():
        raise ValueError("RDKit rdFreeSASA returned nonfinite values")
    return isolated_values, complex_values, options


def _query_atom_indices(mol, query):
    if query is None:
        return set()
    result = set()
    for match in mol.GetSubstructMatches(query, uniquify=True):
        result.update(match)
    return result


def _chemical_features(mol):
    from rdkit.Chem import Lipinski
    return ({i for match in Lipinski._HDonors(mol) for i in match},
            {i for match in Lipinski._HAcceptors(mol) for i in match})


def _total_hydrogens(atom):
    try:
        return int(atom.GetTotalNumHs(includeNeighbors=True))
    except Exception:
        return int(atom.GetNumExplicitHs() + atom.GetNumImplicitHs())


def _attachment_chemistry(atom):
    hydrogens = _total_hydrogens(atom)
    symbol = atom.GetSymbol()
    if symbol in {"N", "O"} and hydrogens > 0:
        return {
            "chemical_feasibility": "feasible_NH_or_OH_substitution",
            "hydrogen_count": hydrogens,
            "suggested_handles": [
                f"{symbol}-linked substituent replacing an existing {symbol}-H hydrogen"
            ],
            "caveat": "Reaction compatibility, regioselectivity, protonation, and retained activity are not established.",
        }
    if symbol == "C" and hydrogens > 0:
        return {
            "chemical_feasibility": "C-H_derivatization_required",
            "hydrogen_count": hydrogens,
            "suggested_handles": [],
            "caveat": "A C-H bond exists, but no generally feasible attachment chemistry is inferred.",
        }
    return {
        "chemical_feasibility": "unknown",
        "hydrogen_count": hydrogens,
        "suggested_handles": [],
        "caveat": "No NH/OH attachment hydrogen was identified; feasibility requires chemistry-specific review.",
    }


def _contact_record(record, distance):
    return {
        "label_asym_id": record["label_asym_id"],
        "label_comp_id": record["label_comp_id"],
        "auth_seq_id": record["auth_seq_id"],
        "label_atom_id": record["label_atom_id"],
        "type_symbol": record["type_symbol"],
        "distance_A": round(float(distance), 4),
        "interpretation": "heavy-atom proximity only; not proof of a hydrogen bond or binding contribution",
    }


def _parse_sar(sar, valid_maps):
    if sar is None:
        sar = {}
    if not isinstance(sar, dict):
        raise TypeError("sar must be a dictionary keyed by string atom-map numbers")
    parsed = {}
    for raw_map, entry in sar.items():
        if not isinstance(raw_map, str):
            raise ValueError("sar keys must be string atom-map numbers")
        atom_map = _normalize_map(raw_map, "SAR atom map")
        if atom_map not in valid_maps:
            raise ValueError(f"SAR references absent ligand atom map {atom_map}")
        if not isinstance(entry, dict):
            raise TypeError(f"SAR entry {atom_map} must be a dictionary")
        status = entry.get("status", "unknown")
        if status not in _ALLOWED_SAR_STATUS:
            raise ValueError(f"Unsupported SAR status for atom map {atom_map}: {status}")
        allowed = entry.get("allowed_transformations")
        if allowed is None:
            allowed = []
        if isinstance(allowed, str):
            allowed = [allowed]
        elif not isinstance(allowed, (list, tuple)):
            raise TypeError(f"allowed_transformations for atom map {atom_map} must be a list")
        parsed[atom_map] = {
            "status": status,
            "source": entry.get("source"),
            "locator": entry.get("locator"),
            "allowed_transformations": list(allowed),
        }
    return parsed


def analyze_sites(mol, protein, sar, protected_maps, pose_kind="co_crystal"):
    ligand_coordinates = _require_one_conformer(mol, "mol")
    map_to_index = _atom_maps(mol, require_all=True)
    protected_curated = _normalize_map_set(protected_maps, "protected_maps")
    absent_protected = sorted(protected_curated - set(map_to_index))
    if absent_protected:
        raise ValueError(f"protected_maps references absent ligand maps: {absent_protected}")

    protein_records, protein_coordinates, protein_elements = _protein_arrays(protein)
    parsed_sar = _parse_sar(sar, set(map_to_index))

    heavy_indices = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1
    ]
    heavy_elements = [mol.GetAtomWithIdx(index).GetSymbol() for index in heavy_indices]
    heavy_coordinates = ligand_coordinates[heavy_indices]
    heavy_position = {atom_index: position for position, atom_index in enumerate(heavy_indices)}

    isolated_sasa, complex_sasa, sasa_options = _sasa_values(
        heavy_elements,
        heavy_coordinates,
        protein_elements,
        protein_coordinates,
    )
    hbd_indices, hba_indices = _chemical_features(mol)

    atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        index = atom.GetIdx()
        atom_map = atom.GetAtomMapNum()
        position = heavy_position[index]
        coordinate = ligand_coordinates[index]
        distances = np.linalg.norm(protein_coordinates - coordinate, axis=1)
        if not np.isfinite(distances).all():
            raise ValueError("Nonfinite ligand-protein distance")
        nearest = float(np.min(distances))
        contact_positions = np.where(distances <= 4.0)[0]
        contacts = [
            _contact_record(protein_records[int(i)], distances[int(i)])
            for i in sorted(contact_positions, key=lambda i: float(distances[int(i)]))
        ]

        isolated = float(isolated_sasa[position])
        complex_value = float(complex_sasa[position])
        relative_exposure = None if isolated == 0.0 else complex_value / isolated
        burial_fraction = None if isolated == 0.0 else 1.0 - relative_exposure
        sar_entry = parsed_sar.get(
            atom_map,
            {
                "status": "unknown",
                "source": None,
                "locator": None,
                "allowed_transformations": [],
            },
        )
        evidence_citation_complete = bool(sar_entry["source"] and sar_entry["locator"])
        known_precedent = (
            sar_entry["status"] == "known_attachment_precedent"
            and evidence_citation_complete
        )
        tolerated = (
            sar_entry["status"] == "substitution_tolerated"
            and evidence_citation_complete
        )

        protection_reasons = []
        if atom_map in protected_curated:
            protection_reasons.append("caller_curated_core")
        if sar_entry["status"] == "essential":
            protection_reasons.append("SAR_status_essential")
        if nearest < 3.5 and not known_precedent:
            protection_reasons.append(
                "conservative_heavy_atom_proximity_under_3.5_A_without_cited_attachment_precedent"
            )

        exposure_pass = (
            complex_value > 2.0
            and relative_exposure is not None
            and relative_exposure >= 0.25
        )
        distance_pass = nearest >= 3.5
        sar_pass = known_precedent or tolerated
        chemistry = _attachment_chemistry(atom)

        if protection_reasons:
            state = "PROTECTED"
            why_selected = []
            why_excluded = [
                "Protected status overrides exposure or SAR-based modifiability.",
                *protection_reasons,
            ]
        elif exposure_pass and distance_pass and sar_pass:
            state = "MODIFIABLE"
            why_selected = [
                "Complex SASA exceeds 2.0 A^2.",
                "Relative exposure is at least 0.25.",
                "Nearest protein heavy atom is at least 3.5 A away.",
                "Cited atom-specific SAR reports attachment precedent or tolerated substitution.",
            ]
            why_excluded = []
        else:
            state = "UNKNOWN"
            why_selected = []
            why_excluded = []
            if complex_value <= 2.0:
                why_excluded.append("Complex SASA does not exceed 2.0 A^2.")
            if relative_exposure is None:
                why_excluded.append("Relative exposure is undefined because isolated SASA is zero.")
            elif relative_exposure < 0.25:
                why_excluded.append("Relative exposure is below 0.25.")
            if nearest < 3.5:
                why_excluded.append("Nearest protein heavy atom is below 3.5 A.")
            if not sar_pass:
                if sar_entry["status"] in {
                    "known_attachment_precedent",
                    "substitution_tolerated",
                }:
                    why_excluded.append(
                        "SAR status lacks both a source and locator and is not qualifying evidence."
                    )
                else:
                    why_excluded.append(
                        "No cited attachment precedent or substitution-tolerance evidence was supplied."
                    )

        if known_precedent and pose_kind == "co_crystal":
            confidence = "medium"
            confidence_basis = "cited attachment precedent plus caller-labeled co-crystal pose"
        elif atom_map not in parsed_sar or not evidence_citation_complete:
            confidence = "low"
            confidence_basis = "SAR is missing, unknown, or lacks source/locator"
        else:
            confidence = "low"
            confidence_basis = "developer heuristic classification without qualifying co-crystal precedent"

        atoms.append(
            {
                "atom_map": atom_map,
                "rdkit_index_zero_based": index,
                "element": atom.GetSymbol(),
                "formal_charge": atom.GetFormalCharge(),
                "state": state,
                "protected": state == "PROTECTED",
                "protection_reasons": protection_reasons,
                "SASA_complex_A2": complex_value,
                "SASA_isolated_same_pose_A2": isolated,
                "relative_exposure": relative_exposure,
                "burial_fraction": burial_fraction,
                "nearest_protein_heavy_atom_A": nearest,
                "contacts_within_4_A": contacts,
                "chemical_HBD": index in hbd_indices,
                "chemical_HBA": index in hba_indices,
                "attachment_chemistry": chemistry,
                "SAR": sar_entry,
                "SAR_evidence_has_source_and_locator": evidence_citation_complete,
                "why_selected": why_selected,
                "why_excluded": why_excluded,
                "confidence": confidence,
                "confidence_basis": confidence_basis,
                "evidence": {
                    "SASA": {
                        "complex_A2": complex_value,
                        "isolated_same_pose_A2": isolated,
                        "relative_exposure": relative_exposure,
                        "burial_fraction": burial_fraction,
                    },
                    "contacts": {
                        "nearest_heavy_atom_A": nearest,
                        "within_4_A": contacts,
                        "claim_scope": "proximity only",
                    },
                    "SAR": sar_entry,
                },
                "limitations": [
                    "SASA and distances describe only the supplied coordinates and selected protein atoms.",
                    "Heavy-atom proximity is not a hydrogen-bond assignment and does not prove binding contribution.",
                    "Chemical HBD/HBA labels are topology-based and independent of observed geometry.",
                    "Modifiability thresholds are versioned developer heuristics, not expert scientific cutoffs.",
                    "No potency, efficacy, permeability, selectivity, or synthetic success is inferred.",
                ],
            }
        )

    return {
        "atoms": atoms,
        "method": {
            "kind": "evidence_bound_descriptive_exit_vector_analysis",
            "version": _METHOD_VERSION,
            "pose_kind": pose_kind,
            "pose_kind_provenance": "caller supplied; not independently verified",
            "SASA_tool": "RDKit rdFreeSASA",
            "SASA_algorithm": str(sasa_options.algorithm),
            "probe_radius_A": 1.4,
            "radii": "explicit RDKit periodic-table van der Waals radius for every included atom",
            "included_atoms": "ligand and protein heavy atoms only",
            "contact_reporting_radius_A": 4.0,
            "modifiable_heuristics": {
                "complex_SASA_strictly_greater_than_A2": 2.0,
                "relative_exposure_minimum": 0.25,
                "nearest_protein_heavy_atom_minimum_A": 3.5,
                "SAR_requirement": "known_attachment_precedent or substitution_tolerated with source and locator",
            },
            "conservative_protection_heuristics": {
                "proximity_strictly_less_than_A": 3.5,
                "proximity_exception": "cited known_attachment_precedent only",
            },
            "protected_precedence": True,
            "hydrogen_bond_policy": "No hydrogen bond is claimed from heavy-atom proximity.",
            "threshold_status": "versioned developer heuristics; not expert scientific cutoffs",
        },
    }


def _record_tier(record):
    raw = record.get("evidence_tier", record.get("tier"))
    if raw is None:
        return "unknown", True
    normalized = str(raw).strip().lower().replace("-", "_").replace(" ", "_")
    compact = normalized.replace("_", "")
    if normalized in _EVIDENCE_TIERS:
        key = normalized
    elif compact in _EVIDENCE_TIERS:
        key = compact
    else:
        key = "unknown"
    canonical = {
        0: "co_crystal",
        1: "experimentalbinding",
        2: "knownPROTAC",
        3: "dockingsupported",
        4: "unknown",
    }[_EVIDENCE_TIERS[key]]
    return canonical, key == "unknown"


def funnel(records, limit=10):
    if not isinstance(records, (list, tuple)):
        raise TypeError("records must be a list or tuple")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")

    ranked = []
    for input_order, record in enumerate(records):
        if not isinstance(record, dict):
            raise TypeError(f"record {input_order} must be a dictionary")
        tier, tier_missing_or_unknown = _record_tier(record)
        output = dict(record)
        missing_states = []
        if record.get("state") in (None, ""):
            missing_states.append("state")
        if record.get("evidence_tier", record.get("tier")) in (None, ""):
            missing_states.append("evidence_tier")
        output["funnel_evidence_tier"] = tier
        output["funnel_missing_states"] = missing_states
        output["funnel_requires_review"] = bool(missing_states or tier_missing_or_unknown)
        output["funnel_note"] = (
            "Ranked only by categorical evidence tier; Kd, Ki, and IC50 values were not compared."
        )
        output["_funnel_input_order"] = input_order
        ranked.append(output)

    ranked.sort(
        key=lambda item: (
            _EVIDENCE_TIERS[_record_tier(item)[0].lower().replace("_", "")],
            item["_funnel_input_order"],
        )
    )
    result = []
    for rank, item in enumerate(ranked[:limit], start=1):
        clean = dict(item)
        clean.pop("_funnel_input_order", None)
        clean["funnel_rank"] = rank
        result.append(clean)
    return result


def _direct_rmsd(reference, moving):
    reference = np.asarray(reference, dtype=float)
    moving = np.asarray(moving, dtype=float)
    if reference.shape != moving.shape or reference.ndim != 2 or reference.shape[1] != 3:
        raise ValueError("RMSD coordinate arrays must have matching Nx3 shapes")
    if len(reference) == 0:
        return None
    if not np.isfinite(reference).all() or not np.isfinite(moving).all():
        raise ValueError("Nonfinite RMSD coordinates")
    return float(np.sqrt(np.mean(np.sum((moving - reference) ** 2, axis=1))))


def _residue_key(record):
    return (
        str(record["label_asym_id"]),
        str(record["auth_seq_id"]),
        str(record["label_comp_id"]),
    )


def _residue_groups(records, coordinates):
    groups = {}
    for record, coordinate in zip(records, coordinates):
        groups.setdefault(_residue_key(record), []).append(coordinate)
    return {
        key: np.asarray(values, dtype=float)
        for key, values in groups.items()
    }


def _clash_diagnostics(mol, ligand_coordinates, protein_coordinates, protein_elements):
    ligand_indices = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1
    ]
    ligand_xyz = ligand_coordinates[ligand_indices]
    ligand_radii = np.asarray(
        [_vdw_radius(mol.GetAtomWithIdx(index).GetSymbol()) for index in ligand_indices],
        dtype=float,
    )
    protein_radii = np.asarray([_vdw_radius(symbol) for symbol in protein_elements], dtype=float)
    difference = ligand_xyz[:, None, :] - protein_coordinates[None, :, :]
    distances = np.linalg.norm(difference, axis=2)
    if not np.isfinite(distances).all():
        raise ValueError("Nonfinite clash distance")
    overlaps = ligand_radii[:, None] + protein_radii[None, :] - distances
    maximum = float(max(0.0, np.max(overlaps)))
    return {
        "heavy_atom_pairs_with_vdw_overlap_over_0.6_A": int(np.sum(overlaps > 0.6)),
        "maximum_heavy_atom_vdw_overlap_A": maximum,
        "minimum_ligand_protein_heavy_atom_distance_A": float(np.min(distances)),
    }


def pose_preservation(parent, poses, protein, protected_maps, policy=None):
    parent_coordinates = _require_one_conformer(parent, "parent")
    parent_maps = _atom_maps(parent, require_all=True)
    protected = _normalize_map_set(protected_maps, "protected_maps")
    protein_records, protein_coordinates, protein_elements = _protein_arrays(protein)
    residue_groups = _residue_groups(protein_records, protein_coordinates)

    defaults = {
        "core_rmsd_max_A": 2.0,
        "clash_max_overlap_A": 1.0,
        "clash_reporting_overlap_A": 0.6,
        "proximity_retention_minimum": 0.75,
        "protected_proximity_cutoff_A": 4.0,
        "all_common_rmsd_max_A": None,
    }
    if policy is not None:
        if not isinstance(policy, dict):
            raise TypeError("policy must be a dictionary")
        unknown = set(policy) - set(defaults)
        if unknown:
            raise ValueError(f"Unknown policy fields: {sorted(unknown)}")
        defaults.update(policy)
    evaluated_policy = dict(defaults)
    for key, value in list(evaluated_policy.items()):
        if value is None and key == "all_common_rmsd_max_A":
            continue
        evaluated_policy[key] = _finite_float(value, f"policy {key}")
    if evaluated_policy["core_rmsd_max_A"] < 0:
        raise ValueError("core_rmsd_max_A must be nonnegative")
    if evaluated_policy["clash_max_overlap_A"] < 0:
        raise ValueError("clash_max_overlap_A must be nonnegative")
    if evaluated_policy["protected_proximity_cutoff_A"] <= 0:
        raise ValueError("protected_proximity_cutoff_A must be positive")
    if not 0 <= evaluated_policy["proximity_retention_minimum"] <= 1:
        raise ValueError("proximity_retention_minimum must be between zero and one")

    if isinstance(poses, Chem.Mol):
        pose_entries = [poses]
    elif isinstance(poses, (list, tuple)):
        pose_entries = list(poses)
    else:
        raise TypeError("poses must be an RDKit Mol or a list/tuple of poses")

    baseline_pairs = []
    cutoff = evaluated_policy["protected_proximity_cutoff_A"]
    for atom_map in sorted(protected & set(parent_maps)):
        coordinate = parent_coordinates[parent_maps[atom_map]]
        for residue, residue_xyz in residue_groups.items():
            minimum = float(np.min(np.linalg.norm(residue_xyz - coordinate, axis=1)))
            if minimum <= cutoff:
                baseline_pairs.append((atom_map, residue, minimum))

    parent_hbd, parent_hba = _chemical_features(parent)
    parent_feature_maps = {}
    for atom_map, atom_index in parent_maps.items():
        feature_types = []
        if atom_index in parent_hbd:
            feature_types.append("HBD")
        if atom_index in parent_hba:
            feature_types.append("HBA")
        if feature_types:
            parent_feature_maps[atom_map] = feature_types

    diagnostics = []
    any_pass = False
    any_review = False

    for pose_number, entry in enumerate(pose_entries):
        if isinstance(entry, Chem.Mol):
            pose = entry
            score = None
            supplied_metadata = {}
        elif isinstance(entry, dict) and isinstance(entry.get("mol"), Chem.Mol):
            pose = entry["mol"]
            score = entry.get("score")
            supplied_metadata = {
                key: value for key, value in entry.items() if key not in {"mol", "score"}
            }
        else:
            raise TypeError(f"pose {pose_number} must be a Mol or a dictionary containing 'mol'")

        pose_coordinates = _require_one_conformer(pose, f"pose {pose_number}")
        pose_maps = _atom_maps(pose, require_all=True)
        common_maps = sorted(m for m in set(parent_maps) & set(pose_maps)
            if parent.GetAtomWithIdx(parent_maps[m]).GetAtomicNum()>1 and pose.GetAtomWithIdx(pose_maps[m]).GetAtomicNum()>1)
        core_maps = sorted(protected & set(parent_maps) & set(pose_maps))
        missing_protected_parent = sorted(protected - set(parent_maps))
        missing_protected_pose = sorted((protected & set(parent_maps)) - set(pose_maps))

        all_common_rmsd = _direct_rmsd(
            [parent_coordinates[parent_maps[value]] for value in common_maps],
            [pose_coordinates[pose_maps[value]] for value in common_maps],
        )
        core_rmsd = _direct_rmsd(
            [parent_coordinates[parent_maps[value]] for value in core_maps],
            [pose_coordinates[pose_maps[value]] for value in core_maps],
        )

        clash = _clash_diagnostics(
            pose, pose_coordinates, protein_coordinates, protein_elements
        )

        retained_pairs = []
        lost_pairs = []
        evaluable_pairs = 0
        for atom_map, residue, parent_distance in baseline_pairs:
            if atom_map not in pose_maps:
                lost_pairs.append(
                    {
                        "atom_map": atom_map,
                        "residue": {
                            "label_asym_id": residue[0],
                            "auth_seq_id": residue[1],
                            "label_comp_id": residue[2],
                        },
                        "parent_minimum_distance_A": parent_distance,
                        "pose_minimum_distance_A": None,
                    }
                )
                continue
            evaluable_pairs += 1
            pose_distance = float(
                np.min(
                    np.linalg.norm(
                        residue_groups[residue] - pose_coordinates[pose_maps[atom_map]],
                        axis=1,
                    )
                )
            )
            detail = {
                "atom_map": atom_map,
                "residue": {
                    "label_asym_id": residue[0],
                    "auth_seq_id": residue[1],
                    "label_comp_id": residue[2],
                },
                "parent_minimum_distance_A": parent_distance,
                "pose_minimum_distance_A": pose_distance,
                "interpretation": "heavy-atom residue proximity proxy; not a hydrogen-bond assignment",
            }
            if pose_distance <= cutoff:
                retained_pairs.append(detail)
            else:
                lost_pairs.append(detail)

        if baseline_pairs:
            proximity_retention = len(retained_pairs) / len(baseline_pairs)
        else:
            proximity_retention = None

        pose_hbd, pose_hba = _chemical_features(pose)
        feature_positions = []
        for atom_map in sorted(set(parent_feature_maps) & set(pose_maps)):
            parent_index = parent_maps[atom_map]
            pose_index = pose_maps[atom_map]
            pose_types = []
            if pose_index in pose_hbd:
                pose_types.append("HBD")
            if pose_index in pose_hba:
                pose_types.append("HBA")
            displacement = float(
                np.linalg.norm(
                    pose_coordinates[pose_index] - parent_coordinates[parent_index]
                )
            )
            feature_positions.append(
                {
                    "atom_map": atom_map,
                    "parent_feature_types": parent_feature_maps[atom_map],
                    "pose_feature_types": pose_types,
                    "position_displacement_A_in_receptor_frame": displacement,
                }
            )

        review_reasons = []
        failure_reasons = []
        for mapping in core_maps:
            a,b=parent.GetAtomWithIdx(parent_maps[mapping]),pose.GetAtomWithIdx(pose_maps[mapping])
            if (a.GetAtomicNum(),a.GetIsotope(),a.GetFormalCharge(),a.GetIsAromatic()) != (b.GetAtomicNum(),b.GetIsotope(),b.GetFormalCharge(),b.GetIsAromatic()):
                failure_reasons.append('Protected atom chemistry mismatch at map '+str(mapping))
        for bond in parent.GetBonds():
            a,b=bond.GetBeginAtom().GetAtomMapNum(),bond.GetEndAtom().GetAtomMapNum()
            if a in core_maps and b in core_maps:
                other=pose.GetBondBetweenAtoms(pose_maps[a],pose_maps[b])
                if other is None or other.GetBondType()!=bond.GetBondType():failure_reasons.append('Protected bond mismatch')
        if len(common_maps) < 3:
            review_reasons.append("Fewer than three common mapped atoms.")
        if len(core_maps) < 3:
            review_reasons.append("Fewer than three common protected/core mapped atoms.")
        if missing_protected_parent:
            review_reasons.append("Some protected maps are absent from the parent.")
        if missing_protected_pose:
            review_reasons.append("Some parent protected maps are absent from the pose.")

        if core_rmsd is not None and core_rmsd > evaluated_policy["core_rmsd_max_A"]:
            failure_reasons.append("Core RMSD exceeds policy.")
        common_limit = evaluated_policy["all_common_rmsd_max_A"]
        if (
            common_limit is not None
            and all_common_rmsd is not None
            and all_common_rmsd > common_limit
        ):
            failure_reasons.append("All-common-atom RMSD exceeds policy.")
        if clash["maximum_heavy_atom_vdw_overlap_A"] >= evaluated_policy["clash_max_overlap_A"]:
            failure_reasons.append("Maximum protein-ligand VdW overlap is not below policy.")
        if (
            proximity_retention is not None
            and proximity_retention < evaluated_policy["proximity_retention_minimum"]
        ):
            failure_reasons.append("Protected-atom residue-proximity retention is below policy.")

        if failure_reasons:
            pose_status = "fail"
            pose_preserved = False
        elif review_reasons:
            pose_status = "review"
            pose_preserved = "review"
            any_review = True
        else:
            pose_status = "pass"
            pose_preserved = True
            any_pass = True

        diagnostics.append(
            {
                "pose_index_zero_based": pose_number,
                "caller_supplied_score": score,
                "score_used_for_pass": False,
                "caller_supplied_metadata": supplied_metadata,
                "pose_origin": "caller supplied and not inferred or verified",
                "common_mapped_heavy_atom_count": len(
                    [
                        value
                        for value in common_maps
                        if parent.GetAtomWithIdx(parent_maps[value]).GetAtomicNum() > 1
                        and pose.GetAtomWithIdx(pose_maps[value]).GetAtomicNum() > 1
                    ]
                ),
                "core_mapped_atom_count": len(core_maps),
                "common_atom_maps": common_maps,
                "core_atom_maps": core_maps,
                "missing_protected_maps_in_parent": missing_protected_parent,
                "missing_protected_maps_in_pose": missing_protected_pose,
                "all_common_atom_RMSD_A_in_receptor_frame": all_common_rmsd,
                "core_RMSD_A_in_receptor_frame": core_rmsd,
                "alignment_applied": False,
                "clashes": clash,
                "protected_atom_near_residue_preservation": {
                    "cutoff_A": cutoff,
                    "baseline_pair_count": len(baseline_pairs),
                    "evaluable_pair_count": evaluable_pairs,
                    "retained_pair_count": len(retained_pairs),
                    "retention_fraction": proximity_retention,
                    "retained": retained_pairs,
                    "lost_or_missing": lost_pairs,
                    "interpretation": "heavy-atom proximity proxy only; not proof of hydrogen bonding or binding",
                },
                "feature_positions": feature_positions,
                "docking_pose_preserved": pose_preserved,
                "status": pose_status,
                "review_reasons": review_reasons,
                "failure_reasons": failure_reasons,
            }
        )

    if not diagnostics:
        overall = "review"
        status = "review"
        reason = "No poses were supplied."
    elif any_pass:
        overall = True
        status = "pass"
        reason = "At least one pose passed geometry policy; caller scores did not select the result."
    elif any_review:
        overall = "review"
        status = "review"
        reason = "No pose passed and at least one pose lacked sufficient mapped/core evidence."
    else:
        overall = False
        status = "fail"
        reason = "All supplied poses failed one or more geometry policies."

    return {
        "docking_pose_preserved": overall,
        "status": status,
        "reason": reason,
        "all_poses_diagnostics": diagnostics,
        "policy": {
            **evaluated_policy,
            "version": _METHOD_VERSION,
            "threshold_status": "developer heuristics; not expert scientific cutoffs",
            "pass_selection": "geometry policy only; pose score never drives pass",
        },
        "method": {
            "coordinate_frame": "supplied receptor frame",
            "ligand_only_alignment": False,
            "RMSD": "direct mapped-coordinate RMSD without superposition",
            "clash_definition": "sum of explicit RDKit VdW radii minus heavy-atom distance",
            "reported_clash_pair_threshold_A": 0.6,
            "protein_interaction_claim": "heavy-atom proximity proxy only",
            "pose_origin": "not inferred; only the caller can establish whether poses are actual docking outputs",
        },
        "limitations": [
            "A preserved geometry does not establish affinity, efficacy, selectivity, or a correct binding mode.",
            "Protected-atom residue retention uses heavy-atom proximity and is not a hydrogen-bond analysis.",
            "Missing maps or fewer than three common core atoms force review rather than a true result.",
            "The best caller-supplied score is reported but never used to determine preservation.",
        ],
    }
