"""Build and load hash-verified, receptor-aligned literature reference parents.

The builder is deliberately offline. It consumes an immutable collector-v2 catalog,
validates every referenced input hash, globally aligns deposited protein sequences,
and transforms bound ligand coordinates into the supplied 6HAZ chain-A frame.
Unknown or technically inadequate alignments are reported as not ready rather than
being replaced by author-number offsets or guessed coordinates.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path, PurePosixPath

import gemmi
import numpy as np
from rdkit import Chem


_FORMAT_VERSION = "reference-parents-v1"
_PARENT_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")

_MIN_PAIRED_CA = 60
_MIN_SEQUENCE_IDENTITY = 0.8
_MAX_RMSD_A = 3.0

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def _json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode("utf-8")


def _sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def _sha256_path(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_relative(root, relative_path):
    if not isinstance(relative_path, str) or not relative_path:
        raise ValueError("Catalog relative_path must be a nonempty string")
    if "\\" in relative_path or ":" in relative_path or relative_path.startswith("//"):
        raise ValueError(f"Unsafe catalog relative_path: {relative_path!r}")
    raw_parts = relative_path.split("/")
    relative = PurePosixPath(relative_path)
    if (
        relative.is_absolute()
        or any(part in {"", ".", ".."} for part in raw_parts)
        or any(part in {".", ".."} for part in relative.parts)
    ):
        raise ValueError(f"Unsafe catalog relative_path: {relative_path!r}")

    root_path = Path(root).absolute()
    if root_path.is_symlink():
        raise ValueError(f"Catalog root must not be a symbolic link: {root_path}")
    candidate = root_path
    for part in relative.parts:
        candidate = candidate / part
        if candidate.is_symlink():
            raise ValueError(f"Catalog path contains a symbolic link: {relative_path!r}")

    resolved_root = root_path.resolve()
    resolved_candidate = candidate.resolve()
    try:
        resolved_candidate.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(f"Catalog path escapes its root: {relative_path!r}") from exc
    if not resolved_candidate.is_file():
        raise FileNotFoundError(resolved_candidate)
    return resolved_candidate


def _checked_catalog_file(root, descriptor):
    if not isinstance(descriptor, dict):
        raise TypeError("Catalog file descriptor must be a dictionary")
    expected = str(descriptor.get("sha256", "")).lower()
    if not _SHA256.fullmatch(expected):
        raise ValueError("Catalog file descriptor has an invalid SHA-256")
    path = _safe_relative(root, descriptor.get("relative_path"))
    actual = _sha256_path(path)
    if actual != expected:
        raise ValueError(f"Hash mismatch for catalog input {descriptor['relative_path']}")
    return path, actual


def _rows(block, category):
    values = block.get_mmcif_category(category)
    return [dict(zip(values, row)) for row in zip(*values.values())] if values else []


def _clean_cif(value):
    return None if value in (None, False, "", ".", "?") else str(value)


def kabsch_transform(reference_xyz, moving_xyz):
    """Return the proper row-vector rigid transform moving -> reference."""
    reference = np.asarray(reference_xyz, dtype=float)
    moving = np.asarray(moving_xyz, dtype=float)
    if reference.shape != moving.shape or reference.ndim != 2 or reference.shape[1] != 3:
        raise ValueError("Kabsch arrays must have matching Nx3 shapes")
    if len(reference) < 3:
        raise ValueError("At least three coordinate pairs are required")
    if not np.isfinite(reference).all() or not np.isfinite(moving).all():
        raise ValueError("Kabsch coordinates must be finite")
    x = moving - moving.mean(axis=0)
    y = reference - reference.mean(axis=0)
    if np.linalg.matrix_rank(x) < 2 or np.linalg.matrix_rank(y) < 2:
        raise ValueError("Collinear coordinates do not define an unambiguous transform")
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.eye(3)
    correction[-1, -1] = -1.0 if np.linalg.det(u @ vt) < 0 else 1.0
    rotation = u @ correction @ vt
    translation = reference.mean(axis=0) - moving.mean(axis=0) @ rotation
    fitted = moving @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - reference) ** 2, axis=1))))
    if not math.isfinite(rmsd) or not np.isfinite(rotation).all() or not np.isfinite(translation).all():
        raise ValueError("Kabsch calculation produced nonfinite values")
    return {
        "rotation_row_vectors": rotation.tolist(),
        "translation_A": translation.tolist(),
        "rmsd_A": rmsd,
    }


def apply_transform(xyz, transform):
    coordinates = np.asarray(xyz, dtype=float)
    rotation = np.asarray(transform["rotation_row_vectors"], dtype=float)
    translation = np.asarray(transform["translation_A"], dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("Coordinates must be Nx3")
    if rotation.shape != (3, 3) or translation.shape != (3,):
        raise ValueError("Malformed rigid transform")
    result = coordinates @ rotation + translation
    if not np.isfinite(result).all():
        raise ValueError("Transform produced nonfinite coordinates")
    return result


def global_sequence_alignment(reference_sequence, moving_sequence):
    """Needleman-Wunsch alignment retaining only pairs stable across all optima."""
    reference = str(reference_sequence)
    moving = str(moving_sequence)
    allowed = set("ACDEFGHIKLMNPQRSTVWY")
    if not reference or not moving or set(reference) - allowed or set(moving) - allowed:
        raise ValueError("Protein sequences must be nonempty unambiguous amino-acid strings")

    n, m = len(reference), len(moving)
    forward = np.empty((n + 1, m + 1), dtype=np.int32)
    reverse = np.empty((n + 1, m + 1), dtype=np.int32)
    paths = np.zeros((n + 1, m + 1), dtype=np.uint8)
    forward[:, 0] = -2 * np.arange(n + 1)
    forward[0, :] = -2 * np.arange(m + 1)
    paths[:, 0] = 1
    paths[0, :] = 1

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = forward[i - 1, j - 1] + (2 if reference[i - 1] == moving[j - 1] else -1)
            up = forward[i - 1, j] - 2
            left = forward[i, j - 1] - 2
            best = max(diagonal, up, left)
            forward[i, j] = best
            count = 0
            if diagonal == best:
                count += int(paths[i - 1, j - 1])
            if up == best:
                count += int(paths[i - 1, j])
            if left == best:
                count += int(paths[i, j - 1])
            paths[i, j] = min(2, count)

    reverse[n, :] = -2 * np.arange(m, -1, -1)
    reverse[:, m] = -2 * np.arange(n, -1, -1)
    for i in range(n - 1, -1, -1):
        for j in range(m - 1, -1, -1):
            reverse[i, j] = max(
                reverse[i + 1, j + 1] + (2 if reference[i] == moving[j] else -1),
                reverse[i + 1, j] - 2,
                reverse[i, j + 1] - 2,
            )

    optimum = int(forward[n, m])
    reference_diagonals = {i: set() for i in range(1, n + 1)}
    moving_diagonals = {j: set() for j in range(1, m + 1)}
    reference_can_gap = {i: False for i in range(1, n + 1)}
    moving_can_gap = {j: False for j in range(1, m + 1)}

    for i in range(n + 1):
        for j in range(m + 1):
            if i < n and j < m:
                edge = 2 if reference[i] == moving[j] else -1
                if int(forward[i, j]) + edge + int(reverse[i + 1, j + 1]) == optimum:
                    reference_diagonals[i + 1].add(j + 1)
                    moving_diagonals[j + 1].add(i + 1)
            if i < n and int(forward[i, j]) - 2 + int(reverse[i + 1, j]) == optimum:
                reference_can_gap[i + 1] = True
            if j < m and int(forward[i, j]) - 2 + int(reverse[i, j + 1]) == optimum:
                moving_can_gap[j + 1] = True

    stable_pairs = []
    for ri in range(1, n + 1):
        possible = reference_diagonals[ri]
        if len(possible) != 1 or reference_can_gap[ri]:
            continue
        mi = next(iter(possible))
        if moving_diagonals[mi] != {ri} or moving_can_gap[mi]:
            continue
        stable_pairs.append((ri, mi, reference[ri - 1], moving[mi - 1]))

    stable_reference = {ri for ri, _, _, _ in stable_pairs}
    stable_moving = {mi for _, mi, _, _ in stable_pairs}
    matches = sum(a == b for _, _, a, b in stable_pairs)
    identity = matches / len(stable_pairs) if stable_pairs else 0.0
    excluded_reference = [
        {
            "reference_label_seq_id": i,
            "possible_moving_label_seq_ids": sorted(reference_diagonals[i]),
            "can_be_aligned_to_gap_in_an_optimal_path": reference_can_gap[i],
        }
        for i in range(1, n + 1) if i not in stable_reference
    ]
    excluded_moving = [
        {
            "moving_label_seq_id": j,
            "possible_reference_label_seq_ids": sorted(moving_diagonals[j]),
            "can_be_aligned_to_gap_in_an_optimal_path": moving_can_gap[j],
        }
        for j in range(1, m + 1) if j not in stable_moving
    ]
    return {
        "algorithm": "Needleman-Wunsch global alignment with forward/reverse optimal-edge analysis",
        "scoring": {"match": 2, "mismatch": -1, "gap": -2},
        "score": optimum,
        "reference_length": n,
        "moving_length": m,
        "aligned_residue_pair_count": len(stable_pairs),
        "identical_residue_pair_count": matches,
        "sequence_identity_over_aligned_residue_pairs": identity,
        "optimal_alignment_ambiguous": bool(paths[n, m] > 1),
        "stable_pair_policy": "A pair is retained only when each residue is paired exclusively to the other and neither residue can be gapped in any optimal alignment.",
        "ambiguity_outside_used_core": bool(paths[n, m] > 1 and (excluded_reference or excluded_moving)),
        "excluded_reference_positions": excluded_reference,
        "excluded_moving_positions": excluded_moving,
        "pairs": [
            {
                "reference_label_seq_id": ri,
                "moving_label_seq_id": mi,
                "reference_one_letter": ra,
                "moving_one_letter": ma,
                "identity_match": ra == ma,
                "stable_across_all_optimal_alignments": True,
            }
            for ri, mi, ra, ma in stable_pairs
        ],
    }


def _cif_chain(path, chain):
    block = gemmi.cif.read_file(str(path)).sole_block()
    asym_to_entity = {
        str(row["id"]): str(row["entity_id"])
        for row in _rows(block, "_struct_asym.")
    }
    if chain not in asym_to_entity:
        raise ValueError(f"mmCIF has no label chain {chain!r}")
    entity = asym_to_entity[chain]
    sequence_rows = [
        row for row in _rows(block, "_entity_poly_seq.")
        if str(row["entity_id"]) == entity
    ]
    if not sequence_rows:
        raise ValueError(f"Chain {chain!r} has no _entity_poly_seq records")
    sequence_by_id = {}
    for row in sequence_rows:
        number = int(row["num"])
        residue = _THREE_TO_ONE.get(str(row["mon_id"]).upper())
        if residue is None:
            raise ValueError(f"Unsupported sequence residue {row['mon_id']!r}")
        if number in sequence_by_id and sequence_by_id[number] != residue:
            raise ValueError("Conflicting _entity_poly_seq records")
        sequence_by_id[number] = residue
    expected = list(range(1, max(sequence_by_id) + 1))
    if sorted(sequence_by_id) != expected:
        raise ValueError("_entity_poly_seq numbering must be complete from one")
    sequence = "".join(sequence_by_id[i] for i in expected)

    ca = {}
    ca_records = {}
    for row in _rows(block, "_atom_site."):
        if str(row.get("label_asym_id")) != chain:
            continue
        if _clean_cif(row.get("pdbx_PDB_model_num")) not in (None, "1"):
            continue
        if str(row.get("label_atom_id")) != "CA" or str(row.get("type_symbol")).upper() != "C":
            continue
        alt = _clean_cif(row.get("label_alt_id"))
        if alt not in (None, "A") or float(row.get("occupancy", 1.0)) <= 0:
            continue
        raw_seq = _clean_cif(row.get("label_seq_id"))
        if raw_seq is None:
            continue
        seq_id = int(raw_seq)
        if seq_id in ca:
            raise ValueError(f"Ambiguous CA for chain {chain} label_seq_id {seq_id}")
        if seq_id not in sequence_by_id:
            raise ValueError("CA label_seq_id is absent from entity sequence")
        comp = _THREE_TO_ONE.get(str(row.get("label_comp_id", "")).upper())
        if comp != sequence_by_id[seq_id]:
            raise ValueError("CA residue identity disagrees with _entity_poly_seq")
        xyz = np.asarray([float(row[f"Cartn_{axis}"]) for axis in "xyz"], dtype=float)
        if xyz.shape != (3,) or not np.isfinite(xyz).all():
            raise ValueError("Malformed CA coordinate")
        ca[seq_id] = xyz
        ca_records[seq_id] = {
            "label_seq_id": seq_id,
            "label_comp_id": str(row["label_comp_id"]),
            "auth_seq_id": _clean_cif(row.get("auth_seq_id")),
            "pdbx_PDB_ins_code": _clean_cif(row.get("pdbx_PDB_ins_code")),
            "atom_site_id": str(row.get("id")),
        }
    if not ca:
        raise ValueError(f"Chain {chain!r} has no eligible CA coordinates")
    return {"sequence": sequence, "ca": ca, "ca_records": ca_records, "entity_id": entity}


def align_cif_chains(reference_cif, reference_chain, moving_cif, moving_chain):
    reference = _cif_chain(reference_cif, reference_chain)
    moving = _cif_chain(moving_cif, moving_chain)
    alignment = global_sequence_alignment(reference["sequence"], moving["sequence"])

    reference_xyz = []
    moving_xyz = []
    ca_pairs = []
    for pair in alignment["pairs"]:
        if not pair["identity_match"]:
            continue
        ri = pair["reference_label_seq_id"]
        mi = pair["moving_label_seq_id"]
        if ri not in reference["ca"] or mi not in moving["ca"]:
            continue
        reference_xyz.append(reference["ca"][ri])
        moving_xyz.append(moving["ca"][mi])
        ca_pairs.append({
            "reference": reference["ca_records"][ri],
            "moving": moving["ca_records"][mi],
            "one_letter": pair["reference_one_letter"],
        })

    result = {
        "reference_chain": reference_chain,
        "source_chain": moving_chain,
        "reference_entity_id": reference["entity_id"],
        "source_entity_id": moving["entity_id"],
        "sequence_alignment": alignment,
        "paired_identical_CA_count": len(ca_pairs),
        "paired_CA": ca_pairs,
        "technical_limits": {
            "minimum_paired_CA": _MIN_PAIRED_CA,
            "minimum_sequence_identity": _MIN_SEQUENCE_IDENTITY,
            "maximum_CA_RMSD_A": _MAX_RMSD_A,
            "status": "developer review limits, not biological similarity cutoffs",
        },
    }
    reasons = []
    if alignment["sequence_identity_over_aligned_residue_pairs"] < _MIN_SEQUENCE_IDENTITY:
        reasons.append("Sequence identity is below 0.8.")
    if len(ca_pairs) < _MIN_PAIRED_CA:
        reasons.append("Fewer than 60 residue-identity-matched CA pairs are modeled.")

    transform = None
    if len(ca_pairs) >= 3:
        try:
            transform = kabsch_transform(reference_xyz, moving_xyz)
        except ValueError as exc:
            reasons.append(str(exc))
        else:
            if transform["rmsd_A"] > _MAX_RMSD_A:
                reasons.append("Aligned CA RMSD exceeds 3 A.")
    else:
        reasons.append("Fewer than three CA pairs are available for rigid-transform metrics.")

    result.update({"ready": not reasons, "reasons": reasons, "transform": transform})
    return result


def _canonical(mol):
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    return Chem.MolToSmiles(copy, isomericSmiles=True)


def _read_single_sdf(path):
    with Path(path).open("rb") as stream:
        supplier = Chem.ForwardSDMolSupplier(stream, sanitize=True, removeHs=False)
        molecules = [mol for mol in supplier if mol is not None]
    if len(molecules) != 1:
        raise ValueError(f"Bound SDF must contain exactly one readable molecule: {path}")
    mol = molecules[0]
    if mol.GetNumConformers() != 1:
        raise ValueError("Bound ligand must contain exactly one conformer")
    coordinates = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    if coordinates.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(coordinates).all():
        raise ValueError("Bound ligand has malformed coordinates")
    maps = [a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    if any(value <= 0 for value in maps) or len(maps) != len(set(maps)):
        raise ValueError("Every bound ligand heavy atom must retain a unique positive atom map")
    return mol


def _transform_molecule(mol, transform):
    result = Chem.Mol(mol)
    before = _canonical(result)
    coordinates = apply_transform(result.GetConformer().GetPositions(), transform)
    conformer = result.GetConformer()
    for index, xyz in enumerate(coordinates):
        conformer.SetAtomPosition(index, tuple(float(value) for value in xyz))
    if _canonical(result) != before:
        raise ValueError("Rigid coordinate transformation changed ligand graph identity")
    return result


def _sdf_bytes(mol):
    from io import StringIO
    stream = StringIO()
    writer = Chem.SDWriter(stream)
    writer.write(mol)
    writer.close()
    return stream.getvalue().replace("\r\n", "\n").encode("utf-8")


def _find_catalog(catalog_dir):
    root = Path(catalog_dir).resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)
    path = root / "catalog.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TypeError("collector catalog.json must contain a dictionary")
    if not isinstance(value.get("cards"), list) or not isinstance(value.get("sources"), list):
        raise ValueError("collector catalog.json must contain cards and sources lists")
    if value.get("source_root", "source_snapshots") != "source_snapshots":
        raise ValueError("collector catalog source_root must be source_snapshots")
    return root, path.resolve(), value


def _source_descriptor(catalog, expected_hash, name):
    expected_hash = str(expected_hash).lower()
    if not _SHA256.fullmatch(expected_hash):
        raise ValueError("Representative record has an invalid source SHA-256")
    if not isinstance(name, str) or not name or "/" in name or "\\" in name or ":" in name:
        raise ValueError("Source snapshot name must be a plain file name")
    expected_relative = f"{expected_hash}/{name}"
    matches = [
        source for source in catalog["sources"]
        if isinstance(source, dict)
        and str(source.get("sha256", "")).lower() == expected_hash
        and source.get("relative_path") == expected_relative
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one catalog source at {expected_relative!r}, found {len(matches)}"
        )
    return matches[0]


def _representative_record(card):
    versions = card.get("experimental_versions")
    if not isinstance(versions, list) or not versions:
        raise ValueError("Card has no experimental_versions")
    requested = card.get("representative_record_id")
    matches = [record for record in versions if record.get("record_id") == requested]
    if len(matches) != 1:
        raise ValueError("Card representative_record_id does not identify exactly one version")
    return matches[0]


def _choose_source_chain(record):
    rows = record.get("contacts", {}).get("by_exact_target_chain", [])
    choices = []
    for row in rows:
        chain = row.get("target_label_asym_id")
        count = row.get("contact_pair_count")
        if not isinstance(chain, str) or not chain or isinstance(count, bool):
            continue
        choices.append((int(count), chain))
    if not choices:
        raise ValueError("Representative record has no exact target-chain contact counts")
    maximum = max(count for count, _ in choices)
    tied = sorted({chain for count, chain in choices if count == maximum})
    return tied[0], {
        "criterion": "maximum deposited ligand-protein heavy-atom contact_pair_count",
        "maximum_contact_pair_count": maximum,
        "tied_source_chains": tied,
        "ambiguous_tie": len(tied) > 1,
        "tie_break": "lexicographically first label_asym_id",
        "selected": tied[0],
    }


def _parent_id(card, record):
    if str(record.get("pdb", "")).upper() == "6HAZ" and str(record.get("ccd", "")).upper() == "FX5":
        return "SMARCA2-FX5"
    value = f"SMARCA2-{str(record.get('pdb', '')).upper()}-{str(record.get('ccd', '')).upper()}"
    if not _PARENT_ID.fullmatch(value):
        raise ValueError(f"Cannot construct safe parent id from record {record.get('record_id')!r}")
    return value


def _curated_map_references(value):
    references = set()

    def visit(item):
        if isinstance(item, dict):
            for key, child in item.items():
                if key == "atom_map" and isinstance(child, int) and not isinstance(child, bool):
                    references.add(child)
                visit(child)
        elif isinstance(item, list):
            for child in item:
                visit(child)

    visit(value)
    return references


def _manifest_sha256_values(value):
    hashes = set()
    if isinstance(value, dict):
        for key, child in value.items():
            if "sha256" in str(key).lower() and isinstance(child, str) and _SHA256.fullmatch(child.lower()):
                hashes.add(child.lower())
            hashes.update(_manifest_sha256_values(child))
    elif isinstance(value, list):
        for child in value:
            hashes.update(_manifest_sha256_values(child))
    return hashes


def _mapping_by_ccd_name(rows, label):
    if not isinstance(rows, list) or not rows:
        raise ValueError(f"{label} atom_mapping must be a nonempty list")
    result = {}
    used_maps = set()
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError(f"{label} atom_mapping rows must be dictionaries")
        name = next(
            (row.get(key) for key in ("ccd_atom_id", "ccd_atom_name", "atom_id") if isinstance(row.get(key), str) and row.get(key)),
            None,
        )
        atom_map = row.get("atom_map")
        if name is None or not isinstance(atom_map, int) or isinstance(atom_map, bool) or atom_map <= 0:
            raise ValueError(f"{label} atom_mapping requires CCD atom names and positive integer atom maps")
        if name in result or atom_map in used_maps:
            raise ValueError(f"{label} atom_mapping must be one-to-one")
        result[name] = atom_map
        used_maps.add(atom_map)
    return result


def _remap_curated_atom_maps(value, old_to_new):
    if isinstance(value, dict):
        return {
            key: (
                old_to_new[child]
                if key == "atom_map" and isinstance(child, int) and child in old_to_new
                else _remap_curated_atom_maps(child, old_to_new)
            )
            for key, child in value.items()
        }
    if isinstance(value, list):
        return [_remap_curated_atom_maps(child, old_to_new) for child in value]
    return value


def _curated_fx5(reference_cif, parent_id, record, ligand, hashes):
    if parent_id != "SMARCA2-FX5":
        return None
    path = Path(reference_cif).with_name("warhead.json")
    manifest_path = Path(reference_cif).with_name("manifest.json")
    if not path.is_file():
        raise FileNotFoundError("Curated FX5 metadata warhead.json is required beside the reference CIF")
    if not manifest_path.is_file():
        raise FileNotFoundError("Curated FX5 inputs must be covered by the adjacent manifest.json")
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("id") != parent_id or value.get("pdb") != "6HAZ" or value.get("ccd") != "FX5":
        raise ValueError("Curated FX5 metadata does not identify SMARCA2-FX5/6HAZ/FX5")

    curated_hash = _sha256_path(path)
    manifest_hashes = _manifest_sha256_values(json.loads(manifest_path.read_text(encoding="utf-8")))
    required_hashes = {
        hashes["reference_6HAZ_mmcif_sha256"],
        curated_hash,
    }
    if not required_hashes.issubset(manifest_hashes):
        raise ValueError("Design-source manifest does not cover the verified 6HAZ CIF and warhead.json")

    declared_identity = value.get("canonical_isomeric_smiles")
    if declared_identity is None and isinstance(value.get("parent_identity"), dict):
        declared_identity = value["parent_identity"].get("canonical_isomeric_smiles")
    if declared_identity is None and isinstance(value.get("identity"), dict):
        declared_identity = value["identity"].get("canonical_isomeric_smiles")
    declared_mol = Chem.MolFromSmiles(declared_identity) if isinstance(declared_identity, str) else None
    if declared_mol is None or _canonical(declared_mol) != _canonical(ligand):
        raise ValueError("Curated FX5 canonical identity differs from the verified source ligand")

    actual_maps = {
        atom.GetAtomMapNum() for atom in ligand.GetAtoms() if atom.GetAtomicNum() > 1
    }
    source_by_name = _mapping_by_ccd_name(record.get("atom_mapping"), "Representative-record")
    curated_by_name = _mapping_by_ccd_name(value.get("atom_mapping"), "Curated FX5")
    if set(source_by_name) != set(curated_by_name):
        raise ValueError("Curated FX5 and source mappings do not cover the same CCD atom names")
    if set(source_by_name.values()) != actual_maps:
        raise ValueError("FX5 source atom maps disagree with representative-record CCD atom_mapping")

    old_to_new = {curated_by_name[name]: source_by_name[name] for name in source_by_name}
    if len(old_to_new) != len(curated_by_name):
        raise ValueError("Curated FX5 map-number translation is not one-to-one")
    remapped = _remap_curated_atom_maps(value, old_to_new)
    protected_original = value.get("protected_maps", [])
    if (
        not isinstance(protected_original, list)
        or any(not isinstance(item, int) or isinstance(item, bool) for item in protected_original)
        or len(protected_original) != len(set(protected_original))
        or any(item not in old_to_new for item in protected_original)
    ):
        raise ValueError("Curated FX5 protected_maps must be unique mapped integer atom maps")
    remapped["protected_maps"] = [old_to_new[item] for item in protected_original]
    referenced = set(remapped["protected_maps"]) | _curated_map_references(remapped.get("sar", {}))
    if not referenced.issubset(actual_maps):
        raise ValueError("Curated FX5 metadata references atom maps absent from the source ligand")
    return remapped, path, curated_hash


def _metadata(card, record, parent_id, alignment, hashes, chain_choice, curated, observed_charge):
    if curated is not None:
        curated_value, curated_path, curated_hash = curated
        sar = curated_value.get("sar", {})
        protected_maps = curated_value.get("protected_maps", [])
        curated_source = {"path_name": curated_path.name, "sha256": curated_hash}
        limitations = list(curated_value.get("limitations", []))
    else:
        sar = {}
        protected_maps = []
        curated_source = None
        limitations = [
            "No atom-specific SAR was curated; strict design has zero approved modifiable sites and requires expert review.",
            "UNKNOWN sites remain UNKNOWN and may only be explored under an explicitly exploratory policy.",
            "A co-crystal pose and receptor-frame comparison do not establish affinity, selectivity, or attachment tolerance.",
        ]
    has_qualifying_sar = bool(sar)
    return {
        "format_version": _FORMAT_VERSION,
        "id": parent_id,
        "target": record.get("target", "SMARCA2"),
        "pdb": record.get("pdb"),
        "ccd": record.get("ccd"),
        "evidence_tier": "co_crystal",
        "source_record_id": record.get("record_id"),
        "source_group_id": card.get("group_id"),
        "source_chain": alignment["source_chain"],
        "reference_chain": alignment["reference_chain"],
        "source_chain_selection": chain_choice,
        "source_input_hashes": hashes,
        "identity": record.get("identity"),
        "canonical_isomeric_smiles": record.get("identity", {}).get(
            "canonical_isomeric_smiles", card.get("canonical_isomeric_smiles")
        ),
        "atom_mapping": record.get("atom_mapping", []),
        "sar": sar,
        "protected_maps": protected_maps,
        "unknown_site_policy": "No missing or UNKNOWN SAR state is promoted to MODIFIABLE.",
        "chemical_state_policy": "This aligned co-crystal record preserves the deposited ligand graph and charge. It does not replace the platform's separately defined explicitly neutral original default.",
        "observed_source_chemical_state": {
            "label": "deposited source state",
            "formal_charge": observed_charge,
            "coordinate_provenance": "source bound-ligand coordinates after receptor-frame rigid transformation",
            "neutral_parent_handling": "Handled separately by the platform; this builder does not claim or emit a neutralized state.",
        },
        "alignment": alignment,
        "curated_metadata_source": curated_source,
        "measured_evidence": {
            "Kd": record.get("measured_Kd"),
            "Ki": record.get("measured_Ki"),
            "IC50": record.get("measured_IC50"),
            "join_policy": "Stored separately; no cross-assay potency ordering or medchem graph join is performed here.",
        },
        "suitability": {
            "fixed_receptor_pose_comparison": True,
            "strict_design": has_qualifying_sar,
            "exploratory_design": True,
            "strict_modifiable_site_explanation": (
                "Curated atom-specific SAR is present; downstream evidence rules still control modifiability."
                if has_qualifying_sar else
                "No SAR-qualified modifiable site is supplied; expert curation is required rather than treating this as an error."
            ),
        },
        "limitations": limitations + [
            "Coordinates are a rigidly transformed deposited ligand pose, not a docking result.",
            "Alignment limits are technical developer-review limits, not biological acceptance criteria.",
            "The fixed 6HAZ receptor enables pose comparison and hypothesis generation only.",
            "The requested two SAR-qualified sites, each with 20 valid graphs, remain unmet unless explicit evidence is curated later.",
        ],
    }


def build_reference_parents(catalog_dir, reference_cif, output):
    """Build a new immutable parent folder and return its index dictionary."""
    catalog_root, catalog_path, catalog = _find_catalog(catalog_dir)
    reference_cif = Path(reference_cif).resolve()
    if not reference_cif.is_file():
        raise FileNotFoundError(reference_cif)
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Output already exists; destructive overwrite is forbidden: {output}")

    reference_hash = _sha256_path(reference_cif)
    prepared = []
    not_ready = []
    ids = set()

    for card in catalog["cards"]:
        try:
            if card.get("domain_partition") != "bromodomain":
                not_ready.append({"group_id": card.get("group_id"), "reason": "Card is not a bromodomain card."})
                continue
            record = _representative_record(card)
            parent_id = _parent_id(card, record)
            if parent_id in ids:
                raise ValueError(f"Duplicate generated parent id: {parent_id}")
            ids.add(parent_id)

            bound_path, bound_hash = _checked_catalog_file(catalog_root, card.get("bound_sdf"))
            source_hash = str(record.get("source_hashes", {}).get("pdb_mmcif_sha256", "")).lower()
            source_descriptor = _source_descriptor(catalog, source_hash, f"{record['pdb']}.cif")
            source_cif, verified_source_hash = _checked_catalog_file(
                catalog_root / "source_snapshots", source_descriptor
            )
            source_chain, chain_choice = _choose_source_chain(record)

            ligand = _read_single_sdf(bound_path)
            expected_graph = record.get("identity", {}).get(
                "canonical_isomeric_smiles", card.get("canonical_isomeric_smiles")
            )
            expected_mol = Chem.MolFromSmiles(expected_graph) if expected_graph else None
            if expected_mol is None or _canonical(ligand) != _canonical(expected_mol):
                raise ValueError("Bound SDF graph differs from the representative record; no automatic graph approval")
            observed_charge = int(Chem.GetFormalCharge(ligand))
            declared_charge = record.get("identity", {}).get("formal_charge")
            if isinstance(declared_charge, bool) or not isinstance(declared_charge, int):
                raise ValueError("Representative record identity must declare an integer formal_charge")
            if observed_charge != declared_charge:
                raise ValueError("Bound SDF formal charge differs from the representative record")
            alignment = align_cif_chains(reference_cif, "A", source_cif, source_chain)
            if not alignment["ready"]:
                not_ready.append({
                    "id": parent_id,
                    "group_id": card.get("group_id"),
                    "record_id": record.get("record_id"),
                    "reason": "Protein alignment is not ready.",
                    "alignment": alignment,
                })
                continue
            transformed = _transform_molecule(ligand, alignment["transform"])
            if _canonical(transformed) != _canonical(ligand):
                raise ValueError("Transformed ligand graph is not identical to source ligand")

            hashes = {
                "catalog_json_sha256": _sha256_path(catalog_path),
                "bound_sdf_sha256": bound_hash,
                "source_pdb_mmcif_sha256": verified_source_hash,
                "reference_6HAZ_mmcif_sha256": reference_hash,
            }
            curated = _curated_fx5(reference_cif, parent_id, record, ligand, hashes)
            metadata = _metadata(
                card, record, parent_id, alignment, hashes, chain_choice, curated,
                observed_charge,
            )
            sdf_data = _sdf_bytes(transformed)
            metadata_data = _json_bytes(metadata)
            prepared.append({
                "id": parent_id,
                "mol": transformed,
                "sdf": sdf_data,
                "metadata": metadata,
                "metadata_bytes": metadata_data,
            })
        except (FileNotFoundError, ValueError, TypeError, KeyError) as exc:
            # Integrity and path failures are operational failures and must not be
            # silently converted into scientific ineligibility.
            if isinstance(exc, FileNotFoundError) or "Hash mismatch" in str(exc) or "Unsafe" in str(exc) or "escapes" in str(exc):
                raise
            not_ready.append({
                "group_id": card.get("group_id") if isinstance(card, dict) else None,
                "reason": str(exc),
            })

    index = {
        "format_version": _FORMAT_VERSION,
        "default_parent_id": "SMARCA2-FX5" if any(p["id"] == "SMARCA2-FX5" for p in prepared) else None,
        "reference": {"pdb": "6HAZ", "chain": "A", "sha256": reference_hash},
        "parents": [
            {
                "id": parent["id"],
                "sdf": f"parents/{parent['id']}.sdf",
                "metadata": f"parents/{parent['id']}.metadata.json",
                "ready": True,
                "suitability": parent["metadata"]["suitability"],
                "limitations": parent["metadata"]["limitations"],
            }
            for parent in sorted(prepared, key=lambda item: item["id"])
        ],
        "not_ready": not_ready,
        "limitations": [
            "This index records deposited co-crystal ligands after validated rigid receptor alignment.",
            "It does not invent SAR, docking scores, binding measurements, or automatic attachment approvals.",
        ],
    }

    index_data = _json_bytes(index)
    files = {"index.json": index_data}
    for parent in prepared:
        files[f"parents/{parent['id']}.sdf"] = parent["sdf"]
        files[f"parents/{parent['id']}.metadata.json"] = parent["metadata_bytes"]
    manifest = {
        "format_version": _FORMAT_VERSION,
        "hash_algorithm": "sha256",
        "files": [
            {"path": path, "bytes": len(data), "sha256": _sha256_bytes(data)}
            for path, data in sorted(files.items())
        ],
        "manifest_self_hash_policy": "manifest.json is excluded because a file cannot contain its own stable cryptographic hash",
    }
    files["manifest.json"] = _json_bytes(manifest)

    output.mkdir(parents=True)
    (output / "parents").mkdir()
    for relative, data in sorted(files.items()):
        destination = output / relative
        with destination.open("xb") as stream:
            stream.write(data)
    return index


def _load_verified_manifest(source_root):
    root = Path(source_root).resolve()
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("format_version") != _FORMAT_VERSION or manifest.get("hash_algorithm") != "sha256":
        raise ValueError("Unsupported reference-parent manifest")
    verified = {}
    entries = manifest.get("files", [])
    if not isinstance(entries, list) or not entries:
        raise ValueError("Reference-parent manifest must contain a nonempty files list")
    for entry in entries:
        if not isinstance(entry, dict) or entry.get("path") in verified:
            raise ValueError("Reference-parent manifest contains an invalid or duplicate file entry")
        path = _safe_relative(root, entry.get("path"))
        expected = str(entry.get("sha256", "")).lower()
        if not _SHA256.fullmatch(expected):
            raise ValueError("Manifest contains an invalid SHA-256")
        if path.stat().st_size != int(entry.get("bytes", -1)) or _sha256_path(path) != expected:
            raise ValueError(f"Manifest verification failed for {entry.get('path')}")
        verified[entry["path"]] = path
    if "index.json" not in verified:
        raise ValueError("Manifest does not cover index.json")
    return root, verified


def reference_source_binding(source_root):
    """Return hashes for every manifest-covered parent input, or an explicit absence."""
    root = Path(source_root)
    if not root.exists():
        return {
            "status": "not_configured",
            "format_version": _FORMAT_VERSION,
            "files": {},
        }
    verified_root, verified = _load_verified_manifest(root)
    manifest_path = verified_root / "manifest.json"
    files = {
        "reference_parents/manifest.json": {
            "status": "configured_hash_verified",
            "sha256": _sha256_path(manifest_path),
        }
    }
    for relative, path in sorted(verified.items()):
        files["reference_parents/" + relative] = {
            "status": "configured_hash_verified",
            "sha256": _sha256_path(path),
        }
    return {
        "status": "configured_hash_verified",
        "format_version": _FORMAT_VERSION,
        "files": files,
    }


def reference_parent_catalog(source_root):
    """Expose only technically verified ready parents as selectable inputs."""
    default = {
        "id": "SMARCA2-FX5",
        "label": "SMARCA2-FX5 · original neutral design parent",
        "ready": True,
        "source": "cases/design_sources/SMARCA2-neutral-design.sdf",
        "technical_input_ready": True,
        "scientific_approval": False,
        "limitations": [
            "This is the existing explicitly neutral platform parent and retains the existing curated SAR.",
            "Technical readiness is not scientific approval.",
        ],
    }
    binding = reference_source_binding(source_root)
    result = {
        "status": binding["status"],
        "default_parent_id": default["id"],
        "available": [default],
        "not_ready": [],
        "excluded": [],
        "limitations": [
            "Ready means hash-verified technical input, not affinity, SAR, attachment, or scientific approval.",
            "Literature SAR is parent-specific and is never transferred automatically between parent graphs.",
            "All literature poses use a validated rigid transform into the fixed 6HAZ chain-A receptor frame.",
        ],
        "binding": binding,
    }
    if binding["status"] != "configured_hash_verified":
        result["limitations"].append(
            "cases/reference_parents is not configured; only the original neutral SMARCA2-FX5 parent is selectable."
        )
        return result
    root, verified = _load_verified_manifest(source_root)
    index = json.loads(verified["index.json"].read_text(encoding="utf-8"))
    for entry in index.get("parents", []):
        if not isinstance(entry, dict) or entry.get("ready") is not True:
            continue
        parent_id = entry.get("id")
        if parent_id == default["id"]:
            result["excluded"].append({
                "id": parent_id,
                "reason": "The deposited charged reference is retained as provenance but is not a separate selectable replacement for the original neutral default.",
            })
            continue
        expected = {
            f"parents/{parent_id}.sdf",
            f"parents/{parent_id}.metadata.json",
        }
        if not isinstance(parent_id, str) or not _PARENT_ID.fullmatch(parent_id) or not expected.issubset(verified):
            raise ValueError("Reference-parent index contains an unsafe or uncovered ready parent")
        result["available"].append({
            "id": parent_id,
            "label": parent_id,
            "ready": True,
            "source": "hash_verified_aligned_literature_parent",
            "technical_input_ready": True,
            "scientific_approval": False,
            "suitability": entry.get("suitability", {}),
            "limitations": entry.get("limitations", []),
        })
    result["not_ready"] = index.get("not_ready", []) if isinstance(index.get("not_ready", []), list) else []
    return result


def load_reference_parent(parent_id, source_root):
    """Load a parent selected only through a fully hash-verified index."""
    if not isinstance(parent_id, str) or not _PARENT_ID.fullmatch(parent_id):
        raise ValueError("Invalid parent id")
    root, verified = _load_verified_manifest(source_root)
    index = json.loads(verified["index.json"].read_text(encoding="utf-8"))
    matches = [entry for entry in index.get("parents", []) if entry.get("id") == parent_id]
    if len(matches) != 1:
        raise KeyError(f"Unknown reference parent id: {parent_id}")
    entry = matches[0]
    expected_sdf = f"parents/{parent_id}.sdf"
    expected_metadata = f"parents/{parent_id}.metadata.json"
    if entry.get("sdf") != expected_sdf or entry.get("metadata") != expected_metadata:
        raise ValueError("Index parent paths are not the fixed safe support paths")
    if expected_sdf not in verified or expected_metadata not in verified:
        raise ValueError("Manifest does not cover all selected parent files")
    mol = _read_single_sdf(verified[expected_sdf])
    meta = json.loads(verified[expected_metadata].read_text(encoding="utf-8"))
    if meta.get("id") != parent_id:
        raise ValueError("Parent metadata id mismatch")
    if _canonical(mol) != _canonical(Chem.MolFromSmiles(meta["canonical_isomeric_smiles"])):
        raise ValueError("Loaded parent graph disagrees with metadata identity")
    return mol, meta


__all__ = [
    "apply_transform",
    "align_cif_chains",
    "build_reference_parents",
    "global_sequence_alignment",
    "kabsch_transform",
    "load_reference_parent",
    "reference_parent_catalog",
    "reference_source_binding",
]
