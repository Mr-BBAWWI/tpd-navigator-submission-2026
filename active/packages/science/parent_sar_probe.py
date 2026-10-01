"""Offline, source-bound audit and exact observed-pair SAR probe.

The module verifies immutable parent exports before reading their results.  It only
constructs the observed SMI-6080 to SMI-6085 chlorine-to-bromine pair and never
promotes the transformation into general design policy or scientific acceptance.
"""
from __future__ import annotations

import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem

FORMAT = "parent-sar-probe/1.0"
EXPECTED_PARENT_EXPORTS = 9
PAIR_ID = "SMI-6080__SMI-6085"
PARENT_ID = "SMARCA2-9D12-A1A1P"


class EvidenceError(ValueError):
    """Raised when immutable source evidence is absent, inconsistent, or tampered."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise EvidenceError(f"INVALID_JSON:{path}") from exc
    if not isinstance(value, dict):
        raise EvidenceError(f"JSON_OBJECT_REQUIRED:{path}")
    return value


_ARTIFACT_ID = re.compile(r"a-[a-f0-9]+\Z")
_SHA256 = re.compile(r"[a-f0-9]{64}\Z")


def _manifest_paths(root: Path) -> list[Path]:
    return sorted(path for path in root.glob("*/manifest.json") if path.is_file())


def _contained_file(path: Path, export_dir: Path, label: str) -> Path:
    try:
        resolved = path.resolve(strict=True)
        resolved.relative_to(export_dir.resolve(strict=True))
    except (OSError, ValueError) as exc:
        raise EvidenceError(f"{label}_OUTSIDE_EXPORT:{path}") from exc
    if not resolved.is_file():
        raise EvidenceError(f"{label}_NOT_FILE:{path}")
    return resolved


def _verified_blob(export_dir: Path, artifact_id: Any, digest: Any) -> Path:
    if not isinstance(artifact_id, str) or _ARTIFACT_ID.fullmatch(artifact_id) is None:
        raise EvidenceError(f"INVALID_ARTIFACT_ID:{artifact_id}")
    if not isinstance(digest, str) or _SHA256.fullmatch(digest) is None:
        raise EvidenceError("INVALID_SHA256")
    path = export_dir / "blobs" / artifact_id
    resolved = _contained_file(path, export_dir, "BLOB")
    if resolved.name != artifact_id or resolved.parent != (export_dir / "blobs").resolve():
        raise EvidenceError(f"INVALID_BLOB_LOCATION:{path}")
    if sha256_file(resolved) != digest:
        raise EvidenceError(f"TAMPERED_BLOB:{path}")
    return resolved


def _artifact_refs(value: Any) -> Iterable[dict[str, Any]]:
    if isinstance(value, dict):
        keys = {"artifact_id", "version", "sha256"}
        if "artifact_id" in value:
            if not keys.issubset(value):
                raise EvidenceError("INCOMPLETE_ARTIFACT_REFERENCE")
            yield value
        for child in value.values():
            yield from _artifact_refs(child)
    elif isinstance(value, list):
        for child in value:
            yield from _artifact_refs(child)


def _candidate_id(record: dict[str, Any]) -> str | None:
    for key in ("candidate_id", "warhead_id", "analog_id", "compound_id", "id", "title"):
        value = record.get(key)
        if isinstance(value, str) and value:
            return value
    smiles = record.get("mapped_smiles") or record.get("canonical_smiles")
    if isinstance(smiles, str) and smiles:
        return "graph-sha256:" + hashlib.sha256(smiles.encode("utf-8")).hexdigest()
    return None


def _family(value: Any) -> str | None:
    if not isinstance(value, str) or not value:
        return None
    if value in {"ring_expansion", "ring_contraction"}:
        return "ring_modification"
    return value


def _walk(value: Any, path: tuple[str, ...] = ()) -> Iterable[tuple[tuple[str, ...], Any]]:
    yield path, value
    if isinstance(value, dict):
        for key, child in value.items():
            yield from _walk(child, path + (str(key),))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from _walk(child, path + (str(index),))


def _row_family(row: dict[str, Any]) -> str | None:
    return _family(row.get("broad_family") or row.get("transformation_class"))


def _catalog_families(result: dict[str, Any]) -> set[str]:
    catalog = result.get("rule_catalog")
    values = catalog.values() if isinstance(catalog, dict) else catalog if isinstance(catalog, list) else []
    families: set[str] = set()
    for value in values:
        family = _row_family(value) if isinstance(value, dict) else _family(value)
        if family:
            families.add(family)
    return families


def _records(result: dict[str, Any]) -> list[dict[str, Any]]:
    analogs = result.get("analogs")
    return [row for row in analogs if isinstance(row, dict)] if isinstance(analogs, list) else []


def _explicit_cheap_pass(record: dict[str, Any]) -> bool:
    value = record.get("cheap_filter")
    return isinstance(value, dict) and value.get("valid") is True


def _explicit_pose_qualified(record: dict[str, Any]) -> bool:
    docking = record.get("docking")
    return (
        record.get("pipeline_status") == "qualified"
        and isinstance(docking, dict)
        and docking.get("status") == "completed_with_limits"
        and docking.get("pose_preserved") is True
        and record.get("parent_redocking_supported") is True
    )


def _summary_counts(result: dict[str, Any]) -> dict[str, int]:
    found: dict[str, int] = {}
    summary = result.get("summary")
    for path, value in _walk(summary if isinstance(summary, dict) else {}):
        if path and path[-1].endswith("count") and type(value) is int:
            found["/".join(path)] = value
    return found


def _summarise_result(result: dict[str, Any]) -> dict[str, Any]:
    records = _records(result)
    cheap = [row for row in records if _explicit_cheap_pass(row)]
    strict = [row for row in cheap if row.get("qualified_for_counts") is True]
    exploratory = [row for row in cheap if row.get("qualified_for_counts") is not True]
    pose = [row for row in records if _explicit_pose_qualified(row)]
    selected = [row for row in pose if row.get("selected") is True]
    strict_families = sorted(filter(None, {_row_family(row) for row in strict}))
    rejections = result.get("rejections") if isinstance(result.get("rejections"), list) else []
    attachment = [row for row in records if row.get("pipeline_status") == "rejected_no_supported_attachment"]
    errors = [
        {"candidate_id": _candidate_id(row), "status": (row.get("docking") or {}).get("status"), "error": (row.get("docking") or {}).get("error")}
        for row in records
        if isinstance(row.get("docking"), dict)
        and ((row["docking"].get("error") is not None) or str(row["docking"].get("status", "")).startswith("failed"))
    ]
    actual_counts = {
        "cheap_passed_count": len(cheap),
        "strict_qualified_count": len(strict),
        "pose_qualified_count": len(pose),
        "selected_count": len(selected),
        "attachment_exclusion_count": len(attachment),
    }
    stored = _summary_counts(result)
    aliases = {key: value for key, value in actual_counts.items()}
    mismatches = [
        {"path": path, "stored": value, "actual": aliases[path.rsplit("/", 1)[-1]]}
        for path, value in stored.items()
        if path.rsplit("/", 1)[-1] in aliases and value != aliases[path.rsplit("/", 1)[-1]]
    ]
    return {
        "catalog_families": sorted(_catalog_families(result)),
        "cheap_passed_families": sorted(filter(None, {_row_family(row) for row in cheap})),
        "strict_qualified_families": strict_families,
        "exploratory_cheap_passed_families": sorted(filter(None, {_row_family(row) for row in exploratory})),
        "pose_qualified_families": sorted(filter(None, {_row_family(row) for row in pose})),
        "cheap_passed_candidate_ids": sorted(filter(None, (_candidate_id(row) for row in cheap))),
        "strict_qualified_candidate_ids": sorted(filter(None, (_candidate_id(row) for row in strict))),
        "exploratory_candidate_ids": sorted(filter(None, (_candidate_id(row) for row in exploratory))),
        "pose_qualified_candidate_ids": sorted(filter(None, (_candidate_id(row) for row in pose))),
        "selected_candidate_ids": sorted(filter(None, (_candidate_id(row) for row in selected))),
        "qualified_for_counts_policy": "Potential family coverage only: literal cheap_filter.valid=true and qualified_for_counts=true are required; this is not authenticated policy or scientific approval.",
        "cheap_passed_count": len(cheap),
        "strict_qualified_family_count": len(strict_families),
        "per_parent_at_least_6_gate": len(strict_families) >= 6,
        "stored_summary_counts": stored,
        "stored_summary_count_mismatches": mismatches,
        "retained_docking_errors": errors,
        "retained_docking_error_count": len(errors),
        "rejections": rejections,
        "attachment_exclusions": [{"candidate_id": _candidate_id(row), "pipeline_status": row.get("pipeline_status")} for row in attachment],
        "attachment_exclusion_count": len(attachment),
    }


def verify_parent_export(manifest_path: Path, root: Path) -> dict[str, Any]:
    export_dir = manifest_path.parent.resolve(strict=True)
    if manifest_path.name != "manifest.json" or export_dir.parent != root.resolve(strict=True):
        raise EvidenceError("MANIFEST_NOT_IMMEDIATE_EXPORT_CHILD")
    manifest = _load_json(_contained_file(manifest_path, export_dir, "MANIFEST"))
    receipt_digest = manifest.get("raw_receipt_sha256")
    if not isinstance(receipt_digest, str) or _SHA256.fullmatch(receipt_digest) is None:
        raise EvidenceError("RAW_RECEIPT_SHA256_REQUIRED")
    receipt_path = _contained_file(export_dir / "job-receipt.json", export_dir, "RECEIPT")
    if sha256_file(receipt_path) != receipt_digest:
        raise EvidenceError("JOB_RECEIPT_HASH_MISMATCH")
    receipt = _load_json(receipt_path)
    if receipt.get("state") != "completed" or receipt.get("stage") != "finished":
        raise EvidenceError("PARENT_JOB_NOT_COMPLETED")
    job_id = manifest.get("job_id")
    parent_id = manifest.get("parent_id")
    if not isinstance(job_id, str) or receipt.get("id") != job_id:
        raise EvidenceError("JOB_ID_MISMATCH")
    if not isinstance(parent_id, str) or (receipt.get("parameters") or {}).get("parent_id") != parent_id:
        raise EvidenceError("PARENT_ID_MISMATCH")

    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise EvidenceError("ARTIFACTS_REQUIRED")
    declarations: dict[str, tuple[Any, str]] = {}
    closure: list[dict[str, Any]] = []
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise EvidenceError("INVALID_ARTIFACT_RECORD")
        artifact_id, version, digest = artifact.get("artifact_id"), artifact.get("version"), artifact.get("sha256")
        if artifact_id in declarations:
            raise EvidenceError(f"DUPLICATE_ARTIFACT_ID:{artifact_id}")
        path = _verified_blob(export_dir, artifact_id, digest)
        declarations[artifact_id] = (version, digest)
        closure.append({"artifact_id": artifact_id, "version": version, "sha256": digest, "blob_path": str(path), "bytes": path.stat().st_size})
    for ref in _artifact_refs({"runtime": receipt.get("runtime"), "result": receipt.get("result")}):
        declared = declarations.get(ref.get("artifact_id"))
        if declared is None or declared != (ref.get("version"), ref.get("sha256")):
            raise EvidenceError(f"UNDECLARED_OR_MISMATCHED_ARTIFACT_REF:{ref.get('artifact_id')}")
    result = receipt.get("result")
    if not isinstance(result, dict):
        raise EvidenceError("RECEIPT_RESULT_OBJECT_REQUIRED")
    result_parent = (result.get("parameters") or {}).get("parent_id") or (result.get("parent_scope") or {}).get("actual_design_parent_id")
    if result_parent != parent_id:
        raise EvidenceError("RESULT_PARENT_ID_MISMATCH")
    return {
        "manifest_path": str(manifest_path.resolve()),
        "manifest_sha256": sha256_file(manifest_path),
        "receipt_path": str(receipt_path),
        "raw_receipt_sha256": receipt_digest,
        "job_id": job_id,
        "parent_id": parent_id,
        "artifact_blob_closure": closure,
        "result": result,
    }


def summarise_exports(export_root: str | Path, expected_count: int = EXPECTED_PARENT_EXPORTS) -> dict[str, Any]:
    """Verify and summarize unique parent export directories.

    By default exactly nine exports are required.  Set ``expected_count`` explicitly
    only for isolated tests or a deliberately scoped audit.
    """
    root = Path(export_root).resolve()
    if not root.is_dir():
        raise EvidenceError("EXPORT_ROOT_NOT_DIRECTORY")
    manifests = _manifest_paths(root)
    by_directory: dict[Path, Path] = {}
    for path in manifests:
        directory = path.parent.resolve()
        if directory in by_directory:
            raise EvidenceError(f"MULTIPLE_MANIFESTS_FOR_PARENT_DIRECTORY:{directory}")
        by_directory[directory] = path
    if len(by_directory) != expected_count:
        raise EvidenceError(f"EXPECTED_{expected_count}_UNIQUE_PARENT_EXPORTS_FOUND_{len(by_directory)}")

    parents: list[dict[str, Any]] = []
    parent_ids: set[str] = set()
    for path in sorted(by_directory.values()):
        verified = verify_parent_export(path, root)
        parent_id = verified["parent_id"]
        if parent_id in parent_ids:
            raise EvidenceError(f"DUPLICATE_PARENT_ID:{parent_id}")
        parent_ids.add(parent_id)
        summary = _summarise_result(verified.pop("result"))
        parents.append({**verified, **summary})
    pooled = sorted({family for parent in parents for family in parent["cheap_passed_families"]})
    return {
        "format": FORMAT,
        "status": "verified",
        "parent_count": len(parents),
        "parents": parents,
        "pooled_cheap_passed_family_union": pooled,
        "pooled_union_satisfies_per_parent_gate": False,
        "per_parent_gate_policy": "A pooled union is informational and never satisfies a parent's independent >=6 distinct strict broad-family potential-coverage gate; it is not scientific approval.",
        "all_parents_meet_at_least_6_gate": all(parent["per_parent_at_least_6_gate"] for parent in parents),
        "scientific_approval": False,
        "changes_design_policy": False,
    }


def _source_data(medchem_source: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    if not isinstance(medchem_source, dict) or medchem_source.get("status") != "configured_hash_verified":
        raise EvidenceError("MEDCHEM_SOURCE_NOT_HASH_VERIFIED")
    data = medchem_source.get("data")
    if not isinstance(data, dict):
        raise EvidenceError("MEDCHEM_SOURCE_DATA_MISSING")
    reference = {
        "relative_path": medchem_source.get("relative_path"),
        "sha256": medchem_source.get("sha256"),
    }
    if not isinstance(reference["sha256"], str) or _SHA256.fullmatch(reference["sha256"]) is None:
        raise EvidenceError("MEDCHEM_SOURCE_SHA256_MISSING")
    return data, reference


def _record(data: dict[str, Any], compound_id: str) -> dict[str, Any]:
    records = ((data.get("parent_catalog") or {}).get("records") or [])
    matches = [row for row in records if isinstance(row, dict) and row.get("compound_id") == compound_id]
    if len(matches) != 1:
        raise EvidenceError(f"SOURCE_RECORD_COUNT_{compound_id}:{len(matches)}")
    return matches[0]


def _clean_graph(mol: Chem.Mol) -> str:
    clean = Chem.RemoveHs(Chem.Mol(mol), sanitize=True)
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in clean.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
        bond.SetBondDir(Chem.BondDir.NONE)
    Chem.SanitizeMol(clean)
    return Chem.MolToSmiles(clean, canonical=True, isomericSmiles=False)


def unresolved_potential_stereo(mol: Chem.Mol) -> list[dict[str, Any]]:
    clean = Chem.RemoveHs(Chem.Mol(mol), sanitize=True)
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(clean, cleanIt=True, force=True)
    return [
        {"type": str(info.type), "centered_on": int(info.centeredOn)}
        for info in Chem.FindPotentialStereo(clean)
        if str(info.specified) == "Unspecified"
    ]


def build_exact_halogen_probe(parent: Chem.Mol, linkage: dict[str, Any], medchem_source: dict[str, Any]) -> dict[str, Any]:
    """Build only the exact observed SMI-6080/SMI-6085 Cl-to-Br graph."""
    if not isinstance(parent, Chem.Mol):
        raise TypeError("parent must be an RDKit Mol")
    if linkage.get("status") != "exact_measured_parent_join" or linkage.get("selected_parent_id") != PARENT_ID:
        raise EvidenceError("SELECTED_PARENT_LINKAGE_REQUIRED")
    matched = linkage.get("matched_record") or {}
    if matched.get("compound_id") != "SMI-6080":
        raise EvidenceError("LINKAGE_PARENT_MUST_BE_SMI_6080")
    linked_pairs = [
        row for row in linkage.get("source_sar", [])
        if isinstance(row, dict) and row.get("pair_id") == PAIR_ID
        and row.get("requested_change_label") == "chlorine_to_bromine"
        and row.get("other_compound_id") == "SMI-6085"
    ]
    if len(linked_pairs) != 1 or linked_pairs[0].get("selected_parent_affected_atom_maps") != [9]:
        raise EvidenceError("EXACT_LINKAGE_PAIR_OR_MAP9_MISSING")
    data, source_reference = _source_data(medchem_source)
    left = _record(data, "SMI-6080")
    right = _record(data, "SMI-6085")
    pairs = [row for row in data.get("sar_pair_evidence", []) if isinstance(row, dict) and row.get("pair_id") == PAIR_ID]
    if len(pairs) != 1 or pairs[0].get("requested_change_label") != "chlorine_to_bromine":
        raise EvidenceError("SOURCE_PAIR_MISSING")

    left_mol = Chem.MolFromSmiles(str((left.get("identity") or {}).get("source_smiles")))
    right_mol = Chem.MolFromSmiles(str((right.get("identity") or {}).get("source_smiles")))
    if left_mol is None or right_mol is None:
        raise EvidenceError("SOURCE_SMILES_INVALID")
    parent_copy = Chem.Mol(parent)
    Chem.SanitizeMol(parent_copy)
    if _clean_graph(parent_copy) != _clean_graph(left_mol) or Chem.GetFormalCharge(parent_copy) != Chem.GetFormalCharge(left_mol):
        raise EvidenceError("PARENT_GRAPH_OR_CHARGE_NOT_SMI_6080")

    by_map: dict[int, Chem.Atom] = {}
    for atom in parent_copy.GetAtoms():
        number = atom.GetAtomMapNum()
        if number <= 0 or number in by_map:
            raise EvidenceError("UNIQUE_POSITIVE_PARENT_MAPS_REQUIRED")
        by_map[number] = atom
    anchor = by_map.get(9)
    if anchor is None or anchor.GetAtomicNum() != 6 or not anchor.GetIsAromatic():
        raise EvidenceError("MAP9_MUST_BE_AROMATIC_CARBON")
    terminal_chlorines = [atom for atom in anchor.GetNeighbors() if atom.GetAtomicNum() == 17 and atom.GetDegree() == 1]
    if len(terminal_chlorines) != 1:
        raise EvidenceError("MAP9_REQUIRES_ONE_TERMINAL_CHLORINE_NEIGHBOR")
    halogen = terminal_chlorines[0]
    if halogen.GetAtomMapNum() != 8:
        raise EvidenceError("DISCOVERED_TERMINAL_CHLORINE_MUST_BE_MAP8")

    candidate = Chem.RWMol(parent_copy)
    candidate.GetAtomWithIdx(halogen.GetIdx()).SetAtomicNum(35)
    product = candidate.GetMol()
    Chem.SanitizeMol(product)
    if Chem.GetFormalCharge(product) != Chem.GetFormalCharge(right_mol):
        raise EvidenceError("PRODUCT_CHARGE_MISMATCH")
    if _clean_graph(product) != _clean_graph(right_mol):
        raise EvidenceError("UNEXPECTED_PRODUCT_GRAPH")
    for number, before in by_map.items():
        after = next((atom for atom in product.GetAtoms() if atom.GetAtomMapNum() == number), None)
        if after is None:
            raise EvidenceError("ATOM_MAP_LOST")
        expected_atomic_number = 35 if number == 8 else before.GetAtomicNum()
        if after.GetAtomicNum() != expected_atomic_number:
            raise EvidenceError("UNDECLARED_ELEMENT_CHANGE")

    protected_maps = set()
    for row in linkage.get("sites", []):
        if isinstance(row, dict) and row.get("state") == "PROTECTED":
            protected_maps.add(int(row["atom_map"]))
    protected_maps.update(int(value) for value in linkage.get("protected_atom_maps", []) if not isinstance(value, bool))
    conflict = 8 in protected_maps
    source_pair = pairs[0]
    return {
        "format": FORMAT,
        "proposal_id": PAIR_ID + "__exact_observed_probe",
        "selected_parent_id": PARENT_ID,
        "source_pair_id": PAIR_ID,
        "requested_change_label": "chlorine_to_bromine",
        "left_compound_id": "SMI-6080",
        "right_compound_id": "SMI-6085",
        "affected_anchor_atom_map": 9,
        "changed_atom_map": 8,
        "changed_element": {"from": "Cl", "to": "Br"},
        "mapped_smiles": Chem.MolToSmiles(product, canonical=True, isomericSmiles=True),
        "canonical_map_free_nonstereo_graph": _clean_graph(product),
        "formal_charge": int(Chem.GetFormalCharge(product)),
        "graph_verified_against_source_SMI_6085": True,
        "source_row_locators": source_pair.get("source_row_locators", linked_pairs[0].get("source_row_locators", [])),
        "source_reference": source_reference,
        "observed_source_records": {"SMI-6080": copy.deepcopy(left), "SMI-6085": copy.deepcopy(right)},
        "measured_potency_status": "observed_source_records_only_no_transfer",
        "original_source_reference": {
            "path": ((data.get("freeze_provenance") or {}).get("original_ref")),
            "sha256": ((data.get("freeze_provenance") or {}).get("original_content_sha256")),
        },
        "protected_map_conflict": conflict,
        "protected_conflict_atom_maps": [8] if conflict else [],
        "parent_diagnostic_mask": sorted(protected_maps),
        "candidate_diagnostic_mask": sorted(protected_maps - {8}),
        "candidate_mask_exclusion": {"atom_map": 8, "reason": "only actually changed atom"},
        "strict_gate_comparable": False,
        "qualified_for_counts": False,
        "exploratory": True,
        "changes_design_policy": False,
        "scientific_approval": False,
        "protected_atom_activation": False,
        "measurement_transfer_to_enumerated_stereo_or_docked_candidate": False,
        "stereochemistry_caveat": "Parent stereochemistry is model/reference-assigned; the source CSV is stereochemically unspecified. Measured potency is not transferred to an enumerated stereo form or docked candidate.",
        "parent_unresolved_potential_stereo": unresolved_potential_stereo(parent_copy),
    }


def validate_predeclared_seeds(requested: Iterable[int], results: Iterable[dict[str, Any]], ligand_labels: Iterable[str] = ("parent", "observed_br")) -> None:
    seeds = list(requested)
    labels = list(ligand_labels)
    if not seeds or any(type(seed) is not int or not 0 <= seed <= 0xFFFFFFFF for seed in seeds) or len(seeds) != len(set(seeds)):
        raise ValueError("SEEDS_MUST_BE_UNIQUE_NONNEGATIVE_32BIT_INTEGERS")
    if not labels or any(type(label) is not str or not label for label in labels) or len(labels) != len(set(labels)):
        raise ValueError("LIGAND_LABELS_MUST_BE_UNIQUE_NONEMPTY_STRINGS")
    rows = list(results)
    if any(not isinstance(row, dict) or type(row.get("seed")) is not int or type(row.get("ligand")) is not str for row in rows):
        raise EvidenceError("INVALID_RETAINED_SEED_RECORD")
    pairs = [(row["ligand"], row["seed"]) for row in rows]
    if len(pairs) != len(set(pairs)):
        raise EvidenceError("DUPLICATE_RETAINED_LIGAND_SEED_PAIR")
    expected = {(label, seed) for label in labels for seed in seeds}
    if set(pairs) != expected:
        raise EvidenceError("NOT_EXACTLY_ALL_PREDECLARED_SEEDS_AND_LABELS_RETAINED")


__all__ = [
    "EvidenceError", "FORMAT", "PAIR_ID", "PARENT_ID", "build_exact_halogen_probe",
    "sha256_file", "summarise_exports", "unresolved_potential_stereo",
    "validate_predeclared_seeds", "verify_parent_export",
]
