"""Source-bound, narrowly scoped handling of the relayed FX5 opinion.

This module does not authenticate a reviewer and cannot create scientific or
formal approval. It operates only on an archived result that the service has
already verified through its normal immutable-result path.
"""
from __future__ import annotations

import hashlib
import io
import json
import zipfile
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

from rdkit import Chem

from packages.contracts import ContractError

VERSION = "source-bound-expert-opinion/20260930.2"
ROOT = Path(__file__).resolve().parents[2]
MANIFEST_RELATIVE_PATH = "cases/expert_opinions/20260930/opinion.json"
SOURCE_RELATIVE_PATH = "cases/expert_opinions/20260930/reply.docx"
SOURCE_SHA256 = "42a78a5540f5549347bc8ee2af1fd9758cec43056c6e808f89ff9aedc6d521cb"
SOURCE_BYTES = 46983
DECISION_LOCATOR = "/body/7/row/4"
CONDITIONS_LOCATOR_CANDIDATE = "/body/7/row/7"

FX5_SOURCE_FILES = {
    "design_sources/6HAZ.cif": "7752dda8c54b7e56ab45c84e03ca544e383d834fb30e4c0849fc3544278e2c9d",
    "design_sources/SMARCA2-neutral-design.sdf": "50e689f68d455083f2d7432b5928a93f6f7f971c4f2fc3cfda911bc0164ac6ee",
    "design_sources/SMARCA2-parent.sdf": "7b39c8e3f36f1c7a94c02a4ddda6c1cda1d672d3746094bb8572b81cd298e9b9",
    "design_sources/warhead.json": "0398358d87a7f6e354d51877fefb61a0ebc9145081e02d3dfda296c266c1629f",
}

EXPECTED_DESIGN_IDENTITY = {
    "canonical_isomeric_smiles": "Nc1nnc(-c2ccccc2O)cc1N1CCNCC1",
    "inchikey": "SZKHGLTYXIDOFH-UHFFFAOYSA-N",
    "formal_charge": 0,
    "heavy_atom_count": 20,
}
EXPECTED_PARENT_IDENTITY = {
    "canonical_isomeric_smiles": "Nc1nnc(-c2ccccc2O)cc1N1CC[NH2+]CC1",
    "inchikey": "SZKHGLTYXIDOFH-UHFFFAOYSA-O",
    "formal_charge": 1,
    "heavy_atom_count": 20,
}
EXPECTED_COMMON_MAPS = list(range(1, 21))
EXPECTED_CORE_MAPS = [1, 2, 3, 4, 5, 6, 7, 9, 10, 13, 14, 15, 16, 17, 18, 20]


def _expected_manifest() -> dict[str, Any]:
    return {
        "manifest_type": "source_bound_relayed_expert_opinion",
        "manifest_version": "20260930.2",
        "attribution": "Reviewer A",
        "origin": "user-relayed-document",
        "review_date": "2026-09-30",
        "authority": {
            "authenticated_reviewer": False,
            "direct_decision": False,
            "digital_signature": False,
            "pki_verified": False,
            "formal_approval": False,
        },
        "source": {
            "path": SOURCE_RELATIVE_PATH,
            "sha256": SOURCE_SHA256,
            "bytes": SOURCE_BYTES,
            "decision_locator": DECISION_LOCATOR,
            "conditions_locator_candidate": CONDITIONS_LOCATOR_CANDIDATE,
        },
        "fx5_scope": {
            "parent_id": "SMARCA2-FX5",
            "receptor_frame": "6HAZ chain A",
            "pdb": "6HAZ",
            "target_chain": "A",
            "only_modifiable_atom_map": 19,
            "must_remain_unknown_atom_maps": [8, 11, 12],
            "modifiable_sites_min": 1,
            "distinct_graphs_per_site_min": 30,
        },
        "numerical_defaults": {
            "calibration_seeds_min": 5,
            "novel_ternary_repeats_min": 5,
            "source_kind": "relayed_expert_default_minimum",
            "not_success_threshold": True,
        },
        "execution_plan": {
            "candidate_analogs": [
                "W-c2afc5e73c1a",
                "W-4c0a639c0a41",
                "W-80f8f4a11b5d",
            ],
            "e3_branches": ["VHL", "CRBN"],
            "representative_linker": "alkyl_c6",
            "representative_linker_selected_by": "developer",
            "seeds": [23, 41, 61, 79, 97],
            "planned_new_ternary_computations": 30,
            "execution_status": "pending",
            "preserve_bad_seed_41": True,
            "planned_new_input_mode": "single_sequence",
            "same_single_sequence_input_across_novel_seeds": True,
            "matched_known_case_msa": False,
            "matched_msa_followup_required_for_comparison": True,
            "scientific_threshold_created": False,
        },
        "unresolved": [
            "parent_funnel_5_to_10",
            "broad_families_6",
            "qualified_panel_10_to_20",
            "final_interaction_fingerprint",
            "final_microstate",
            "exact_synthesis_route",
            "calibration_geometry",
            "formal_authenticated_expert_decision",
        ],
    }


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _safe_file(root: Path, relative: str, error_code: str) -> Path:
    if not isinstance(relative, str) or not relative:
        raise ContractError(error_code)
    relative_path = Path(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise ContractError(error_code)
    try:
        strict_root = root.resolve(strict=True)
        candidate = root / relative_path
        current = candidate
        while current != root and current != current.parent:
            if current.exists() and current.is_symlink():
                raise ContractError(error_code)
            current = current.parent
        path = candidate.resolve(strict=True)
    except ContractError:
        raise
    except OSError as error:
        raise ContractError(error_code) from error
    if path.is_symlink() or not path.is_file() or not (path == strict_root or strict_root in path.parents):
        raise ContractError(error_code)
    return path


def _extract_docx_locator(path: Path, locator: str) -> str:
    parts = locator.strip("/").split("/")
    if len(parts) != 4 or parts[0] != "body" or parts[2] != "row":
        raise ContractError("EXPERT_OPINION_SOURCE_LOCATOR")
    try:
        body_index = int(parts[1])
        row_index = int(parts[3])
    except ValueError as error:
        raise ContractError("EXPERT_OPINION_SOURCE_LOCATOR") from error
    try:
        with zipfile.ZipFile(path) as archive:
            raw = archive.read("word/document.xml")
        root = ElementTree.fromstring(raw)
        namespace = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
        body = root.find(f"{namespace}body")
        if body is None:
            raise ContractError("EXPERT_OPINION_SOURCE_LOCATOR")
        block = list(body)[body_index]
        rows = block.findall(f"{namespace}tr")
        row = rows[row_index]
        text = "".join(node.text or "" for node in row.iter(f"{namespace}t")).strip()
    except ContractError:
        raise
    except (IndexError, KeyError, OSError, zipfile.BadZipFile, ElementTree.ParseError) as error:
        raise ContractError("EXPERT_OPINION_SOURCE_LOCATOR") from error
    if not text:
        raise ContractError("EXPERT_OPINION_SOURCE_LOCATOR")
    return text


def load_bundled_opinion(root: Path = ROOT) -> dict[str, Any]:
    root = Path(root)
    manifest_path = _safe_file(root, MANIFEST_RELATIVE_PATH, "EXPERT_OPINION_MANIFEST_FILE")
    source_path = _safe_file(root, SOURCE_RELATIVE_PATH, "EXPERT_OPINION_SOURCE_FILE")
    raw = manifest_path.read_bytes()
    try:
        manifest = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as error:
        raise ContractError("EXPERT_OPINION_MANIFEST_JSON") from error
    expected = _expected_manifest()
    if type(manifest) is not dict or manifest != expected:
        raise ContractError("EXPERT_OPINION_MANIFEST_MISMATCH")
    if raw != _canonical(expected) + b"\n":
        raise ContractError("EXPERT_OPINION_MANIFEST_BYTES")
    source_raw = source_path.read_bytes()
    if len(source_raw) != SOURCE_BYTES or _sha(source_raw) != SOURCE_SHA256:
        raise ContractError("EXPERT_OPINION_SOURCE_HASH")
    decision_text = _extract_docx_locator(source_path, DECISION_LOCATOR)
    conditions_text = _extract_docx_locator(source_path, CONDITIONS_LOCATOR_CANDIDATE)
    return {
        "manifest": manifest,
        "manifest_sha256": _sha(raw),
        "source_path": SOURCE_RELATIVE_PATH,
        "source_sha256": SOURCE_SHA256,
        "decision_locator": DECISION_LOCATOR,
        "decision_extract_sha256": _sha(decision_text.encode("utf-8")),
        "conditions_locator": CONDITIONS_LOCATOR_CANDIDATE,
        "conditions_extract_sha256": _sha(conditions_text.encode("utf-8")),
    }


def _verify_bound_sources(result: dict[str, Any], root: Path) -> tuple[list[str], list[str]]:
    failures: list[str] = []
    verified: list[str] = []
    binding = result.get("input_binding")
    files = binding.get("source_files") if isinstance(binding, dict) else None
    if not isinstance(files, dict) or not set(FX5_SOURCE_FILES).issubset(files):
        return ["INPUT_BINDING_SOURCE_SET_MISMATCH"], []
    for relative, expected_hash in FX5_SOURCE_FILES.items():
        record = files.get(relative)
        if (
            not isinstance(record, dict)
            or not {"status", "sha256"}.issubset(record)
            or record.get("status") != "configured_hash_verified"
            or record.get("sha256") != expected_hash
        ):
            failures.append("INPUT_BINDING_SOURCE_RECORD_MISMATCH")
            continue
        try:
            path = _safe_file(root / "cases", relative, "EXPERT_OPINION_BOUND_SOURCE_FILE")
        except ContractError:
            failures.append("BOUND_SOURCE_CURRENT_FILE_UNAVAILABLE")
            continue
        if _sha(path.read_bytes()) != expected_hash:
            failures.append("BOUND_SOURCE_CURRENT_HASH_MISMATCH")
            continue
        verified.append(relative)
    return list(dict.fromkeys(failures)), verified


def _single_sdf_molecule(path: Path) -> Chem.Mol:
    try:
        raw = path.read_bytes()
        supplier = Chem.ForwardSDMolSupplier(
            io.BytesIO(raw), removeHs=False, sanitize=True
        )
        molecules: list[Chem.Mol] = []
        invalid_record = False
        for molecule in supplier:
            if molecule is None:
                invalid_record = True
            else:
                molecules.append(molecule)
    except OSError as error:
        raise ContractError("EXPERT_OPINION_SDF_SINGLE_MOLECULE_REQUIRED") from error
    if invalid_record or len(molecules) != 1:
        raise ContractError("EXPERT_OPINION_SDF_SINGLE_MOLECULE_REQUIRED")
    molecule = Chem.Mol(molecules[0])
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return molecule


def _identity_matches(molecule: Chem.Mol, identity: Any, expected: dict[str, Any]) -> bool:
    if not isinstance(identity, dict):
        return False
    observed = {
        "canonical_isomeric_smiles": Chem.MolToSmiles(
            Chem.RemoveHs(molecule), canonical=True, isomericSmiles=True
        ),
        "inchikey": Chem.MolToInchiKey(Chem.RemoveHs(molecule)),
        "formal_charge": Chem.GetFormalCharge(molecule),
        "heavy_atom_count": molecule.GetNumHeavyAtoms(),
    }
    supplied = {key: identity.get(key) for key in expected}
    return observed == expected and supplied == expected


def _graph_and_mapping_failures(result: dict[str, Any], root: Path) -> list[str]:
    failures: list[str] = []
    warhead = result.get("warhead")
    if not isinstance(warhead, dict):
        return ["WARHEAD_MISSING"]
    try:
        design = _single_sdf_molecule(
            _safe_file(root / "cases", "design_sources/SMARCA2-neutral-design.sdf", "EXPERT_OPINION_DESIGN_SDF")
        )
        parent = _single_sdf_molecule(
            _safe_file(root / "cases", "design_sources/SMARCA2-parent.sdf", "EXPERT_OPINION_PARENT_SDF")
        )
    except ContractError:
        return ["PARENT_CANONICAL_GRAPH_UNVERIFIED"]
    if not _identity_matches(design, warhead.get("design_identity"), EXPECTED_DESIGN_IDENTITY):
        failures.append("DESIGN_IDENTITY_GRAPH_MISMATCH")
    if not _identity_matches(parent, warhead.get("parent_identity"), EXPECTED_PARENT_IDENTITY):
        failures.append("PARENT_IDENTITY_GRAPH_MISMATCH")

    mapping = warhead.get("atom_mapping")
    if not isinstance(mapping, list) or len(mapping) != 20:
        failures.append("ATOM_MAPPING_CARDINALITY_MISMATCH")
    else:
        maps, indices = [], []
        for item in mapping:
            if not isinstance(item, dict) or type(item.get("atom_map")) is not int or type(item.get("rdkit_index_zero_based")) is not int:
                failures.append("ATOM_MAPPING_INVALID")
                break
            maps.append(item["atom_map"])
            indices.append(item["rdkit_index_zero_based"])
        if len(set(maps)) != len(maps) or sorted(maps) != EXPECTED_COMMON_MAPS:
            failures.append("ATOM_MAPPING_DUPLICATE_OR_INCOMPLETE")
        if len(set(indices)) != len(indices) or sorted(indices) != list(range(20)):
            failures.append("RDKIT_INDEX_MAPPING_DUPLICATE_OR_INCOMPLETE")
    return failures


def _site_failures(result: dict[str, Any]) -> tuple[list[str], dict[str, str]]:
    failures: list[str] = []
    atoms = result.get("sites", {}).get("atoms") if isinstance(result.get("sites"), dict) else None
    states: dict[int, str] = {}
    if not isinstance(atoms, list) or len(atoms) != 20:
        return ["SITE_MAP_CARDINALITY_MISMATCH"], states
    for atom in atoms:
        if not isinstance(atom, dict) or type(atom.get("atom_map")) is not int:
            failures.append("SITE_MAP_INVALID")
            continue
        atom_map = atom["atom_map"]
        if atom_map in states:
            failures.append("SITE_MAP_DUPLICATE")
        states[atom_map] = atom.get("state")
    if sorted(states) != EXPECTED_COMMON_MAPS:
        failures.append("SITE_MAP_INCOMPLETE")
    if states.get(19) != "MODIFIABLE":
        failures.append("MAP19_NOT_MODIFIABLE")
    if any(state == "MODIFIABLE" and atom_map != 19 for atom_map, state in states.items()):
        failures.append("MAP19_NOT_ONLY_MODIFIABLE_SITE")
    if any(states.get(atom_map) != "UNKNOWN" for atom_map in (8, 11, 12)):
        failures.append("MAP8_11_12_NOT_ALL_UNKNOWN")
    return list(dict.fromkeys(failures)), {str(key): states[key] for key in sorted(states)}


def _redocking_failures(result: dict[str, Any]) -> list[str]:
    failures: list[str] = []
    parent_scope = result.get("parent_scope")
    redocking = result.get("parent_redocking")
    if not isinstance(parent_scope, dict) or type(parent_scope.get("parent_redock_required_for_qualification")) is not bool or parent_scope.get("parent_redock_required_for_qualification") is not True:
        failures.append("PARENT_REDOCK_SCOPE_BOOLEAN_REQUIRED")
    if not isinstance(redocking, dict):
        return failures + ["PARENT_REDOCKING_MISSING"]
    if type(redocking.get("real_docking")) is not bool or redocking.get("real_docking") is not True:
        failures.append("REAL_DOCKING_BOOLEAN_TRUE_REQUIRED")
    if type(redocking.get("parent_redocking")) is not bool or redocking.get("parent_redocking") is not True:
        failures.append("PARENT_REDOCKING_BOOLEAN_TRUE_REQUIRED")
    if redocking.get("status") != "completed_with_limits":
        failures.append("PARENT_REDOCKING_STATUS_MISMATCH")
    results = redocking.get("results")
    preservation = results.get("pose_preservation") if isinstance(results, dict) else None
    if (
        not isinstance(preservation, dict)
        or type(preservation.get("docking_pose_preserved")) is not bool
        or preservation.get("docking_pose_preserved") is not True
        or preservation.get("status") != "pass"
    ):
        failures.append("PARENT_POSE_PRESERVATION_PASS_REQUIRED")
        return failures
    diagnostics = preservation.get("all_poses_diagnostics")
    if not isinstance(diagnostics, list) or not diagnostics:
        return failures + ["PARENT_POSE_DIAGNOSTICS_REQUIRED"]
    indices: list[int] = []
    passing_exact_pose = False
    for diagnostic in diagnostics:
        if not isinstance(diagnostic, dict) or type(diagnostic.get("pose_index_zero_based")) is not int:
            failures.append("PARENT_POSE_DIAGNOSTIC_INVALID")
            continue
        indices.append(diagnostic["pose_index_zero_based"])
        exact_maps = (
            diagnostic.get("common_mapped_heavy_atom_count") == 20
            and diagnostic.get("core_mapped_atom_count") == 16
            and diagnostic.get("common_atom_maps") == EXPECTED_COMMON_MAPS
            and diagnostic.get("core_atom_maps") == EXPECTED_CORE_MAPS
        )
        if not exact_maps:
            failures.append("PARENT_POSE_MAPPING_COUNTS_MISMATCH")
        if exact_maps and diagnostic.get("docking_pose_preserved") is True and diagnostic.get("status") == "pass":
            passing_exact_pose = True
    if len(indices) != len(set(indices)):
        failures.append("PARENT_POSE_INDEX_DUPLICATE")
    if not passing_exact_pose:
        failures.append("EXACT_MAPPED_PASSING_POSE_REQUIRED")
    return list(dict.fromkeys(failures))


def _direct_numerical_ids(active: list[dict[str, Any]], key: str) -> list[str]:
    identifiers: list[str] = []
    for decision in active:
        if (
            isinstance(decision, dict)
            and decision.get("action") == "accept"
            and decision.get("kind") == "numerical_criteria"
            and isinstance(decision.get("data"), dict)
            and isinstance(decision["data"].get("values"), dict)
            and key in decision["data"]["values"]
            and isinstance(decision.get("id"), str)
        ):
            identifiers.append(decision["id"])
    return identifiers


def _analysis_provenance(result: dict[str, Any]) -> dict[str, Any]:
    method = result.get("sites", {}).get("method") if isinstance(result.get("sites"), dict) else None
    if not isinstance(method, dict):
        return {"status": "pending_metadata"}
    version = method.get("SASA_tool_version")
    tool = method.get("SASA_tool")
    algorithm = method.get("SASA_algorithm")
    if not all(isinstance(value, str) and value for value in (version, tool, algorithm)):
        return {"status": "pending_metadata"}
    return {
        "status": "result_provenance_present",
        "SASA_tool": tool,
        "SASA_tool_version": version,
        "SASA_algorithm": algorithm,
    }


def apply_bundled_opinion(
    policy: dict[str, Any],
    result: dict[str, Any],
    active_direct_decisions: list[dict[str, Any]],
    *,
    archived_result_verified: bool,
    root: Path = ROOT,
) -> dict[str, Any]:
    if type(archived_result_verified) is not bool or archived_result_verified is not True:
        raise ContractError("EXPERT_OPINION_ARCHIVED_RESULT_NOT_VERIFIED")
    bundle = load_bundled_opinion(root)
    numerical = policy.get("numerical_criteria")
    if not isinstance(numerical, dict):
        raise ContractError("EXPERT_OPINION_POLICY_NUMERICAL")
    for key in ("modifiable_sites_min", "distinct_graphs_per_site_min", "calibration_seeds_min", "novel_ternary_repeats_min"):
        if type(numerical.get(key)) is not int:
            raise ContractError("EXPERT_OPINION_POLICY_NUMERICAL_TYPE")

    parent_scope = result.get("parent_scope")
    warhead = result.get("warhead")
    scope_failures: list[str] = []
    if not isinstance(parent_scope, dict) or parent_scope.get("actual_design_parent_id") != "SMARCA2-FX5":
        scope_failures.append("PARENT_NOT_SMARCA2_FX5")
    if not isinstance(parent_scope, dict) or parent_scope.get("receptor_frame") != "6HAZ chain A":
        scope_failures.append("RECEPTOR_FRAME_NOT_EXACT_6HAZ_CHAIN_A")
    if (
        not isinstance(warhead, dict)
        or warhead.get("id") != "SMARCA2-FX5"
        or warhead.get("pdb") != "6HAZ"
        or warhead.get("target_chain") != "A"
    ):
        scope_failures.append("WARHEAD_SCOPE_MISMATCH")

    source_failures, verified_sources = _verify_bound_sources(result, Path(root))
    identity_failures = _graph_and_mapping_failures(result, Path(root))
    base_failures = list(dict.fromkeys(scope_failures + source_failures + identity_failures))
    fx5_bound_scope = not base_failures

    site_failures, preserved_states = _site_failures(result)
    redocking_failures = _redocking_failures(result)
    waiver_failures = list(dict.fromkeys(base_failures + redocking_failures + site_failures))
    if numerical["distinct_graphs_per_site_min"] != 30:
        waiver_failures.append("DISTINCT_GRAPH_MINIMUM_NOT_30")

    site_direct_ids = _direct_numerical_ids(active_direct_decisions, "modifiable_sites_min")
    if site_direct_ids:
        waiver_failures.append("AUTHENTICATED_DIRECT_NUMERICAL_DECISION_HAS_PRIORITY")
    site_waiver_applied = not waiver_failures
    original_site_minimum = numerical["modifiable_sites_min"]
    if site_waiver_applied:
        numerical["modifiable_sites_min"] = 1

    floors: dict[str, Any] = {}
    for key in ("calibration_seeds_min", "novel_ternary_repeats_min"):
        original = numerical[key]
        direct_ids = _direct_numerical_ids(active_direct_decisions, key)
        applied = fx5_bound_scope and not direct_ids
        if applied:
            numerical[key] = max(original, 5)
        floors[key] = {
            "scope_applicable": fx5_bound_scope,
            "applied": applied,
            "original": original,
            "effective": numerical[key],
            "minimum": 5,
            "source_kind": "relayed_expert_default_minimum",
            "scientific_success_threshold": False,
            "direct_decision_priority": bool(direct_ids),
            "direct_decision_ids": direct_ids,
        }

    manifest = bundle["manifest"]
    policy_material = {
        "manifest_sha256": bundle["manifest_sha256"],
        "source_sha256": bundle["source_sha256"],
        "verified_input_sources": {path: FX5_SOURCE_FILES[path] for path in verified_sources},
        "fx5_scope": manifest["fx5_scope"],
        "numerical_defaults": manifest["numerical_defaults"],
    }
    policy["opinion"] = {
        "format": VERSION,
        "attribution": manifest["attribution"],
        "origin": manifest["origin"],
        "authority": manifest["authority"],
        "source": {
            "path": bundle["source_path"],
            "sha256": bundle["source_sha256"],
            "decision_locator": bundle["decision_locator"],
            "decision_extract_sha256": bundle["decision_extract_sha256"],
            "conditions_locator": bundle["conditions_locator"],
            "conditions_extract_sha256": bundle["conditions_extract_sha256"],
            "manifest_sha256": bundle["manifest_sha256"],
        },
        "policy_digest": _sha(_canonical(policy_material)),
        "verified_input_sources": {path: FX5_SOURCE_FILES[path] for path in verified_sources},
        "scope": {
            "fx5_bound_scope": fx5_bound_scope,
            "condition_failures": base_failures,
        },
        "site_waiver": {
            "applied": site_waiver_applied,
            "condition_failures": waiver_failures,
            "original_modifiable_sites_min": original_site_minimum,
            "effective_modifiable_sites_min": numerical["modifiable_sites_min"],
            "distinct_graphs_per_site_min": numerical["distinct_graphs_per_site_min"],
            "only_modifiable_atom_map": 19,
            "unknown_atom_maps_required_unchanged": [8, 11, 12],
            "site_states_observed_without_mutation": preserved_states,
            "direct_decision_ids": site_direct_ids,
        },
        "seed_floors": floors,
        "analysis_provenance": _analysis_provenance(result),
        "planning": manifest["execution_plan"],
        "unresolved": manifest["unresolved"],
        "formal_approval": False,
        "scientific_accepted": False,
    }
    return policy


__all__ = [
    "VERSION",
    "MANIFEST_RELATIVE_PATH",
    "SOURCE_RELATIVE_PATH",
    "SOURCE_SHA256",
    "FX5_SOURCE_FILES",
    "load_bundled_opinion",
    "apply_bundled_opinion",
]
