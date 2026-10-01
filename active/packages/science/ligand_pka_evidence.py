"""Auditable MolGpKa evidence generation with lazy optional dependencies."""
from __future__ import annotations

import argparse
import contextlib
import copy
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterator

ANALOG_IDS = ("W-c2afc5e73c1a", "W-4c0a639c0a41", "W-80f8f4a11b5d")
EXPECTED_MAP_ELEMENTS = {
    "W-c2afc5e73c1a": {19: "N", 5001: "O"},
    "W-4c0a639c0a41": {19: "N", 5001: "C"},
    "W-80f8f4a11b5d": {19: "N", 5001: "N"},
}
EXPECTED_PKA_MAPS = {
    "W-c2afc5e73c1a": [19],
    "W-4c0a639c0a41": [19],
    "W-80f8f4a11b5d": [19, 5001],
}
UPSTREAM_FILES = frozenset({
    "benchmark_delta_pka/README.md",
    "LICENSE.md",
    "models/weight_acid.pth",
    "models/weight_base.pth",
    "molgpka_api_example.py",
    "README.md",
    "src/__init__.py",
    "src/baseline/prepare_dataset_ap.py",
    "src/baseline/train_ap.py",
    "src/datasets/README.md",
    "src/predict_pka.py",
    "src/prepare_dataset_graph.py",
    "src/protonate.py",
    "src/train_graph.py",
    "src/utils/__init__.py",
    "src/utils/descriptor.py",
    "src/utils/gcn_conv.py",
    "src/utils/inits.py",
    "src/utils/ionization_group.py",
    "src/utils/net.py",
    "src/utils/README.md",
    "src/utils/smarts_pattern.tsv",
})
RUNTIME_UPSTREAM_FILES = frozenset({
    "LICENSE.md",
    "models/weight_acid.pth",
    "models/weight_base.pth",
    "src/__init__.py",
    "src/predict_pka.py",
    "src/utils/__init__.py",
    "src/utils/descriptor.py",
    "src/utils/gcn_conv.py",
    "src/utils/inits.py",
    "src/utils/ionization_group.py",
    "src/utils/net.py",
    "src/utils/smarts_pattern.tsv",
})
ARTIFACT_ID = re.compile(r"a-[0-9a-f]{32}\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
FORMAT = "ligand-pka-evidence/20261001.2"
UPSTREAM_URL = "https://github.com/Xundrug/MolGpKa"
CITATION = "https://pubs.acs.org/doi/10.1021/acs.jcim.1c00075"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _reject_constant(value: str) -> None:
    raise ValueError("JSON_NONFINITE:" + value)


def _pairs(pairs: list[tuple[str, Any]]) -> dict:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("JSON_DUPLICATE_KEY:" + key)
        result[key] = value
    return result


def _finite_json(value: Any) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("JSON_NONFINITE_NUMBER")
    if isinstance(value, dict):
        for item in value.values():
            _finite_json(item)
    elif isinstance(value, list):
        for item in value:
            _finite_json(item)


def parse_json_bytes(data: bytes) -> Any:
    value = json.loads(
        data.decode("utf-8"), object_pairs_hook=_pairs,
        parse_constant=_reject_constant,
    )
    _finite_json(value)
    return value


def _json_bytes(value: Any) -> bytes:
    _finite_json(value)
    return json.dumps(
        value, sort_keys=True, ensure_ascii=False, allow_nan=False,
        separators=(",", ":"),
    ).encode("utf-8")


def canonical_digest(value: Any) -> str:
    return sha256_bytes(_json_bytes(value))


def _write_json(path: Path, value: Any) -> None:
    path.write_bytes(_json_bytes(value) + b"\n")


def _strict_number(value: Any, code: str) -> float:
    if type(value) not in (int, float) or not math.isfinite(float(value)):
        raise ValueError(code)
    return float(value)


def validate_ph(value: Any) -> float:
    ph = _strict_number(value, "PH_MUST_BE_FINITE_NUMBER")
    if not 0.0 <= ph <= 14.0:
        raise ValueError("PH_OUT_OF_RANGE")
    return ph


def site_fraction(ph: float, pka: float, kind: str) -> float:
    ph = _strict_number(ph, "PH_MUST_BE_FINITE_NUMBER")
    pka = _strict_number(pka, "PKA_MUST_BE_FINITE_NUMBER")
    if kind not in {"base", "acid"}:
        raise ValueError("SITE_KIND_INVALID")
    exponent = ph - pka if kind == "base" else pka - ph
    if exponent > 308.0:
        return 0.0
    if exponent < -308.0:
        return 1.0
    if exponent >= 0.0:
        return 1.0 / (1.0 + 10.0 ** exponent)
    power = 10.0 ** (-exponent)
    return power / (1.0 + power)


def sensitivity(ph: float, pka: float, kind: str) -> dict:
    return {
        "classification": "uncoupled_henderson_hasselbalch_sensitivity",
        "not_a_confidence_interval": True,
        "not_a_coupled_microstate_population": True,
        "fraction_semantics": (
            "protonated_fraction" if kind == "base" else "deprotonated_fraction"
        ),
        "pka_minus_1": site_fraction(ph, pka - 1.0, kind),
        "predicted_pka": site_fraction(ph, pka, kind),
        "pka_plus_1": site_fraction(ph, pka + 1.0, kind),
    }


def _safe_child(root: Path, relative: str, must_exist: bool = True) -> Path:
    if not relative or Path(relative).is_absolute() or ".." in Path(relative).parts:
        raise ValueError("PATH_ESCAPE:" + relative)
    candidate = root / relative
    if candidate.is_symlink():
        raise ValueError("SYMLINK_FORBIDDEN:" + relative)
    resolved = candidate.resolve(strict=must_exist)
    if resolved != root and root not in resolved.parents:
        raise ValueError("PATH_ESCAPE:" + relative)
    return resolved


def verify_upstream(root: Path, manifest_path: Path) -> dict:
    root = root.resolve(strict=True)
    manifest_path = manifest_path.resolve(strict=True)
    document = parse_json_bytes(manifest_path.read_bytes())
    if not isinstance(document, dict):
        raise ValueError("UPSTREAM_MANIFEST_OBJECT_REQUIRED")
    if not isinstance(document.get("source"), str) or not document["source"]:
        raise ValueError("UPSTREAM_MANIFEST_SOURCE_REQUIRED")
    if not isinstance(document.get("commit"), str) or not document["commit"]:
        raise ValueError("UPSTREAM_MANIFEST_COMMIT_REQUIRED")
    if not isinstance(document.get("license"), str) or not document["license"]:
        raise ValueError("UPSTREAM_MANIFEST_LICENSE_REQUIRED")
    rows = document.get("files")
    if not isinstance(rows, list):
        raise ValueError("UPSTREAM_MANIFEST_FILES_LIST_REQUIRED")
    entries: dict[str, dict] = {}
    for row in rows:
        if not isinstance(row, dict) or set(row) != {"path", "sha256", "bytes"}:
            raise ValueError("UPSTREAM_MANIFEST_FILE_SHAPE")
        relative, digest, size = row["path"], row["sha256"], row["bytes"]
        if not isinstance(relative, str) or relative in entries:
            raise ValueError("UPSTREAM_MANIFEST_PATH_INVALID")
        if not isinstance(digest, str) or SHA256.fullmatch(digest) is None:
            raise ValueError("UPSTREAM_MANIFEST_SHA256_INVALID")
        if type(size) is not int or size < 0:
            raise ValueError("UPSTREAM_MANIFEST_BYTES_INVALID")
        entries[relative] = row
    if set(entries) != UPSTREAM_FILES:
        raise ValueError("UPSTREAM_MANIFEST_EXACT_22_FILES_REQUIRED")
    if not RUNTIME_UPSTREAM_FILES.issubset(entries):
        raise ValueError("UPSTREAM_RUNTIME_FILES_MISSING")
    for relative, row in entries.items():
        candidate = _safe_child(root, relative)
        if not candidate.is_file():
            raise ValueError("UPSTREAM_FILE_REQUIRED:" + relative)
        if candidate.stat().st_size != row["bytes"]:
            raise ValueError("UPSTREAM_SIZE_MISMATCH:" + relative)
        if sha256_file(candidate) != row["sha256"]:
            raise ValueError("UPSTREAM_HASH_MISMATCH:" + relative)
    return document


def _manifest_file_map(document: dict) -> dict[str, dict]:
    rows = document.get("files")
    if not isinstance(rows, list):
        rows = document.get("artifacts")
    if not isinstance(rows, list):
        raise ValueError("CORE_MANIFEST_FILES_LIST_REQUIRED")
    result = {}
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("CORE_MANIFEST_ENTRY_INVALID")
        name = row.get("path")
        digest = row.get("sha256")
        size = row.get("bytes")
        if not isinstance(name, str) or name in result:
            raise ValueError("CORE_MANIFEST_PATH_INVALID")
        if not isinstance(digest, str) or SHA256.fullmatch(digest) is None:
            raise ValueError("CORE_MANIFEST_SHA_INVALID")
        if type(size) is not int or size < 0:
            raise ValueError("CORE_MANIFEST_BYTES_INVALID")
        result[name] = row
    return result


def _verify_core(core: Path) -> tuple[dict, bytes, dict, Path]:
    core = core.resolve(strict=True)
    manifest_path = _safe_child(core, "manifest.json")
    manifest = parse_json_bytes(manifest_path.read_bytes())
    if not isinstance(manifest, dict):
        raise ValueError("CORE_MANIFEST_OBJECT_REQUIRED")
    files = _manifest_file_map(manifest)
    for name in ("summary.json", "all-states.sdf"):
        if name not in files:
            raise ValueError("CORE_REQUIRED_FILE_MISSING:" + name)
    for name, row in files.items():
        path = _safe_child(core, name)
        if not path.is_file() or path.stat().st_size != row["bytes"]:
            raise ValueError("CORE_ARTIFACT_SIZE_MISMATCH:" + name)
        if sha256_file(path) != row["sha256"]:
            raise ValueError("CORE_ARTIFACT_HASH_MISMATCH:" + name)
    summary = parse_json_bytes(_safe_child(core, "summary.json").read_bytes())
    if not isinstance(summary, dict):
        raise ValueError("CORE_SUMMARY_OBJECT_REQUIRED")
    return summary, _safe_child(core, "all-states.sdf").read_bytes(), manifest, manifest_path


def split_sdf_bytes(data: bytes) -> list[bytes]:
    blocks: list[bytes] = []
    start = 0
    while True:
        index = data.find(b"$$$$", start)
        if index < 0:
            if data[start:].strip():
                raise ValueError("SDF_TRAILING_DATA")
            return blocks
        end = index + 4
        if data[end:end + 2] == b"\r\n":
            end += 2
        elif data[end:end + 1] in (b"\r", b"\n"):
            end += 1
        blocks.append(data[start:end])
        start = end


def validate_core_records(summary: dict, sdf: bytes) -> tuple[list[dict], list[bytes]]:
    records = summary.get("rows")
    if not isinstance(records, list) or not all(isinstance(x, dict) for x in records):
        raise ValueError("CORE_ROWS_LIST_REQUIRED")
    blocks = split_sdf_bytes(sdf)
    if len(records) != 40 or len(blocks) != 40:
        raise ValueError("CORE_REQUIRES_EXACTLY_40_PAIRED_RECORDS")
    for index, (record, block) in enumerate(zip(records, blocks)):
        declared = record.get("sdf_sha256")
        if declared is not None and declared != sha256_bytes(block):
            raise ValueError("CORE_SDF_RECORD_HASH_MISMATCH:" + str(index))
        if type(record.get("required_contact_evidence")) is not bool:
            raise ValueError("CORE_REQUIRED_CONTACT_EVIDENCE_BOOL:" + str(index))
    true_count = sum(record["required_contact_evidence"] for record in records)
    if true_count != 16 or len(records) - true_count != 24:
        raise ValueError("CORE_CONTACT_COUNTS_MUST_BE_16_TRUE_24_FALSE")
    selected = [
        i for i, record in enumerate(records)
        if record.get("analog_id") in ANALOG_IDS
        and record.get("pose_index") == 0
        and record.get("state_id") == "source"
    ]
    if len(selected) != 3 or {records[i]["analog_id"] for i in selected} != set(ANALOG_IDS):
        raise ValueError("CORE_SOURCE_POSE0_SELECTION_INVALID")
    return records, blocks


def contact_summary(records: list[dict]) -> tuple[list[dict], dict]:
    rows = []
    groups: dict[tuple[str, int], list[tuple[bool, Any]]] = {}
    for index, record in enumerate(records):
        analog_id = record.get("analog_id")
        pose = record.get("pose_index")
        contact = record["required_contact_evidence"]
        if analog_id in ANALOG_IDS and type(pose) is int and pose in range(5):
            groups.setdefault((analog_id, pose), []).append((contact, record.get("failure")))
        rows.append({
            "core_index": index,
            "analog_id": analog_id,
            "pose_index": pose,
            "state_index": record.get("state_index"),
            "state_id": record.get("state_id"),
            "status": record.get("status"),
            "failure": record.get("failure"),
            "required_contact_evidence": contact,
            "core_rmsd_auxiliary_only": record.get("auxiliary_core_RMSD_A"),
        })
    supported = {}
    for analog_id in ANALOG_IDS:
        supported[analog_id] = {}
        for pose in range(5):
            states = groups.get((analog_id, pose), [])
            expected_states = 4 if analog_id == "W-80f8f4a11b5d" else 2
            supported[analog_id][str(pose)] = (
                len(states) == expected_states
                and all(contact and failure is None for contact, failure in states)
            )
    return rows, supported


def _artifact_reference(value: Any) -> bool:
    return isinstance(value, dict) and isinstance(value.get("artifact_id"), str) and isinstance(value.get("sha256"), str)


def _artifact_blob(connection: sqlite3.Connection, blob_root: Path, project: str, ref: dict) -> bytes:
    artifact_id = ref.get("artifact_id")
    if not isinstance(artifact_id, str) or ARTIFACT_ID.fullmatch(artifact_id) is None:
        raise ValueError("ARTIFACT_ID_INVALID")
    if not isinstance(ref.get("sha256"), str) or SHA256.fullmatch(ref["sha256"]) is None:
        raise ValueError("ARTIFACT_SHA256_INVALID")
    row = connection.execute(
        "SELECT metadata FROM artifacts WHERE id=? AND project=?", (artifact_id, project)
    ).fetchone()
    if row is None:
        raise ValueError("ARTIFACT_NOT_PROJECT_SCOPED")
    metadata = parse_json_bytes(row[0].encode("utf-8") if isinstance(row[0], str) else row[0])
    if metadata != ref:
        raise ValueError("ARTIFACT_METADATA_REFERENCE_MISMATCH")
    blob_root = blob_root.resolve(strict=True)
    blob_path = blob_root / artifact_id
    if blob_path.is_symlink() or blob_path.name != artifact_id:
        raise ValueError("ARTIFACT_BLOB_SYMLINK_OR_SUFFIX")
    resolved = blob_path.resolve(strict=True)
    if resolved.parent != blob_root or not resolved.is_file():
        raise ValueError("ARTIFACT_BLOB_PATH_ESCAPE")
    data = resolved.read_bytes()
    if sha256_bytes(data) != ref["sha256"]:
        raise ValueError("ARTIFACT_BLOB_HASH_MISMATCH")
    return data


def _count_equal(values: Any, target: dict) -> int:
    if not isinstance(values, list):
        return 0
    return sum(item.get("ref") == target for item in values if isinstance(item, dict))


def read_bound_design(store: Path, project: str, job_id: str) -> tuple[dict, dict]:
    store = store.resolve(strict=True)
    db_path = store if store.is_file() else _safe_child(store, "index.sqlite3")
    blob_root = db_path.parent / "blobs"
    connection = sqlite3.connect("file:" + str(db_path) + "?mode=ro", uri=True)
    try:
        row = connection.execute(
            "SELECT body,state FROM workbench_jobs WHERE id=? AND project=?",
            (job_id, project),
        ).fetchone()
        if row is None:
            raise ValueError("WORKBENCH_JOB_NOT_FOUND")
        body = parse_json_bytes(row[0].encode("utf-8") if isinstance(row[0], str) else row[0])
        if not isinstance(body, dict):
            raise ValueError("WORKBENCH_JOB_BODY_OBJECT_REQUIRED")
        if row[1] != "completed" or body.get("state") != "completed":
            raise ValueError("WORKBENCH_JOB_NOT_COMPLETED")
        if body.get("id") != job_id or body.get("project_id") != project:
            raise ValueError("WORKBENCH_JOB_BODY_SCOPE_MISMATCH")
        result = body.get("result")
        if not isinstance(result, dict):
            raise ValueError("WORKBENCH_JOB_RESULT_REQUIRED")
        files = result.get("files")
        if not isinstance(files, dict):
            raise ValueError("WORKBENCH_JOB_RESULT_FILES_REQUIRED")
        json_ref, report_ref = files.get("json"), files.get("report")
        if not _artifact_reference(json_ref) or not _artifact_reference(report_ref):
            raise ValueError("WORKBENCH_JOB_JSON_REPORT_REFS_REQUIRED")
        outputs = body.get("outputs")
        if not isinstance(outputs, list) or not all(
            isinstance(item, dict)
            and set(item) == {"name", "ref"}
            and isinstance(item["name"], str)
            and item["name"]
            and _artifact_reference(item["ref"])
            for item in outputs
        ):
            raise ValueError("WORKBENCH_JOB_OUTPUT_WRAPPER_SHAPE_REQUIRED")
        if _count_equal(outputs, json_ref) != 1 or _count_equal(outputs, report_ref) != 1:
            raise ValueError("WORKBENCH_JOB_OUTPUT_REFS_EXACTLY_ONCE")
        json_blob = _artifact_blob(connection, blob_root, project, json_ref)
        _artifact_blob(connection, blob_root, project, report_ref)
        archived = copy.deepcopy(result)
        archived_files = archived.get("files")
        archived_files.pop("json", None)
        archived_files.pop("report", None)
        parsed_archive = parse_json_bytes(json_blob)
        if parsed_archive != archived:
            raise ValueError("ARCHIVED_DESIGN_JSON_NOT_EXACT_RESULT")
        body_digest = body.get("input_digest")
        binding_digest = body.get("binding", {}).get("digest") if isinstance(body.get("binding"), dict) else None
        result_digest = result.get("input_binding", {}).get("digest") if isinstance(result.get("input_binding"), dict) else None
        if not isinstance(body_digest, str) or SHA256.fullmatch(body_digest) is None:
            raise ValueError("INPUT_DIGEST_INVALID")
        if body_digest != binding_digest or body_digest != result_digest:
            raise ValueError("INPUT_BINDING_DIGEST_MISMATCH")
        parent_scope = result.get("parent_scope")
        parent_id = parent_scope.get("actual_design_parent_id") if isinstance(parent_scope, dict) else None
        if not isinstance(parent_id, str) or not parent_id:
            raise ValueError("ACTUAL_DESIGN_PARENT_ID_REQUIRED")
        binding = {
            "project": project,
            "job_id": job_id,
            "input_sha256": body_digest,
            "result_sha256": json_ref["sha256"],
            "result_binding_kind": "exact_archived_design_json_bytes",
            "parent_id": parent_id,
        }
        return parsed_archive, binding
    finally:
        connection.close()


def select_analogs(design: dict) -> dict[str, dict]:
    analogs = design.get("analogs")
    if not isinstance(analogs, list):
        raise ValueError("DESIGN_ANALOGS_LIST_REQUIRED")
    selected = {}
    for analog_id in ANALOG_IDS:
        matches = [x for x in analogs if isinstance(x, dict) and x.get("id") == analog_id]
        if len(matches) != 1:
            raise ValueError("SUPPORTED_ANALOG_NOT_UNIQUE:" + analog_id)
        analog = matches[0]
        if not isinstance(analog.get("mapped_smiles"), str) or not analog["mapped_smiles"]:
            raise ValueError("ANALOG_MAPPED_SMILES_REQUIRED:" + analog_id)
        selected[analog_id] = analog
    return selected


def mapped_graph(mol: Any) -> str:
    from rdkit import Chem
    seen = set()
    value = Chem.RemoveHs(Chem.Mol(mol))
    for atom in value.GetAtoms():
        atom_map = atom.GetAtomMapNum()
        if atom_map <= 0 or atom_map in seen:
            raise ValueError("POSITIVE_UNIQUE_HEAVY_ATOM_MAPS_REQUIRED")
        seen.add(atom_map)
    return Chem.MolToSmiles(value, canonical=True, isomericSmiles=True)


def _unmapped_canonical(mol: Any) -> str:
    from rdkit import Chem
    value = Chem.Mol(mol)
    for atom in value.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(value, canonical=True, isomericSmiles=True)


def chemistry_identity(sdf_block: bytes, mapped_smiles: str, analog_id: str) -> dict:
    from rdkit import Chem
    source = Chem.MolFromMolBlock(sdf_block.decode("utf-8"), sanitize=True, removeHs=False)
    result = Chem.MolFromSmiles(mapped_smiles)
    if source is None or result is None:
        raise ValueError("MAPPED_MOLECULE_PARSE_FAILED")
    source_graph = mapped_graph(Chem.RemoveHs(source))
    result_graph = mapped_graph(Chem.RemoveHs(result))
    if source_graph != result_graph:
        raise ValueError("CORE_AND_RESULT_EXACT_MAPPED_GRAPH_MISMATCH")
    elements = {atom.GetAtomMapNum(): atom.GetSymbol() for atom in result.GetAtoms() if atom.GetAtomicNum() != 1}
    for atom_map, element in EXPECTED_MAP_ELEMENTS[analog_id].items():
        if elements.get(atom_map) != element:
            raise ValueError(f"EXPECTED_MAP_ELEMENT_MISMATCH:{analog_id}:{atom_map}")
    canonical = Chem.MolToSmiles(result, canonical=True, isomericSmiles=True)
    roundtrip = Chem.MolFromSmiles(canonical)
    if roundtrip is None or mapped_graph(roundtrip) != result_graph:
        raise ValueError("CANONICAL_MAPPED_ISOMERIC_GRAPH_ROUNDTRIP_FAILED")
    return {
        "mapped_smiles": mapped_smiles,
        "canonical_mapped_isomeric_smiles": canonical,
        "canonical_unmapped_isomeric_graph": _unmapped_canonical(result),
        "mapped_graph_sha256": canonical_digest(result_graph),
    }


@contextlib.contextmanager
def _isolated_upstream(root: Path, dependency: Path | None) -> Iterator[Any]:
    src = _safe_child(root.resolve(strict=True), "src")
    old_cwd, old_path = Path.cwd(), list(sys.path)
    prefixes = ("utils", "_trusted_molgpka_predict")
    old_modules = {k: v for k, v in sys.modules.items() if k == "utils" or k.startswith("utils.") or k == prefixes[1]}
    for key in list(old_modules):
        sys.modules.pop(key, None)
    try:
        if dependency is not None:
            sys.path.insert(0, str(dependency.resolve(strict=True)))
        sys.path.insert(0, str(src))
        os.chdir(src)
        spec = importlib.util.spec_from_file_location("_trusted_molgpka_predict", src / "predict_pka.py")
        if spec is None or spec.loader is None:
            raise RuntimeError("UPSTREAM_IMPORT_SPEC_FAILED")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        yield module
    finally:
        os.chdir(old_cwd)
        sys.path[:] = old_path
        for key in list(sys.modules):
            if key == "utils" or key.startswith("utils.") or key == prefixes[1]:
                sys.modules.pop(key, None)
        sys.modules.update(old_modules)


def _bind_model_site(molecule, site_kind: str, model_target_atom_index: int) -> dict:
    """Bind an official MolGpKa target to its mapped heavy-atom site."""
    if site_kind not in {"acid", "base"}:
        raise ValueError("MODEL_SITE_KIND_INVALID")
    target_index = int(model_target_atom_index)
    target = molecule.GetAtomWithIdx(target_index)

    if site_kind == "base":
        if target.GetAtomicNum() == 1 or target.GetAtomMapNum() <= 0:
            raise ValueError("BASE_MODEL_SITE_NOT_MAPPED_HEAVY_ATOM")
        parent = target
    else:
        if target.GetAtomicNum() != 1:
            raise ValueError("ACID_MODEL_SITE_NOT_HYDROGEN")
        heavy_neighbors = [neighbor for neighbor in target.GetNeighbors() if neighbor.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            raise ValueError("ACID_MODEL_HYDROGEN_REQUIRES_ONE_HEAVY_NEIGHBOR")
        parent = heavy_neighbors[0]
        parent_map = parent.GetAtomMapNum()
        if parent_map <= 0:
            raise ValueError("ACID_MODEL_HYDROGEN_PARENT_NOT_POSITIVELY_MAPPED")
        mapped_parents = [
            atom for atom in molecule.GetAtoms()
            if atom.GetAtomicNum() != 1 and atom.GetAtomMapNum() == parent_map
        ]
        if len(mapped_parents) != 1:
            raise ValueError("ACID_MODEL_HYDROGEN_PARENT_MAP_NOT_UNIQUE")

    atom_map = parent.GetAtomMapNum()
    site = {
        "computed_site_id": f"{site_kind}/{target_index}/{atom_map}",
        "site_kind": site_kind,
        "atom_map": atom_map,
        "element": parent.GetSymbol(),
        "prepared_atom_index": target_index,
        "raw_prepared_atom_index": target_index,
        "model_target_atom_index": target_index,
        "model_target_element": target.GetSymbol(),
        "input_formal_charge": parent.GetFormalCharge(),
        "parent_input_formal_charge": parent.GetFormalCharge(),
        "model_target_formal_charge": target.GetFormalCharge(),
    }
    if site_kind == "acid":
        site["acid_proton_attached_to_atom_map"] = atom_map
    return site


def worker(request_path: Path, result_path: Path, root: Path, manifest: Path, dependency: Path | None) -> int:
    import sys

    if dependency is not None:
        sys.path.insert(0, str(dependency.resolve(strict=True)))

    import numpy
    import pandas
    import torch
    import torch_geometric
    import torch_scatter
    from rdkit import Chem
    version_match = re.match(r"^(\d+)\.(\d+)", str(torch.__version__))
    if version_match is None or tuple(map(int, version_match.groups())) < (2, 6):
        raise RuntimeError("TORCH_2_6_OR_NEWER_REQUIRED_FOR_SAFE_OFFICIAL_MODEL_LOAD")
    request = parse_json_bytes(request_path.read_bytes())
    upstream = verify_upstream(root, manifest)
    hashes = {row["path"]: row["sha256"] for row in upstream["files"]}
    torch.set_num_threads(1)
    results = {}
    with _isolated_upstream(root, dependency) as module:
        acid_model = module.load_model(root / "models/weight_acid.pth", "cpu")
        base_model = module.load_model(root / "models/weight_base.pth", "cpu")
        for analog_id in ANALOG_IDS:
            mapped_smiles = request.get("analogs", {}).get(analog_id)
            try:
                molecule = Chem.MolFromSmiles(mapped_smiles)
                if molecule is None:
                    raise ValueError("MAPPED_SMILES_PARSE_FAILED")
                before = mapped_graph(molecule)
                prepared = Chem.AddHs(molecule)
                if mapped_graph(prepared) != before:
                    raise ValueError("ADD_HS_CHANGED_HEAVY_MAPPED_GRAPH")
                sites = []
                for kind, model in (("base", base_model), ("acid", acid_model)):
                    atom_ids = module.get_ionization_aid(prepared, acid_or_base=kind)
                    for atom_index in atom_ids:
                        site = _bind_model_site(prepared, kind, int(atom_index))
                        # Official acid targets are hydrogen indices. Preserve that
                        # index as the model input; the mapped heavy atom is metadata.
                        pka = float(module.model_pred(prepared, atom_index, model, "cpu"))
                        if not math.isfinite(pka):
                            raise ValueError("MODEL_PKA_NOT_FINITE")
                        site["pKa"] = pka
                        sites.append(site)
                if not sites:
                    raise ValueError("NO_ACTUAL_FINITE_PKA_SITES")
                result = {"status": "success", "sites": sites}
                if analog_id == "W-80f8f4a11b5d":
                    expected_maps = [19, 5001]
                    predicted_maps = sorted({
                        site["atom_map"] for site in sites
                        if site["atom_map"] in expected_maps
                    })
                    missing_maps = [atom_map for atom_map in expected_maps if atom_map not in predicted_maps]
                    result["model_coverage"] = {
                        "expected_atom_maps": expected_maps,
                        "predicted_expected_atom_maps": predicted_maps,
                        "missing_expected_atom_maps": missing_maps,
                        "complete": not missing_maps,
                    }
                    if 5001 in missing_maps:
                        result["limitations"] = [
                            "Official MolGpKa get_ionization_aid did not detect the expected terminal amine at atom map 5001; no model score was produced for that site."
                        ]
                results[analog_id] = result
            except Exception as error:
                results[analog_id] = {
                    "status": "failed", "sites": [],
                    "failure": f"{type(error).__name__}:{error}",
                }
    _write_json(result_path, {
        "results": results,
        "source_hashes": hashes,
        "source_commit": upstream["commit"],
        "worker_code_sha256": sha256_file(Path(__file__).resolve()),
        "upstream_module_sha256": hashes["src/predict_pka.py"],
        "versions": {
            "python": sys.version.split()[0],
            "torch": str(torch.__version__),
            "torch_geometric": str(torch_geometric.__version__),
            "torch_scatter": str(torch_scatter.__version__),
            "numpy": str(numpy.__version__),
            "pandas": str(pandas.__version__),
            "rdkit": getattr(sys.modules.get("rdkit"), "__version__", "unknown"),
        },
    })
    return 0


def _run_worker(request: dict, root: Path, manifest: Path, dependency: Path | None) -> tuple[dict | None, dict]:
    with tempfile.TemporaryDirectory(prefix="molgpka-worker-") as temporary:
        request_path = Path(temporary) / "request.json"
        result_path = Path(temporary) / "result.json"
        _write_json(request_path, request)
        repository_root = Path(__file__).resolve().parents[2]
        command = [
            sys.executable, "-X", "utf8", "-B", "-m",
            "packages.science.ligand_pka_evidence",
            "--worker", str(request_path), str(result_path),
            "--upstream-root", str(root), "--upstream-manifest", str(manifest),
        ]
        if dependency is not None:
            command.extend(("--dependency-directory", str(dependency)))
        environment = dict(os.environ)
        environment["CUDA_VISIBLE_DEVICES"] = ""
        environment["PYTHONUTF8"] = "1"
        try:
            completed = subprocess.run(
                command, cwd=repository_root, env=environment, capture_output=True,
                text=True, encoding="utf-8", errors="replace",
                timeout=600, check=False,
            )
            receipt = {
                "command": command,
                "cwd": str(repository_root),
                "returncode": completed.returncode,
                "stdout": completed.stdout,
                "stderr": completed.stderr,
                "cuda_visible_devices": "",
            }
        except subprocess.TimeoutExpired as error:
            receipt = {
                "command": command,
                "cwd": str(repository_root),
                "returncode": None,
                "timeout_seconds": 600,
                "stdout": error.stdout.decode(errors="replace") if isinstance(error.stdout, bytes) else (error.stdout or ""),
                "stderr": error.stderr.decode(errors="replace") if isinstance(error.stderr, bytes) else (error.stderr or ""),
                "cuda_visible_devices": "",
            }
            return None, receipt
        if completed.returncode != 0 or not result_path.is_file():
            return None, receipt
        return parse_json_bytes(result_path.read_bytes()), receipt


def _interaction_reports(store: Path, project: str, job_id: str, result_sha256: str) -> dict[str, tuple[dict, dict]]:
    store = store.resolve(strict=True)
    db_path = store if store.is_file() else _safe_child(store, "index.sqlite3")
    connection = sqlite3.connect("file:" + str(db_path) + "?mode=ro", uri=True)
    reports = {}
    try:
        rows = connection.execute(
            "SELECT binding,ref FROM scientific_evidence WHERE project=? AND job_id=? AND kind='interaction' AND superseded=0 ORDER BY rowid DESC",
            (project, job_id),
        ).fetchall()
        for binding_raw, ref_raw in rows:
            binding = parse_json_bytes(binding_raw.encode() if isinstance(binding_raw, str) else binding_raw)
            ref = parse_json_bytes(ref_raw.encode() if isinstance(ref_raw, str) else ref_raw)
            if not isinstance(binding, dict) or binding.get("result_sha256") != result_sha256:
                continue
            blob = _artifact_blob(connection, db_path.parent / "blobs", project, ref)
            report = parse_json_bytes(blob)
            if not isinstance(report, dict) or sha256_bytes(_json_bytes(report) + b"\n") != binding.get("evidence_sha256"):
                raise ValueError("INTERACTION_EVIDENCE_DIGEST_MISMATCH")
            analog_id = report.get("analog_id")
            report_id = report.get("report_id")
            if analog_id in ANALOG_IDS and isinstance(report_id, str) and report_id.startswith("interaction-W-") and analog_id not in reports:
                reports[analog_id] = (ref, report)
        return reports
    finally:
        connection.close()


def quantitative_source(project: str, job_id: str, binding: dict, predictions: dict, identities: dict, reports: dict, ph: float) -> dict:
    records = []
    for analog_id in ANALOG_IDS:
        prediction = predictions[analog_id]
        if prediction["status"] != "success":
            continue
        if analog_id not in reports:
            raise ValueError("CURRENT_INTERACTION_REPORT_MISSING:" + analog_id)
        interaction_ref, report = reports[analog_id]
        report_ph = _strict_number(report.get("pH_context"), "INTERACTION_REPORT_PH_CONTEXT_INVALID")
        if report_ph != ph:
            raise ValueError("INTERACTION_REPORT_PH_CONTEXT_MISMATCH")
        alternatives = report.get("state_alternatives")
        if not isinstance(alternatives, list) or not alternatives or not isinstance(alternatives[0], dict):
            raise ValueError("INTERACTION_REPORT_STATE_ZERO_REQUIRED")
        state = alternatives[0]
        if state.get("status") != "computed" or not isinstance(state.get("mapped_smiles"), str):
            raise ValueError("INTERACTION_REPORT_STATE_ZERO_NOT_COMPUTED")
        state_ph = _strict_number(state.get("pH_context"), "INTERACTION_STATE_PH_CONTEXT_INVALID")
        if state_ph != ph:
            raise ValueError("INTERACTION_STATE_PH_CONTEXT_MISMATCH")
        state_identity = chemistry_identity(
            prediction["source_sdf_block"], state["mapped_smiles"], analog_id
        )
        expected_graph = identities[analog_id]["canonical_unmapped_isomeric_graph"]
        if state_identity["canonical_unmapped_isomeric_graph"] != expected_graph:
            raise ValueError("INTERACTION_STATE_NOT_ACTUAL_PREDICTION_INPUT")
        for site in prediction["sites"]:
            records.append({
                "project_id": project,
                "job_id": job_id,
                "result_sha256": binding["result_sha256"],
                "analog_id": analog_id,
                "selected_state_index": 0,
                "analog_canonical_isomeric_graph": expected_graph,
                "state_canonical_isomeric_graph": expected_graph,
                "pH": ph,
                "interaction_ref": interaction_ref,
                "report_id": report["report_id"],
                "method": "MolGpKa graph-neural-network site pKa prediction",
                "protocol": "Chem.AddHs followed by official get_ionization_aid and model_pred calls",
                "conditions": f"Isolated mapped analog; CPU inference; pH context {ph}",
                "uncertainty_description": "Model prediction without calibrated confidence interval; ±1 pKa HH values are sensitivity only",
                "interpretation": "Prediction input state, not a final state selection or population estimate",
                "source_identifier": "MolGpKa:" + binding["result_sha256"],
                "source_url": UPSTREAM_URL,
                "source_kind": "computed",
                "pKa": site["pKa"],
                "atom_map": site["atom_map"],
                "site_kind": site["site_kind"],
                "computed_site_id": site["computed_site_id"],
                "prepared_atom_index": site["prepared_atom_index"],
                "raw_prepared_atom_index": site["raw_prepared_atom_index"],
                "model_target_atom_index": site["model_target_atom_index"],
                "model_target_element": site["model_target_element"],
                "parent_element": site["element"],
                "parent_input_formal_charge": site["parent_input_formal_charge"],
                "model_target_formal_charge": site["model_target_formal_charge"],
                "acid_proton_attached_to_atom_map": site.get("acid_proton_attached_to_atom_map"),
                "model_commit": prediction["model_commit"],
                "weights_sha256": prediction["weights_sha256"],
                "input_graph": identities[analog_id]["canonical_mapped_isomeric_smiles"],
                "source_input_graph_sha256": identities[analog_id]["mapped_graph_sha256"],
                "source_code_sha256": prediction["source_code_sha256"],
            })
    if not records:
        raise ValueError("QUANTITATIVE_SOURCE_REQUIRES_ACTUAL_PKA_RECORD")
    return {
        "type": "computed_quantitative_source",
        "scientific_approved": False,
        "state_selection_performed": False,
        "binding_role": "prediction_input_state_not_final_selection",
        "records": records,
    }


def run(args: argparse.Namespace) -> int:
    ph = validate_ph(args.ph)
    output = Path(args.output)
    if output.exists():
        raise ValueError("OUTPUT_DIRECTORY_MUST_NOT_EXIST")
    upstream_root = Path(args.upstream_root).resolve(strict=True)
    upstream_manifest_path = Path(args.upstream_manifest).resolve(strict=True)
    dependency = Path(args.dependency_directory).resolve(strict=True) if args.dependency_directory else None
    upstream = verify_upstream(upstream_root, upstream_manifest_path)
    summary, sdf, core_manifest, core_manifest_path = _verify_core(Path(args.core_directory))
    records, blocks = validate_core_records(summary, sdf)
    design, binding = read_bound_design(Path(args.store), args.project, args.job_id)
    analogs = select_analogs(design)
    output.mkdir()
    (output / "inputs").mkdir()
    shutil.copyfile(upstream_manifest_path, output / "upstream-manifest.json")
    shutil.copyfile(upstream_root / "LICENSE.md", output / "LICENSE.md")
    shutil.copyfile(core_manifest_path, output / "core-manifest.json")

    identities, source_blocks = {}, {}
    for analog_id in ANALOG_IDS:
        indexes = [i for i, row in enumerate(records) if row.get("analog_id") == analog_id and row.get("pose_index") == 0 and row.get("state_id") == "source"]
        source_blocks[analog_id] = blocks[indexes[0]]
        identities[analog_id] = chemistry_identity(blocks[indexes[0]], analogs[analog_id]["mapped_smiles"], analog_id)
    combined = b"".join(source_blocks[x] for x in ANALOG_IDS)
    (output / "inputs" / "source-pose0.sdf").write_bytes(combined)

    request = {"analogs": {x: identities[x]["mapped_smiles"] for x in ANALOG_IDS}}
    worker_result, receipt = _run_worker(request, upstream_root, upstream_manifest_path, dependency)
    _write_json(output / "worker-receipt.json", receipt)
    if worker_result is None:
        failure = {
            "format": FORMAT,
            "status": "worker_failure",
            "scientific_approved": False,
            "population_prediction_performed": False,
            "final_state_selection_performed": False,
            "receipt": receipt,
            "source": upstream["source"],
            "commit": upstream["commit"],
            "source_files": upstream["files"],
            "worker_code_sha256": sha256_file(Path(__file__).resolve()),
            "source_input_graph_sha256": {x: identities[x]["mapped_graph_sha256"] for x in ANALOG_IDS},
        }
        _write_json(output / "raw-predictions.json", failure)
        _write_json(output / "model-provenance.json", failure)
        raise RuntimeError("MOLGPKA_WORKER_FAILED")

    weight_hashes = {
        name: next(row["sha256"] for row in upstream["files"] if row["path"] == name)
        for name in ("models/weight_acid.pth", "models/weight_base.pth")
    }
    predictions = {}
    failures = 0
    for analog_id in ANALOG_IDS:
        prediction = copy.deepcopy(worker_result.get("results", {}).get(analog_id, {"status": "failed", "sites": [], "failure": "WORKER_RESULT_MISSING"}))
        valid_sites = prediction.get("sites") if isinstance(prediction.get("sites"), list) else []
        for site in valid_sites:
            _strict_number(site.get("pKa"), "MODEL_PKA_NOT_FINITE")
            site["pH"] = ph
            site["uncoupled_titration_sensitivity"] = sensitivity(ph, site["pKa"], site["site_kind"])
        if prediction.get("status") != "success" or not valid_sites:
            prediction["status"] = "failed"
            failures += 1
        detected = {site.get("atom_map") for site in valid_sites}
        expected = EXPECTED_PKA_MAPS[analog_id]
        prediction.update({
            "actual_finite_site_prediction_performed": prediction["status"] == "success" and bool(valid_sites),
            "population_prediction_performed": False,
            "final_state_selection_performed": False,
            "scientific_approved": False,
            "expected_maps": expected,
            "missing_required_maps": [x for x in expected if x not in detected],
            "model_commit": upstream["commit"],
            "weights_sha256": weight_hashes,
            "source_code_sha256": worker_result.get("upstream_module_sha256"),
            "worker_code_sha256": worker_result.get("worker_code_sha256"),
            "source_input_graph_sha256": identities[analog_id]["mapped_graph_sha256"],
            "source_sdf_block": source_blocks[analog_id],
        })
        predictions[analog_id] = prediction

    core_rows, supported_poses = contact_summary(records)
    raw_predictions = copy.deepcopy(predictions)
    for prediction in raw_predictions.values():
        prediction.pop("source_sdf_block", None)
    raw = {
        "format": FORMAT,
        "job_binding": binding,
        "pH": ph,
        "predictions": raw_predictions,
        "core_records": core_rows,
        "supported_poses": supported_poses,
        "required_contact_definition": {
            "polar": "map 17 to ASN1464 OD1",
            "hydrophobic": ["VAL1408", "PHE1409", "ILE1470"],
        },
        "scientific_approved": False,
        "population_prediction_performed": False,
        "final_state_selection_performed": False,
        "limitations": [
            "Predictions are computed for isolated warhead analogues, not assembled PROTACs or bound-protein pKa.",
            "Effects of conjugation, binding, and multiple ionization have not been validated.",
            "HH ±1 pKa values are sensitivity calculations, not confidence intervals.",
            "No coupled microstate population was calculated.",
            "Core RMSD is auxiliary and does not establish required-contact support.",
            "State index zero is prediction input state, not a human or final selection.",
        ],
    }
    _write_json(output / "raw-predictions.json", raw)

    reports = _interaction_reports(Path(args.store), args.project, args.job_id, binding["result_sha256"])
    quantitative = quantitative_source(args.project, args.job_id, binding, predictions, identities, reports, ph)
    _write_json(output / "quantitative-source.json", quantitative)

    compact = {
        "format": FORMAT,
        "status": "partial_failure" if failures else "completed_with_limits",
        "pH": ph,
        "analog_results": {
            analog_id: {
                "status": predictions[analog_id]["status"],
                "sites": [{
                    "atom_map": site["atom_map"],
                    "site_kind": site["site_kind"],
                    "pKa": site["pKa"],
                    "uncoupled_titration_sensitivity": site["uncoupled_titration_sensitivity"],
                    "source_input_graph_sha256": identities[analog_id]["mapped_graph_sha256"],
                } for site in predictions[analog_id]["sites"]],
            } for analog_id in ANALOG_IDS
        },
        "missing_required_maps": {x: predictions[x]["missing_required_maps"] for x in ANALOG_IDS},
        "required_contact_true_count": sum(row["required_contact_evidence"] for row in core_rows),
        "required_contact_false_count": sum(not row["required_contact_evidence"] for row in core_rows),
        "supported_poses": supported_poses,
        "input_sdf_sha256": sha256_bytes(combined),
        "result_sha256": binding["result_sha256"],
        "upstream_manifest_sha256": sha256_file(upstream_manifest_path),
        "core_manifest_sha256": sha256_file(core_manifest_path),
        "model": {"name": "MolGpKa", "commit": upstream["commit"], "weights_sha256": weight_hashes},
        "scientific_approved": False,
        "population_prediction_performed": False,
        "final_state_selection_performed": False,
        "limitations": raw["limitations"],
    }
    _write_json(output / "compact-summary.json", compact)

    provenance = {
        "format": FORMAT,
        "source": upstream["source"],
        "commit": upstream["commit"],
        "license": upstream["license"],
        "source_files": upstream["files"],
        "weights_sha256": weight_hashes,
        "runtime_versions": worker_result.get("versions", {}),
        "core_manifest": core_manifest,
        "input_sdf_sha256": sha256_bytes(combined),
        "official_loader_used": True,
        "safe_torch_default_preserved": True,
        "unsafe_pickle_fallback_used": False,
        "official_get_ionization_aid_used": True,
        "official_model_pred_used": True,
        "hidden_neutralization_performed": False,
        "scientific_approved": False,
        "population_prediction_performed": False,
        "final_state_selection_performed": False,
        "citations": [CITATION, UPSTREAM_URL],
    }
    _write_json(output / "model-provenance.json", provenance)

    readme = """# 리간드 pKa 계산 근거

이 결과는 분리된 warhead 매핑 아날로그에 대한 MolGpKa 계산 예측이며, 조립된 PROTAC 또는 단백질 결합 상태의 pKa가 아닙니다. 접합, 결합 및 다중 이온화 효과는 검증되지 않았습니다. 사람의 과학적 승인, 최종 상태 선택, 결합 미세상태 집단 예측을 수행하지 않았습니다. Henderson–Hasselbalch ±1 pKa 결과는 민감도 분석이며 신뢰구간이 아닙니다.

## 일반 휴대형 검증

```sh
python -m unittest tests.test_ligand_pka_evidence
python scripts/run_ligand_pka_evidence.py --help
```

## 선택적 pKa 실행 환경

MolGpKa 소스는 upstream-manifest.json의 정확한 22개 파일, 지정 커밋, LICENSE.md 및 해시가 일치해야 합니다. torch, torch-geometric, torch-scatter, RDKit은 별도 선택 환경에 설치해야 하며 결과 팩에 환경이나 가중치를 포함하지 않습니다.

```sh
python scripts/run_ligand_pka_evidence.py --upstream-root <MOLGPKA_ROOT> --upstream-manifest <UPSTREAM_MANIFEST_JSON> --dependency-directory <OPTIONAL_PKA_ENV_SITE_PACKAGES> --core-directory <CORE_DIRECTORY> --store <READONLY_STORE> --project <PROJECT_ID> --job-id <JOB_ID> --output <NEW_OUTPUT_DIRECTORY> --ph 7.4
```

모델 소스, 커밋, 가중치 파일 SHA-256 및 런타임 버전은 model-provenance.json에 기록됩니다.
"""
    (output / "README.md").write_text(readme, encoding="utf-8")

    artifacts = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            artifacts.append({"path": path.relative_to(output).as_posix(), "sha256": sha256_file(path), "bytes": path.stat().st_size})
    _write_json(output / "manifest.json", {
        "format": FORMAT,
        "files": artifacts,
        "scientific_approved": False,
        "population_prediction_performed": False,
        "final_state_selection_performed": False,
    })
    return 1 if failures else 0


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser()
    value.add_argument("--upstream-root", required=True)
    value.add_argument("--upstream-manifest", required=True)
    value.add_argument("--dependency-directory")
    value.add_argument("--core-directory")
    value.add_argument("--store")
    value.add_argument("--project")
    value.add_argument("--job-id")
    value.add_argument("--output")
    value.add_argument("--ph", type=float, default=7.4)
    value.add_argument("--worker", nargs=2, metavar=("REQUEST", "RESULT"))
    return value


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.worker:
        return worker(Path(args.worker[0]), Path(args.worker[1]), Path(args.upstream_root), Path(args.upstream_manifest), Path(args.dependency_directory) if args.dependency_directory else None)
    for field in ("core_directory", "store", "project", "job_id", "output"):
        if not getattr(args, field):
            parser().error("--" + field.replace("_", "-") + " is required")
    try:
        validate_ph(args.ph)
    except ValueError as error:
        parser().error(str(error))
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
