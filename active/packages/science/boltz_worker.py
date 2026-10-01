"""독립 Boltz2 CRBN 보정 실행과 보수적 구조 비교.

참조 구조는 입력 템플릿이나 좌표 제약으로 전달하지 않는다. 결과 지표는 알려진
구조에 대한 기술적 재현 비교이며 효능, 승인 또는 E3 우월성을 뜻하지 않는다.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path, PurePosixPath
from typing import Any, Callable

import gemmi
import numpy as np
import yaml
from rdkit import Chem

from .molecules import rows
from .structures import aligned_rmsd

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}
REQUIRED_HELP_OPTIONS = (
    "--out_dir", "--cache", "--checkpoint", "--accelerator", "--recycling_steps",
    "--sampling_steps", "--diffusion_samples", "--max_parallel_samples", "--num_workers",
    "--preprocessing-threads", "--no_kernels", "--use_msa_server", "--max_msa_seqs",
)
FROZEN_CRBN_INPUT_SHA256 = "c1fddb2c8fa57055e57af21b81d8730d10ad6219087d433087654296a67d6661"
FROZEN_CRBN_CHECKPOINT_SHA256 = "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1"
MAX_FROZEN_CACHE_COPY_BYTES = 2 * 1024 * 1024 * 1024
MAX_FROZEN_RECEIPT_BYTES = 16 * 1024 * 1024


class BoltzWorkerError(ValueError):
    """검증 또는 실행 준비 실패."""


def _check(condition: bool, code: str) -> None:
    if not condition:
        raise BoltzWorkerError(code)


def sha256_file(path: Path) -> str:
    path = Path(path)
    _check(path.is_file() and not path.is_symlink(), "REGULAR_FILE_REQUIRED")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sequence(value: Any) -> str:
    return "".join(str(value or "").split()).upper()


def _block(path: Path):
    path = Path(path)
    _check(path.is_file() and not path.is_symlink(), "REFERENCE_REGULAR_FILE_REQUIRED")
    return gemmi.cif.read_file(str(path)).sole_block()


def _entity_material(block) -> tuple[dict[str, dict], dict[str, str], dict[str, str]]:
    entities = {str(r["id"]): r for r in rows(block, "_entity.")}
    polymers = {str(r["entity_id"]): r for r in rows(block, "_entity_poly.")}
    asym = {str(r["id"]): str(r["entity_id"]) for r in rows(block, "_struct_asym.")}
    return entities, polymers, asym


def _eligible_atom_rows(block, chain: str) -> list[dict]:
    result, seen = [], set()
    for row in rows(block, "_atom_site."):
        if str(row.get("label_asym_id")) != chain:
            continue
        model = row.get("pdbx_PDB_model_num")
        if model not in (None, False, "", "1", 1):
            continue
        if str(row.get("type_symbol", "")).upper() in {"H", "D"}:
            continue
        alt = row.get("label_alt_id")
        if alt not in (None, False, "", ".", "?", "A"):
            continue
        try:
            occupancy = float(row.get("occupancy", 1.0))
        except (TypeError, ValueError):
            raise BoltzWorkerError("ATOM_OCCUPANCY_INVALID")
        if occupancy <= 0:
            continue
        key = (str(row.get("label_seq_id")), str(row.get("label_atom_id")))
        _check(key not in seen, "AMBIGUOUS_ATOM_SITE")
        seen.add(key)
        try:
            row = dict(row)
            row["xyz"] = [float(row["Cartn_" + axis]) for axis in "xyz"]
        except (KeyError, TypeError, ValueError):
            raise BoltzWorkerError("ATOM_COORDINATE_INVALID")
        _check(np.isfinite(row["xyz"]).all(), "ATOM_COORDINATE_NONFINITE")
        result.append(row)
    return result


def prepare_input(reference: Path, metadata: Path, destination: Path, msa_mode: str) -> dict:
    """6BOY의 검증된 deposited construct로 좌표 없는 Boltz YAML을 만든다."""
    _check(msa_mode in {"server", "single_sequence"}, "MSA_MODE_REQUIRED")
    reference, metadata, destination = Path(reference), Path(metadata), Path(destination)
    meta = json.loads(metadata.read_text(encoding="utf-8"))
    required = {
        "pdb": "6BOY", "ccd": "RN6", "reference_target_chain": "C",
        "reference_e3_chain": "B", "reference_ligand_chain": "E", "e3_type": "CRBN",
    }
    for key, expected in required.items():
        _check(meta.get(key) == expected, "REFERENCE_METADATA_MISMATCH")

    block = _block(reference)
    entities, polymers, asym = _entity_material(block)
    selected = {
        "target": str(meta["reference_target_chain"]),
        "e3": str(meta["reference_e3_chain"]),
    }
    role_terms = {
        "target": ("bromodomain-containing protein 4", "brd4"),
        "e3": ("cereblon",),
    }
    proteins, chain_mapping = {}, {}
    for role, chain in selected.items():
        _check(chain in asym and asym[chain] in polymers, "REFERENCE_POLYMER_CHAIN_MISSING")
        entity_id = asym[chain]
        poly = polymers[entity_id]
        _check(str(poly.get("type", "")).lower() == "polypeptide(l)", "REFERENCE_NOT_L_POLYPEPTIDE")
        sequence = _sequence(poly.get("pdbx_seq_one_letter_code_can"))
        _check(bool(sequence) and not set(sequence) - AMINO_ACIDS, "REFERENCE_SEQUENCE_UNSUPPORTED")
        description = str(entities.get(entity_id, {}).get("pdbx_description", ""))
        _check(any(term in description.lower() for term in role_terms[role]), "REFERENCE_ROLE_NOT_VERIFIED")
        atoms = _eligible_atom_rows(block, chain)
        _check(bool(atoms), "REFERENCE_CHAIN_HAS_NO_ATOMS")
        auth = sorted({str(r.get("auth_asym_id")) for r in atoms})
        _check(len(auth) == 1 and auth[0] not in {"", "None", "False"}, "REFERENCE_AUTH_CHAIN_AMBIGUOUS")
        proteins[role] = sequence
        chain_mapping[role] = {
            "label_asym_id": chain, "auth_asym_id": auth[0], "entity_id": entity_id,
            "description": description, "sequence_length": len(sequence),
            "sequence_sha256": _sha_bytes(sequence.encode("ascii")),
            "sequence_source": "_entity_poly.pdbx_seq_one_letter_code_can (deposited construct 포함)",
        }

    excluded = []
    ddb1_found = False
    for chain, entity_id in asym.items():
        if chain in selected.values() or entity_id not in polymers:
            continue
        description = str(entities.get(entity_id, {}).get("pdbx_description", ""))
        if "dna damage-binding protein 1" in description.lower() or "ddb1" in description.lower():
            ddb1_found = True
        excluded.append({"label_asym_id": chain, "entity_id": entity_id, "description": description})
    _check(ddb1_found, "REFERENCE_DDB1_EXCLUSION_NOT_VERIFIED")

    ligand_chain = str(meta["reference_ligand_chain"])
    ligand_rows = _eligible_atom_rows(block, ligand_chain)
    _check(bool(ligand_rows), "REFERENCE_LIGAND_MISSING")
    _check({str(r.get("label_comp_id")) for r in ligand_rows} == {"RN6"}, "REFERENCE_LIGAND_CCD_MISMATCH")
    atom_names = [str(r.get("label_atom_id")) for r in ligand_rows]
    _check(len(atom_names) == len(set(atom_names)), "REFERENCE_LIGAND_ATOM_NAMES_AMBIGUOUS")

    protein_extra = {"msa": "empty"} if msa_mode == "single_sequence" else {}
    document = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "T", "sequence": proteins["target"], **protein_extra}},
            {"protein": {"id": "E", "sequence": proteins["e3"], **protein_extra}},
            {"ligand": {"id": "L", "ccd": "RN6"}},
        ],
    }
    encoded = yaml.safe_dump(document, sort_keys=False, allow_unicode=False).encode("utf-8")
    _check(not destination.exists(), "INPUT_EXISTS")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("xb") as stream:
        stream.write(encoded)
    return {
        "input_path": str(destination), "input_sha256": _sha_bytes(encoded), "msa_mode": msa_mode,
        "protein_sequences": proteins, "reference_chain_mapping": chain_mapping,
        "ligand": {"input_chain": "L", "ccd": "RN6", "reference_label_asym_id": ligand_chain,
                   "reference_heavy_atom_names": sorted(atom_names)},
        "input_chain_roles": {"target": "T", "e3": "E", "ligand": "L"},
        "explicit_exclusions": {"DDB1": "입력에서 제외됨", "other_reference_polymers": excluded},
        "omissions": ["참조 좌표", "template", "constraint", "potential", "DDB1"],
    }


def inspect_help(executable: Path, timeout: float = 30.0) -> dict:
    executable = Path(executable)
    _check(executable.is_file() and not executable.is_symlink(), "BOLTZ_EXECUTABLE_REQUIRED")
    try:
        proc = subprocess.run([str(executable), "predict", "--help"], capture_output=True,
                              text=True, encoding="utf-8", errors="replace", timeout=timeout,
                              shell=False, check=False,
                              creationflags=subprocess.CREATE_NO_WINDOW if hasattr(subprocess, "CREATE_NO_WINDOW") else 0)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise BoltzWorkerError("BOLTZ_HELP_FAILED") from error
    text = (proc.stdout or "") + (proc.stderr or "")
    _check(proc.returncode == 0 and bool(text.strip()), "BOLTZ_HELP_FAILED")
    missing = [option for option in REQUIRED_HELP_OPTIONS if option not in text]
    _check(not missing, "BOLTZ_OFFICIAL_OPTIONS_MISSING:" + ",".join(missing))
    _check("--use_potentials" not in REQUIRED_HELP_OPTIONS, "POTENTIALS_FORBIDDEN")
    has_seed = re.search(r"(?<![\w-])--seed(?:[\s=,]|$)", text) is not None
    return {"text": text, "sha256": _sha_bytes(text.encode("utf-8")), "has_seed": has_seed}


def build_command(*, executable: Path, input_yaml: Path, out_dir: Path, cache: Path,
                  checkpoint: Path, accelerator: str, seed: int, msa_mode: str,
                  max_msa_seqs: int, help_has_seed: bool) -> list[str]:
    _check(type(seed) is int and seed >= 0, "SEED_MUST_BE_NONNEGATIVE_INTEGER")
    _check(accelerator in {"gpu", "cpu"}, "ACCELERATOR_INVALID")
    _check(msa_mode in {"server", "single_sequence"}, "MSA_MODE_REQUIRED")
    _check(type(max_msa_seqs) is int and max_msa_seqs > 0, "MAX_MSA_SEQS_INVALID")
    command = [
        str(executable), "predict", str(input_yaml), "--out_dir", str(out_dir),
        "--cache", str(cache), "--checkpoint", str(checkpoint), "--accelerator", accelerator,
        "--recycling_steps", "3", "--sampling_steps", "200", "--diffusion_samples", "1",
        "--max_parallel_samples", "1", "--num_workers", "0", "--preprocessing-threads", "1",
        "--no_kernels", "--max_msa_seqs", str(max_msa_seqs),
    ]
    if msa_mode == "server":
        command.append("--use_msa_server")
    if help_has_seed:
        command.extend(["--seed", str(seed)])
    _check("--use_potentials" not in command, "POTENTIALS_FORBIDDEN")
    return command


def _wrapper_python(executable: Path) -> Path:
    candidates = [
        executable.parent / "python", executable.parent / "python.exe",
        executable.parent.parent / "python.exe",
    ]
    for candidate in candidates:
        if candidate.is_file() and not candidate.is_symlink():
            return candidate
    raise BoltzWorkerError("BOLTZ_SEED_WRAPPER_PYTHON_NOT_FOUND")


def _write_seed_wrapper(path: Path) -> None:
    source = '''import importlib.metadata, os, random, sys
seed=int(sys.argv[1]); args=sys.argv[2:]
os.environ["PYTHONHASHSEED"]=str(seed)
random.seed(seed)
import numpy as np
np.random.seed(seed % (2**32))
import torch
torch.manual_seed(seed)
if torch.cuda.is_available(): torch.cuda.manual_seed_all(seed)
sys.argv=["boltz", *args]
points=[p for p in importlib.metadata.entry_points(group="console_scripts") if p.name=="boltz"]
if len(points)!=1: raise RuntimeError("BOLTZ_CONSOLE_ENTRYPOINT_AMBIGUOUS")
result=points[0].load()()
raise SystemExit(0 if result is None else result)
'''
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(source)


def _run_process(command: list[str], log: Path, timeout: float,
                 cancel: Any = None) -> tuple[int | None, str, float]:
    _check(timeout > 0 and math.isfinite(timeout), "TIMEOUT_INVALID")
    started = time.monotonic()
    with log.open("xb") as output:
        try:
            proc = subprocess.Popen(
                command, stdout=output, stderr=subprocess.STDOUT, shell=False,
                creationflags=subprocess.CREATE_NO_WINDOW if hasattr(subprocess, "CREATE_NO_WINDOW") else 0,
            )
        except OSError as error:
            raise BoltzWorkerError("BOLTZ_PROCESS_START_FAILED") from error
        reason = "completed"
        while proc.poll() is None:
            cancelled = cancel() if callable(cancel) else bool(cancel is not None and cancel.is_set())
            if cancelled:
                reason = "cancelled"
            elif time.monotonic() - started > timeout:
                reason = "timeout"
            else:
                time.sleep(0.2)
                continue
            proc.terminate()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=10)
            break
    return proc.returncode, reason, time.monotonic() - started


def _entity_poly_seq_sequences(block) -> dict[str, str]:
    """Build strict canonical sequences from numbered entity_poly_seq rows."""
    grouped: dict[str, dict[int, str]] = {}
    for row in rows(block, "_entity_poly_seq."):
        entity = str(row.get("entity_id"))
        try:
            index = int(str(row.get("num")))
        except (TypeError, ValueError):
            raise BoltzWorkerError("PREDICTION_ENTITY_POLY_SEQ_INDEX_INVALID")
        _check(index > 0, "PREDICTION_ENTITY_POLY_SEQ_INDEX_INVALID")
        hetero = row.get("hetero")
        _check(hetero in (None, False, "", ".", "?", "n", "N", "no", "NO", "No", 0, "0"),
               "PREDICTION_ENTITY_POLY_SEQ_HETEROGENEOUS")
        monomer = str(row.get("mon_id", "")).upper()
        _check(monomer in THREE_TO_ONE, "PREDICTION_ENTITY_POLY_SEQ_UNKNOWN_MONOMER")
        entity_rows = grouped.setdefault(entity, {})
        _check(index not in entity_rows, "PREDICTION_ENTITY_POLY_SEQ_DUPLICATE_INDEX")
        entity_rows[index] = THREE_TO_ONE[monomer]

    result = {}
    for entity, entity_rows in grouped.items():
        indices = sorted(entity_rows)
        _check(indices == list(range(1, len(indices) + 1)),
               "PREDICTION_ENTITY_POLY_SEQ_GAP")
        result[entity] = "".join(entity_rows[index] for index in indices)
    return result


def _prediction_sequences(block) -> dict[str, str]:
    _, polymers, asym = _entity_material(block)
    poly_seq = _entity_poly_seq_sequences(block)
    result = {}
    for chain, entity in asym.items():
        if entity not in polymers:
            continue
        row = polymers[entity]
        if str(row.get("type", "")).lower() != "polypeptide(l)":
            continue
        derived = poly_seq.get(entity)
        _check(bool(derived), "PREDICTION_ENTITY_POLY_SEQ_MISSING")
        canonical = _sequence(row.get("pdbx_seq_one_letter_code_can"))
        # An all-X value carries no identity information. Any valid canonical
        # value remains an independent assertion and must agree exactly.
        if canonical and not set(canonical) - AMINO_ACIDS:
            _check(canonical == derived, "PREDICTION_SEQUENCE_SOURCES_DISAGREE")
        result[chain] = derived
    return result


def exact_sequence_mapping(reference: str, prediction: str) -> list[tuple[int, int]]:
    """전역 정렬의 완전 일치 특수 경우. 삽입·삭제·치환은 보수적으로 거부한다."""
    reference, prediction = _sequence(reference), _sequence(prediction)
    _check(reference == prediction and bool(reference), "PREDICTION_SEQUENCE_NOT_EXACT")
    return [(index, index) for index in range(1, len(reference) + 1)]


def _assign_prediction_chains(block, packet: dict) -> dict[str, str]:
    sequences = _prediction_sequences(block)
    assigned = {}
    for role in ("target", "e3"):
        expected = packet["protein_sequences"][role]
        candidates = [chain for chain, sequence in sequences.items() if sequence == expected]
        _check(len(candidates) == 1, "PREDICTION_CHAIN_MAPPING_AMBIGUOUS")
        exact_sequence_mapping(expected, sequences[candidates[0]])
        assigned[role] = candidates[0]
    _check(assigned["target"] != assigned["e3"], "PREDICTION_CHAIN_REUSED")
    return assigned


def _protein_atoms(block, chain: str, role: str) -> dict[tuple[str, int, str], dict]:
    result = {}
    for row in _eligible_atom_rows(block, chain):
        value = row.get("label_seq_id")
        try:
            sequence_index = int(str(value))
        except (TypeError, ValueError):
            raise BoltzWorkerError("PROTEIN_LABEL_SEQ_ID_REQUIRED")
        atom = str(row.get("label_atom_id"))
        key = (role, sequence_index, atom)
        _check(key not in result, "PROTEIN_ATOM_MAPPING_AMBIGUOUS")
        result[key] = row
    _check(bool(result), "PREDICTION_PROTEIN_ATOMS_MISSING")
    return result


def _ligand_atoms(block, expected_names: set[str]) -> tuple[str, dict[str, dict]]:
    _, polymers, asym = _entity_material(block)
    polymer_chains = {chain for chain, entity in asym.items() if entity in polymers}
    candidates = []
    for chain in asym:
        if chain in polymer_chains:
            continue
        atoms = _eligible_atom_rows(block, chain)
        if not atoms or {str(r.get("label_comp_id")) for r in atoms} != {"RN6"}:
            continue
        named = {str(r.get("label_atom_id")): r for r in atoms}
        if len(named) == len(atoms) and set(named) == expected_names:
            candidates.append((chain, named))
    # 일부 Boltz CIF는 비고분자 struct_asym 행을 생략하므로 atom_site도 확인한다.
    if not candidates:
        chains = sorted({str(r.get("label_asym_id")) for r in rows(block, "_atom_site.")})
        for chain in chains:
            if chain in polymer_chains:
                continue
            atoms = _eligible_atom_rows(block, chain)
            if atoms and {str(r.get("label_comp_id")) for r in atoms} == {"RN6"}:
                named = {str(r.get("label_atom_id")): r for r in atoms}
                if len(named) == len(atoms) and set(named) == expected_names:
                    candidates.append((chain, named))
    _check(len(candidates) == 1, "PREDICTION_LIGAND_MAPPING_AMBIGUOUS")
    return candidates[0]


def _apply_fit(xyz: np.ndarray, fit: dict) -> np.ndarray:
    return xyz @ np.asarray(fit["rotation_row_vectors"], float) + np.asarray(fit["translation_A"], float)


def _rmsd(reference: np.ndarray, prediction_fitted: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((prediction_fitted - reference) ** 2, axis=1))))


def _contact_set(proteins: dict, ligand: dict, cutoff: float = 4.0) -> set[tuple]:
    contacts = set()
    for (role, sequence_index, _atom_name), row in proteins.items():
        p = np.asarray(row["xyz"], float)
        for ligand_name, ligand_row in ligand.items():
            if np.linalg.norm(p - np.asarray(ligand_row["xyz"], float)) <= cutoff:
                contacts.add((role, sequence_index, ligand_name))
    return contacts


def _clashes(proteins: dict, ligand: dict) -> list[dict]:
    periodic = Chem.GetPeriodicTable()
    result = []
    for (role, sequence_index, atom_name), row in proteins.items():
        protein_element = str(row.get("type_symbol", "")).title()
        try:
            protein_radius = periodic.GetRvdw(protein_element)
        except RuntimeError:
            raise BoltzWorkerError("UNKNOWN_PROTEIN_ELEMENT")
        for ligand_name, ligand_row in ligand.items():
            ligand_element = str(ligand_row.get("type_symbol", "")).title()
            try:
                ligand_radius = periodic.GetRvdw(ligand_element)
            except RuntimeError:
                raise BoltzWorkerError("UNKNOWN_LIGAND_ELEMENT")
            distance = float(np.linalg.norm(np.asarray(row["xyz"]) - np.asarray(ligand_row["xyz"])))
            threshold = 0.75 * (protein_radius + ligand_radius)
            if distance < threshold:
                result.append({"role": role, "sequence_index": sequence_index,
                               "protein_atom": atom_name, "ligand_atom": ligand_name,
                               "distance_A": distance, "descriptive_threshold_A": threshold})
    return result


def _load_groups(path: Path | None, ligand_names: set[str]) -> dict[str, list[str]]:
    if path is None:
        return {}
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    _check(isinstance(value, dict) and bool(value), "LIGAND_GROUP_MAPPING_INVALID")
    result = {}
    for name, atoms in value.items():
        _check(isinstance(name, str) and name and isinstance(atoms, list) and atoms,
               "LIGAND_GROUP_MAPPING_INVALID")
        _check(all(isinstance(atom, str) for atom in atoms) and len(atoms) == len(set(atoms)),
               "LIGAND_GROUP_MAPPING_INVALID")
        _check(set(atoms) <= ligand_names, "LIGAND_GROUP_UNKNOWN_ATOM")
        result[name] = atoms
    return result


def compare_structures(reference: Path, prediction: Path, packet: dict,
                       ligand_groups: Path | None = None) -> dict:
    _check(sha256_file(reference) != sha256_file(prediction), "REFERENCE_PREDICTION_IDENTICAL_FILE")
    ref_block, pred_block = _block(reference), _block(prediction)
    pred_chains = _assign_prediction_chains(pred_block, packet)
    ref_chains = {role: packet["reference_chain_mapping"][role]["label_asym_id"] for role in ("target", "e3")}
    ref_proteins, pred_proteins = {}, {}
    for role in ("target", "e3"):
        ref_proteins.update(_protein_atoms(ref_block, ref_chains[role], role))
        pred_proteins.update(_protein_atoms(pred_block, pred_chains[role], role))

    ca_pairs = {}
    for role in ("target", "e3"):
        keys = sorted(set(k for k in ref_proteins if k[0] == role and k[2] == "CA") &
                      set(k for k in pred_proteins if k[0] == role and k[2] == "CA"), key=lambda k: k[1])
        _check(len(keys) >= 3, "INSUFFICIENT_PAIRED_CA_ATOMS")
        ca_pairs[role] = keys
    target_ref = np.asarray([ref_proteins[k]["xyz"] for k in ca_pairs["target"]], float)
    target_pred = np.asarray([pred_proteins[k]["xyz"] for k in ca_pairs["target"]], float)
    fit = aligned_rmsd(target_ref, target_pred)
    e3_ref = np.asarray([ref_proteins[k]["xyz"] for k in ca_pairs["e3"]], float)
    e3_pred = _apply_fit(np.asarray([pred_proteins[k]["xyz"] for k in ca_pairs["e3"]], float), fit)

    ligand_names = set(packet["ligand"]["reference_heavy_atom_names"])
    _, ref_ligand = _ligand_atoms(ref_block, ligand_names)
    pred_ligand_chain, pred_ligand = _ligand_atoms(pred_block, ligand_names)
    ordered_ligand = sorted(ligand_names)
    ligand_ref = np.asarray([ref_ligand[name]["xyz"] for name in ordered_ligand], float)
    ligand_pred = _apply_fit(np.asarray([pred_ligand[name]["xyz"] for name in ordered_ligand], float), fit)

    common_protein = set(ref_proteins) & set(pred_proteins)
    ref_common = {key: ref_proteins[key] for key in common_protein}
    pred_common = {key: pred_proteins[key] for key in common_protein}
    ref_contacts = _contact_set(ref_common, ref_ligand)
    pred_contacts = _contact_set(pred_common, pred_ligand)
    union, intersection = ref_contacts | pred_contacts, ref_contacts & pred_contacts
    groups = _load_groups(ligand_groups, ligand_names)
    group_metrics = {}
    name_index = {name: index for index, name in enumerate(ordered_ligand)}
    for name, atoms in groups.items():
        indices = [name_index[atom] for atom in atoms]
        group_metrics[name] = {"explicit_ccd_atoms": atoms,
                               "RMSD_after_target_alignment_A": _rmsd(ligand_ref[indices], ligand_pred[indices])}

    return {
        "mapping": {
            "reference_chains": ref_chains, "prediction_chains": pred_chains,
            "prediction_ligand_chain": pred_ligand_chain,
            "protein_correspondence": "deposited construct 전체 서열의 완전 일치 전역 대응",
            "ligand_correspondence": "RN6 CCD heavy-atom name 완전 일치",
            "paired_target_CA_count": len(ca_pairs["target"]),
            "paired_e3_CA_count": len(ca_pairs["e3"]),
            "paired_ligand_heavy_atom_count": len(ordered_ligand),
        },
        "metrics": {
            "target_CA_RMSD_A": fit["rmsd_A"],
            "e3_CA_RMSD_after_target_alignment_A": _rmsd(e3_ref, e3_pred),
            "ligand_heavy_atom_RMSD_after_target_alignment_A": _rmsd(ligand_ref, ligand_pred),
            "contact_definition": "공통 관측 protein heavy atoms와 ligand heavy atoms, residue/ligand-atom contact <= 4.0 Å",
            "reference_contact_count": len(ref_contacts), "prediction_contact_count": len(pred_contacts),
            "contact_intersection_count": len(intersection), "contact_union_count": len(union),
            "contact_jaccard": (len(intersection) / len(union)) if union else None,
            "reference_descriptive_clashes": _clashes(ref_common, ref_ligand),
            "prediction_descriptive_clashes": _clashes(pred_common, pred_ligand),
            "clash_definition": "intermolecular distance < 0.75 × van der Waals radii sum; 기술 지표만 해당",
            "explicit_ligand_group_metrics": group_metrics,
        },
        "target_alignment": fit,
        "interpretation_ko": "구조 재현 기술 지표이며 효능, 분해 성공, 승인 또는 E3 우월성 판정 기준이 아니다.",
    }


def _validate_confidence(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise BoltzWorkerError("CONFIDENCE_JSON_INVALID") from error
    _check(isinstance(value, dict), "CONFIDENCE_JSON_INVALID")
    numeric = 0

    def visit(item: Any) -> None:
        nonlocal numeric
        if isinstance(item, bool) or item is None or isinstance(item, str):
            return
        if isinstance(item, (int, float)):
            _check(math.isfinite(float(item)), "CONFIDENCE_NONFINITE")
            numeric += 1
        elif isinstance(item, list):
            for child in item:
                visit(child)
        elif isinstance(item, dict):
            for child in item.values():
                visit(child)
        else:
            raise BoltzWorkerError("CONFIDENCE_VALUE_INVALID")

    visit(value)
    _check(numeric > 0, "CONFIDENCE_NUMERIC_VALUE_MISSING")
    return value


def discover_outputs(out_dir: Path) -> dict:
    out_dir = Path(out_dir)
    if not out_dir.is_dir():
        return {"status": "missing_prediction", "reason": "BOLTZ_OUTPUT_DIRECTORY_MISSING"}
    cifs = sorted(path for path in out_dir.rglob("*.cif") if path.is_file() and not path.is_symlink())
    if not cifs:
        return {"status": "missing_prediction", "reason": "PREDICTION_CIF_MISSING"}
    if len(cifs) != 1:
        return {"status": "malformed_prediction", "reason": "PREDICTION_CIF_AMBIGUOUS",
                "candidates": [str(path) for path in cifs]}
    confidence = sorted(path for path in out_dir.rglob("*.json")
                        if path.is_file() and not path.is_symlink() and "confidence" in path.name.lower())
    if len(confidence) != 1:
        return {"status": "malformed_prediction", "reason": "CONFIDENCE_JSON_MISSING_OR_AMBIGUOUS",
                "prediction": str(cifs[0]), "confidence_candidates": [str(path) for path in confidence]}
    return {"status": "located", "prediction": str(cifs[0]), "confidence": str(confidence[0])}


def assess_outputs(out_dir: Path, reference: Path, packet: dict,
                   ligand_groups: Path | None = None) -> dict:
    found = discover_outputs(out_dir)
    if found["status"] != "located":
        return found
    try:
        confidence = _validate_confidence(Path(found["confidence"]))
        comparison = compare_structures(reference, Path(found["prediction"]), packet, ligand_groups)
    except (BoltzWorkerError, OSError, ValueError, RuntimeError) as error:
        return {"status": "inspection_failed", "reason": str(error),
                "prediction": found["prediction"], "confidence": found["confidence"]}
    return {
        "status": "success", "prediction": found["prediction"], "confidence_file": found["confidence"],
        "model_confidence": confidence,
        "confidence_interpretation_ko": "모델이 출력한 예측 신뢰도이며 측정된 구조 재현도나 효능이 아니다.",
        "comparison": comparison,
    }


def _hash_tree(root: Path) -> dict[str, str]:
    if not root.is_dir():
        return {}
    result = {}
    for path in sorted(root.rglob("*")):
        _check(not path.is_symlink(), "OUTPUT_SYMLINK_FORBIDDEN")
        if path.is_file():
            result[path.relative_to(root).as_posix()] = sha256_file(path)
    return result


def _tool_environment(executable: Path) -> dict:
    """Boltz와 같은 환경의 Python에서 패키지 버전과 실제 가속기 상태를 읽는다."""
    try:
        python = _wrapper_python(executable)
        source = (
            "import importlib.metadata,json,torch;"
            "print(json.dumps({'boltz_version':importlib.metadata.version('boltz'),"
            "'torch_version':torch.__version__,'cuda_available':torch.cuda.is_available(),"
            "'cuda_runtime':torch.version.cuda,'cuda_device_count':torch.cuda.device_count(),"
            "'cuda_devices':[torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]}))"
        )
        proc = subprocess.run(
            [str(python), "-c", source], capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=30, shell=False, check=False,
            creationflags=subprocess.CREATE_NO_WINDOW if hasattr(subprocess, "CREATE_NO_WINDOW") else 0,
        )
        _check(proc.returncode == 0, "BOLTZ_ENVIRONMENT_PROBE_FAILED")
        value = json.loads(proc.stdout)
        _check(isinstance(value, dict) and bool(value.get("boltz_version")),
               "BOLTZ_ENVIRONMENT_PROBE_FAILED")
        value["python_path"] = str(python)
        return value
    except (OSError, subprocess.TimeoutExpired, json.JSONDecodeError, BoltzWorkerError) as error:
        return {"probe_succeeded": False, "error": str(error) or type(error).__name__}


def _validated_frozen_relative_path(value: str) -> PurePosixPath:
    """Validate every receipt path before deciding whether it is reusable."""
    _check(isinstance(value, str) and value and "\\" not in value,
           "FROZEN_CACHE_PATH_INVALID")
    path = PurePosixPath(value)
    _check(not path.is_absolute() and path.as_posix() == value and
           all(part not in {"", ".", ".."} for part in path.parts),
           "FROZEN_CACHE_PATH_TRAVERSAL")
    _check(len(path.parts) >= 2 and path.parts[0] == "boltz_results_6BOY_CRBN",
           "FROZEN_CACHE_PATH_INVALID")
    return path


def _frozen_relative_path(value: str) -> PurePosixPath:
    path = _validated_frozen_relative_path(value)
    _check(len(path.parts) >= 3 and path.parts[1] in {"msa", "processed"},
           "FROZEN_CACHE_SUBTREE_FORBIDDEN")
    return path


def _bounded_bytes(path: Path, limit: int, code: str) -> bytes:
    _check(path.is_file() and not path.is_symlink(), code)
    _check(path.stat().st_size <= limit, code)
    with path.open("rb") as stream:
        value = stream.read(limit + 1)
    _check(len(value) <= limit, code)
    return value


def _assert_existing_parents_not_symlinks(path: Path, code: str) -> None:
    current = Path(path)
    while True:
        if current.exists():
            _check(not current.is_symlink(), code)
        if current.parent == current:
            break
        current = current.parent


def _tree_files(root: Path) -> dict[str, Path]:
    result: dict[str, Path] = {}
    _check(root.is_dir() and not root.is_symlink(), "FROZEN_SOURCE_OUTPUT_DIRECTORY_REQUIRED")
    resolved_root = root.resolve(strict=True)
    for path in sorted(root.rglob("*")):
        _check(not path.is_symlink(), "FROZEN_SOURCE_SYMLINK_FORBIDDEN")
        if not path.is_file():
            continue
        try:
            path.resolve(strict=True).relative_to(resolved_root)
        except ValueError as error:
            raise BoltzWorkerError("FROZEN_CACHE_PATH_TRAVERSAL") from error
        relative = path.relative_to(root).as_posix()
        _validated_frozen_relative_path(relative)
        result[relative] = path
    return result


def _frozen_record_is_bound(path: Path) -> bool:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise BoltzWorkerError("FROZEN_CACHE_RECORD_INVALID") from error
    text = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return path.stem == "6BOY_CRBN" and "6BOY_CRBN" in text and "msa" in text.lower()


def _install_frozen_seed_cache(*, source_directory: Path, destination: Path,
                               input_yaml: Path, reference: Path, checkpoint: Path,
                               executable: Path, environment: dict, settings: dict,
                               help_sha256: str) -> dict:
    source_directory = Path(source_directory)
    _check(source_directory.is_dir() and not source_directory.is_symlink(),
           "FROZEN_SEED_DIRECTORY_REQUIRED")
    receipt_path = source_directory / "receipt.json"
    _check(receipt_path.is_file() and not receipt_path.is_symlink(),
           "FROZEN_SOURCE_RECEIPT_REQUIRED")
    receipt_bytes = _bounded_bytes(receipt_path, MAX_FROZEN_RECEIPT_BYTES,
                                   "FROZEN_SOURCE_RECEIPT_INVALID")
    try:
        original = json.loads(receipt_bytes)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise BoltzWorkerError("FROZEN_SOURCE_RECEIPT_INVALID") from error
    hashes = original.get("hashes")
    _check(isinstance(hashes, dict), "FROZEN_SOURCE_HASHES_REQUIRED")
    expected = hashes.get("processed_output_files_sha256")
    _check(isinstance(expected, dict) and expected, "FROZEN_SOURCE_OUTPUT_HASHES_REQUIRED")
    _check(original.get("exit_code") == 0 and original.get("process_state") == "completed",
           "FROZEN_SOURCE_PROCESS_NOT_COMPLETED")
    _check(hashes.get("input_yaml_sha256") == FROZEN_CRBN_INPUT_SHA256 == sha256_file(input_yaml),
           "FROZEN_SOURCE_INPUT_MISMATCH")
    original_input = source_directory / "6BOY_CRBN.yaml"
    _check(original_input.is_file() and not original_input.is_symlink(),
           "FROZEN_SOURCE_INPUT_FILE_REQUIRED")
    _check(sha256_file(original_input) == hashes.get("input_yaml_sha256"),
           "FROZEN_SOURCE_ORIGINAL_INPUT_MISMATCH")
    _check(hashes.get("checkpoint_sha256") == FROZEN_CRBN_CHECKPOINT_SHA256 == sha256_file(checkpoint),
           "FROZEN_SOURCE_CHECKPOINT_MISMATCH")
    _check(hashes.get("reference_sha256") == sha256_file(reference),
           "FROZEN_SOURCE_REFERENCE_MISMATCH")
    _check(hashes.get("executable_sha256") == sha256_file(executable),
           "FROZEN_SOURCE_EXECUTABLE_MISMATCH")
    _check(original.get("tool_environment") == environment, "FROZEN_SOURCE_ENVIRONMENT_MISMATCH")
    _check(original.get("settings") == settings, "FROZEN_SOURCE_SETTINGS_MISMATCH")
    _check(original.get("tool", {}).get("help_sha256") == help_sha256,
           "FROZEN_SOURCE_HELP_MISMATCH")

    source_output = source_directory / "boltz_output"
    actual_files = _tree_files(source_output)
    validated_expected: dict[str, str] = {}
    for stored, digest in expected.items():
        relative = _validated_frozen_relative_path(stored)
        _check(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
               "FROZEN_SOURCE_HASH_INVALID")
        validated_expected[relative.as_posix()] = digest
    _check(set(actual_files) == set(validated_expected),
           "FROZEN_SOURCE_OUTPUT_FILE_SET_MISMATCH")

    selected: dict[str, tuple[Path, str, int]] = {}
    source_snapshot: dict[str, str] = {}
    total = 0
    for stored, digest in validated_expected.items():
        source = actual_files[stored]
        actual_digest = sha256_file(source)
        _check(actual_digest == digest, "FROZEN_SOURCE_FILE_HASH_MISMATCH:" + stored)
        source_snapshot[stored] = actual_digest
        relative = PurePosixPath(stored)
        if len(relative.parts) < 3 or relative.parts[1] not in {"msa", "processed"}:
            continue
        size = source.stat().st_size
        total += size
        _check(total <= MAX_FROZEN_CACHE_COPY_BYTES, "FROZEN_CACHE_COPY_TOO_LARGE")
        selected[stored] = (source, digest, size)
    records = [name for name in selected if "/processed/records/" in "/" + name]
    processed_msa = [name for name in selected if "/processed/msa/" in "/" + name]
    _check(bool(records) and bool(processed_msa), "FROZEN_CACHE_REQUIRED_FILES_MISSING")
    _check(all(_frozen_record_is_bound(selected[name][0]) for name in records),
           "FROZEN_CACHE_RECORD_BINDING_MISMATCH")

    _assert_existing_parents_not_symlinks(destination.parent,
                                          "FROZEN_CACHE_DESTINATION_PARENT_SYMLINK")
    _check(not destination.exists(), "FROZEN_CACHE_DESTINATION_EXISTS")
    destination_parent = destination.parent.resolve(strict=True)
    copied = {}
    for relative, (source, digest, size) in selected.items():
        target = destination.joinpath(*PurePosixPath(relative).parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        _assert_existing_parents_not_symlinks(target.parent,
                                              "FROZEN_CACHE_DESTINATION_PARENT_SYMLINK")
        try:
            target.parent.resolve(strict=True).relative_to(destination_parent)
        except ValueError as error:
            raise BoltzWorkerError("FROZEN_CACHE_PATH_TRAVERSAL") from error
        copied_bytes = 0
        with source.open("rb") as reader, target.open("xb") as writer:
            while True:
                chunk = reader.read(min(1024 * 1024, size - copied_bytes + 1))
                if not chunk:
                    break
                copied_bytes += len(chunk)
                _check(copied_bytes <= size, "FROZEN_SOURCE_CHANGED_DURING_COPY")
                writer.write(chunk)
        _check(copied_bytes == size and sha256_file(source) == digest,
               "FROZEN_SOURCE_CHANGED_DURING_COPY")
        _check(sha256_file(target) == digest, "FROZEN_CACHE_COPY_HASH_MISMATCH")
        copied[relative] = digest
    return {
        "source_directory": str(source_directory),
        "source_receipt_path": str(receipt_path),
        "source_receipt_sha256": _sha_bytes(receipt_bytes),
        "original_status": original.get("status"),
        "original_inspection": original.get("inspection"),
        "original_failure_reason": original.get("failure_reason"),
        "original_exit_code": original.get("exit_code"),
        "original_process_state": original.get("process_state"),
        "reuse_scope": ["msa", "processed"],
        "predictions_reused": False,
        "copied_files_sha256": copied,
        "source_output_files_sha256": source_snapshot,
        "copied_bytes": total,
    }


def _verify_frozen_post_run(reuse: dict, boltz_out: Path) -> None:
    copied = reuse["copied_files_sha256"]
    current_tree = _hash_tree(boltz_out)
    current_reusable = {name: digest for name, digest in current_tree.items()
                        if len(PurePosixPath(name).parts) >= 3 and
                        PurePosixPath(name).parts[0] == "boltz_results_6BOY_CRBN" and
                        PurePosixPath(name).parts[1] in {"msa", "processed"}}
    _check(set(current_reusable) == set(copied), "FROZEN_CACHE_POST_RUN_FILE_SET_MISMATCH")
    manifest_changes = []
    immutable_changes = []
    for name, expected_digest in copied.items():
        if current_reusable[name] == expected_digest:
            continue
        if name.endswith("/processed/manifest.json"):
            manifest_changes.append({"path": name, "source_sha256": expected_digest,
                                     "post_run_sha256": current_reusable[name],
                                     "change_source": "official Boltz processed-manifest rewrite"})
        else:
            immutable_changes.append(name)
    reuse["post_run_manifest_changes"] = manifest_changes
    reuse["post_run_immutable_changed_files"] = sorted(immutable_changes)
    reuse["post_run_files_unchanged"] = not immutable_changes and not manifest_changes
    _check(not immutable_changes, "FROZEN_CACHE_CHANGED_AFTER_RUN")

    source_directory = Path(reuse["source_directory"])
    source_tree = _tree_files(source_directory / "boltz_output")
    observed_source = {name: sha256_file(path) for name, path in source_tree.items()}
    _check(observed_source == reuse["source_output_files_sha256"],
           "FROZEN_SOURCE_CHANGED_AFTER_RUN")
    reuse["source_cache_unchanged"] = True


def run_seed(*, reference: Path, metadata: Path, output_root: Path, executable: Path,
             checkpoint: Path, cache: Path, seed: int, msa_mode: str,
             accelerator: str = "gpu", max_msa_seqs: int = 8192,
             timeout: float = 21600.0, cancel: Any = None,
             ligand_groups: Path | None = None,
             frozen_seed_directory: Path | None = None) -> dict:
    """한 seed를 새 전용 디렉터리에서 실행하고 성공/실패 receipt를 항상 기록한다."""
    _check(type(seed) is int and seed >= 0, "SEED_MUST_BE_NONNEGATIVE_INTEGER")
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    run_dir = output_root / f"seed-{seed}"
    run_dir.mkdir(exist_ok=False)
    receipt_path = run_dir / "receipt.json"
    input_yaml, boltz_out, log = run_dir / "6BOY_CRBN.yaml", run_dir / "boltz_output", run_dir / "inference.log"
    started_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    receipt: dict[str, Any] = {
        "format": "tpd-boltz-crbn-calibration/1.0", "seed": seed, "started_utc": started_utc,
        "status": "failure", "성공": False, "calibration_complete": False,
        "known_structure_training_overlap_possible": True,
        "single_seed_is_broad_calibration": False,
        "효능_판정": "하지 않음",
    }
    try:
        executable, checkpoint, cache = Path(executable), Path(checkpoint), Path(cache)
        _check(checkpoint.is_file() and not checkpoint.is_symlink(), "CHECKPOINT_REGULAR_FILE_REQUIRED")
        _check(executable.is_file() and not executable.is_symlink(), "BOLTZ_EXECUTABLE_REQUIRED")
        # 입력 모델과 실행 파일의 증거는 어떤 외부 프로세스도 시작하기 전에 고정한다.
        receipt["hashes"] = {
            "checkpoint_sha256": sha256_file(checkpoint),
            "reference_sha256": sha256_file(reference),
            "executable_sha256": sha256_file(executable),
        }
        cache.mkdir(parents=True, exist_ok=True)
        _check(cache.is_dir() and not cache.is_symlink(), "CACHE_DIRECTORY_REQUIRED")
        receipt["tool_environment"] = _tool_environment(executable)
        packet = prepare_input(reference, metadata, input_yaml, msa_mode)
        help_info = inspect_help(executable)
        base_command = build_command(executable=executable, input_yaml=input_yaml, out_dir=boltz_out,
                                     cache=cache, checkpoint=checkpoint, accelerator=accelerator, seed=seed,
                                     msa_mode=msa_mode, max_msa_seqs=max_msa_seqs,
                                     help_has_seed=help_info["has_seed"])
        settings = {"accelerator": accelerator, "recycling_steps": 3, "sampling_steps": 200,
                    "diffusion_samples": 1, "max_parallel_samples": 1, "num_workers": 0,
                    "preprocessing_threads": 1, "no_kernels": True, "use_potentials": False,
                    "msa_mode": msa_mode, "max_msa_seqs": max_msa_seqs}
        if frozen_seed_directory is not None:
            receipt["frozen_cache_reuse"] = _install_frozen_seed_cache(
                source_directory=Path(frozen_seed_directory), destination=boltz_out,
                input_yaml=input_yaml, reference=Path(reference), checkpoint=checkpoint,
                executable=executable, environment=receipt["tool_environment"], settings=settings,
                help_sha256=help_info["sha256"])
        seed_application = {"seed": seed, "method": "official_cli_--seed", "actually_applied": True}
        command = base_command
        if not help_info["has_seed"]:
            wrapper = run_dir / "seeded_boltz_entrypoint.py"
            _write_seed_wrapper(wrapper)
            python = _wrapper_python(executable)
            command = [str(python), str(wrapper), str(seed), *base_command[1:]]
            seed_application = {
                "seed": seed, "method": "same-environment Python console-entrypoint wrapper",
                "actually_applied": True,
                "applied_to": ["random", "numpy", "torch CPU", "torch CUDA when available"],
                "wrapper_sha256": sha256_file(wrapper),
            }
        exit_code, process_state, elapsed = _run_process(command, log, timeout, cancel)
        receipt.update({
            "exit_code": exit_code, "process_state": process_state, "elapsed_seconds": elapsed,
            "settings": settings,
            "seed_application": seed_application, "input": packet,
            "tool": {"help_sha256": help_info["sha256"],
                     "version_source": "importlib.metadata from the Boltz environment"},
        })
        receipt["hashes"]["input_yaml_sha256"] = sha256_file(input_yaml)
        assessment = assess_outputs(boltz_out, reference, packet, ligand_groups)
        receipt["inspection"] = assessment
        reuse = receipt.get("frozen_cache_reuse")
        if isinstance(reuse, dict):
            _verify_frozen_post_run(reuse, boltz_out)
        if exit_code == 0 and process_state == "completed" and assessment.get("status") == "success":
            receipt["status"], receipt["성공"] = "success", True
        else:
            receipt["failure_reason"] = (assessment.get("reason") or
                                         ("BOLTZ_NONZERO_EXIT" if exit_code else process_state))
    except Exception as error:
        if isinstance(error, KeyboardInterrupt):
            raise
        receipt["failure_reason"] = str(error) or type(error).__name__
        receipt["exception_type"] = type(error).__name__
    finally:
        output_hashes = _hash_tree(boltz_out)
        receipt.setdefault("hashes", {})["processed_output_files_sha256"] = output_hashes
        receipt["hashes"]["msa_files_sha256"] = {
            name: digest for name, digest in output_hashes.items()
            if "msa" in name.lower() or Path(name).suffix.lower() in {".a3m", ".sto"}
        }
        receipt["hashes"]["processed_msa_files_sha256"] = {
            name: digest for name, digest in output_hashes.items()
            if "/processed/msa/" in "/" + name
        }
        if log.exists():
            receipt["hashes"]["inference_log_sha256"] = sha256_file(log)
        receipt["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        with receipt_path.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(receipt, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
    return receipt
