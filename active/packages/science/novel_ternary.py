"""Reference-free exploratory Boltz 2 ternary-complex planning and execution."""
from __future__ import annotations

import hashlib
import json
import math
import os
import subprocess
import time
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
import yaml
from rdkit import Chem, rdBase

from . import boltz_worker as worker
from .molecules import canonical, rows

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
PLAN_FORMAT = "tpd-novel-ternary-plan/1.0"
RECEIPT_FORMAT = "tpd-novel-ternary-receipt/1.0"
SUPPORTED_BOLTZ_VERSION = "2.2.1"
ROLE_KEYS = ("warhead_maps", "linker_maps", "recruiter_maps")
SOURCE_RULES = {
    "target": {"pdb": "6HAZ", "chain": "A", "length": 123,
               "terms": ("snf2l2", "smarca2", "probable global transcription activator snf2l2")},
    "CRBN": {"pdb": "6BOY", "chain": "B", "length": 463,
             "terms": ("cereblon", "crbn")},
    "VHL": {"pdb": "6HAY", "chain": "B", "length": 162,
            "terms": ("von hippel-lindau", "vhl")},
}


class NovelTernaryError(ValueError):
    """A frozen plan, source, prediction, or execution invariant failed."""


def _check(condition: bool, code: str) -> None:
    if not condition:
        raise NovelTernaryError(code)


def _json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                      allow_nan=False).encode("utf-8")


def _digest(value: Any) -> str:
    return hashlib.sha256(_json_bytes(value)).hexdigest()


def _sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sequence(value: Any) -> str:
    return "".join(str(value or "").split()).upper()


def _source_block(path: Path):
    path = Path(path)
    _check(path.is_file() and not path.is_symlink(), "SOURCE_REGULAR_FILE_REQUIRED")
    return gemmi.cif.read_file(str(path)).sole_block()


def _entity_tables(block):
    entities = {str(row["id"]): row for row in rows(block, "_entity.")}
    polymers = {str(row["entity_id"]): row for row in rows(block, "_entity_poly.")}
    asym = {str(row["id"]): str(row["entity_id"]) for row in rows(block, "_struct_asym.")}
    return entities, polymers, asym


def _entry_id(block) -> str | None:
    entries = rows(block, "_entry.")
    if not entries:
        return None
    value = str(entries[0].get("id", "")).strip().upper()
    return value or None


def _normalized_description(value: Any) -> str:
    return " ".join(str(value or "").lower().replace("–", "-").replace("—", "-")
                    .replace("-", " ").split())


def _validate_source(source: dict, expected: str) -> dict:
    _check(isinstance(source, dict), "SOURCE_OBJECT_REQUIRED")
    required = {"path", "label_asym_id", "expected_sha256", "role", "description_terms"}
    _check(required <= set(source), "SOURCE_FIELDS_MISSING")
    _check(source["role"] == ("target" if expected == "target" else "e3"), "SOURCE_ROLE_MISMATCH")
    path = Path(source["path"])
    actual_hash = worker.sha256_file(path)
    _check(source["expected_sha256"] == actual_hash, "SOURCE_HASH_MISMATCH")
    terms = source["description_terms"]
    _check(isinstance(terms, list) and terms and
           all(isinstance(value, str) and value.strip() for value in terms),
           "SOURCE_DESCRIPTION_TERMS_INVALID")

    rule = SOURCE_RULES[expected]
    _check(str(source["label_asym_id"]) == rule["chain"], "SOURCE_CHAIN_NOT_APPROVED")
    block = _source_block(path)
    deposit = _entry_id(block)
    if deposit is not None:
        _check(deposit == rule["pdb"], "SOURCE_DEPOSIT_ID_MISMATCH")
    entities, polymers, asym = _entity_tables(block)
    chain = str(source["label_asym_id"])
    _check(chain in asym and asym[chain] in polymers, "SOURCE_POLYMER_CHAIN_MISSING")
    entity_id = asym[chain]
    polymer = polymers[entity_id]
    _check(str(polymer.get("type", "")).lower() == "polypeptide(l)", "SOURCE_NOT_L_POLYPEPTIDE")
    sequence = _sequence(polymer.get("pdbx_seq_one_letter_code_can"))
    _check(len(sequence) == rule["length"], "SOURCE_CONSTRUCT_LENGTH_MISMATCH")
    _check(bool(sequence) and not set(sequence) - AMINO_ACIDS,
           "SOURCE_CANONICAL_SEQUENCE_UNSUPPORTED")
    description = str(entities.get(entity_id, {}).get("pdbx_description", ""))
    normalized = _normalized_description(description)
    _check(any(_normalized_description(term) in normalized for term in terms),
           "SOURCE_DECLARED_ROLE_TERMS_ABSENT")
    _check(any(_normalized_description(term) in normalized for term in rule["terms"]),
           "SOURCE_APPROVED_ROLE_NOT_VERIFIED")

    atoms = worker._eligible_atom_rows(block, chain)
    _check(bool(atoms), "SOURCE_CHAIN_HAS_NO_ATOMS")
    auth = sorted({str(row.get("auth_asym_id")) for row in atoms})
    _check(len(auth) == 1 and auth[0] not in {"", "None", "False", ".", "?"},
           "SOURCE_AUTH_CHAIN_AMBIGUOUS")

    excluded = []
    for other_chain, other_entity in sorted(asym.items()):
        if other_chain == chain or other_entity not in polymers:
            continue
        other_description = str(entities.get(other_entity, {}).get("pdbx_description", ""))
        low = _normalized_description(other_description)
        if "ddb1" in low or "dna damage binding protein 1" in low:
            category = "DDB1"
        elif "elongin c" in low:
            category = "elongin_C"
        elif "elongin b" in low:
            category = "elongin_B"
        else:
            category = "other_protein"
        excluded.append({"label_asym_id": other_chain, "entity_id": other_entity,
                         "description": other_description, "category": category})

    return {
        "path": str(path), "expected_sha256": actual_hash, "label_asym_id": chain,
        "auth_asym_id": auth[0], "entity_id": entity_id, "role": source["role"],
        "description": description, "description_terms": list(terms),
        "canonical_sequence": sequence, "sequence_length": len(sequence),
        "sequence_sha256": _sha_bytes(sequence.encode("ascii")),
        "sequence_source": "_entity_poly.pdbx_seq_one_letter_code_can",
        "expected_deposit": rule["pdb"], "observed_entry_id": deposit,
        "excluded_proteins": excluded,
    }


def _registered_ref(value: Any) -> bool:
    if not isinstance(value, dict):
        return False
    required = {"artifact_id", "version", "sha256", "media_type", "schema_id", "provenance"}
    return (required <= set(value) and
            isinstance(value.get("artifact_id"), str) and bool(value["artifact_id"]) and
            type(value.get("version")) is int and value["version"] > 0 and
            isinstance(value.get("sha256"), str) and len(value["sha256"]) == 64 and
            all(character in "0123456789abcdef" for character in value["sha256"]) and
            isinstance(value.get("media_type"), str) and bool(value["media_type"]) and
            isinstance(value.get("schema_id"), str) and bool(value["schema_id"]) and
            isinstance(value.get("provenance"), str) and bool(value["provenance"]))


def _qualification(candidate: dict) -> dict:
    _check(candidate.get("assembly_mode") == "pose_supported_hypothesis",
           "CANDIDATE_ASSEMBLY_MODE_NOT_POSE_SUPPORTED_HYPOTHESIS")
    analog_id = candidate.get("warhead_analog_id")
    _check(isinstance(analog_id, str) and analog_id.strip(), "WARHEAD_ANALOG_ID_REQUIRED")
    qualification = candidate.get("analog_qualification")
    _check(isinstance(qualification, dict), "ANALOG_QUALIFICATION_REQUIRED")
    expected = {"id", "selected", "pipeline_status", "assembly_eligible",
                "parent_redocking_supported", "docking"}
    _check(set(qualification) == expected, "ANALOG_QUALIFICATION_FIELDS_INVALID")
    _check(qualification.get("id") == analog_id, "ANALOG_QUALIFICATION_ID_MISMATCH")
    _check(qualification.get("selected") is True, "ANALOG_NOT_SELECTED")
    _check(qualification.get("pipeline_status") == "qualified", "ANALOG_PIPELINE_NOT_QUALIFIED")
    _check(qualification.get("assembly_eligible") is True, "ANALOG_NOT_ASSEMBLY_ELIGIBLE")
    _check(qualification.get("parent_redocking_supported") is True,
           "ANALOG_PARENT_REDOCKING_NOT_SUPPORTED")
    docking = qualification.get("docking")
    _check(isinstance(docking, dict), "ANALOG_DOCKING_REQUIRED")
    _check(docking.get("status") == "completed_with_limits", "ANALOG_DOCKING_STATUS_INVALID")
    _check(docking.get("pose_preserved") is True, "ANALOG_DOCKING_POSE_NOT_PRESERVED")
    passing = docking.get("passing_pose_count")
    _check(type(passing) is int and passing > 0, "ANALOG_DOCKING_NO_PASSING_POSE")
    files = docking.get("files")
    _check(isinstance(files, dict) and _registered_ref(files.get("poses_sdf")),
           "ANALOG_DOCKING_POSES_SDF_REF_INVALID")
    return json.loads(json.dumps(qualification, allow_nan=False))


def _candidate_smiles(candidate: dict) -> tuple[str, str]:
    mapped = candidate.get("actualmapped_smiles", candidate.get("actual_mapped_smiles",
                                                                  candidate.get("mapped_smiles")))
    plain = candidate.get("canonical_smiles")
    _check(isinstance(mapped, str) and mapped and isinstance(plain, str) and plain,
           "CANDIDATE_SMILES_REQUIRED")
    return mapped, plain


def _mapped_graph(candidate: dict) -> dict:
    qualification = _qualification(candidate)
    mapped, declared_plain = _candidate_smiles(candidate)
    mol = Chem.MolFromSmiles(mapped)
    _check(mol is not None and len(Chem.GetMolFrags(mol)) == 1, "CANDIDATE_MAPPED_GRAPH_INVALID")
    _check(all(atom.GetAtomicNum() > 0 for atom in mol.GetAtoms()), "CANDIDATE_DUMMY_ATOM_FORBIDDEN")
    maps = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
    _check(all(type(value) is int and value > 0 for value in maps) and len(maps) == len(set(maps)),
           "CANDIDATE_ATOM_MAPS_NOT_UNIQUE_COMPLETE")

    parsed_plain = canonical(mol)
    declared_mol = Chem.MolFromSmiles(declared_plain)
    _check(declared_mol is not None and len(Chem.GetMolFrags(declared_mol)) == 1,
           "CANDIDATE_CANONICAL_SMILES_INVALID")
    declared_canonical = canonical(declared_mol)
    _check(parsed_plain == declared_canonical, "CANDIDATE_CANONICAL_GRAPH_MISMATCH")
    identity = candidate.get("identity")
    if isinstance(identity, dict) and identity.get("canonical_isomeric_smiles"):
        identity_mol = Chem.MolFromSmiles(identity["canonical_isomeric_smiles"])
        _check(identity_mol is not None and parsed_plain == canonical(identity_mol),
               "CANDIDATE_IDENTITY_GRAPH_MISMATCH")

    atom_roles = candidate.get("atom_roles")
    _check(isinstance(atom_roles, dict) and set(atom_roles) == set(ROLE_KEYS), "ATOM_ROLES_INVALID")
    role_sets = {}
    for key in ROLE_KEYS:
        values = atom_roles[key]
        _check(isinstance(values, list) and values and
               all(type(value) is int and value > 0 for value in values), "ATOM_ROLE_MAPS_INVALID")
        _check(len(values) == len(set(values)), "ATOM_ROLE_MAP_DUPLICATE")
        role_sets[key] = set(values)
    _check(all(not role_sets[a] & role_sets[b] for a, b in
               ((ROLE_KEYS[0], ROLE_KEYS[1]), (ROLE_KEYS[0], ROLE_KEYS[2]),
                (ROLE_KEYS[1], ROLE_KEYS[2]))), "ATOM_ROLES_NOT_DISJOINT")
    _check(set().union(*role_sets.values()) == set(maps), "ATOM_ROLES_NOT_EXACT_GRAPH_PARTITION")

    attachment = candidate.get("attachment_metadata")
    _check(isinstance(attachment, dict), "ATTACHMENT_METADATA_REQUIRED")
    endpoints = {}
    for key in ("warhead_linker_bond", "recruiter_linker_bond"):
        value = attachment.get(key)
        _check(isinstance(value, dict) and value.get("bond_type") == "SINGLE",
               "ATTACHMENT_BOND_INVALID")
        a, b = value.get("attachment_atom_map"), value.get("partner_atom_map")
        _check(type(a) is int and type(b) is int, "ATTACHMENT_MAP_INVALID")
        atom_a = next((atom for atom in mol.GetAtoms() if atom.GetAtomMapNum() == a), None)
        atom_b = next((atom for atom in mol.GetAtoms() if atom.GetAtomMapNum() == b), None)
        _check(atom_a is not None and atom_b is not None, "ATTACHMENT_MAP_NOT_IN_GRAPH")
        actual_bond = mol.GetBondBetweenAtoms(atom_a.GetIdx(), atom_b.GetIdx())
        _check(actual_bond is not None, "ATTACHMENT_BOND_ABSENT_FROM_GRAPH")
        _check(actual_bond.GetBondType() == Chem.BondType.SINGLE,
               "ATTACHMENT_BOND_TYPE_GRAPH_MISMATCH")
        endpoints[key] = (a, b)
    _check(endpoints["warhead_linker_bond"][0] in role_sets["warhead_maps"] and
           endpoints["warhead_linker_bond"][1] in role_sets["linker_maps"] and
           endpoints["recruiter_linker_bond"][0] in role_sets["recruiter_maps"] and
           endpoints["recruiter_linker_bond"][1] in role_sets["linker_maps"],
           "ATTACHMENT_ROLE_MISMATCH")

    with_h = Chem.AddHs(Chem.MolFromSmiles(mapped))
    ranks = list(Chem.CanonicalRankAtoms(with_h))
    Chem.AssignStereochemistry(with_h, cleanIt=True, force=True)
    predicted_names, elements = {}, {}
    for atom in with_h.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        atom_map = atom.GetAtomMapNum()
        predicted_names[str(atom_map)] = atom.GetSymbol().upper() + str(ranks[atom.GetIdx()] + 1)
        elements[str(atom_map)] = atom.GetSymbol().upper()
    _check(len(set(predicted_names.values())) == len(predicted_names), "BOLTZ_ATOM_NAMES_NOT_UNIQUE")

    order_names = {Chem.BondType.SINGLE: "SING", Chem.BondType.DOUBLE: "DOUB",
                   Chem.BondType.TRIPLE: "TRIP", Chem.BondType.AROMATIC: "AROM"}
    bonds = []
    for bond in mol.GetBonds():
        a = mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetAtomMapNum()
        b = mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetAtomMapNum()
        bonds.append({"maps": sorted([a, b]), "order": order_names.get(bond.GetBondType())})
    _check(all(item["order"] for item in bonds), "CANDIDATE_BOND_TYPE_UNSUPPORTED")
    bonds.sort(key=lambda item: (item["maps"], item["order"]))

    normalized = {
        "candidate_id": candidate.get("candidate_id"), "e3_type": candidate.get("e3_type"),
        "assembly_mode": candidate.get("assembly_mode"),
        "canonical_full_graph": {
            "derived": True,
            "method": "RDKit canonical graph equality of mapped_smiles and canonical_smiles",
            "mapped_graph_canonical_smiles": parsed_plain,
            "declared_graph_canonical_smiles": declared_canonical,
        },
        "warhead_analog_id": qualification["id"],
        "analog_qualification": qualification, "actualmapped_smiles": mapped,
        "canonical_smiles": declared_plain, "parsed_canonical_smiles": parsed_plain,
        "atom_roles": {key: list(atom_roles[key]) for key in ROLE_KEYS},
        "attachment_metadata": attachment, "atom_map_to_predicted_atom_name": predicted_names,
        "atom_map_to_element": elements, "heavy_bonds": bonds,
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "stereochemistry": {"source": "input_mapped_smiles", "coordinate_inference": False},
    }
    normalized["graph_sha256"] = _digest({
        "mapped_smiles": mapped, "parsed_canonical_smiles": parsed_plain,
        "atom_roles": normalized["atom_roles"], "attachment_metadata": attachment,
        "stereochemistry_source": "input_mapped_smiles",
    })
    return normalized


def _validate_binding(binding: dict, kind: str) -> dict:
    _check(isinstance(binding, dict), f"{kind.upper()}_BINDING_REQUIRED")
    required = ({"project", "job_id", "input_sha256", "result_sha256", "candidate_graph_sha256"}
                if kind == "job" else {"module", "revision", "digest"})
    _check(set(binding) == required, f"{kind.upper()}_BINDING_FIELDS_INVALID")
    _check(all(isinstance(binding[key], str) and binding[key] for key in required),
           f"{kind.upper()}_BINDING_VALUE_INVALID")
    for key in required & {"input_sha256", "result_sha256", "candidate_graph_sha256", "digest"}:
        value = binding[key]
        _check(len(value) == 64 and all(character in "0123456789abcdef" for character in value),
               f"{kind.upper()}_BINDING_DIGEST_INVALID")
    return dict(binding)


def _document(target: str, e3: str, smiles: str, msa_mode: str) -> dict:
    extra = {"msa": "empty"} if msa_mode == "single_sequence" else {}
    return {"version": 1, "sequences": [
        {"protein": {"id": "T", "sequence": target, **extra}},
        {"protein": {"id": "E", "sequence": e3, **extra}},
        {"ligand": {"id": "L", "smiles": smiles}},
    ]}


def prepare_plan(candidate: dict, *, target_source: dict, e3_source: dict,
                 job_binding: dict, policy_binding: dict, seeds: list[int],
                 msa_mode: str = "single_sequence") -> dict:
    """Validate and freeze a coordinate-free novel ternary prediction plan."""
    _check(msa_mode in {"single_sequence", "server"}, "MSA_MODE_INVALID")
    _check(isinstance(candidate, dict), "CANDIDATE_REQUIRED")
    graph = _mapped_graph(candidate)
    _check(graph["e3_type"] in {"CRBN", "VHL"}, "E3_TYPE_INVALID")
    _check(isinstance(graph["candidate_id"], str) and graph["candidate_id"], "CANDIDATE_ID_REQUIRED")
    _check(isinstance(seeds, list) and 0 < len(seeds) <= 32 and
           all(type(seed) is int and 0 <= seed < 2**32 for seed in seeds), "SEEDS_INVALID")
    _check(len(seeds) == len(set(seeds)), "SEEDS_NOT_UNIQUE")

    target = _validate_source(target_source, "target")
    e3 = _validate_source(e3_source, graph["e3_type"])
    job = _validate_binding(job_binding, "job")
    policy = _validate_binding(policy_binding, "policy")
    _check(job["candidate_graph_sha256"] == graph["graph_sha256"],
           "JOB_CANDIDATE_GRAPH_BINDING_MISMATCH")
    _check(job["result_sha256"] != job["candidate_graph_sha256"],
           "JOB_RESULT_AND_GRAPH_DIGEST_NOT_DISTINCT")

    document = _document(target["canonical_sequence"], e3["canonical_sequence"],
                         graph["actualmapped_smiles"], msa_mode)
    yaml_text = yaml.safe_dump(document, sort_keys=False, allow_unicode=False)
    plan = {
        "format": PLAN_FORMAT,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "candidate_graph": graph, "sources": {"target": target, "e3": e3},
        "bindings": {"job": job, "policy": policy}, "seeds": list(seeds), "msa_mode": msa_mode,
        "boltz_input": {"yaml": yaml_text, "sha256": _sha_bytes(yaml_text.encode("utf-8")),
                        "chain_roles": {"target": "T", "e3": "E", "ligand": "L"},
                        "templates": [], "constraints": [], "potentials": False,
                        "ligand_input": "mapped SMILES; no SDF or input 3D coordinates"},
        "rdkit": {"planning_version": rdBase.rdkitVersion,
                  "required_boltz_environment_version": rdBase.rdkitVersion,
                  "atom_naming_algorithm": "Boltz 2.2.1 mapped-SMILES atom naming"},
        "run_authorization": {"classification": "exploratory_not_expert",
                              "authorized_accelerators": ["cpu", "gpu"],
                              "gpu_authorization_source": "operator CLI/API invocation only"},
        "scientific_scope": {"status_vocabulary": "computed_hypothesis",
                             "experimental_ternary_reference": None,
                             "reference_rmsd_permitted": False, "reference_free": True,
                             "affinity_model": False,
                             "cofactors_excluded": ["DDB1", "elongin C", "elongin B",
                                                    "all other source proteins"],
                             "statement": "Simplified reference-free exploratory ternary model."},
    }
    plan["plan_digest"] = _digest(plan)
    return plan


def verify_plan(plan: dict) -> dict:
    _check(isinstance(plan, dict) and plan.get("format") == PLAN_FORMAT, "PLAN_FORMAT_INVALID")
    supplied = plan.get("plan_digest")
    unsigned = dict(plan)
    unsigned.pop("plan_digest", None)
    _check(isinstance(supplied, str) and supplied == _digest(unsigned), "PLAN_DIGEST_INVALID")
    graph = _mapped_graph(plan.get("candidate_graph", {}))
    _check(graph == plan.get("candidate_graph"), "PLAN_CANDIDATE_GRAPH_CORRUPTED")
    target = _validate_source(plan["sources"]["target"], "target")
    e3 = _validate_source(plan["sources"]["e3"], graph["e3_type"])
    _check(target["canonical_sequence"] == plan["sources"]["target"].get("canonical_sequence"),
           "PLAN_TARGET_SEQUENCE_CORRUPTED")
    _check(e3["canonical_sequence"] == plan["sources"]["e3"].get("canonical_sequence"),
           "PLAN_E3_SEQUENCE_CORRUPTED")
    job = _validate_binding(plan["bindings"]["job"], "job")
    _validate_binding(plan["bindings"]["policy"], "policy")
    _check(job["candidate_graph_sha256"] == graph["graph_sha256"],
           "PLAN_JOB_GRAPH_BINDING_INVALID")
    _check(job["result_sha256"] != job["candidate_graph_sha256"],
           "PLAN_JOB_RESULT_GRAPH_DIGEST_COLLISION")
    _check(plan.get("msa_mode") in {"single_sequence", "server"}, "PLAN_MSA_MODE_INVALID")
    seeds = plan.get("seeds")
    _check(isinstance(seeds, list) and 0 < len(seeds) <= 32 and len(seeds) == len(set(seeds)) and
           all(type(seed) is int and 0 <= seed < 2**32 for seed in seeds), "PLAN_SEEDS_INVALID")
    expected_yaml = yaml.safe_dump(_document(target["canonical_sequence"], e3["canonical_sequence"],
                                              graph["actualmapped_smiles"], plan["msa_mode"]),
                                   sort_keys=False, allow_unicode=False)
    boltz = plan.get("boltz_input", {})
    _check(boltz.get("yaml") == expected_yaml and
           boltz.get("sha256") == _sha_bytes(expected_yaml.encode("utf-8")),
           "PLAN_BOLTZ_INPUT_CORRUPTED")
    _check(boltz.get("templates") == [] and boltz.get("constraints") == [] and
           boltz.get("potentials") is False, "PLAN_UNSAFE_INPUT_FEATURE")
    _check(plan.get("scientific_scope", {}).get("reference_free") is True and
           plan["scientific_scope"].get("experimental_ternary_reference") is None,
           "PLAN_REFERENCE_FREE_SCOPE_INVALID")
    return plan


def write_plan(plan: dict, destination: Path) -> dict:
    verify_plan(plan)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(plan, ensure_ascii=False, indent=2, sort_keys=True,
                         allow_nan=False).encode("utf-8") + b"\n"
    try:
        with destination.open("xb") as stream:
            stream.write(encoded)
    except FileExistsError as error:
        raise NovelTernaryError("PLAN_DESTINATION_EXISTS") from error
    return {"path": str(destination), "sha256": _sha_bytes(encoded),
            "plan_digest": plan["plan_digest"]}


def _probe_boltz_environment(executable: Path) -> dict:
    python = worker._wrapper_python(executable)
    source = ("import importlib.metadata,json; import rdkit; import torch;"
              "print(json.dumps({'boltz_version':importlib.metadata.version('boltz'),"
              "'rdkit_version':rdkit.__version__,'torch_version':torch.__version__,"
              "'cuda_available':torch.cuda.is_available(),'cuda_runtime':torch.version.cuda,"
              "'cuda_device_count':torch.cuda.device_count(),"
              "'cuda_devices':[torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]}))")
    flags = subprocess.CREATE_NO_WINDOW if hasattr(subprocess, "CREATE_NO_WINDOW") else 0
    try:
        result = subprocess.run([str(python), "-c", source], capture_output=True, text=True,
                                encoding="utf-8", errors="replace", timeout=30, shell=False,
                                check=False, creationflags=flags)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise NovelTernaryError("BOLTZ_ENVIRONMENT_PROBE_FAILED") from error
    _check(result.returncode == 0, "BOLTZ_ENVIRONMENT_PROBE_FAILED")
    try:
        value = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        raise NovelTernaryError("BOLTZ_ENVIRONMENT_PROBE_INVALID") from error
    _check(value.get("boltz_version") == SUPPORTED_BOLTZ_VERSION, "BOLTZ_VERSION_NOT_2_2_1")
    value["python_path"] = str(python)
    return value


def _prediction_chain_mapping(block, packet: dict) -> dict[str, str]:
    sequences = worker._prediction_sequences(block)
    result = {}
    for role in ("target", "e3"):
        expected = packet["sources"][role]["canonical_sequence"]
        matches = [chain for chain, sequence in sequences.items() if sequence == expected]
        _check(len(matches) == 1, "PREDICTION_SEQUENCE_ROLE_MAPPING_AMBIGUOUS")
        result[role] = matches[0]
    _check(result["target"] != result["e3"], "PREDICTION_PROTEIN_CHAIN_REUSED")
    return result


def _ligand_atoms(block, graph: dict) -> tuple[str, str, dict[str, dict]]:
    expected = set(graph["atom_map_to_predicted_atom_name"].values())
    expected_elements = {graph["atom_map_to_predicted_atom_name"][key]: value
                         for key, value in graph["atom_map_to_element"].items()}
    _, polymers, asym = _entity_tables(block)
    polymer_chains = {chain for chain, entity in asym.items() if entity in polymers}
    atom_chains = sorted({str(row.get("label_asym_id")) for row in rows(block, "_atom_site.")})
    candidates = []
    for chain in atom_chains:
        if chain in polymer_chains:
            continue
        atoms = worker._eligible_atom_rows(block, chain)
        named = {str(row.get("label_atom_id")): row for row in atoms}
        comp_ids = {str(row.get("label_comp_id")) for row in atoms}
        if len(named) != len(atoms) or set(named) != expected or len(comp_ids) != 1:
            continue
        if any(str(row.get("type_symbol", "")).upper() != expected_elements[name]
               for name, row in named.items()):
            continue
        candidates.append((chain, next(iter(comp_ids)), named))
    _check(len(candidates) == 1, "PREDICTION_LIGAND_NAME_ELEMENT_MAPPING_FAILED")
    return candidates[0]


def _bond_semantics(mol: Chem.Mol) -> dict[tuple[int, int], str]:
    result = {}
    for bond in mol.GetBonds():
        edge = tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
        result[edge] = ("AROM" if bond.GetIsAromatic() else
                        {Chem.BondType.SINGLE: "SING", Chem.BondType.DOUBLE: "DOUB",
                         Chem.BondType.TRIPLE: "TRIP"}.get(bond.GetBondType(), "UNSUPPORTED"))
    return result


def _validate_chem_comp_bonds(block, graph: dict, comp_id: str) -> str:
    bond_rows = rows(block, "_chem_comp_bond.")
    if not bond_rows:
        return "not_available"
    expected_names = set(graph["atom_map_to_predicted_atom_name"].values())
    relevant = [row for row in bond_rows
                if str(row.get("comp_id")) == comp_id and
                str(row.get("atom_id_1")) in expected_names and
                str(row.get("atom_id_2")) in expected_names]
    if not relevant:
        return "not_available"
    name_by_map = graph["atom_map_to_predicted_atom_name"]
    expected_edges = {tuple(sorted((name_by_map[str(item["maps"][0])],
                                    name_by_map[str(item["maps"][1])]))): item["order"]
                      for item in graph["heavy_bonds"]}
    aliases = {"SINGLE": "SING", "DOUBLE": "DOUB", "TRIPLE": "TRIP", "AROMATIC": "AROM"}
    observed = {}
    for row in relevant:
        edge = tuple(sorted((str(row["atom_id_1"]), str(row["atom_id_2"]))))
        order = aliases.get(str(row.get("value_order", "")).upper(),
                            str(row.get("value_order", "")).upper())
        _check(edge not in observed, "PREDICTION_CHEM_COMP_DUPLICATE_BOND")
        _check(order in {"SING", "DOUB", "TRIP", "AROM"},
               "PREDICTION_CHEM_COMP_BOND_ORDER_UNSUPPORTED")
        observed[edge] = order
    _check(set(observed) == set(expected_edges), "PREDICTION_CHEM_COMP_BOND_GRAPH_MISMATCH")

    source = Chem.MolFromSmiles(graph["actualmapped_smiles"])
    _check(source is not None, "PREDICTION_EXPECTED_GRAPH_INVALID")
    source_by_name = {name_by_map[str(atom.GetAtomMapNum())]: atom for atom in source.GetAtoms()}
    ordered_names = sorted(expected_names)
    indices = {name: index for index, name in enumerate(ordered_names)}
    editable = Chem.RWMol()
    for name in ordered_names:
        original = source_by_name[name]
        atom = Chem.Atom(original.GetAtomicNum())
        atom.SetFormalCharge(original.GetFormalCharge())
        atom.SetIsotope(original.GetIsotope())
        atom.SetNoImplicit(original.GetNoImplicit())
        atom.SetNumExplicitHs(original.GetNumExplicitHs())
        editable.AddAtom(atom)
    types = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
             "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for edge, order in observed.items():
        editable.AddBond(indices[edge[0]], indices[edge[1]], types[order])
        if order == "AROM":
            bond = editable.GetBondBetweenAtoms(indices[edge[0]], indices[edge[1]])
            bond.SetIsAromatic(True)
            editable.GetAtomWithIdx(indices[edge[0]]).SetIsAromatic(True)
            editable.GetAtomWithIdx(indices[edge[1]]).SetIsAromatic(True)
    observed_mol = editable.GetMol()
    try:
        Chem.SanitizeMol(observed_mol)
    except Exception as error:
        raise NovelTernaryError("PREDICTION_CHEM_COMP_GRAPH_NOT_SANITIZABLE") from error

    expected_mol = Chem.RWMol()
    for name in ordered_names:
        original = source_by_name[name]
        atom = Chem.Atom(original.GetAtomicNum())
        atom.SetFormalCharge(original.GetFormalCharge())
        atom.SetIsotope(original.GetIsotope())
        atom.SetNoImplicit(original.GetNoImplicit())
        atom.SetNumExplicitHs(original.GetNumExplicitHs())
        expected_mol.AddAtom(atom)
    for edge, order in expected_edges.items():
        expected_mol.AddBond(indices[edge[0]], indices[edge[1]], types[order])
        if order == "AROM":
            bond = expected_mol.GetBondBetweenAtoms(indices[edge[0]], indices[edge[1]])
            bond.SetIsAromatic(True)
            expected_mol.GetAtomWithIdx(indices[edge[0]]).SetIsAromatic(True)
            expected_mol.GetAtomWithIdx(indices[edge[1]]).SetIsAromatic(True)
    expected_mol = expected_mol.GetMol()
    try:
        Chem.SanitizeMol(expected_mol)
    except Exception as error:
        raise NovelTernaryError("PREDICTION_EXPECTED_GRAPH_INVALID") from error
    _check(_bond_semantics(observed_mol) == _bond_semantics(expected_mol),
           "PREDICTION_CHEM_COMP_BOND_ORDER_MISMATCH")
    return "semantic_match"


def _contacts(proteins: dict, ligand: dict, names_by_role: dict[str, set[str]], cutoff=4.0) -> dict:
    result = {role: {group: 0 for group in names_by_role} for role in ("target", "e3")}
    seen = set()
    for (role, residue, _), atom in proteins.items():
        point = np.asarray(atom["xyz"], float)
        for group, names in names_by_role.items():
            for name in names:
                if np.linalg.norm(point - np.asarray(ligand[name]["xyz"], float)) <= cutoff:
                    seen.add((role, group, residue, name))
    for role, group, _, _ in seen:
        result[role][group] += 1
    return result


def _interface_clashes(proteins: dict, limit: int = 256) -> dict:
    target = [((residue, atom), row) for (role, residue, atom), row in proteins.items()
              if role == "target"]
    e3 = [((residue, atom), row) for (role, residue, atom), row in proteins.items()
          if role == "e3"]
    periodic = Chem.GetPeriodicTable()

    def radius(row: dict) -> float:
        symbol = str(row.get("element", row.get("type_symbol", "C"))).strip().title()
        try:
            atomic_number = periodic.GetAtomicNumber(symbol)
            return float(periodic.GetRvdw(atomic_number)) if atomic_number else 1.7
        except Exception:
            return 1.7

    target_xyz = np.asarray([row["xyz"] for _, row in target], dtype=float)
    e3_xyz = np.asarray([row["xyz"] for _, row in e3], dtype=float)
    target_radii = np.asarray([radius(row) for _, row in target], dtype=float)
    e3_radii = np.asarray([radius(row) for _, row in e3], dtype=float)
    count, examples = 0, []
    chunk_size = 512
    for start in range(0, len(target), chunk_size):
        stop = min(start + chunk_size, len(target))
        delta = target_xyz[start:stop, None, :] - e3_xyz[None, :, :]
        distances = np.sqrt(np.einsum("ijk,ijk->ij", delta, delta))
        thresholds = 0.75 * (target_radii[start:stop, None] + e3_radii[None, :])
        hits = np.argwhere(distances < thresholds)
        count += int(hits.shape[0])
        remaining = limit - len(examples)
        if remaining > 0:
            for local_target, e3_index in hits[:remaining]:
                target_index = start + int(local_target)
                examples.append({
                    "target": list(target[target_index][0]), "e3": list(e3[int(e3_index)][0]),
                    "distance_A": float(distances[local_target, e3_index]),
                    "threshold_A": float(thresholds[local_target, e3_index]),
                })
    return {"count": count, "examples": examples, "examples_truncated": count > len(examples),
            "example_limit": limit}


def assess_novel_outputs(out_dir: Path, packet: dict) -> dict:
    """Inspect a prediction without comparing it to deposited coordinates."""
    found = worker.discover_outputs(Path(out_dir))
    if found["status"] != "located":
        return found
    try:
        confidence = worker._validate_confidence(Path(found["confidence"]))
        block = _source_block(Path(found["prediction"]))
        chains = _prediction_chain_mapping(block, packet)
        graph = packet["candidate_graph"]
        ligand_chain, ligand_comp_id, ligand = _ligand_atoms(block, graph)
        chem_graph = _validate_chem_comp_bonds(block, graph, ligand_comp_id)
        proteins = {}
        for role in ("target", "e3"):
            proteins.update(worker._protein_atoms(block, chains[role], role))
        groups = {key.removesuffix("_maps"):
                  {graph["atom_map_to_predicted_atom_name"][str(value)]
                   for value in graph["atom_roles"][key]} for key in ROLE_KEYS}
        attachment = graph["attachment_metadata"]
        endpoint_maps = [attachment["warhead_linker_bond"]["partner_atom_map"],
                         attachment["recruiter_linker_bond"]["partner_atom_map"]]
        endpoint_names = [graph["atom_map_to_predicted_atom_name"][str(value)]
                          for value in endpoint_maps]
        endpoint_distance = float(np.linalg.norm(np.asarray(ligand[endpoint_names[0]]["xyz"], float) -
                                                 np.asarray(ligand[endpoint_names[1]]["xyz"], float)))
        ca = {}
        for role in ("target", "e3"):
            ca[role] = np.asarray([row["xyz"] for (item_role, _, atom), row in proteins.items()
                                   if item_role == role and atom == "CA"], float)
            _check(len(ca[role]) > 0, "PREDICTION_CA_ATOMS_MISSING")
        distances = np.linalg.norm(ca["target"][:, None, :] - ca["e3"][None, :, :], axis=2)
        result = {
            "status": "computed_hypothesis", "reference_free": True,
            "actual_computation": True, "prediction": found["prediction"],
            "confidence_file": found["confidence"], "model_confidence": confidence,
            "mapping": {"protein_chains": chains, "ligand_chain": ligand_chain,
                        "ligand_comp_id": ligand_comp_id,
                        "protein_identity": "exact canonical sequence from _entity_poly_seq",
                        "ligand_identity": "exact Boltz atom names and elements",
                        "chem_comp_bond_graph": chem_graph,
                        "stereochemistry_source": "input_mapped_smiles_only",
                        "near_position_never_used_for_identity": True},
            "descriptive_metrics": {
                "contact_definition": "protein-heavy-atom/ligand-heavy-atom distance <= 4.0 A",
                "contacts_by_protein_role_and_ligand_group": _contacts(proteins, ligand, groups),
                "linker_endpoint_atom_names": endpoint_names,
                "linker_endpoint_distance_A": endpoint_distance,
                "protein_ligand_clashes": worker._clashes(proteins, ligand),
                "target_e3_heavy_atom_clashes": _interface_clashes(proteins),
                "clash_definition": "distance < 0.75 times van der Waals radii sum; descriptive",
                "target_e3_interface": {
                    "minimum_target_e3_CA_separation_A": float(np.min(distances)),
                    "target_e3_CA_pairs_within_8A": int(np.sum(distances <= 8.0)),
                    "CA_centroid_separation_A": float(np.linalg.norm(
                        ca["target"].mean(axis=0) - ca["e3"].mean(axis=0))),
                }},
            "interpretation": "Reference-free exploratory prediction; not evidence of binding, efficacy, degradation, approval, or ubiquitination.",
        }
        _json_bytes(result)
        return result
    except (NovelTernaryError, worker.BoltzWorkerError, OSError, ValueError, RuntimeError) as error:
        return {"status": "inspection_failed", "reference_free": True,
                "actual_computation": False, "reason": str(error) or type(error).__name__,
                "prediction": found.get("prediction"), "confidence": found.get("confidence")}


def _cancelled(cancel: Any) -> bool:
    return bool(cancel() if callable(cancel) else cancel is not None and cancel.is_set())


def _atomic_progress(path: Path, item: dict) -> None:
    previous = path.read_bytes() if path.exists() else b""
    temporary = path.with_suffix(".tmp")
    temporary.write_bytes(previous + _json_bytes(item) + b"\n")
    os.replace(temporary, path)


def _pairwise_consistency(valid: list[dict], packet: dict) -> list[dict]:
    loaded = []
    for receipt in valid:
        block = _source_block(Path(receipt["inspection"]["prediction"]))
        chains = _prediction_chain_mapping(block, packet)
        _, _, ligand = _ligand_atoms(block, packet["candidate_graph"])
        loaded.append((receipt["seed"], worker._protein_atoms(block, chains["target"], "target"),
                       worker._protein_atoms(block, chains["e3"], "e3"), ligand))
    names = sorted(packet["candidate_graph"]["atom_map_to_predicted_atom_name"].values())
    results = []
    for index, first in enumerate(loaded):
        for second in loaded[index + 1:]:
            common = [key for key in sorted(set(first[1]) & set(second[1]),
                                             key=lambda key: (key[1], key[2])) if key[2] == "CA"]
            if len(common) < 3:
                continue
            fit = worker.aligned_rmsd(np.asarray([first[1][key]["xyz"] for key in common]),
                                      np.asarray([second[1][key]["xyz"] for key in common]))
            rotation = np.asarray(fit["rotation_row_vectors"])
            translation = np.asarray(fit["translation_A"])
            ligand_a = np.asarray([first[3][name]["xyz"] for name in names])
            ligand_b = np.asarray([second[3][name]["xyz"] for name in names]) @ rotation + translation
            common_e3 = [key for key in sorted(set(first[2]) & set(second[2]),
                                                key=lambda key: (key[1], key[2])) if key[2] == "CA"]
            e3_rmsd = None
            if common_e3:
                a = np.asarray([first[2][key]["xyz"] for key in common_e3])
                b = np.asarray([second[2][key]["xyz"] for key in common_e3]) @ rotation + translation
                e3_rmsd = float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))
            results.append({"seed_a": first[0], "seed_b": second[0],
                            "target_CA_alignment_RMSD_A": fit["rmsd_A"],
                            "ligand_RMSD_after_target_CA_alignment_A": float(np.sqrt(
                                np.mean(np.sum((ligand_a - ligand_b) ** 2, axis=1)))),
                            "e3_CA_RMSD_after_target_CA_alignment_A": e3_rmsd,
                            "interpretation": "prediction consistency, not reference accuracy"})
    return results


def _safe_hash(path: Path, tree: bool = False):
    try:
        return worker._hash_tree(path) if tree else worker.sha256_file(path)
    except Exception as error:
        return {"hash_error": str(error) or type(error).__name__}


def run_plan(plan: dict, *, output_root: Path, executable: Path, checkpoint: Path,
             cache: Path, timeout: float = 1800, cancel: Any = None,
             accelerator: str = "gpu", max_msa_seqs: int = 256) -> dict:
    """Run all frozen seeds serially in a newly reserved output directory."""
    verify_plan(plan)
    _check(accelerator in {"cpu", "gpu"}, "ACCELERATOR_INVALID")
    _check(type(max_msa_seqs) is int and max_msa_seqs > 0, "MAX_MSA_SEQS_INVALID")
    _check(isinstance(timeout, (int, float)) and timeout > 0 and math.isfinite(timeout),
           "TIMEOUT_INVALID")
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=False)
    progress, aggregate_path = output_root / "run-progress.jsonl", output_root / "receipt.json"
    receipts, run_environment = [], None
    help_info = checkpoint_hash = executable_hash = None
    setup_error = None
    try:
        executable, checkpoint, cache = Path(executable), Path(checkpoint), Path(cache)
        _check(executable.is_file() and not executable.is_symlink(), "BOLTZ_EXECUTABLE_REQUIRED")
        _check(checkpoint.is_file() and not checkpoint.is_symlink(), "CHECKPOINT_REGULAR_FILE_REQUIRED")
        executable_hash, checkpoint_hash = worker.sha256_file(executable), worker.sha256_file(checkpoint)
        cache.mkdir(parents=True, exist_ok=True)
        _check(cache.is_dir() and not cache.is_symlink(), "CACHE_DIRECTORY_REQUIRED")
        run_environment = _probe_boltz_environment(executable)
        _check(run_environment.get("rdkit_version") == plan["rdkit"]["required_boltz_environment_version"],
               "BOLTZ_RDKIT_VERSION_MISMATCH")
        help_info = worker.inspect_help(executable)
    except Exception as error:
        setup_error = str(error) or type(error).__name__

    binding_summary = {"project": plan["bindings"]["job"]["project"],
                       "job_id": plan["bindings"]["job"]["job_id"],
                       "policy_digest": plan["bindings"]["policy"]["digest"]}
    for seed in plan["seeds"]:
        seed_dir = output_root / f"seed-{seed}"
        seed_dir.mkdir(exist_ok=False)
        receipt_path, input_path = seed_dir / "receipt.json", seed_dir / "novel_ternary.yaml"
        boltz_out, log = seed_dir / "boltz_output", seed_dir / "inference.log"
        receipt = {
            "format": RECEIPT_FORMAT, "seed": seed, "e3_type": plan["candidate_graph"]["e3_type"],
            "candidate_id": plan["candidate_graph"]["candidate_id"], "status": "failed",
            "scientific_status": "not_computed", "actual_computation": False,
            "reference_free": True, "bindings": binding_summary,
            "plan_digest": plan["plan_digest"],
            "seed_application": {"method": None, "actually_applied": False},
            "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "artifacts": {"seed_directory": str(seed_dir), "receipt": str(receipt_path),
                          "input_yaml": str(input_path), "output_directory": str(boltz_out),
                          "log": str(log)},
        }
        try:
            if setup_error:
                raise NovelTernaryError(setup_error)
            if _cancelled(cancel):
                receipt.update({"status": "cancelled", "process_state": "cancelled_before_start"})
                raise NovelTernaryError("CANCELLED_BEFORE_SEED_START")
            input_path.write_text(plan["boltz_input"]["yaml"], encoding="utf-8", newline="\n")
            _check(worker.sha256_file(input_path) == plan["boltz_input"]["sha256"],
                   "WRITTEN_INPUT_HASH_MISMATCH")
            command = worker.build_command(executable=executable, input_yaml=input_path,
                                           out_dir=boltz_out, cache=cache, checkpoint=checkpoint,
                                           accelerator=accelerator, seed=seed, msa_mode=plan["msa_mode"],
                                           max_msa_seqs=max_msa_seqs,
                                           help_has_seed=help_info["has_seed"])
            method = "official_cli_--seed"
            wrapper_hash = None
            if not help_info["has_seed"]:
                wrapper = seed_dir / "seeded_boltz_entrypoint.py"
                worker._write_seed_wrapper(wrapper)
                command = [str(worker._wrapper_python(executable)), str(wrapper), str(seed), *command[1:]]
                method, wrapper_hash = "same-environment seed wrapper", worker.sha256_file(wrapper)
            receipt["seed_application"] = {"method": method, "actually_applied": False}
            if wrapper_hash:
                receipt["seed_application"]["wrapper_sha256"] = wrapper_hash
            exit_code, state, elapsed = worker._run_process(command, log, timeout, cancel)
            receipt["seed_application"]["actually_applied"] = True
            inspection = assess_novel_outputs(boltz_out, plan)
            receipt.update({"exit_code": exit_code, "process_state": state,
                            "elapsed_seconds": elapsed, "inspection": inspection,
                            "settings": {"accelerator": accelerator, "recycling_steps": 3,
                                         "sampling_steps": 200, "diffusion_samples": 1,
                                         "serial": True, "no_kernels": True,
                                         "use_potentials": False, "affinity_model": False,
                                         "msa_mode": plan["msa_mode"],
                                         "max_msa_seqs": max_msa_seqs}})
            if state == "cancelled":
                receipt["status"] = "cancelled"
                receipt["failure_reason"] = "CANCELLED"
            elif state == "timeout":
                receipt["status"] = "timeout"
                receipt["failure_reason"] = "TIMEOUT"
            elif exit_code == 0 and state == "completed" and inspection.get("status") == "computed_hypothesis":
                receipt.update({"status": "completed", "scientific_status": "computed_hypothesis",
                                "actual_computation": True})
            else:
                text = log.read_text(encoding="utf-8", errors="replace").lower() if log.exists() else ""
                receipt["failure_reason"] = ("OOM" if "out of memory" in text or "cuda oom" in text
                                             else inspection.get("reason") or state or "BOLTZ_NONZERO_EXIT")
        except Exception as error:
            receipt.setdefault("failure_reason", str(error) or type(error).__name__)
            receipt["exception_type"] = type(error).__name__
        finally:
            receipt["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            try:
                receipt["hashes"] = {
                    "checkpoint_sha256": checkpoint_hash, "executable_sha256": executable_hash,
                    "input_yaml_sha256": _safe_hash(input_path) if input_path.exists() else None,
                    "inference_log_sha256": _safe_hash(log) if log.exists() else None,
                    "output_files_sha256": _safe_hash(boltz_out, tree=True),
                }
            except Exception as error:
                receipt["hashing_failure"] = str(error) or type(error).__name__
            receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True,
                                               allow_nan=False) + "\n", encoding="utf-8", newline="\n")
            receipts.append(receipt)
            _atomic_progress(progress, {"seed": seed, "status": receipt["status"],
                                        "receipt": str(receipt_path),
                                        "finished_utc": receipt["finished_utc"]})

    valid = [receipt for receipt in receipts if receipt["status"] == "completed"]
    consistency, consistency_failure = [], None
    if len(valid) > 1:
        try:
            consistency = _pairwise_consistency(valid, plan)
        except Exception as error:
            consistency_failure = {"type": type(error).__name__,
                                   "reason": str(error) or type(error).__name__}
    aggregate = {
        "format": RECEIPT_FORMAT, "plan_digest": plan["plan_digest"],
        "candidate_id": plan["candidate_graph"]["candidate_id"],
        "e3_type": plan["candidate_graph"]["e3_type"], "bindings": binding_summary,
        "reference_free": True, "actual_computation": bool(valid),
        "completed_count": len(valid),
        "cancelled_count": sum(item["status"] == "cancelled" for item in receipts),
        "timeout_count": sum(item["status"] == "timeout" for item in receipts),
        "failed_count": sum(item["status"] == "failed" for item in receipts),
        "all_seed_receipts": [item["artifacts"]["receipt"] for item in receipts],
        "display_anchor_seed": valid[0]["seed"] if valid else None,
        "display_anchor_policy": "first valid seed; not a best-seed selection",
        "pairwise_prediction_consistency": consistency,
        "consistency_failure": consistency_failure, "tool_environment": run_environment,
        "scientific_status": "computed_hypothesis" if valid else "not_computed",
        "interpretation": "No experimental ternary reference or reference-accuracy claim.",
    }
    aggregate_path.write_text(json.dumps(aggregate, ensure_ascii=False, indent=2, sort_keys=True,
                                         allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    return aggregate
