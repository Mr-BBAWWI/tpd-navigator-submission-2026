"""Conservative offline checks for the existing Boltz YAML draft, not a runner."""
from __future__ import annotations

import importlib.metadata

import yaml
from rdkit import Chem

from .handoff import Snapshot, build_material, check, sha
from .molecules import canonical

REFERENCE = {
    "release": "v2.2.1", "commit": "cb04aeccdd480fd4db707f0bbafde538397fa2ac",
    "documentation": "https://github.com/jwohlwend/boltz/blob/cb04aeccdd480fd4db707f0bbafde538397fa2ac/docs/prediction.md",
    "documentation_sha256": "b1c0296954ce2e4f2a7a8c64f02c910c9211f205fab83059ac092d8231dfba24",
    "version_meaning": "Input-format reference only; installed model/runtime version is not selected or validated.",
}


class UniqueLoader(yaml.SafeLoader):
    pass


def unique_mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        check(isinstance(key, str) and key not in result, "YAML_DUPLICATE_OR_NONSTRING_KEY")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique_mapping)


def check_input(data: bytes, candidate: dict) -> dict:
    """Only the current unconstrained, SMILES-based, no-MSA draft is supported."""
    value = yaml.load(data, Loader=UniqueLoader)
    check(isinstance(value, dict) and set(value) == {"version", "sequences"}
          and type(value["version"]) is int and value["version"] == 1, "BOLTZ_INPUT_SCOPE")
    check(isinstance(value["sequences"], list), "BOLTZ_SEQUENCES_REQUIRED")
    proteins = {p["label_asym_id"]: p["sequence"] for p in candidate["proteins"]}
    check(len(proteins) == len(candidate["proteins"]) and "L" not in proteins, "BOLTZ_EXPECTED_CHAIN_CONFLICT")
    seen, ligand_seen = set(), False
    for row in value["sequences"]:
        check(isinstance(row, dict) and len(row) == 1, "BOLTZ_ENTITY_INVALID")
        kind, entity = next(iter(row.items()))
        check(isinstance(entity, dict) and isinstance(entity.get("id"), str), "BOLTZ_CHAIN_ID_INVALID")
        chain = entity["id"]
        check(chain not in seen, "BOLTZ_DUPLICATE_CHAIN")
        seen.add(chain)
        if kind == "protein":
            check(set(entity) == {"id", "sequence"}, "BOLTZ_MSA_OR_MODIFICATION_NOT_REVIEWED")
            check(chain in proteins and entity["sequence"] == proteins[chain], "BOLTZ_PROTEIN_SEQUENCE_MISMATCH")
            check(bool(entity["sequence"]) and not set(entity["sequence"]) - set("ACDEFGHIKLMNPQRSTVWY"), "BOLTZ_PROTEIN_ALPHABET")
        elif kind == "ligand":
            check(set(entity) == {"id", "smiles"} and chain == "L" and not ligand_seen, "BOLTZ_LIGAND_SCOPE")
            check(isinstance(entity["smiles"], str), "BOLTZ_SMILES_REQUIRED")
            mol = Chem.MolFromSmiles(entity["smiles"])
            check(mol is not None and len(Chem.GetMolFrags(mol)) == 1, "BOLTZ_LIGAND_INVALID")
            check(canonical(mol) == candidate["identity"]["canonical_isomeric_smiles"], "BOLTZ_MOLECULE_MISMATCH")
            ligand_seen = True
        else:
            check(False, "BOLTZ_ENTITY_UNSUPPORTED")
    check(ligand_seen and seen == set(proteins) | {"L"}, "BOLTZ_CHAIN_SET_MISMATCH")
    return {"compound_id": candidate["id"], "molecule_id": candidate["molecule_id"],
            "input_sha256": sha(data), "input_validation": "passed_offline_subset",
            "protein_chains": sorted(proteins), "sequence_sha256": {k: sha(v.encode()) for k,v in proteins.items()},
            "ligand_chain": "L", "msa_status": "not_prepared", "templates": "absent", "constraints": "absent",
            "affinity_requested": False, "prediction_status": "not_run"}


def preflight(snapshot: Snapshot) -> dict:
    material = build_material(snapshot)
    candidates = []
    for c in snapshot.handoff["candidates"]:
        name = c["prediction"]["input"]
        check(name == c["id"] + "-boltz-input.yaml" and name in snapshot.files, "BOLTZ_INPUT_REFERENCE_MISMATCH")
        candidates.append(check_input(snapshot.files[name], c))
    return {"format": "tpd-boltz-preflight/0.1.0-draft", "case_id": snapshot.handoff["case_id"],
            "material_digest": material["material_digest"], "official_reference": REFERENCE,
            "scope": "Offline YAML/sequence/ligand identity checks, not the official Boltz parser or an inference run.",
            "candidates": candidates, "execution_readiness": "blocked",
            "blockers": ["G1_AND_M2_EXECUTION_SCOPE_NOT_BOUND", "MSA_NOT_PREPARED",
                         "BOLTZ_RUNTIME_MODEL_AND_OFFICIAL_PARSER_NOT_VALIDATED", "GPU_RESOURCE_AND_BUDGET_NOT_BOUND"],
            "official_parser_validation": "not_run", "gpu_used": False, "network_used": False,
            "tools": {p: importlib.metadata.version(p) for p in ("rdkit", "PyYAML")}}
