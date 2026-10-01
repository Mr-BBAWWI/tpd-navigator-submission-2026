"""CPU-only official Boltz parser adapter. This module never dispatches inference."""
from __future__ import annotations

import importlib.metadata
from pathlib import Path

from rdkit import Chem

from .boltz_preflight import REFERENCE, check_input
from .handoff import ACTIVE, Snapshot, build_material, check, encoded, parse, sha
from .molecules import canonical

FORMAT = "tpd-boltz-parser/0.1.0-draft"
ARCHIVE_SHA256 = "39e076d96dbec6b4e86982bbda16f3a53a2a60c9bdc17828d88f6f9a0c7d1fd7"
CCD_REVISION = "6fdef46d763fee7fbb83ca5501ccceff43b85607"
CCD_URL = f"https://huggingface.co/boltz-community/boltz-2/resolve/{CCD_REVISION}/mols.tar"
PARSER_SOURCE_SHA256 = {
    "data/parse/yaml.py": "19115132d0c3e0f4383c8a40fd2d638f4bd1684bb62e565dfe0e2b375cd4acda",
    "data/parse/schema.py": "c0dae4c89c4c4a22175a433b5a0c68860329d2a292c76224328bdde19b73b743",
}


def seal(value):
    return {**value, "digest": sha(encoded(value))}


def validate_receipt(receipt, snapshot):
    from jsonschema import Draft202012Validator
    schema = parse((ACTIVE / "contracts/drafts/b_boltz_parser.schema.json").read_bytes())
    check(not list(Draft202012Validator(schema).iter_errors(receipt)), "BOLTZ_RECEIPT_SCHEMA")
    check(receipt["digest"] == sha(encoded({k: v for k, v in receipt.items() if k != "digest"})),
          "BOLTZ_RECEIPT_DIGEST")
    check(receipt["material_digest"] == build_material(snapshot)["material_digest"], "STALE_BOLTZ_RECEIPT")
    check(receipt["case_id"] == snapshot.handoff["case_id"], "BOLTZ_RECEIPT_CASE")
    check(all(receipt["installed_source_sha256"].get(n) == h for n, h in PARSER_SOURCE_SHA256.items()),
          "BOLTZ_PARSER_SOURCE_VERSION")
    expected = {c["id"]: c for c in snapshot.handoff["candidates"]}
    check(len(receipt["candidates"]) == len(expected), "BOLTZ_RECEIPT_CANDIDATES")
    check({c["compound_id"] for c in receipt["candidates"]} == set(expected), "BOLTZ_RECEIPT_CANDIDATES")
    for c in receipt["candidates"]:
        source = expected[c["compound_id"]]
        check(c["molecule_id"] == source["molecule_id"] and
              c["input_sha256"] == sha(snapshot.files[source["prediction"]["input"]]), "BOLTZ_RECEIPT_INPUT")
        check(c["proteins"] == {p["label_asym_id"]: p["sequence"] for p in source["proteins"]},
              "BOLTZ_RECEIPT_SEQUENCES")
        check(c["chain_index_to_id"] == {str(i): ch for i, ch in enumerate([*c["proteins"], "L"])},
              "BOLTZ_RECEIPT_CHAINS")
        check(c["mapping_count"] == len(c["ligand_mappings"]), "BOLTZ_RECEIPT_MAP_COUNT")
        ref_atoms = snapshot.json(source["id"] + "-atom-map.json")["reference_atoms"]
        expected_atoms = {(a["ccd_atom_id"], a["atom_map"]) for a in ref_atoms}
        first = c["ligand_mappings"][0]
        names = {a["boltz_atom_name"] for a in first}
        check(len(names) == len(expected_atoms), "BOLTZ_RECEIPT_ATOM_NAMES")
        for mapping in c["ligand_mappings"]:
            check(len(mapping) == len(expected_atoms) and
                  {(a["ccd_atom_id"], a["atom_map"]) for a in mapping} == expected_atoms and
                  {a["boltz_atom_name"] for a in mapping} == names, "BOLTZ_RECEIPT_ATOM_MAP")
    return receipt


def ligand_mappings(parsed_mol, sdf_bytes):
    """Enumerate complete stereo-preserving graph isomorphisms, retaining symmetry."""
    import io
    source = next(Chem.ForwardSDMolSupplier(io.BytesIO(sdf_bytes), removeHs=True))
    check(source is not None and canonical(source) == canonical(parsed_mol), "BOLTZ_PARSED_IDENTITY")
    ccd_names = source.GetProp("ccd_atom_names_in_sdf_order").split(",")
    query = Chem.Mol(source)
    for atom in query.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(query, cleanIt=True, force=True)
    matches = parsed_mol.GetSubstructMatches(query, useChirality=True, uniquify=False, maxMatches=257)
    check(0 < len(matches) < 257, "BOLTZ_MAPPING_MISSING_OR_ENUMERATION_LIMIT")
    result = []
    for match in sorted(set(matches)):
        check(len(set(match)) == parsed_mol.GetNumAtoms() == source.GetNumAtoms(), "BOLTZ_MAPPING_INCOMPLETE")
        result.append([{"boltz_atom_name": parsed_mol.GetAtomWithIdx(j).GetProp("name"),
                        "ccd_atom_id": ccd_names[i], "atom_map": source.GetAtomWithIdx(i).GetAtomMapNum(),
                        "element": source.GetAtomWithIdx(i).GetSymbol()}
                       for i, j in enumerate(match)])
    return result


def prepare_receipt(snapshot: Snapshot, mol_dir: Path):
    """Run unmodified parse_yaml(..., boltz2=True) using verified official CCD data."""
    check(importlib.metadata.version("boltz") == "2.2.1", "BOLTZ_VERSION_NOT_PINNED")
    from boltz.data import const
    from boltz.data.mol import load_molecules
    from boltz.data.parse.yaml import parse_yaml
    import boltz

    root = Path(boltz.__file__).parent
    check(all(sha((root / n).read_bytes()) == h for n, h in PARSER_SOURCE_SHA256.items()),
          "BOLTZ_PARSER_SOURCE_VERSION")
    manifest = parse((mol_dir / "manifest.json").read_bytes())
    check(manifest["archive_sha256"] == ARCHIVE_SHA256 and manifest["url"] == CCD_URL,
          "BOLTZ_CCD_PROVENANCE")
    needed = sorted({const.prot_letter_to_token[a] for c in snapshot.handoff["candidates"]
                     for p in c["proteins"] for a in p["sequence"]})
    for name in needed:
        check(sha((mol_dir / (name + ".pkl")).read_bytes()) == manifest["files"][name + ".pkl"],
              "BOLTZ_CCD_HASH_MISMATCH")
    ccd = load_molecules(mol_dir, needed)
    candidates, targets = [], {}
    for c in snapshot.handoff["candidates"]:
        name = c["prediction"]["input"]
        preflight = check_input(snapshot.files[name], c)
        target = parse_yaml(snapshot.directory / name, ccd, mol_dir, boltz2=True)
        check(sha((snapshot.directory / name).read_bytes()) == preflight["input_sha256"], "BOLTZ_INPUT_CHANGED_DURING_PARSE")
        chains = {str(int(ch["asym_id"])): str(ch["name"]) for ch in target.structure.chains}
        check(list(chains.values()) == [p["label_asym_id"] for p in c["proteins"]] + ["L"],
              "BOLTZ_PARSED_CHAIN_ORDER")
        check(len(target.extra_mols) == 1, "BOLTZ_PARSED_LIGAND_COUNT")
        component, mol = next(iter(target.extra_mols.items()))
        mappings = ligand_mappings(mol, snapshot.files[c["id"] + ".sdf"])
        candidates.append({"compound_id": c["id"], "molecule_id": c["molecule_id"],
                           "input_name": name, "input_sha256": preflight["input_sha256"],
                           "record_id": target.record.id, "chain_index_to_id": chains,
                           "proteins": {p["label_asym_id"]: p["sequence"] for p in c["proteins"]},
                           "ligand_chain": "L", "ligand_component": component,
                           "ligand_mappings": mappings, "mapping_count": len(mappings),
                           "parsed_atom_count": len(target.structure.atoms),
                           "parsed_residue_count": len(target.structure.residues)})
        targets[c["id"]] = target
    source_names = ["data/parse/yaml.py", "data/parse/schema.py", "data/types.py",
                    "data/write/mmcif.py", "data/write/writer.py"]
    receipt = seal({"format": FORMAT, "case_id": snapshot.handoff["case_id"],
                    "material_digest": build_material(snapshot)["material_digest"],
                    "official_release": REFERENCE["release"], "official_commit": REFERENCE["commit"],
                    "installed_source_sha256": {n: sha((root / n).read_bytes()) for n in source_names},
                    "ccd": manifest,
                    "runtime": {p: importlib.metadata.version(p) for p in
                                ("boltz", "torch", "rdkit", "numpy", "scipy", "gemmi", "PyYAML", "modelcif")},
                    "candidates": candidates, "parser_validation": "passed",
                    "prediction_status": "not_run", "execution_readiness": "blocked",
                    "blockers": ["G1_M2_SCOPE_NOT_BOUND", "MSA_NOT_PREPARED", "INFERENCE_ENVIRONMENT_WEIGHTS_AND_GPU_NOT_VALIDATED"],
                    "symmetry_policy": "All enumerated stereo-preserving mappings retained; no unique CCD-name claim.",
                    "coordinate_kind": "parser_reference_conformer_not_prediction",
                    "gpu_used": False, "network_used_during_parse": False})
    validate_receipt(receipt, snapshot)
    return receipt, targets
