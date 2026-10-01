"""Read one explicitly bound Boltz result and compare with its experimental reference.

No inference, ranking, approvals or efficacy thresholds. Synthetic inputs remain synthetic.
"""
from __future__ import annotations

import io
import math
from pathlib import Path

import gemmi
import numpy as np
from jsonschema import Draft202012Validator
from rdkit import Chem

from .boltz_parser import validate_receipt
from .handoff import ACTIVE, Snapshot, check, local_file, parse, sha
from .molecules import canonical, rows
from .structures import aligned_rmsd

CONFIDENCE_KEYS = ("confidence_score", "ptm", "iptm", "ligand_iptm", "protein_iptm",
                   "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde")


def read_confidence(data, chains):
    if data is None:
        return {"status": "missing", "raw": None}
    value = parse(data)
    check(isinstance(value, dict) and set(CONFIDENCE_KEYS) | {"chains_ptm", "pair_chains_iptm"} <= set(value),
          "BOLTZ_CONFIDENCE_FIELDS")

    def number(x, distance=False):
        check(type(x) in (float, int) and math.isfinite(x) and x >= 0 and (distance or x <= 1),
              "BOLTZ_CONFIDENCE_RANGE")

    for key in CONFIDENCE_KEYS:
        number(value[key], key in ("complex_pde", "complex_ipde"))
    check(isinstance(value["chains_ptm"], dict) and isinstance(value["pair_chains_iptm"], dict) and
          set(value["chains_ptm"]) == set(chains) == set(value["pair_chains_iptm"]), "BOLTZ_CONFIDENCE_CHAINS")
    for chain in chains:
        number(value["chains_ptm"][chain])
        pairs = value["pair_chains_iptm"][chain]
        check(isinstance(pairs, dict) and set(pairs) == set(chains), "BOLTZ_CONFIDENCE_PAIR_CHAINS")
        for x in pairs.values():
            number(x)
    return {"status": "present", "raw": value, "chain_index_to_id": chains,
            "meaning": "Model confidence, not measured binding/degradation or a TPD efficacy score.",
            "units": {"confidence_and_tm_fields": "dimensionless_0_to_1", "pde_fields": "angstrom"}}


def read_sites(data):
    block = gemmi.cif.read_string(data.decode("utf-8")).sole_block()
    result, seen = [], set()
    for row in rows(block, "_atom_site."):
        if row["type_symbol"].upper() in ("H", "D"):
            continue
        check(str(row.get("pdbx_PDB_model_num", "1")) == "1", "BOLTZ_MULTIPLE_MODELS_UNSUPPORTED")
        check(row.get("label_alt_id") in (False, None, "", "A"), "BOLTZ_ALTLOC_UNSUPPORTED")
        xyz = [float(row["Cartn_" + axis]) for axis in "xyz"]
        check(all(math.isfinite(v) for v in xyz), "BOLTZ_NONFINITE_COORDINATES")
        occupancy = float(row["occupancy"])
        check(math.isfinite(occupancy) and 0 < occupancy <= 1, "BOLTZ_OCCUPANCY_INVALID")
        # label_seq_id is null for the single non-polymer residue emitted by modelcif.
        seq = row.get("label_seq_id")
        key = (row["label_asym_id"], seq, row["label_atom_id"])
        check(key not in seen, "BOLTZ_DUPLICATE_ATOM")
        seen.add(key)
        result.append({**row, "xyz": xyz})
    check(bool(result), "BOLTZ_EMPTY_STRUCTURE")
    return block, result


def ca_coordinates(block, sites, chain, sequence, *, complete):
    asyms = rows(block, "_struct_asym.")
    entities = [r["entity_id"] for r in asyms if r["id"] == chain]
    check(len(entities) == 1, "BOLTZ_CHAIN_ENTITY_MISMATCH")
    seq_rows = [r for r in rows(block, "_entity_poly_seq.") if r["entity_id"] == entities[0]]
    check(len(seq_rows) == len(sequence) and {int(r["num"]) for r in seq_rows} == set(range(1, len(sequence) + 1)),
          "BOLTZ_DECLARED_SEQUENCE_LENGTH")
    declared = {int(r["num"]): r["mon_id"] for r in seq_rows}
    actual = "".join(gemmi.find_tabulated_residue(declared[i]).one_letter_code for i in range(1, len(sequence) + 1))
    check(actual == sequence, "BOLTZ_DECLARED_SEQUENCE_MISMATCH")
    result = {}
    for row in sites:
        if row["label_asym_id"] != chain:
            continue
        pos = int(row["label_seq_id"])
        check(pos in declared and row["label_comp_id"] == declared[pos], "BOLTZ_OBSERVED_SEQUENCE_MISMATCH")
        if row["label_atom_id"] == "CA":
            check(row["type_symbol"].upper() == "C", "BOLTZ_CA_ELEMENT")
            result[pos] = row
    check(len(result) >= 3, "BOLTZ_TOO_FEW_CA_PAIRS")
    if complete:
        check(set(result) == set(declared), "BOLTZ_MISSING_PREDICTED_CA")
    return result


def role_chain(candidate, accession):
    matches = [p["label_asym_id"] for p in candidate["proteins"]
               if any(r["pdbx_db_accession"] == accession for r in p["database_references"])]
    check(len(matches) == 1, "BOLTZ_REFERENCE_ROLE_AMBIGUOUS")
    return matches[0]


def compare_structure(data, receipt_candidate, candidate, snapshot):
    predicted_block, predicted_sites = read_sites(data)
    chains = set(receipt_candidate["proteins"]) | {receipt_candidate["ligand_chain"]}
    check({r["label_asym_id"] for r in predicted_sites} == chains, "BOLTZ_OUTPUT_CHAIN_SET")
    pdb_name = candidate["pdb"] + ".cif"
    # Experimental sources can contain alternate sites outside the chosen assembly.
    reference_block = gemmi.cif.read_string(snapshot.raw[pdb_name].decode("utf-8")).sole_block()
    reference_sites = []
    for row in rows(reference_block, "_atom_site."):
        if (row["label_asym_id"] not in set(receipt_candidate["proteins"]) | {candidate["ligand_label_asym_id"]}
                or row["pdbx_PDB_model_num"] != "1" or row["label_alt_id"] not in (False, None, "", "A")
                or row["type_symbol"].upper() in ("H", "D") or float(row["occupancy"]) <= 0):
            continue
        row["xyz"] = [float(row["Cartn_" + axis]) for axis in "xyz"]
        check(all(math.isfinite(v) for v in row["xyz"]), "BOLTZ_REFERENCE_COORDINATES")
        reference_sites.append(row)
    ref_ca, pred_ca = {}, {}
    for chain, sequence in receipt_candidate["proteins"].items():
        ref_ca[chain] = ca_coordinates(reference_block, reference_sites, chain, sequence, complete=False)
        pred_ca[chain] = ca_coordinates(predicted_block, predicted_sites, chain, sequence, complete=True)
    target = role_chain(candidate, snapshot.config["target"]["uniprot"])
    vhl = role_chain(candidate, "P40337")
    paired = sorted(ref_ca[target])
    fit = aligned_rmsd([ref_ca[target][i]["xyz"] for i in paired], [pred_ca[target][i]["xyz"] for i in paired])
    rotation, translation = np.array(fit["rotation_row_vectors"]), np.array(fit["translation_A"])

    def distances(ref, pred):
        return np.linalg.norm(np.asarray(pred) @ rotation + translation - np.asarray(ref), axis=1)

    def rmsd(values):
        return float(np.sqrt(np.mean(np.square(values))))

    protein_comparison = {}
    for chain in receipt_candidate["proteins"]:
        positions = sorted(ref_ca[chain])
        dist = distances([ref_ca[chain][i]["xyz"] for i in positions], [pred_ca[chain][i]["xyz"] for i in positions])
        protein_comparison[chain] = {
            "rmsd_after_target_alignment_A": rmsd(dist), "paired_ca_count": len(positions),
            "declared_sequence_length": len(receipt_candidate["proteins"][chain]),
            "missing_reference_label_seq_ids": sorted(set(pred_ca[chain]) - set(positions)),
            "pairs": [{"label_seq_id": i, "reference_auth_seq_id": ref_ca[chain][i]["auth_seq_id"],
                       "reference_atom_site_id": ref_ca[chain][i]["id"],
                       "predicted_atom_site_id": pred_ca[chain][i]["id"], "distance_A": float(d)}
                      for i, d in zip(positions, dist)]}
    ligand_sites = [r for r in predicted_sites if r["label_asym_id"] == receipt_candidate["ligand_chain"]]
    check({r["label_comp_id"] for r in ligand_sites} == {receipt_candidate["ligand_component"]}, "BOLTZ_LIGAND_COMPONENT")
    pred_ligand = {r["label_atom_id"]: r for r in ligand_sites}
    ref_ligand = {r["label_atom_id"]: r for r in reference_sites if r["label_asym_id"] == candidate["ligand_label_asym_id"]}
    check(len(pred_ligand) == len(ligand_sites) and set(pred_ligand) ==
          {a["boltz_atom_name"] for a in receipt_candidate["ligand_mappings"][0]}, "BOLTZ_LIGAND_ATOM_SET")
    parts = snapshot.json(candidate["id"] + "-parts.json")["parts"]
    part_maps = {p["role"]: {a.GetAtomMapNum() for a in Chem.MolFromSmiles(p["mapped_smiles"]).GetAtoms()
                             if a.GetAtomicNum() > 0} for p in parts}
    source = next(Chem.ForwardSDMolSupplier(io.BytesIO(snapshot.files[candidate["id"] + ".sdf"])))
    comparisons = []
    for index, mapping in enumerate(receipt_candidate["ligand_mappings"]):
        check(set(ref_ligand) == {a["ccd_atom_id"] for a in mapping}, "BOLTZ_REFERENCE_LIGAND_ATOMS")
        check(all(pred_ligand[a["boltz_atom_name"]]["type_symbol"].upper() == a["element"].upper() ==
                  ref_ligand[a["ccd_atom_id"]]["type_symbol"].upper() for a in mapping), "BOLTZ_LIGAND_ELEMENTS")
        dist = distances([ref_ligand[a["ccd_atom_id"]]["xyz"] for a in mapping],
                         [pred_ligand[a["boltz_atom_name"]]["xyz"] for a in mapping])
        conformer = Chem.Conformer(source.GetNumAtoms())
        by_map = {a["atom_map"]: pred_ligand[a["boltz_atom_name"]]["xyz"] for a in mapping}
        for atom in source.GetAtoms():
            conformer.SetAtomPosition(atom.GetIdx(), by_map[atom.GetAtomMapNum()])
        mol = Chem.Mol(source)
        mol.RemoveAllConformers()
        mol.AddConformer(conformer)
        Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
        lengths = [conformer.GetAtomPosition(b.GetBeginAtomIdx()).Distance(conformer.GetAtomPosition(b.GetEndAtomIdx())) for b in mol.GetBonds()]
        comparisons.append({"mapping_index": index, "whole_ligand_rmsd_A": rmsd(dist),
                            "parts_rmsd_A": {role: rmsd([d for a, d in zip(mapping, dist) if a["atom_map"] in maps])
                                             for role, maps in part_maps.items()},
                            "coordinates_match_input_stereochemistry": canonical(mol) == canonical(source),
                            "bond_length_range_A": [min(lengths), max(lengths)],
                            "pairs": [{**a, "distance_A": float(d)} for a, d in zip(mapping, dist)]})
    values = [r["whole_ligand_rmsd_A"] for r in comparisons]
    return {"reference_pdb": candidate["pdb"], "reference_sha256": sha(snapshot.raw[pdb_name]),
            "target_chain": target, "vhl_chain": vhl, "target_alignment": fit,
            "method": "Kabsch on observed target C-alpha atoms paired by deposited label_seq_id; apply ONE transform to all chains and ligand. No ligand/VHL refit.",
            "reference_selection": "assembly protein chains and configured ligand; model 1, blank/A altloc, positive occupancy, heavy atoms",
            "proteins": protein_comparison, "ligand": {
                "mapping_count": len(comparisons), "all_mappings": comparisons,
                "whole_ligand_rmsd_range_A": [min(values), max(values)],
                "lowest_whole_ligand_rmsd_mapping_index": values.index(min(values)),
                "selection_meaning": "Geometric minimum across equivalent atom labels; not candidate ranking or unique atom correspondence."},
            "interpretation": "Descriptive structural comparison only. Known complexes may be in model training data; this is not prospective generalization or degradation efficacy validation."}


def read_result(directory: Path, receipt, snapshot: Snapshot, *, allow_synthetic=False):
    validate_receipt(receipt, snapshot)
    manifest_bytes = local_file(directory, "manifest.json").read_bytes()
    manifest = parse(manifest_bytes)
    schema = parse((ACTIVE / "contracts/drafts/b_boltz_output_manifest.schema.json").read_bytes())
    check(not list(Draft202012Validator(schema).iter_errors(manifest)), "BOLTZ_OUTPUT_MANIFEST_SCHEMA")
    check(manifest["origin"] != "synthetic_test" or allow_synthetic, "BOLTZ_SYNTHETIC_NOT_ALLOWED")
    check(manifest["receipt_digest"] == receipt["digest"], "BOLTZ_OUTPUT_STALE_RECEIPT")
    candidates = [c for c in receipt["candidates"] if c["compound_id"] == manifest["compound_id"]]
    check(len(candidates) == 1, "BOLTZ_OUTPUT_CANDIDATE")
    c = candidates[0]
    check(all(manifest[k] == c[k] for k in ("molecule_id", "input_sha256", "record_id")), "BOLTZ_OUTPUT_INPUT_BINDING")
    names = {"input": c["input_name"], "structure": f'{c["record_id"]}_model_{manifest["model_rank"]}.cif',
             "confidence": f'confidence_{c["record_id"]}_model_{manifest["model_rank"]}.json'}
    files = {}
    for kind, info in manifest["files"].items():
        if info is None:
            files[kind] = None
            continue
        check(info["name"] == names[kind], "BOLTZ_OUTPUT_FILENAME_BINDING")
        data = local_file(directory, info["name"]).read_bytes()
        check(sha(data) == info["sha256"], "BOLTZ_OUTPUT_FILE_HASH")
        files[kind] = data
    check(sha(files["input"]) == c["input_sha256"], "BOLTZ_OUTPUT_INPUT_HASH")
    confidence = read_confidence(files["confidence"], c["chain_index_to_id"])
    candidate = next(x for x in snapshot.handoff["candidates"] if x["id"] == c["compound_id"])
    comparison = compare_structure(files["structure"], c, candidate, snapshot)
    return {"format": "tpd-boltz-comparison/0.1.0-draft", "origin": manifest["origin"],
            "compound_id": c["compound_id"], "molecule_id": c["molecule_id"],
            "material_digest": receipt["material_digest"], "receipt_digest": receipt["digest"],
            "manifest_sha256": sha(manifest_bytes), "source_manifest": manifest,
            "execution_provenance": "caller_reported_not_authenticated_by_this_reader",
            "confidence": confidence, "comparison": comparison,
            "human_review": "pending", "approval_record_created": False, "dispatch_authorized": False,
            "efficacy_claim": "not_established"}
