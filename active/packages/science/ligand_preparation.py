"""Bounded, offline hydrogen preparation with unchanged input chemistry/heavy atoms.

Produces a separate ligand SDF for review. It does not choose a protonation state,
optimize a complex, assign Probe atom types, or change existing quality findings.
Scientific imports are deferred so archive integrity can be checked in core.
"""
from __future__ import annotations

import copy
import io
import math
from pathlib import Path

from .evidence_common import file_index, safe_relative, schema_check, seal, verify_seal, write_new_directory
from .handoff import check, encoded, sha
from .structure_quality import strict_json

REQUEST = "tpd-ligand-preparation-request/0.1.0-draft"
REPORT = "tpd-ligand-preparation-report/0.1.0-draft"
ARCHIVE = "tpd-ligand-preparation-archive/0.1.0-draft"
LIMITS = ["PROTONATION_NOT_PREDICTED", "HYDROGEN_ORIENTATION_NOT_OPTIMIZED",
          "PROBE_ATOM_TYPING_UNVERIFIED", "CONTACT_VALIDATION_NOT_RUN",
          "WHOLE_COMPLEX_NOT_PREPARED", "NOT_A_SCIENTIFIC_APPROVAL"]
INPUTS = ("source_structure", "source_molecule", "model_input")


def validate_request(request, files, *, allow_synthetic=False):
    schema_check(request, "b_ligand_preparation.schema.json")
    check(request["data_mode"] != "synthetic_test" or allow_synthetic, "SYNTHETIC_LIGAND_PREPARATION_NOT_ALLOWED")
    names = [request[k]["path"] for k in INPUTS]
    check(len(set(names)) == len(names) and all(n.startswith("raw/") for n in names), "LIGAND_INPUT_PATHS")
    for key in INPUTS:
        ref = request[key]
        check(sha(files[ref["path"]]) == ref["sha256"], "LIGAND_INPUT_HASH")
    check(request["source_structure"]["sha256"] == request["binding"]["structure_sha256"]
          and request["model_input"]["sha256"] == request["binding"]["input_sha256"], "LIGAND_SOURCE_BINDING")


def _selection(data, selection):
    import gemmi
    structure = gemmi.make_structure_from_block(gemmi.cif.read_string(data.decode("utf-8")).sole_block())
    check(len(structure) == 1, "LIGAND_SINGLE_MODEL_REQUIRED")
    residues = [r for c in structure[0] if c.name == selection["auth_chain_id"] for r in c
                if r.seqid.num == selection["auth_seq_id"] and r.seqid.icode.strip() == selection["insertion_code"]
                and r.name == selection["component_id"]]
    check(len(residues) == 1, "LIGAND_SELECTION_AMBIGUOUS_OR_MISSING")
    atoms = {}
    for a in residues[0]:
        check(a.element.name not in ("H", "D"), "LIGAND_SOURCE_HYDROGENS_REQUIRE_EXPLICIT_POLICY")
        check(not a.altloc.replace("\x00", "").strip(), "LIGAND_ALTERNATE_CONFORMER_UNSUPPORTED")
        check(math.isfinite(a.occ) and abs(a.occ - 1.0) < 1e-6, "LIGAND_PARTIAL_OCCUPANCY_UNSUPPORTED")
        xyz = (a.pos.x, a.pos.y, a.pos.z)
        check(all(math.isfinite(v) for v in xyz), "LIGAND_NONFINITE_COORDINATES")
        check(a.name not in atoms, "LIGAND_DUPLICATE_SOURCE_ATOM")
        atoms[a.name] = {"element": a.element.name, "xyz": xyz}
    check(bool(atoms) and len({v["xyz"] for v in atoms.values()}) == len(atoms), "LIGAND_EMPTY_OR_COINCIDENT_ATOMS")
    return atoms


def _prepare(request, files):
    import gemmi
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    from .molecules import canonical

    check(rdkit.__version__ == "2026.03.6" and gemmi.__version__ == "0.7.5", "LIGAND_RUNTIME_VERSION_UNSUPPORTED")
    source = list(Chem.ForwardSDMolSupplier(io.BytesIO(files[request["source_molecule"]["path"]]), removeHs=False))
    check(len(source) == 1 and source[0] is not None, "LIGAND_SINGLE_VALID_SDF_REQUIRED")
    mol = source[0]
    check(len(Chem.GetMolFrags(mol)) == 1, "LIGAND_DISCONNECTED_INPUT_UNSUPPORTED")
    supported = {5, 6, 7, 8, 9, 15, 16, 17, 35, 53}
    check(all(a.GetAtomicNum() in supported and not a.GetIsotope() and not a.GetNumRadicalElectrons()
              and not a.HasQuery() for a in mol.GetAtoms()), "LIGAND_ELEMENT_ISOTOPE_RADICAL_OR_QUERY_UNSUPPORTED")
    maps = [a.GetAtomMapNum() for a in mol.GetAtoms()]
    check(all(m > 0 for m in maps) and len(set(maps)) == len(maps), "LIGAND_UNIQUE_HEAVY_MAPS_REQUIRED")
    reference = Chem.MolFromSmiles(request["expected_isomeric_smiles"])
    check(reference is not None and canonical(reference) == canonical(mol), "LIGAND_GRAPH_IDENTITY_MISMATCH")
    identity = canonical(mol)
    check(request["binding"]["molecule_id"] == request["binding"]["compound_id"] + ":" + sha(identity.encode()),
          "LIGAND_MOLECULE_VERSION_MISMATCH")
    stereo = Chem.MolFromSmiles(identity)
    check(not any(str(s.specified) == "Unspecified" for s in Chem.FindPotentialStereo(stereo)),
          "LIGAND_UNSPECIFIED_STEREOCHEMISTRY")
    mapping = request["atom_mapping"]
    check(len(mapping) == len(maps) and {a["atom_map"] for a in mapping} == set(maps)
          and len({a["source_atom_name"] for a in mapping}) == len(maps), "LIGAND_ATOM_MAPPING_NOT_BIJECTIVE")
    names = {a["atom_map"]: a["source_atom_name"] for a in mapping}
    atoms = _selection(files[request["source_structure"]["path"]], request["selection"])
    check(set(atoms) == set(names.values()), "LIGAND_HEAVY_ATOM_SET_MISMATCH")
    conf = Chem.Conformer(len(maps)); conf.Set3D(True)
    for a in mol.GetAtoms():
        site = atoms[names[a.GetAtomMapNum()]]
        check(site["element"] == a.GetSymbol(), "LIGAND_ELEMENT_MAPPING_MISMATCH")
        conf.SetAtomPosition(a.GetIdx(), site["xyz"])
    mol.RemoveAllConformers(); mol.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    check(canonical(mol) == identity, "LIGAND_COORDINATE_STEREOCHEMISTRY_MISMATCH")
    expected_h = {a.GetAtomMapNum(): a.GetTotalNumHs() for a in mol.GetAtoms()}
    prepared = Chem.AddHs(mol, addCoords=True)
    check(canonical(Chem.RemoveHs(prepared)) == identity, "LIGAND_H_ADDITION_CHANGED_IDENTITY")
    check(Chem.GetFormalCharge(prepared) == Chem.GetFormalCharge(mol), "LIGAND_FORMAL_CHARGE_CHANGED")
    # Original ligand coordinates are never optimized or replaced by SDF depiction coordinates.
    for a in mol.GetAtoms():
        before = conf.GetAtomPosition(a.GetIdx()); after = prepared.GetConformer().GetAtomPosition(a.GetIdx())
        check(before.Distance(after) < 1e-12, "LIGAND_HEAVY_COORDINATES_CHANGED")
    prepared.SetProp("_Name", "RDKit H geometry on unchanged heavy atoms; review only")
    text = Chem.MolToMolBlock(prepared, forceV3000=True) + "\n$$$$\n"
    reread = list(Chem.ForwardSDMolSupplier(io.BytesIO(text.encode()), removeHs=False))
    check(len(reread) == 1 and reread[0] is not None, "LIGAND_PREPARED_SDF_UNREADABLE")
    reread = reread[0]
    check(canonical(Chem.RemoveHs(reread)) == identity, "LIGAND_SERIALIZED_IDENTITY_CHANGED")
    coordinate_stereo = Chem.RemoveHs(reread)
    Chem.AssignStereochemistryFrom3D(coordinate_stereo, replaceExistingTags=True)
    check(canonical(coordinate_stereo) == identity, "LIGAND_SERIALIZED_STEREOCHEMISTRY_CHANGED")
    check(reread.GetNumAtoms() == prepared.GetNumAtoms(), "LIGAND_SERIALIZED_ATOM_COUNT_CHANGED")
    inventory = []; counts = {m: 0 for m in maps}; max_delta = 0.0; h_lengths = []
    for a in reread.GetAtoms():
        p = reread.GetConformer().GetAtomPosition(a.GetIdx()); xyz = [p.x, p.y, p.z]
        check(all(math.isfinite(v) for v in xyz), "LIGAND_H_COORDINATES_NONFINITE")
        if a.GetAtomicNum() == 1:
            check(a.GetDegree() == 1 and a.GetFormalCharge() == 0, "LIGAND_INVALID_H_CONNECTIVITY")
            parent = a.GetNeighbors()[0]; parent_map = parent.GetAtomMapNum()
            check(parent_map in counts, "LIGAND_INVALID_H_PARENT")
            counts[parent_map] += 1
            name = f"H{parent_map}_{counts[parent_map]}"
            length = p.Distance(reread.GetConformer().GetAtomPosition(parent.GetIdx()))
            check(math.isfinite(length) and length > 0, "LIGAND_INVALID_H_BOND_LENGTH")
            h_lengths.append(length)
        else:
            parent_map = None; name = names[a.GetAtomMapNum()]
            delta = math.dist(xyz, atoms[name]["xyz"]); max_delta = max(max_delta, delta)
            check(delta < 1e-5, "LIGAND_SERIALIZED_HEAVY_COORDINATES_CHANGED")
        inventory.append({"sdf_atom_index": a.GetIdx(), "element": a.GetSymbol(), "atom_map": a.GetAtomMapNum(),
                          "source_or_generated_name": name, "formal_charge": a.GetFormalCharge(),
                          "parent_heavy_atom_map": parent_map, "xyz_A": xyz})
    check(counts == expected_h, "LIGAND_H_COUNTS_MISMATCH")
    details = {"canonical_isomeric_smiles": identity, "formula": rdMolDescriptors.CalcMolFormula(reread),
               "formal_charge": Chem.GetFormalCharge(reread), "heavy_atom_count": len(maps),
               "expected_hydrogen_count": sum(expected_h.values()), "generated_hydrogen_count": len(h_lengths),
               "hydrogens_by_heavy_atom_map": {str(k): v for k, v in counts.items()},
               "max_serialized_heavy_displacement_A": max_delta,
               "hydrogen_bond_length_range_A": [min(h_lengths), max(h_lengths)] if h_lengths else None,
               "atom_inventory": inventory}
    return details, text.encode(), {"rdkit": rdkit.__version__, "gemmi": gemmi.__version__}


def _report(request, outcome, files):
    check(set(outcome) == {"status", "error", "chemistry", "runtime"}, "LIGAND_OUTCOME_FIELDS")
    prepared = outcome["status"] == "prepared_for_review"
    check(outcome["status"] in ("prepared_for_review", "not_prepared"), "LIGAND_OUTCOME_STATUS")
    check(prepared == (outcome["error"] is None) == (outcome["chemistry"] is not None), "LIGAND_OUTCOME_CONSISTENCY")
    check(prepared == ("prepared/ligand.sdf" in files), "LIGAND_OUTPUT_STATUS_MISMATCH")
    if prepared:
        c = outcome["chemistry"]; atoms = c["atom_inventory"]
        heavy = [a for a in atoms if a["element"] != "H"]; hydrogens = [a for a in atoms if a["element"] == "H"]
        check(len(heavy) == c["heavy_atom_count"] and len(hydrogens) == c["expected_hydrogen_count"]
              == c["generated_hydrogen_count"], "LIGAND_INVENTORY_COUNTS")
        check([a["sdf_atom_index"] for a in atoms] == list(range(len(atoms))), "LIGAND_INVENTORY_ORDER")
        check(sum(a["formal_charge"] for a in atoms) == c["formal_charge"], "LIGAND_INVENTORY_CHARGE")
        check({a["atom_map"] for a in heavy} == {a["atom_map"] for a in request["atom_mapping"]}, "LIGAND_INVENTORY_MAPS")
        check(len({a["atom_map"] for a in heavy}) == len(heavy), "LIGAND_INVENTORY_DUPLICATE_MAPS")
        names = {a["atom_map"]: a["source_atom_name"] for a in request["atom_mapping"]}
        check(all(a["source_or_generated_name"] == names[a["atom_map"]] for a in heavy), "LIGAND_INVENTORY_NAMES")
        counts = {str(a["atom_map"]): 0 for a in heavy}
        for a in hydrogens:
            check(str(a["parent_heavy_atom_map"]) in counts, "LIGAND_INVENTORY_H_PARENT")
            counts[str(a["parent_heavy_atom_map"])] += 1
        check(counts == c["hydrogens_by_heavy_atom_map"], "LIGAND_INVENTORY_H_COUNTS")
    report = seal({"format": REPORT, "data_mode": request["data_mode"], "binding": copy.deepcopy(request["binding"]),
                   "request_sha256": sha(encoded(request)), "status": outcome["status"], "error": outcome["error"],
                   "chemistry": outcome["chemistry"], "runtime": outcome["runtime"],
                   "prepared_structure": {"path": "prepared/ligand.sdf", "sha256": sha(files["prepared/ligand.sdf"])} if prepared else None,
                   "method": "RDKit AddHs(addCoords=True), unchanged input chemical state and heavy coordinates",
                   "hydrogen_convention": "rdkit_geometry_not_xray_or_neutron_calibrated",
                   "limitations": LIMITS, "contact_quality_status": "not_assessed", "probe_atom_typing": "unverified",
                   "human_review": "pending", "automatic_acceptance": False})
    schema_check(report, "b_ligand_preparation_report.schema.json")
    return report


def prepare_ligand(request_path, destination, *, allow_synthetic=False):
    request_path = Path(request_path)
    check(not Path(destination).exists(), "OUTPUT_EXISTS")
    request = strict_json(request_path.read_bytes())
    schema_check(request, "b_ligand_preparation.schema.json")
    files = {request[k]["path"]: safe_relative(request_path.parent, request[k]["path"]).read_bytes() for k in INPUTS}
    validate_request(request, files, allow_synthetic=allow_synthetic)
    try:
        chemistry, sdf, runtime = _prepare(request, files)
        files["prepared/ligand.sdf"] = sdf
        outcome = {"status": "prepared_for_review", "error": None, "chemistry": chemistry, "runtime": runtime}
    except ValueError as error:
        # Only declared preparation rejections become not_prepared; programming errors propagate.
        if not str(error).startswith("LIGAND_"):
            raise
        outcome = {"status": "not_prepared", "error": str(error), "chemistry": None, "runtime": {}}
    report = _report(request, outcome, files)
    files.update({"request.json": encoded(request), "outcome.json": encoded(outcome), "report.json": encoded(report)})
    files["manifest.json"] = encoded(seal({"format": ARCHIVE, "files": file_index(files), "report_digest": report["digest"]}))
    write_new_directory(destination, files)
    return report


def read_preparation(directory, *, allow_synthetic=False):
    """Core-only integrity/projection check; does not rerun RDKit chemistry."""
    directory = Path(directory)
    manifest = strict_json((directory / "manifest.json").read_bytes()); verify_seal(manifest)
    check(set(manifest) == {"format", "files", "report_digest", "digest"}, "LIGAND_ARCHIVE_FIELDS")
    check(manifest["format"] == ARCHIVE, "LIGAND_ARCHIVE_FORMAT")
    files = {}
    for name, ref in manifest["files"].items():
        path = safe_relative(directory, name)
        check(not path.is_symlink(), "LIGAND_ARCHIVE_SYMLINK")
        data = path.read_bytes()
        check(sha(data) == ref["sha256"] and len(data) == ref["bytes"], "LIGAND_ARCHIVE_HASH")
        files[name] = data
    check(not any(p.is_symlink() for p in directory.rglob("*")), "LIGAND_ARCHIVE_SYMLINK")
    check({p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file()}
          == set(files) | {"manifest.json"}, "LIGAND_ARCHIVE_UNINDEXED_FILES")
    request = strict_json(files["request.json"]); validate_request(request, files, allow_synthetic=allow_synthetic)
    expected = {request[k]["path"] for k in INPUTS} | {"request.json", "outcome.json", "report.json"}
    if "prepared/ligand.sdf" in files:
        expected.add("prepared/ligand.sdf")
    check(set(files) == expected, "LIGAND_ARCHIVE_UNEXPECTED_FILES")
    report = _report(request, strict_json(files["outcome.json"]), files)
    check(encoded(report) == files["report.json"] and manifest["report_digest"] == report["digest"], "LIGAND_REPORT_PROJECTION")
    return report
