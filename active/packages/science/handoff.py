"""Read-only H1/H2 material from a B CPU snapshot; no approval or dispatch."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from jsonschema import Draft202012Validator
from referencing import Registry, Resource

ACTIVE = Path(__file__).resolve().parents[2]
SCHEMA = ACTIVE / "contracts/drafts/b_review_material.schema.json"


class HandoffError(ValueError):
    """Invalid or inconsistent handoff input; never a scientific verdict."""


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def encoded(value) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True,
                       separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def parse(data: bytes):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise HandoffError("DUPLICATE_JSON_KEY")
            result[key] = value
        return result
    def invalid(value):
        raise HandoffError("NONFINITE_JSON")
    return json.loads(data, object_pairs_hook=pairs, parse_constant=invalid)


def check(condition, code):
    if not condition:
        raise HandoffError(code)


def local_file(directory: Path, name: str) -> Path:
    check(isinstance(name, str) and name not in {"", ".", ".."}
          and not any(c in name for c in '/\\:'), "INVALID_ARTIFACT_PATH")
    path = (directory / name).resolve()
    check(path.parent == directory.resolve(), "ARTIFACT_OUTSIDE_DIRECTORY")
    return path


def artifact_ref(data: bytes, media="application/json", schema="urn:tpd-navigator:raw:1"):
    digest = sha(data)
    return {"artifact_id": "b-material-" + digest, "version": 1, "sha256": digest,
            "media_type": media, "schema_id": schema, "provenance": "computed"}


class Snapshot:
    """Load bytes once, verify both manifests, and retain their exact identity."""
    def __init__(self, directory: Path, raw_directory: Path):
        self.directory = Path(directory).resolve()
        self.handoff_bytes = local_file(self.directory, "handoff.json").read_bytes()
        self.handoff = parse(self.handoff_bytes)
        schema = parse((ACTIVE / "contracts/drafts/b_case_preparation.schema.json").read_bytes())
        check(not list(Draft202012Validator(schema).iter_errors(self.handoff)), "B_HANDOFF_SCHEMA")
        self.files = {}
        for row in self.handoff["artifacts"]:
            name = row["path"]
            check(name not in self.files and name != "handoff.json", "DUPLICATE_ARTIFACT")
            data = local_file(self.directory, name).read_bytes()
            check(sha(data) == row["sha256"] and len(data) == row["bytes"], "ARTIFACT_HASH_MISMATCH")
            self.files[name] = data
        self.config = self.json("case-input.json")
        check(sha(self.files["case-input.json"]) == self.handoff["input_sha256"], "CASE_INPUT_MISMATCH")
        check(self.config["case_id"] == self.handoff["case_id"]
              and self.config["version"] == self.handoff["case_version"]
              and self.config["target"] == self.handoff["target"], "CASE_IDENTITY_MISMATCH")
        self.raw = {}
        rows = self.json("sources.json")
        declared = self.config["sources"]
        check(len({r["file"] for r in declared}) == len(declared)
              and len(rows) == len(declared), "SOURCE_MANIFEST_MISMATCH")
        expected = {r["file"]: r for r in declared}
        for row in rows:
            name = row["file"]
            check(name in expected and name not in self.raw, "SOURCE_MANIFEST_MISMATCH")
            check(all(row[k] == expected[name][k] for k in ("file", "url", "sha256")),
                  "SOURCE_MANIFEST_MISMATCH")
            data = local_file(Path(raw_directory), name).read_bytes()
            check(sha(data) == row["sha256"] and len(data) == row["bytes"], "SOURCE_HASH_MISMATCH")
            self.raw[name] = data
        self.literature = self.json("literature-evidence.json")
        ev = self.literature
        check(ev["doi"] == self.config["doi"] and ev["source_file"] in self.raw
              and ev["source_sha256"] == sha(self.raw[ev["source_file"]]), "LITERATURE_SOURCE_MISMATCH")
        check(all(ev[k] == v for k, v in self.config["literature_annotations"].items()),
              "LITERATURE_ANNOTATION_MISMATCH")
        self.article = ET.fromstring(self.raw[ev["source_file"]])
        check(ev["doi"] in [e.text for e in self.article.findall(".//article-id[@pub-id-type='doi']")],
              "ARTICLE_DOI_MISMATCH")

    def json(self, name):
        check(name in self.files, "REQUIRED_ARTIFACT_MISSING: " + name)
        return parse(self.files[name])

    def file_ref(self, name):
        check(name in self.files, "REQUIRED_ARTIFACT_MISSING: " + name)
        media = {".sdf": "chemical/x-mdl-sdfile", ".svg": "image/svg+xml",
                 ".smi": "chemical/x-daylight-smiles", ".yaml": "application/yaml"}.get(Path(name).suffix, "application/json")
        return {"path": name, "bytes": len(self.files[name]),
                "ref": artifact_ref(self.files[name], media)}

    def citation(self, locator):
        node = self.article.find(locator)
        check(node is not None, "SOURCE_LOCATOR_MISSING")
        # Preserve a locator/text digest, not another unbounded copy of the article.
        text = "".join(node.itertext())
        return {"doi": self.config["doi"], "source_file": self.literature["source_file"],
                "source_sha256": self.literature["source_sha256"], "locator": locator,
                "text_sha256": sha(text.encode("utf-8"))}


def atom_table(rows, count, chain):
    check(len(rows) == count and len({r["ccd_atom_id"] for r in rows}) == count
          and len({r["atom_map"] for r in rows}) == count
          and {r["rdkit_index_zero_based"] for r in rows} == set(range(count))
          and all(r["label_asym_id"] == chain for r in rows), "ATOM_TABLE_MISMATCH")
    return {r["ccd_atom_id"]: r for r in rows}


def molecule_row(snapshot, cid, c, table, stem, map_name):
    identity = copy.deepcopy(c["identity"])
    version = sha(identity["canonical_isomeric_smiles"].encode("utf-8"))
    molecule_id = cid + ":" + version
    check(c.get("molecule_id", molecule_id) == molecule_id, "MOLECULE_VERSION_MISMATCH")
    attachment = c["attachment_atom"]
    check(attachment in table, "ATTACHMENT_ATOM_MISSING")
    check(snapshot.files.get(stem + ".smi", b"").decode("utf-8").strip()
          == identity["canonical_isomeric_smiles"], "SMILES_ARTIFACT_MISMATCH")
    return {"compound_id": cid, "role": "starting_ligand" if cid == "START" else "known_reference_protac",
            "paper": {"doi": snapshot.config["doi"], "name": c["paper_name"],
                      "compound_number": c["paper_compound_number"]},
            "molecule_id": molecule_id, "molecule_version": version, "identity": identity,
            "structure": {"pdb": c["pdb"], "ccd": c["ccd"], "ligand_label_asym_id": c["ligand_label_asym_id"]},
            "attachment_atom": copy.deepcopy(table[attachment]),
            "files": {"sdf_2d": snapshot.file_ref(stem + ".sdf"),
                      "depiction_2d": snapshot.file_ref(stem + ".svg"),
                      "atom_map": snapshot.file_ref(map_name)},
            "review_status": "pending", "prediction_status": "not_run"}


def build_material(snapshot: Snapshot) -> dict:
    h, config = snapshot.handoff, snapshot.config
    start = h["starting_ligand"]
    check(all(start[k] == v for k, v in config["starting_ligand"].items()), "START_IDENTITY_MISMATCH")
    start_atoms = atom_table(snapshot.json("starting-atom-map.json"), start["identity"]["heavy_atom_count"],
                             start["ligand_label_asym_id"])
    compounds = [molecule_row(snapshot, "START", start, start_atoms, "starting-ligand", "starting-atom-map.json")]
    expected = {c["id"]: c for c in config["candidates"]}
    check(len(expected) == len(config["candidates"])
          and len(h["candidates"]) == len(expected)
          and {c["id"] for c in h["candidates"]} == set(expected), "CANDIDATE_SET_MISMATCH")
    reconstructions = []
    for c in h["candidates"]:
        cid = c["id"]
        check(all(c[k] == v for k, v in expected[cid].items()), "CANDIDATE_IDENTITY_MISMATCH")
        maps = snapshot.json(cid + "-atom-map.json")
        check(maps["molecule_id"] == c["molecule_id"], "ATOM_MOLECULE_MISMATCH")
        atoms = atom_table(maps["reference_atoms"], c["identity"]["heavy_atom_count"], c["ligand_label_asym_id"])
        order = maps["assembled_sdf_atom_order"]
        check(len(order) == len(atoms) and {r["sdf_index_one_based"] for r in order} == set(range(1, len(atoms)+1))
              and {r["ccd_atom_id"] for r in order} == set(atoms)
              and all(atoms[r["ccd_atom_id"]]["atom_map"] == r["atom_map"] for r in order), "SDF_ATOM_MAP_MISMATCH")
        mapping = maps["warhead_mapping"]
        pairs = mapping["atoms"]
        check(len(pairs) == len(start_atoms) and {r["source_ccd_atom"] for r in pairs} == set(start_atoms)
              and len({r["candidate_ccd_atom"] for r in pairs}) == len(pairs)
              and all(r["candidate_ccd_atom"] in atoms and atoms[r["candidate_ccd_atom"]]["atom_map"] == r["candidate_map"] for r in pairs),
              "WARHEAD_MAP_MISMATCH")
        check(any(r["source_ccd_atom"] == start["attachment_atom"] and r["candidate_ccd_atom"] == c["attachment_atom"] for r in pairs),
              "ATTACHMENT_MAPPING_MISMATCH")
        part = snapshot.json(cid + "-parts.json")
        check(part["reference_ccd"] == c["ccd"] and part["cuts"] == c["cuts"]
              and part["reference_identity_match"] is True, "PARTS_IDENTITY_MISMATCH")
        row = molecule_row(snapshot, cid, c, atoms, cid, cid + "-atom-map.json")
        row["attachment_atom"]["assembled_sdf_index_one_based"] = next(r["sdf_index_one_based"] for r in order if r["ccd_atom_id"] == c["attachment_atom"])
        compounds.append(row)
        reconstructions.append({"compound_id": cid, "molecule_id": c["molecule_id"],
            "warhead_mapping": copy.deepcopy(mapping), "boundary_bonds": copy.deepcopy(c["cuts"]),
            "parts": snapshot.file_ref(cid + "-parts.json"),
            "meaning": "Known reference cut/rejoin; matching-only charge normalization does not modify the source ligand or establish synthesis.",
            "human_review": "pending"})
    observations = []
    for index, obs in enumerate(snapshot.literature["observations"]):
        c = expected.get(obs["candidate_id"])
        check(c is not None and all(obs[k] == c[k] for k in ("paper_name", "paper_compound_number")), "OBSERVATION_COMPOUND_MISMATCH")
        check(obs["review_status"] == "pending" and ((obs["kind"] == "missing" and obs["value"] is None)
              or (obs["kind"] == "literature_observation" and obs["value"] is not None and obs["locator"])), "OBSERVATION_KIND_MISMATCH")
        observations.append({"observation_id": f'{obs["candidate_id"]}:observation:{index+1}',
                             "original": copy.deepcopy(obs),
                             "citation": snapshot.citation(obs["locator"]) if obs["locator"] else None})
    geometry = snapshot.json("attachment-geometry.json")
    check(geometry["attachment_ccd_atom"] == start["attachment_atom"]
          and geometry["ligand_label_asym_id"] == start["ligand_label_asym_id"]
          and geometry["protein_label_asym_id"] == start["protein_label_asym_id"], "GEOMETRY_ATOM_MISMATCH")
    material = {"format": "tpd-b-review-material/0.1.0-draft", "case_id": h["case_id"],
        "case_version": h["case_version"], "data_mode": "real",
        "source_snapshot": {"handoff_sha256": sha(snapshot.handoff_bytes), "case_input_sha256": h["input_sha256"],
                            "b_run_id": h["run_id"], "raw_sha256": {k: sha(v) for k,v in sorted(snapshot.raw.items())}},
        "h1": {"target": copy.deepcopy(h["target"]), "compounds": compounds, "observations": observations,
               "provenance": snapshot.literature["authoring"], "coverage": snapshot.literature["coverage"]},
        "h2": {"hypothesis_id": "START:" + start["attachment_atom"], "subject_molecule_id": compounds[0]["molecule_id"],
               "attachment_literature": {"summary": snapshot.literature["attachment"]["summary"],
                    "citation": snapshot.citation(snapshot.literature["attachment"]["locator"]), "review_status": "pending"},
               "descriptive_geometry": copy.deepcopy(geometry), "reference_reconstructions": reconstructions,
               "proposed_scope": {"reference_compound_ids": sorted(expected), "new_analog_design": False,
                                  "recruiter": "VHL-binding reference fragments", "scope_status": "proposal_pending_G1"},
               "unresolved": ["named_pharmaceutical_review", "charge_normalization_review", "symmetry_mapping_review",
                              "literature_observation_conditions", "supplementary_synthesis_review",
                              "A_consumer_confirmation", "M2_G1_integration"]},
        "authority": {"human_review": "pending", "approval_record_created": False, "dispatch_authorized": False,
                      "scientific_validation": "not_established_by_handoff_checks"}}
    material["material_digest"] = sha(encoded(material))
    validate_material(material)
    return material


def validate_material(material):
    schema = parse(SCHEMA.read_bytes())
    boundary = parse((ACTIVE / "contracts/v0.1.0/module_boundary.schema.json").read_bytes())
    registry = Registry().with_resource("urn:tpd-navigator:module-boundary:0.1.0", Resource.from_contents(boundary))
    check(not list(Draft202012Validator(schema, registry=registry).iter_errors(material)), "REVIEW_MATERIAL_SCHEMA")
    unsigned = {k:v for k,v in material.items() if k != "material_digest"}
    check(sha(encoded(unsigned)) == material["material_digest"], "MATERIAL_DIGEST_MISMATCH")
    compounds = material["h1"]["compounds"]
    check(len({c["compound_id"] for c in compounds}) == len(compounds), "DUPLICATE_COMPOUND")
    check(sum(c["role"] == "starting_ligand" for c in compounds) == 1, "STARTING_LIGAND_REQUIRED")
    for c in compounds:
        version = sha(c["identity"]["canonical_isomeric_smiles"].encode("utf-8"))
        check(c["molecule_version"] == version and c["molecule_id"] == c["compound_id"] + ":" + version,
              "MOLECULE_VERSION_MISMATCH")
    by_id = {c["compound_id"]: c for c in compounds}
    start = next(c for c in compounds if c["role"] == "starting_ligand")
    h2 = material["h2"]
    check(h2["subject_molecule_id"] == start["molecule_id"]
          and h2["descriptive_geometry"]["attachment_ccd_atom"] == start["attachment_atom"]["ccd_atom_id"],
          "HYPOTHESIS_SUBJECT_MISMATCH")
    expected = {c["compound_id"] for c in compounds if c["role"] == "known_reference_protac"}
    refs = h2["reference_reconstructions"]
    check({r["compound_id"] for r in refs} == expected and len(refs) == len(expected)
          and set(h2["proposed_scope"]["reference_compound_ids"]) == expected, "REFERENCE_SCOPE_MISMATCH")
    for r in refs:
        check(r["molecule_id"] == by_id[r["compound_id"]]["molecule_id"], "REFERENCE_VERSION_MISMATCH")
    obs_ids = set()
    for r in material["h1"]["observations"]:
        o = r["original"]
        c = by_id.get(o["candidate_id"])
        check(c is not None and c["paper"]["name"] == o["paper_name"]
              and c["paper"]["compound_number"] == o["paper_compound_number"], "OBSERVATION_COMPOUND_MISMATCH")
        check(r["observation_id"] not in obs_ids, "DUPLICATE_OBSERVATION")
        obs_ids.add(r["observation_id"])
        if o["kind"] == "missing":
            check(o["value"] is None and o["relation"] is None, "MISSING_OBSERVATION_HAS_VALUE")
        else:
            check(o["value"] is not None and r["citation"] is not None, "OBSERVATION_SOURCE_REQUIRED")
        if r["citation"] is not None:
            check(r["citation"]["doi"] == c["paper"]["doi"]
                  and r["citation"]["locator"] == o["locator"]
                  and material["source_snapshot"]["raw_sha256"].get(r["citation"]["source_file"]) == r["citation"]["source_sha256"],
                  "OBSERVATION_CITATION_MISMATCH")
    return material


def match_compound(material, *, doi, paper_name=None, compound_number=None, expected_digest=None):
    """Exact article-scoped identity only. Never fuzzy-match a PROTAC number."""
    validate_material(material)
    if expected_digest is not None:
        check(expected_digest == material["material_digest"], "STALE_MATERIAL")
    check(isinstance(doi, str) and bool(doi) and (paper_name is not None or compound_number is not None), "LITERATURE_IDENTITY_REQUIRED")
    check(paper_name is None or isinstance(paper_name, str) and bool(paper_name), "LITERATURE_NAME_INVALID")
    check(compound_number is None or isinstance(compound_number, str) and bool(compound_number), "LITERATURE_NUMBER_INVALID")
    pool = [c for c in material["h1"]["compounds"] if c["paper"]["doi"] == doi]
    matches = [c for c in pool if (paper_name is None or c["paper"]["name"] == paper_name)
               and (compound_number is None or c["paper"]["compound_number"] == compound_number)]
    if len(matches) == 1:
        c = matches[0]
        return {"status": "matched", "compound_id": c["compound_id"], "molecule_id": c["molecule_id"],
                "material_digest": material["material_digest"], "human_review": "pending"}
    intersects = any(c["paper"]["name"] == paper_name or c["paper"]["compound_number"] == compound_number for c in pool)
    return {"status": "ambiguous" if len(matches) > 1 else "conflict" if intersects else "not_found",
            "compound_id": None, "molecule_id": None, "material_digest": material["material_digest"], "human_review": "pending"}


def compare_observation(material, query, observation):
    """Return a comparison; never overwrite/normalize the incoming observation."""
    match = match_compound(material, **query)
    result = {"identity_match": match, "incoming": copy.deepcopy(observation),
              "status": "identity_unresolved", "differences": [], "reference_observation_ids": []}
    if match["status"] != "matched":
        return result
    compound = next(c for c in material["h1"]["compounds"] if c["compound_id"] == match["compound_id"])
    for key, expected in (("doi", compound["paper"]["doi"]), ("candidate_id", match["compound_id"]), ("paper_name", compound["paper"]["name"]),
                          ("paper_compound_number", compound["paper"]["compound_number"])):
        if key in observation and expected is not None and observation[key] != expected:
            result["status"] = "identity_conflict"
            result["differences"] = [key]
            return result
    references = [r for r in material["h1"]["observations"] if r["original"]["candidate_id"] == match["compound_id"]
                  and r["original"]["endpoint"] == observation.get("endpoint")]
    result["reference_observation_ids"] = [r["observation_id"] for r in references]
    if len(references) != 1 or references[0]["original"]["kind"] == "missing":
        result["status"] = "reference_missing_requires_review"
        return result
    reference = references[0]["original"]
    fields = ("kind", "value", "unit", "relation", "target", "cell_line", "assay", "biological_replicates", "exposure_time_h", "locator")
    result["differences"] = [k for k in fields if k not in observation or observation[k] != reference.get(k)]
    result["status"] = "difference_requires_review" if result["differences"] else "same_record_pending_review"
    return result


def read_material(reference, read_bytes, *, expected_digest=None):
    """Consume an existing ArtifactRef through an injected reader, without M2 SQL."""
    boundary = parse((ACTIVE / "contracts/v0.1.0/module_boundary.schema.json").read_bytes())
    validator = Draft202012Validator({"$ref": "#/$defs/ArtifactRef", "$defs": boundary["$defs"]})
    check(not list(validator.iter_errors(reference)), "ARTIFACT_REF_SCHEMA")
    check(reference["provenance"] == "computed" and reference["media_type"] == "application/json"
          and reference["schema_id"] == "urn:tpd-navigator:b-review-material:0.1.0-draft", "MATERIAL_REF_KIND")
    data = read_bytes(reference)
    check(sha(data) == reference["sha256"], "MATERIAL_ARTIFACT_HASH_MISMATCH")
    material = validate_material(parse(data))
    if expected_digest is not None:
        check(material["material_digest"] == expected_digest, "STALE_MATERIAL")
    return material


def compare_versions(previous, current):
    """Describe re-review impact. Only M2 can actually invalidate approval records."""
    validate_material(previous)
    validate_material(current)
    fields = ("case_id", "case_version", "source_snapshot", "h1", "h2")
    changes = [k for k in fields if previous[k] != current[k]]
    return {"previous_digest": previous["material_digest"], "current_digest": current["material_digest"],
            "changed_sections": changes, "requires_M2_reassessment": bool(changes),
            "approval_state_modified": False}
