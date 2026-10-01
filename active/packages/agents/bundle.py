"""Turn a verified B result into read-only, candidate-scoped evidence cards."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from jsonschema import Draft202012Validator

from .provider import AgentError

ACTIVE = Path(__file__).resolve().parents[2]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def json_digest(value) -> str:
    return digest(json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode())


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def build_bundle(directory: Path, raw_directory: Path) -> dict:
    """Hashes establish snapshot identity, not source truth or reviewer approval."""
    directory = directory.resolve()
    handoff_bytes = (directory / "handoff.json").read_bytes()
    handoff = json.loads(handoff_bytes)
    schema = read_json(ACTIVE / "contracts/drafts/b_case_preparation.schema.json")
    if list(Draft202012Validator(schema).iter_errors(handoff)):
        raise AgentError("B_HANDOFF_SCHEMA")
    artifacts = {}
    for entry in handoff["artifacts"]:
        name = entry["path"]
        path = (directory / name).resolve()
        if path.parent != directory or name in artifacts:
            raise AgentError("B_ARTIFACT_PATH_OR_DUPLICATE")
        content = path.read_bytes()
        if digest(content) != entry["sha256"] or len(content) != entry["bytes"]:
            raise AgentError("B_ARTIFACT_HASH_MISMATCH")
        artifacts[name] = entry
    required = {"case-input.json", "sources.json", "literature-evidence.json", "attachment-geometry.json"}
    if not required <= artifacts.keys() or artifacts["case-input.json"]["sha256"] != handoff["input_sha256"]:
        raise AgentError("B_INPUT_MANIFEST_MISMATCH")
    literature = read_json(directory / "literature-evidence.json")
    sources_manifest = read_json(directory / "sources.json")
    raw_hashes = {}
    for entry in sources_manifest:
        path = (raw_directory / entry["file"]).resolve()
        if path.parent != raw_directory.resolve() or entry["file"] in raw_hashes:
            raise AgentError("B_SOURCE_PATH_OR_DUPLICATE")
        data = path.read_bytes()
        if digest(data) != entry["sha256"] or len(data) != entry["bytes"]:
            raise AgentError("B_SOURCE_HASH_MISMATCH")
        raw_hashes[entry["file"]] = entry["sha256"]
    if raw_hashes.get(literature["source_file"]) != literature["source_sha256"]:
        raise AgentError("LITERATURE_SOURCE_MISMATCH")
    article = ET.parse(raw_directory / literature["source_file"]).getroot()
    sources = {}
    for figure_id in ("F1", "F2", "F3"):
        locator = f".//fig[@id='{figure_id}']/caption"
        caption = article.find(locator)
        if caption is None:
            raise AgentError("SOURCE_LOCATOR_MISSING")
        text = "".join(caption.itertext())
        sources["paper:" + figure_id] = {
            "kind": "article_caption", "locator": locator,
            "source_sha256": literature["source_sha256"], "doi": literature["doi"],
            "text": text, "text_sha256": digest(text.encode()),
        }
    for name in ("attachment-geometry.json", "literature-evidence.json"):
        sources["artifact:" + name] = {
            "kind": "b_artifact", "path": name, "sha256": artifacts[name]["sha256"],
            "content": read_json(directory / name),
        }
    sources["b:handoff"] = {
        "kind": "b_artifact_excerpt", "sha256": digest(handoff_bytes),
        "run_id": handoff["run_id"], "target": handoff["target"],
        "starting_ligand": handoff["starting_ligand"],
        "candidates": [{k: c[k] for k in ("id", "molecule_id", "paper_name", "paper_compound_number", "pdb", "ccd", "attachment_atom", "identity", "origin", "reconstruction", "prediction", "human_review")} for c in handoff["candidates"]],
        "not_completed": handoff["not_completed"],
    }
    facts = []

    def add(fid, candidate, kind, domain, value, refs):
        facts.append({"id": fid, "candidate_id": candidate, "kind": kind,
                      "domain": domain, "value": value, "source_ids": refs})

    add("START:attachment-literature", "START", "literature_annotation", "literature",
        {"summary": literature["attachment"]["summary"], "review_status": "pending"}, ["paper:F1"])
    add("START:geometry", "START", "computed_descriptive_geometry", "molecule",
        read_json(directory / "attachment-geometry.json"), ["artifact:attachment-geometry.json", "b:handoff"])
    candidate_ids = [c["id"] for c in handoff["candidates"]]
    if len(set(candidate_ids)) != len(candidate_ids):
        raise AgentError("DUPLICATE_CANDIDATE")
    for c in handoff["candidates"]:
        cid = c["id"]
        add(cid + ":identity", cid, "reference_reconstruction", "molecule",
            {k: c[k] for k in ("paper_name", "paper_compound_number", "molecule_id", "pdb", "ccd", "attachment_atom", "identity", "origin")}, ["b:handoff"])
        add(cid + ":prediction-status", cid, "not_run", "molecule", c["prediction"], ["b:handoff"])
    for observation in literature["observations"]:
        if observation["candidate_id"] not in candidate_ids:
            raise AgentError("OBSERVATION_CANDIDATE_UNKNOWN")
        refs = ["artifact:literature-evidence.json"]
        if observation["locator"] is not None:
            matches = [sid for sid, source in sources.items() if source.get("locator") == observation["locator"]]
            if not matches:
                raise AgentError("OBSERVATION_LOCATOR_UNKNOWN")
            refs += matches
        add(observation["candidate_id"] + ":" + observation["endpoint"], observation["candidate_id"],
            observation["kind"], "literature", observation, refs)
    if len({f["id"] for f in facts}) != len(facts):
        raise AgentError("DUPLICATE_FACT")
    bundle = {
        "format": "tpd-agent-evidence/0.1.0-experimental", "case_id": handoff["case_id"],
        "b_run_id": handoff["run_id"], "b_handoff_sha256": digest(handoff_bytes),
        "b_input_sha256": handoff["input_sha256"], "source_hashes": raw_hashes,
        "facts": facts, "sources": sources,
        "coverage": literature["coverage"], "not_completed": handoff["not_completed"],
        "human_review": "pending", "trust": "hash-verified B snapshot; not a scientific approval",
    }
    bundle["input_digest"] = json_digest(bundle)
    return bundle
