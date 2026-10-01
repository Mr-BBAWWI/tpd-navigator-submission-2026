"""Small JSON/file helpers shared by offline evidence tools, without scientific imports."""
from pathlib import Path

from jsonschema import Draft202012Validator
from referencing import Registry, Resource

from .handoff import ACTIVE, check, encoded, parse, sha


def schema_check(value, name):
    schema = parse((ACTIVE / "contracts/drafts" / name).read_bytes())
    boundary = parse((ACTIVE / "contracts/v0.1.0/module_boundary.schema.json").read_bytes())
    registry = Registry().with_resource("urn:tpd-navigator:module-boundary:0.1.0", Resource.from_contents(boundary))
    quality_schemas = ("b_structure_quality.schema.json", "b_structure_quality_report.schema.json", "b_quality_evidence.schema.json")
    related = quality_schemas if name in quality_schemas else ()
    if name == "b_review_evidence.schema.json":
        related = quality_schemas + ("b_candidate_evidence.schema.json", "b_ligand_preparation_report.schema.json", "b_ligand_contacts.schema.json")
    if related:
        for filename in related:
            definition = parse((ACTIVE / "contracts/drafts" / filename).read_bytes())
            registry = registry.with_resource(definition["$id"], Resource.from_contents(definition))
    check(not list(Draft202012Validator(schema, registry=registry).iter_errors(value)), "SCHEMA_INVALID: " + name)


def seal(value):
    return {**value, "digest": sha(encoded(value))}


def verify_seal(value):
    check(value.get("digest") == sha(encoded({k: v for k, v in value.items() if k != "digest"})), "DIGEST_MISMATCH")


def safe_relative(directory, name):
    check(isinstance(name, str) and bool(name) and "\\" not in name and ":" not in name
          and not name.startswith("/") and all(p not in ("", ".", "..") for p in name.split("/")), "INVALID_RELATIVE_PATH")
    root = Path(directory).resolve()
    path = (root / name).resolve()
    check(path.is_relative_to(root) and path != root, "PATH_OUTSIDE_DIRECTORY")
    return path


def write_new_directory(destination, files):
    destination = Path(destination)
    check(not destination.exists(), "OUTPUT_EXISTS")
    for name in files:
        safe_relative(destination, name)
    destination.mkdir(parents=True, exist_ok=False)
    for name, data in files.items():
        path = safe_relative(destination, name)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(data)


def file_index(files):
    return {name: {"sha256": sha(data), "bytes": len(data)} for name, data in sorted(files.items())}
