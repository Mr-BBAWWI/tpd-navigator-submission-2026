#!/usr/bin/env python3
"""Freeze local acceptance artifacts into a new immutable portable bundle."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import stat
import tempfile
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any, Iterable

VERSION = "acceptance-sources/20260930.1"
TOP_FILES = (
    "warhead_catalog.json",
    "medchem_evidence.json",
    "route_evidence.json",
    "preparation.json",
    "crbn_calibration.json",
)
SAFE_EXTENSIONS = {
    ".json", ".yaml", ".yml", ".csv", ".txt", ".log", ".sdf", ".mol",
    ".pdb", ".pqr", ".cif", ".mmcif", ".pdf", ".npz", ".pkl", ".gz",
    ".a3m", ".m8", ".sh", ".pka",
}
FORBIDDEN_PARTS = {
    ".git", ".venv", ".venv-boltz", "venv", "checkpoints", "checkpoint",
    "weights", "models",
}
CHECKPOINT_EXTENSIONS = {".ckpt", ".pt", ".pth", ".safetensors", ".bin"}
CAPTURED_TEXT_BIN_NAMES = {"pdb2pqr.stdout.bin", "pdb2pqr.stderr.bin"}
SCIENTIFIC_EXTENSIONS = {".sdf", ".mol", ".pdb", ".pqr", ".cif", ".mmcif", ".npz", ".a3m", ".m8", ".pka"}
SCIENTIFIC_PATH_KEYS = {
    "path", "path_as_supplied", "relative_path", "source_path", "source_paths",
    "input_path", "output_path", "output_pdb_path", "output_pqr_path",
    "prepared_structure", "before_sdf", "after_sdf", "protein_pdb", "report",
    "prediction", "actual_source", "actual_input", "actual_output",
}
EXECUTION_PATH_KEYS = {
    "executable", "executable_path", "python", "python_path", "checkpoint",
    "checkpoint_path", "checkpoint_file", "weights_path",
}
HUMAN_TEXT_KEYS = {
    "stdout", "stderr", "stdout_tail", "stderr_tail", "text", "message",
    "raw_line", "log_tail",
}
HEX64 = re.compile(r"^[0-9a-f]{64}$")
BACKSLASH = chr(92)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(label + "_INVALID_JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(label + "_OBJECT_REQUIRED")
    return value


def safe_relative(value: str, label: str = "PATH") -> PurePosixPath:
    if not isinstance(value, str) or not value or BACKSLASH in value or ":" in value:
        raise ValueError(label + "_INVALID")
    result = PurePosixPath(value)
    if result.is_absolute() or result.as_posix() in {"", "."} or ".." in result.parts:
        raise ValueError(label + "_TRAVERSAL")
    if result.as_posix() != value or any(part in {"", "."} for part in result.parts):
        raise ValueError(label + "_INVALID")
    return result


def reject_symlink_ancestors(path: Path, label: str) -> None:
    absolute = path.absolute()
    for ancestor in (absolute, *absolute.parents):
        if ancestor.is_symlink():
            raise ValueError(label + "_SYMLINK_FORBIDDEN")


def regular(path: Path, label: str) -> Path:
    reject_symlink_ancestors(path, label)
    if not path.is_file():
        raise ValueError(label + "_REGULAR_FILE_REQUIRED")
    return path


def directory(path: Path, label: str) -> Path:
    reject_symlink_ancestors(path, label)
    if not path.is_dir():
        raise ValueError(label + "_DIRECTORY_REQUIRED")
    return path


def check_source_path(path: Path) -> None:
    reject_symlink_ancestors(path, "SOURCE")
    lowered = {part.lower() for part in path.parts}
    suffix = path.suffix.lower()
    captured_log = path.name.lower() in CAPTURED_TEXT_BIN_NAMES
    if lowered & FORBIDDEN_PARTS or (suffix in CHECKPOINT_EXTENSIONS and not captured_log):
        raise ValueError("FORBIDDEN_SOURCE:" + path.name)


def json_bytes(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


class FreezeBuilder:
    def __init__(self, output: Path) -> None:
        self.output = output.absolute()
        if self.output.exists() or self.output.is_symlink():
            raise FileExistsError("OUTPUT_ALREADY_EXISTS:" + str(self.output))
        self.output.parent.mkdir(parents=True, exist_ok=True)
        reject_symlink_ancestors(self.output.parent, "OUTPUT_WORKSPACE")
        workspace = self.output.parent.resolve(strict=True)
        if self.output.parent != workspace or self.output.name in {"", ".", ".."}:
            raise ValueError("OUTPUT_WORKSPACE_ESCAPE")
        self.scratch = Path(tempfile.mkdtemp(prefix="." + self.output.name + ".failed-", dir=workspace))
        self.destinations: set[str] = set()
        self.path_bindings: list[tuple[str, str]] = []
        self.digest_bindings: dict[str, list[str]] = {}
        self.destination_digests: dict[str, str] = {}

    def destination(self, relative: str) -> Path:
        posix = safe_relative(relative, "DESTINATION")
        normalized = posix.as_posix()
        if normalized in self.destinations:
            raise ValueError("DESTINATION_COLLISION:" + normalized)
        candidate = self.scratch.joinpath(*posix.parts)
        resolved_root = self.scratch.resolve()
        if not candidate.resolve(strict=False).is_relative_to(resolved_root):
            raise ValueError("DESTINATION_TRAVERSAL")
        self.destinations.add(normalized)
        candidate.parent.mkdir(parents=True, exist_ok=True)
        return candidate

    def bind(self, source: Path, relative: str, digest: str, aliases: Iterable[str] = ()) -> None:
        resolved_source = source.resolve(strict=True)
        exact = str(resolved_source).replace(BACKSLASH, "/")
        self.path_bindings.append((exact, relative))
        cwd = Path.cwd().resolve(strict=True)
        if resolved_source.is_relative_to(cwd):
            cwd_relative = resolved_source.relative_to(cwd).as_posix()
            if cwd_relative and cwd_relative != ".":
                self.path_bindings.append((cwd_relative, relative))
        for alias in aliases:
            normalized = alias.replace(BACKSLASH, "/")
            if normalized and normalized != ".":
                self.path_bindings.append((normalized, relative))
        self.digest_bindings.setdefault(digest, []).append(relative)
        previous = self.destination_digests.get(relative)
        if previous is not None and previous != digest:
            raise ValueError("DESTINATION_DIGEST_CONFLICT:" + relative)
        self.destination_digests[relative] = digest

    def copy(self, source: Path, relative: str, aliases: Iterable[str] = ()) -> str:
        source = regular(source, "SOURCE")
        check_source_path(source)
        destination = self.destination(relative)
        with source.open("rb") as reader, destination.open("xb") as writer:
            shutil.copyfileobj(reader, writer, 1024 * 1024)
        digest = sha256_file(destination)
        if digest != sha256_file(source):
            raise ValueError("COPY_HASH_MISMATCH:" + relative)
        self.bind(source, relative, digest, aliases)
        return relative

    def write_json(self, relative: str, value: Any) -> None:
        destination = self.destination(relative)
        destination.write_bytes(json_bytes(value))

    def copy_tree(
        self,
        source: Path,
        prefix: str,
        predicate: Any | None = None,
        binding_root: Path | None = None,
    ) -> list[str]:
        source = directory(source, "SOURCE_TREE")
        binding_root = directory(binding_root, "BINDING_ROOT") if binding_root is not None else source
        copied: list[str] = [];
        for root, dirs, files in os.walk(source, topdown=True, followlinks=False):
            root_path = Path(root)
            for name in dirs:
                child = root_path / name
                check_source_path(child)
                if child.is_symlink():
                    raise ValueError("SOURCE_SYMLINK_FORBIDDEN")
            dirs[:] = sorted(name for name in dirs if name.lower() not in FORBIDDEN_PARTS and not name.startswith("."))
            for name in sorted(files):
                child = root_path / name
                check_source_path(child)
                relative = child.relative_to(source).as_posix()
                safe_relative(relative, "SOURCE_RELATIVE")
                if predicate is not None and not predicate(child, relative):
                    continue
                try:
                    alias = child.relative_to(binding_root).as_posix()
                except ValueError as exc:
                    raise ValueError("SOURCE_BINDING_ROOT_MISMATCH") from exc
                aliases = {alias, relative, child.name}
                if source.name:
                    aliases.add(source.name + "/" + relative)
                copied.append(self.copy(child, prefix.rstrip("/") + "/" + relative, aliases))
        return copied

    def find_by_digest(self, digest: str) -> str | None:
        matches = sorted(set(self.digest_bindings.get(digest, [])))
        return matches[0] if matches else None

    @staticmethod
    def is_absolute_path(value: str) -> bool:
        return Path(value).is_absolute() or PureWindowsPath(value).is_absolute()

    @staticmethod
    def is_private_path(value: str) -> bool:
        normalized = value.replace(BACKSLASH, "/")
        return (
            FreezeBuilder.is_absolute_path(value)
            or ".localdata/" in normalized
            or "/.venv" in normalized
        )

    def bound_path(self, value: str, digest: str | None = None) -> str | None:
        normalized = value.replace(BACKSLASH, "/")
        exact_candidates = {normalized}
        if not self.is_absolute_path(value):
            resolved = (Path.cwd() / Path(normalized)).resolve(strict=False)
            exact_candidates.add(str(resolved).replace(BACKSLASH, "/"))
        matches = {
            target
            for source, target in self.path_bindings
            if source in exact_candidates
            and (digest is None or self.destination_digests.get(target) == digest)
        }
        if not matches and isinstance(digest, str) and HEX64.fullmatch(digest):
            return self.find_by_digest(digest)
        if not matches:
            return None
        digests = {self.destination_digests[target] for target in matches}
        if len(digests) == 1:
            return sorted(matches)[0]
        raise ValueError("AMBIGUOUS_SOURCE_BINDING:" + value)

    @staticmethod
    def digest_for_key(parent: dict[str, Any], key: str) -> str | None:
        lowered = key.lower()
        candidates = [key + "_sha256"]
        if lowered.endswith("_path"):
            candidates.append(key[:-5] + "_sha256")
        if lowered in {"path", "path_as_supplied", "relative_path", "filename", "source_filename", "sourcefilename"}:
            candidates.extend(("sha256", "source_sha256"))
        if lowered == "output_pqr_path":
            candidates.append("output_sha256")
        for candidate in candidates:
            digest = parent.get(candidate)
            if isinstance(digest, str) and HEX64.fullmatch(digest):
                return digest
        return None

    def rewrite_string(self, value: str, digest: str | None = None) -> str:
        if HEX64.fullmatch(value) or value.startswith(("http://", "https://")):
            return value
        binding = self.bound_path(value, digest)
        if binding is not None:
            return binding
        if self.is_private_path(value):
            raise ValueError("UNRESOLVED_PRIVATE_PATH:" + value)
        return value

    def rewrite_human_text(self, value: str) -> str:
        canonical: dict[str, str] = {}
        for _, target in self.path_bindings:
            canonical[target] = target
            canonical[target.replace("/", BACKSLASH)] = target
        variants = dict(canonical)
        for source, target in self.path_bindings:
            variants.setdefault(source, target)
            variants.setdefault(source.replace("/", BACKSLASH), target)
        originals = sorted(
            (source for source in variants if source),
            key=lambda source: (source not in canonical, -len(source), source),
        )
        if originals:
            path_character = r"A-Za-z0-9_.~:/\\-"
            pattern = re.compile(
                rf"(?<![{path_character}])(?:"
                + "|".join(re.escape(source) for source in originals)
                + rf")(?![{path_character}])"
            )
            result = pattern.sub(lambda match: variants[match.group(0)], value)
        else:
            result = value
        normalized = result.replace(BACKSLASH, "/")
        if ".localdata/" in normalized or "/.venv" in normalized:
            return "refer to raw log preserved"
        for index in range(max(0, len(normalized) - 2)):
            if normalized[index:index + 3].endswith(":/") and normalized[index].isalpha():
                return "refer to raw log preserved"
        return result

    @staticmethod
    def external_environment_ref(value: str, digest: str | None = None) -> dict[str, Any]:
        normalized = value.replace(BACKSLASH, "/").strip()
        token = normalized.split()[0] if normalized else "external"
        basename = PureWindowsPath(token).name if PureWindowsPath(token).name else PurePosixPath(token).name
        result: dict[str, Any] = {
            "status": "not_packaged_execution_environment",
            "basename": basename or "external",
            "hash_if_known": None,
        }
        if isinstance(digest, str) and HEX64.fullmatch(digest):
            result["hash_if_known"] = digest
        return result

    def rewrite_command_item(self, value: str, index: int) -> Any:
        if index == 0:
            return self.external_environment_ref(value)
        if "=" in value:
            flag, argument = value.split("=", 1)
            binding = self.bound_path(argument)
            if binding is not None:
                return flag + "=" + binding
            if self.is_private_path(argument):
                raise ValueError("UNRESOLVED_PRIVATE_PATH:" + argument)
            return value
        binding = self.bound_path(value)
        if binding is not None:
            return binding
        if self.is_private_path(value):
            raise ValueError("UNRESOLVED_PRIVATE_PATH:" + value)
        return value

    def rewrite(self, value: Any, key: str | None = None, parent: dict[str, Any] | None = None) -> Any:
        if isinstance(value, dict):
            result: dict[str, Any] = {}
            for child_key, item in value.items():
                lowered = child_key.lower()
                if isinstance(item, str) and HEX64.fullmatch(item):
                    result[child_key] = item
                elif lowered == "command" and isinstance(item, list):
                    result[child_key] = [
                        self.rewrite_command_item(command_item, index)
                        if isinstance(command_item, str)
                        else self.rewrite(command_item)
                        for index, command_item in enumerate(item)
                    ]
                elif isinstance(item, str) and lowered in EXECUTION_PATH_KEYS:
                    if lowered.startswith(("checkpoint", "weights")):
                        digest = value.get("checkpoint_sha256")
                    elif lowered in {"executable", "executable_path"}:
                        digest = value.get("executable_sha256")
                    else:
                        digest = value.get(child_key + "_sha256")
                    result[child_key] = self.external_environment_ref(
                        item, digest if isinstance(digest, str) else None
                    )
                elif isinstance(item, str) and lowered in HUMAN_TEXT_KEYS:
                    result[child_key] = self.rewrite_human_text(item)
                elif isinstance(item, str):
                    result[child_key] = self.rewrite_string(
                        item, self.digest_for_key(value, child_key)
                    )
                else:
                    result[child_key] = self.rewrite(item, child_key, value)
            return result
        if isinstance(value, list):
            return [self.rewrite(item, key, parent) for item in value]
        if isinstance(value, str):
            return self.rewrite_string(value)
        return value

    def assert_scientific_refs_packaged(self, value: Any, key: str | None = None) -> None:
        if isinstance(value, dict):
            if value.get("status") == "not_packaged_execution_environment":
                return
            for child_key, item in value.items():
                self.assert_scientific_refs_packaged(item, child_key)
            return
        if isinstance(value, list):
            for item in value:
                self.assert_scientific_refs_packaged(item, key)
            return
        if not isinstance(value, str) or HEX64.fullmatch(value) or value.startswith(("http://", "https://")):
            return
        candidate_value = value.split("=", 1)[1] if key == "command" and "=" in value else value
        normalized = candidate_value.replace(BACKSLASH, "/")
        if key != "command" and (key is None or key.lower() not in SCIENTIFIC_PATH_KEYS):
            return
        if PurePosixPath(normalized).suffix.lower() not in SCIENTIFIC_EXTENSIONS:
            return
        relative = safe_relative(normalized, "SCIENTIFIC_REFERENCE")
        candidate = self.scratch.joinpath(*relative.parts)
        if (
            normalized not in self.destination_digests
            or not candidate.is_file()
            or candidate.is_symlink()
        ):
            raise ValueError("UNRESOLVED_SCIENTIFIC_DATA_FILE:" + candidate_value)

    def finish(self) -> None:
        files: dict[str, str] = {}
        for path in sorted(self.scratch.rglob("*")):
            if path.is_symlink():
                raise ValueError("SCRATCH_SYMLINK")
            if path.is_file():
                relative = path.relative_to(self.scratch).as_posix()
                if relative != "manifest.json":
                    files[relative] = sha256_file(path)
        missing = sorted(set(TOP_FILES) - set(files))
        if missing:
            raise ValueError("TOP_FILES_MISSING:" + ",".join(missing))
        self.write_json("manifest.json", {"files": files, "version": VERSION})
        for path in sorted(self.scratch.rglob("*"), reverse=True):
            if path.is_file():
                path.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
            elif path.is_dir():
                path.chmod(stat.S_IRUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
        self.scratch.chmod(0o555)
        os.replace(self.scratch, self.output)

    def fail(self, exc: BaseException) -> None:
        try:
            marker = self.scratch / "BUILD_FAILED.txt"
            marker.write_text(f"{exc.__class__.__name__}: {exc}\n", encoding="utf-8")
        except OSError:
            pass


def allowed_support(path: Path, relative: str) -> bool:
    if any(part.startswith(".") for part in PurePosixPath(relative).parts):
        return False
    if path.name.lower() in CAPTURED_TEXT_BIN_NAMES:
        return True
    return path.suffix.lower() in SAFE_EXTENSIONS and path.suffix.lower() not in CHECKPOINT_EXTENSIONS


def locate_json(root: Path, names: Iterable[str], label: str) -> Path:
    for name in names:
        candidate = root / name
        if candidate.is_file() and not candidate.is_symlink():
            return candidate
    raise ValueError(label + "_MISSING")


def freeze_catalog(builder: FreezeBuilder, root: Path) -> dict[str, Any]:
    root = directory(root, "CATALOG_DIR")
    catalog_path = locate_json(root, ("warhead_catalog.json", "catalog.json"), "WARHEAD_CATALOG")
    raw = catalog_path.read_bytes()
    builder.copy(catalog_path, "supporting/catalog/" + catalog_path.name)
    for dirname in ("source_snapshots", "bound_ligands"):
        child = root / dirname
        if child.exists():
            builder.copy_tree(child, "supporting/catalog/" + dirname, allowed_support, binding_root=root)
    for child in sorted(root.iterdir()):
        if child.is_file() and child != catalog_path and allowed_support(child, child.name):
            builder.copy(child, "supporting/catalog/" + child.name, (child.name,))
    document = json.loads(raw.decode("utf-8"))
    if not isinstance(document, dict):
        raise ValueError("WARHEAD_CATALOG_OBJECT_REQUIRED")
    document = builder.rewrite(document)
    document["sources_root"] = "supporting/catalog/source_snapshots"
    for artifact in document.get("artifacts", []):
        if isinstance(artifact, dict) and isinstance(artifact.get("relative_path"), str):
            path = artifact["relative_path"]
            if not path.startswith("supporting/"):
                artifact["relative_path"] = "supporting/catalog/" + safe_relative(path, "CATALOG_ARTIFACT").as_posix()
    for source in document.get("sources", []):
        if isinstance(source, dict) and isinstance(source.get("relative_path"), str):
            path = source["relative_path"]
            if not path.startswith("supporting/"):
                source["relative_path"] = "supporting/catalog/source_snapshots/" + safe_relative(path, "CATALOG_SOURCE").as_posix()
    document["freeze_provenance"] = {
        "original_content_sha256": sha256_bytes(raw),
        "original_ref": "supporting/catalog/" + catalog_path.name,
        "evidence_summary": "Source-bound evidence only; no compound, route, efficacy, degradation, or design approval is granted.",
    }
    return document


def freeze_medchem(builder: FreezeBuilder, root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    root = directory(root, "MEDCHEM_DIR")
    evidence_path = locate_json(root, ("medchem_evidence.json",), "MEDCHEM_EVIDENCE")
    route_path = locate_json(root, ("route_evidence.json",), "ROUTE_EVIDENCE")
    evidence_raw = evidence_path.read_bytes()
    route_raw = route_path.read_bytes()
    builder.copy_tree(root, "supporting/medchem", allowed_support, binding_root=root)
    evidence = builder.rewrite(json.loads(evidence_raw.decode("utf-8")))
    routes = builder.rewrite(json.loads(route_raw.decode("utf-8")))
    evidence["freeze_provenance"] = {
        "original_content_sha256": sha256_bytes(evidence_raw),
        "original_ref": "supporting/medchem/" + evidence_path.name,
        "evidence_summary": "Measured values retain their reported endpoint and units; computed properties, including computed pKa if present, are not represented as measured values.",
    }
    routes["freeze_provenance"] = {
        "original_content_sha256": sha256_bytes(route_raw),
        "original_ref": "supporting/medchem/" + route_path.name,
        "evidence_summary": "Reported source routes are precedents only and do not approve synthesis of a new design.",
    }
    return evidence, routes


def preparation_data_paths(value: Any, key: str | None = None) -> Iterable[str]:
    if isinstance(value, dict):
        for child_key, item in value.items():
            if child_key.lower() == "command" and isinstance(item, list):
                for index, command_item in enumerate(item):
                    if index == 0 or not isinstance(command_item, str):
                        continue
                    yield command_item.split("=", 1)[1] if "=" in command_item else command_item
            else:
                yield from preparation_data_paths(item, child_key)
        return
    if isinstance(value, list):
        for item in value:
            yield from preparation_data_paths(item, key)
        return
    if isinstance(value, str) and key is not None and key.lower() in SCIENTIFIC_PATH_KEYS:
        yield value


def freeze_preparation(builder: FreezeBuilder, path: Path) -> dict[str, Any]:
    path = regular(path, "PREPARATION")
    root = directory(path.parent, "PREPARATION_DIR")
    raw = path.read_bytes()
    builder.copy_tree(root, "supporting/preparation", allowed_support, binding_root=root)
    document = json.loads(raw.decode("utf-8"))
    if not isinstance(document, dict):
        raise ValueError("PREPARATION_OBJECT_REQUIRED")
    for supplied in preparation_data_paths(document):
        normalized = supplied.replace(BACKSLASH, "/")
        if PurePosixPath(normalized).suffix.lower() not in SCIENTIFIC_EXTENSIONS:
            continue
        if builder.bound_path(supplied) is not None:
            continue
        source = Path(supplied)
        if not builder.is_absolute_path(supplied) or not source.is_file():
            continue
        source = regular(source, "PREPARATION_EXTERNAL_DATA")
        digest = sha256_file(source)
        existing = builder.find_by_digest(digest)
        if existing is not None:
            builder.bind(source, existing, digest)
            continue
        builder.copy(
            source,
            "supporting/preparation/external/" + digest + "/" + source.name,
        )
    document = builder.rewrite(document)
    document["freeze_provenance"] = {
        "original_content_sha256": sha256_bytes(raw),
        "original_ref": "supporting/preparation/" + path.name,
        "evidence_summary": "Reference-specific preparation and comparison evidence only; it is not a preparation or approval for another ligand.",
    }
    return document


def find_hashed_file(root: Path, expected: str, suffixes: set[str], label: str) -> Path:
    matches = []
    for child in root.rglob("*"):
        if child.is_symlink():
            raise ValueError("SOURCE_SYMLINK_FORBIDDEN")
        if child.is_file() and child.suffix.lower() in suffixes and sha256_file(child) == expected:
            matches.append(child)
    if len(matches) != 1:
        raise ValueError(label + "_HASHED_FILE_COUNT")
    return matches[0]


def model_version(value: Any) -> tuple[int, int, int]:
    if value != "2.2.1":
        raise ValueError("BOLTZ_VERSION_NOT_CAPTURED_2_2_1")
    return (2, 2, 1)


def verified_seed_application(receipt: dict[str, Any]) -> int:
    seed_value = receipt.get("seed")
    application = receipt.get("seed_application")
    if isinstance(seed_value, dict):
        application = seed_value.get("application")
        seed_value = seed_value.get("requested", seed_value.get("seed"))
    if type(seed_value) is not int or not isinstance(application, dict):
        raise ValueError("SEED_APPLICATION_METADATA_MISSING")
    if application.get("actually_applied") is not True:
        raise ValueError("SEED_NOT_ACTUALLY_APPLIED")
    actual = application.get("actual_seed", application.get("seed"))
    if type(actual) is not int or actual != seed_value:
        raise ValueError("SEED_APPLICATION_MISMATCH")
    return seed_value


def metric_assessment(receipt: dict[str, Any], reinspection: dict[str, Any] | None, receipt_digest: str) -> tuple[dict[str, Any], str | None]:
    if receipt.get("status") == "success":
        inspection = receipt.get("inspection")
        if not isinstance(inspection, dict) or inspection.get("status") != "success":
            raise ValueError("SUCCESS_RECEIPT_INSPECTION_INVALID")
        return inspection, None
    if reinspection is None:
        raise ValueError("FAILED_RECEIPT_REINSPECTION_REQUIRED")
    if reinspection.get("original_receipt_sha256") != receipt_digest:
        raise ValueError("REINSPECTION_ORIGINAL_RECEIPT_HASH_MISMATCH")
    if reinspection.get("actual_execution_success") is not True or reinspection.get("structural_assessment_completed") is not True:
        raise ValueError("REINSPECTION_NOT_VALID")
    assessment = reinspection.get("assessment")
    if not isinstance(assessment, dict) or assessment.get("status") != "success":
        raise ValueError("REINSPECTION_ASSESSMENT_INVALID")
    return assessment, "reinspection.json"


def freeze_seed(builder: FreezeBuilder, seed_dir: Path, reference_digest: str) -> dict[str, Any]:
    seed_dir = directory(seed_dir, "SEED_DIR")
    receipt_path = regular(seed_dir / "receipt.json", "SEED_RECEIPT")
    receipt_bytes = receipt_path.read_bytes()
    receipt_digest = sha256_bytes(receipt_bytes)
    receipt = read_json(receipt_path, "SEED_RECEIPT")
    seed = verified_seed_application(receipt)
    if receipt.get("exit_code") != 0 or receipt.get("process_state") != "completed":
        raise ValueError("SEED_EXECUTION_NOT_SUCCESSFUL")
    environment = receipt.get("tool_environment")
    if not isinstance(environment, dict):
        raise ValueError("TOOL_ENVIRONMENT_MISSING")
    version_text = environment.get("boltz_version")
    model_version(version_text)
    prefix = f"supporting/seeds/seed-{seed}"
    receipt_ref = builder.copy(receipt_path, prefix + "/receipt.json")
    hashes = receipt.get("hashes")
    if not isinstance(hashes, dict):
        raise ValueError("SEED_HASHES_MISSING")
    input_digest = hashes.get("input_yaml_sha256")
    if not isinstance(input_digest, str) or not HEX64.fullmatch(input_digest):
        raise ValueError("INPUT_HASH_INVALID")
    yaml_path = find_hashed_file(seed_dir, input_digest, {".yaml", ".yml"}, "INPUT_YAML")
    builder.copy(yaml_path, prefix + "/" + yaml_path.name)
    if hashes.get("reference_sha256") != reference_digest:
        raise ValueError("REFERENCE_HASH_MISMATCH")
    log_digest = hashes.get("inference_log_sha256")
    if not isinstance(log_digest, str) or not HEX64.fullmatch(log_digest):
        raise ValueError("INFERENCE_LOG_HASH_INVALID")
    log_path = find_hashed_file(seed_dir, log_digest, {".log", ".txt"}, "INFERENCE_LOG")
    builder.copy(log_path, prefix + "/" + log_path.name)

    output_root = directory(seed_dir / "boltz_output", "BOLTZ_OUTPUT")
    expected_outputs = hashes.get("processed_output_files_sha256")
    if not isinstance(expected_outputs, dict) or not expected_outputs:
        raise ValueError("OUTPUT_HASHES_MISSING")
    normalized_outputs: dict[str, str] = {}
    for relative, digest in expected_outputs.items():
        relative_path = safe_relative(relative, "SEED_OUTPUT")
        if not isinstance(digest, str) or not HEX64.fullmatch(digest):
            raise ValueError("SEED_OUTPUT_HASH_INVALID")
        source = output_root.joinpath(*relative_path.parts)
        regular(source, "SEED_OUTPUT")
        if sha256_file(source) != digest:
            raise ValueError("SEED_OUTPUT_HASH_MISMATCH:" + relative)
        normalized_outputs[relative] = digest
        builder.copy(source, prefix + "/boltz_output/" + relative)

    current = set()
    for child in output_root.rglob("*"):
        if child.is_symlink():
            raise ValueError("SOURCE_SYMLINK_FORBIDDEN")
        if child.is_file():
            current.add(child.relative_to(output_root).as_posix())
    if current != set(normalized_outputs):
        raise ValueError("SEED_OUTPUT_FILE_SET_MISMATCH")

    reinspection_path = seed_dir / "reinspection.json"
    reinspection = None
    reinspection_ref = None
    if reinspection_path.exists():
        reinspection = read_json(regular(reinspection_path, "REINSPECTION"), "REINSPECTION")
        verified = reinspection.get("verified_processed_output_files_sha256")
        if verified is not None and verified != normalized_outputs:
            raise ValueError("REINSPECTION_OUTPUT_HASHES_MISMATCH")
        reinspection_ref = builder.copy(reinspection_path, prefix + "/reinspection.json")
    assessment, required_reinspection = metric_assessment(receipt, reinspection, receipt_digest)
    if required_reinspection and reinspection_ref is None:
        raise ValueError("REINSPECTION_REFERENCE_MISSING")
    metrics = assessment.get("comparison", {}).get("metrics")
    if not isinstance(metrics, dict):
        raise ValueError("SEED_METRICS_MISSING")
    required_metrics = (
        "ligand_heavy_atom_RMSD_after_target_alignment_A",
        "e3_CA_RMSD_after_target_alignment_A",
        "contact_jaccard",
    )
    if any(type(metrics.get(name)) not in (int, float) for name in required_metrics):
        raise ValueError("SEED_METRICS_INVALID")
    ligand_value = float(metrics[required_metrics[0]])
    e3_value = float(metrics[required_metrics[1]])
    contact_value = float(metrics[required_metrics[2]])
    if not math.isfinite(ligand_value) or ligand_value < 0 or not math.isfinite(e3_value) or e3_value < 0:
        raise ValueError("SEED_RMSD_INVALID")
    if not math.isfinite(contact_value) or not 0.0 <= contact_value <= 1.0:
        raise ValueError("SEED_CONTACT_JACCARD_INVALID")
    processed_msa = {
        PurePosixPath(path).name: digest
        for path, digest in normalized_outputs.items()
        if "/processed/msa/" in "/" + path
    }
    if not processed_msa:
        raise ValueError("PROCESSED_MSA_HASHES_MISSING")
    return {
        "seed": seed,
        "original_receipt_ref": receipt_ref,
        "original_receipt_sha256": receipt_digest,
        "execution_success": True,
        "seed_application": builder.rewrite(receipt.get("seed_application", receipt.get("seed", {}).get("application") if isinstance(receipt.get("seed"), dict) else None)),
        "inspection": assessment,
        "settings": builder.rewrite(receipt.get("settings")),
        "hashes": hashes,
        "tool_environment": builder.rewrite(environment),
        "reinspection_ref": reinspection_ref,
        "_metric_values": {
            "ligand": ligand_value,
            "e3": e3_value,
            "contact": contact_value,
        },
        "_processed_msa": processed_msa,
        "_msa_mode": receipt.get("input", {}).get("msa_mode", receipt.get("settings", {}).get("msa_mode")),
        "_input_sha256": input_digest,
        "_checkpoint_sha256": hashes.get("checkpoint_sha256"),
        "_boltz_version": version_text,
    }


def aggregate_crbn(rows: list[dict[str, Any]]) -> dict[str, Any]:
    rows.sort(key=lambda row: row["seed"])
    if len(rows) != 3 or len({row["seed"] for row in rows}) != 3:
        raise ValueError("EXACTLY_THREE_DISTINCT_SEEDS_REQUIRED")
    ligand = [row["_metric_values"]["ligand"] for row in rows]
    e3 = [row["_metric_values"]["e3"] for row in rows]
    contacts = [row["_metric_values"]["contact"] for row in rows]
    msa_sets = [row["_processed_msa"] for row in rows]
    modes = [row["_msa_mode"] for row in rows]
    input_hashes = {str(row["seed"]): row["_input_sha256"] for row in rows}
    checkpoint_hashes = {str(row["seed"]): row["_checkpoint_sha256"] for row in rows}
    msa_hashes = {str(row["seed"]): row["_processed_msa"] for row in rows}
    versions = {str(row["seed"]): row["_boltz_version"] for row in rows}
    if any(not isinstance(value, str) or not HEX64.fullmatch(value) for value in checkpoint_hashes.values()):
        raise ValueError("CHECKPOINT_HASH_UNKNOWN")
    if len(set(checkpoint_hashes.values())) != 1:
        raise ValueError("CHECKPOINT_HASH_MISMATCH_ACROSS_SEEDS")
    if len(set(versions.values())) != 1:
        raise ValueError("BOLTZ_VERSION_MISMATCH_ACROSS_SEEDS")
    same_inputs = len(set(input_hashes.values())) == 1
    same_msa = all(item == msa_sets[0] for item in msa_sets[1:])
    controlled_replicates = same_inputs and same_msa
    for row in rows:
        del row["_metric_values"]
        del row["_processed_msa"]
        del row["_msa_mode"]
        del row["_input_sha256"]
        del row["_checkpoint_sha256"]
        del row["_boltz_version"]
    return {
        "seed_receipts": rows,
        "metrics": {
            "seed_count": 3,
            "ligand_rmsd_values_A": ligand,
            "e3_rmsd_values_A": e3,
            "ligand_rmsd_mean_A": sum(ligand) / 3,
            "ligand_rmsd_min_A": min(ligand),
            "ligand_rmsd_max_A": max(ligand),
            "e3_rmsd_mean_A": sum(e3) / 3,
            "e3_rmsd_min_A": min(e3),
            "e3_rmsd_max_A": max(e3),
            "contact_jaccard_values": contacts,
            "contact_jaccard_mean": sum(contacts) / 3,
            "contact_jaccard_min": min(contacts),
            "contact_jaccard_max": max(contacts),
        },
        "reproduction_status": (
            "variable_known_case_reproduction_requires_review"
            if controlled_replicates
            else "uncontrolled_nonreplicate_runs_requires_review"
        ),
        "replicate_claim": controlled_replicates,
        "run_comparability": {
            "boltz_versions_by_seed": versions,
            "checkpoint_sha256_by_seed": checkpoint_hashes,
            "input_yaml_sha256_by_seed": input_hashes,
            "same_input_yaml_sha256": same_inputs,
            "uncontrolled_reason": None if controlled_replicates else "Input YAML or processed MSA hashes differ; these runs are not claimed as replicates.",
        },
        "new_MSA_baseline": {
            "mode": modes[0] if len(set(modes)) == 1 else "mixed",
            "new_vs_VHL": True,
            "same_processed_msa_hashes_across_seeds": same_msa,
            "processed_msa_hashes_by_seed_verified": True,
            "processed_msa_hashes_by_seed": msa_hashes,
        },
        "limitations": [
            "DDB1 excluded from the calibration input.",
            "Known-structure training overlap may exist.",
            "One known case across three seeds is not validation.",
            "These metrics must not be used for cross-E3 or cross-design ranking.",
            "The poor replicate must be retained and reviewed rather than discarded.",
        ],
        "evidence_summary": "Hash-verified known-case reproduction evidence only; no efficacy, degradation, E3 superiority, design, or full acceptance approval is granted.",
    }


def add_extra_sources(builder: FreezeBuilder, sources: list[Path]) -> None:
    names: set[str] = set()
    for index, source in enumerate(sources, 1):
        source = regular(source, "EXTRA_SOURCE")
        if source.name in names:
            raise ValueError("EXTRA_SOURCE_BASENAME_COLLISION:" + source.name)
        names.add(source.name)
        if not allowed_support(source, source.name):
            raise ValueError("EXTRA_SOURCE_TYPE_FORBIDDEN:" + source.name)
        builder.copy(source, f"supporting/extra/{index:02d}_{source.name}", (source.name,))


def build(args: argparse.Namespace) -> None:
    if len(args.seed_dir) != 3:
        raise ValueError("--seed-dir must be supplied exactly three times")
    builder = FreezeBuilder(args.output)
    try:
        reference = regular(args.reference, "REFERENCE")
        reference_digest = sha256_file(reference)
        builder.copy(reference, "supporting/reference/" + reference.name)
        add_extra_sources(builder, args.extra_source)
        catalog = freeze_catalog(builder, args.catalog_dir)
        medchem, routes = freeze_medchem(builder, args.medchem_dir)
        preparation = freeze_preparation(builder, args.preparation)
        seeds = [freeze_seed(builder, path, reference_digest) for path in args.seed_dir]
        crbn = aggregate_crbn(seeds)
        catalog = builder.rewrite(catalog)
        medchem = builder.rewrite(medchem)
        routes = builder.rewrite(routes)
        preparation = builder.rewrite(preparation)
        crbn = builder.rewrite(crbn)
        for document in (catalog, medchem, routes, preparation, crbn):
            builder.assert_scientific_refs_packaged(document)
        builder.write_json("warhead_catalog.json", catalog)
        builder.write_json("medchem_evidence.json", medchem)
        builder.write_json("route_evidence.json", routes)
        builder.write_json("preparation.json", preparation)
        builder.write_json("crbn_calibration.json", crbn)
        builder.finish()
    except Exception as exc:
        builder.fail(exc)
        raise


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--catalog-dir", required=True, type=Path)
    result.add_argument("--medchem-dir", required=True, type=Path)
    result.add_argument("--preparation", required=True, type=Path, help="contact_comparison.json")
    result.add_argument("--seed-dir", required=True, action="append", type=Path)
    result.add_argument("--extra-source", action="append", default=[], type=Path)
    result.add_argument("--reference", required=True, type=Path, help="local 6BOY.cif")
    result.add_argument("--output", type=Path, default=Path("cases/acceptance_sources"))
    return result


def main() -> None:
    build(parser().parse_args())


if __name__ == "__main__":
    main()
