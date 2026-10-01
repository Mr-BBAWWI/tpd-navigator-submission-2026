"""Hash-bound loader for optional acceptance evidence.

The loader exposes source JSON unchanged. It does not infer missing values,
approve designs, or convert reference-specific preparation reports into a
preparation for another ligand.
"""
from __future__ import annotations

import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "cases" / "acceptance_sources"
MANIFEST = SOURCE / "manifest.json"
EXPECTED_FILES = (
    "warhead_catalog.json",
    "medchem_evidence.json",
    "route_evidence.json",
    "preparation.json",
    "crbn_calibration.json",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _valid_hash(value: Any) -> bool:
    return (
        type(value) is str
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _safe_path(root: Path, relative: str) -> Path:
    if type(relative) is not str or not relative:
        raise ValueError("ACCEPTANCE_SOURCE_PATH")
    if (
        "\x00" in relative
        or "\\" in relative
        or ":" in relative
        or relative.startswith(("/", "//"))
        or re.match(r"^[A-Za-z]:", relative)
    ):
        raise ValueError("ACCEPTANCE_SOURCE_PATH")

    candidate = Path(relative)
    if (
        candidate.is_absolute()
        or ".." in candidate.parts
        or "." in candidate.parts
        or candidate.as_posix() != relative
    ):
        raise ValueError("ACCEPTANCE_SOURCE_PATH")

    if root.is_symlink():
        raise ValueError("ACCEPTANCE_SOURCE_SYMLINK")

    path = root / candidate
    cursor = root
    for part in candidate.parts:
        cursor = cursor / part
        if cursor.is_symlink():
            raise ValueError("ACCEPTANCE_SOURCE_SYMLINK")

    resolved_root = root.resolve()
    resolved = path.resolve(strict=False)
    if not resolved.is_relative_to(resolved_root):
        raise ValueError("ACCEPTANCE_SOURCE_PATH")
    return path


def _declared_manifest(manifest: dict[str, Any]) -> dict[str, str]:
    files = manifest["files"]
    supporting = manifest.get("supporting_files", {})
    if supporting is None:
        supporting = {}
    if type(supporting) is not dict or set(files) & set(supporting):
        raise ValueError("ACCEPTANCE_MANIFEST_SCHEMA")
    result: dict[str, str] = {}
    for declarations in (files, supporting):
        for relative, expected_hash in declarations.items():
            _safe_path(SOURCE, relative)
            if not _valid_hash(expected_hash):
                raise ValueError("ACCEPTANCE_MANIFEST_SCHEMA")
            result[relative] = expected_hash
    return result


def load_acceptance_evidence() -> dict[str, Any]:
    """Load the fixed optional catalog, failing closed on configured tampering."""
    blank = {
        name: {
            "status": "not_configured",
            "relative_path": name,
            "sha256": None,
            "data": None,
        }
        for name in EXPECTED_FILES
    }
    if not MANIFEST.exists():
        return {
            "status": "not_configured",
            "manifest": {
                "relative_path": "cases/acceptance_sources/manifest.json",
                "status": "not_configured",
                "version": None,
                "sha256": None,
            },
            "sources": blank,
            "supporting_files": {},
            "approvals": None,
        }
    if MANIFEST.is_symlink():
        raise ValueError("ACCEPTANCE_SOURCE_SYMLINK")
    try:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("ACCEPTANCE_MANIFEST_INVALID") from exc
    if type(manifest) is not dict or type(manifest.get("version")) not in (str, int) or type(manifest.get("files")) is not dict:
        raise ValueError("ACCEPTANCE_MANIFEST_SCHEMA")

    declared_manifest = _declared_manifest(manifest)

    sources = copy.deepcopy(blank)
    supporting_files = {}
    configured = 0
    for relative, expected in sorted(declared_manifest.items()):
        path = _safe_path(SOURCE, relative)
        if not path.is_file():
            raise ValueError("ACCEPTANCE_SOURCE_MISSING")
        actual = _sha256(path)
        if actual != expected:
            raise ValueError("ACCEPTANCE_SOURCE_HASH_MISMATCH")
        if relative in EXPECTED_FILES:
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, UnicodeError, json.JSONDecodeError) as exc:
                raise ValueError("ACCEPTANCE_SOURCE_JSON_INVALID") from exc
            sources[relative] = {
                "status": "configured_hash_verified",
                "relative_path": relative,
                "sha256": actual,
                "data": data,
            }
            configured += 1
        else:
            # Supporting evidence is intentionally not decoded as JSON. Raw
            # RCSB, SDF, receipts, and SI files remain hash-bound bytes.
            supporting_files[relative] = {
                "status": "configured_hash_verified",
                "relative_path": relative,
                "sha256": actual,
                "size_bytes": path.stat().st_size,
            }

    return {
        "status": "configured" if configured == len(EXPECTED_FILES) else "partially_configured",
        "manifest": {
            "relative_path": "cases/acceptance_sources/manifest.json",
            "status": "configured_hash_verified",
            "version": manifest["version"],
            "sha256": _sha256(MANIFEST),
        },
        "sources": sources,
        "supporting_files": supporting_files,
        "approvals": None,
    }


def source_binding() -> dict[str, Any]:
    loaded = load_acceptance_evidence()
    return {
        "status": loaded["status"],
        "manifest": loaded["manifest"],
        "files": {
            name: {
                "status": row["status"],
                "relative_path": row["relative_path"],
                "sha256": row["sha256"],
            }
            for name, row in {
                **loaded["sources"],
                **loaded["supporting_files"],
            }.items()
        },
        "supporting_files": {
            name: {
                "status": row["status"],
                "relative_path": row["relative_path"],
                "sha256": row["sha256"],
            }
            for name, row in loaded["supporting_files"].items()
        },
    }


__all__ = ["EXPECTED_FILES", "load_acceptance_evidence", "source_binding"]
