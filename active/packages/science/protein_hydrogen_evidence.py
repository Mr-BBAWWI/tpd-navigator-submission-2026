"""Hash-bound evidence for fixed-frame PDB2PQR protein preparation.

The deposited mmCIF may contain ligands, waters, and additional chains.  The
source PDB is therefore validated against only the exact matching protein atom
subset of the mmCIF, while the complete source files and selection scope remain
hash-bound.  This module does not infer pKa values or populations and never
constitutes scientific approval.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import os
import re
import shutil
from pathlib import Path
from typing import Any

import gemmi
import numpy as np

from . import chemical_states
from .chemical_states import read_prepared_protein, run_pdb2pqr
from .molecules import rows

VERSION = "protein-hydrogen-evidence/1.1"
_COORDINATE_TOLERANCE_A = 0.0


class EvidenceValidationError(ValueError):
    """Raised when preparation artifacts fail provenance or frame validation."""

    def __init__(self, failures: list[dict[str, Any]]):
        self.failures = failures
        codes = ",".join(str(item.get("code", "validation_failure")) for item in failures)
        super().__init__(f"protein hydrogen evidence rejected: {codes}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _inside(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def _safe_regular(path: Path, root: Path) -> bool:
    return _inside(path, root) and path.is_file() and not path.is_symlink()


def _finite_xyz(values: Any) -> np.ndarray:
    xyz = np.asarray(values, dtype=float)
    if xyz.shape != (3,) or not np.isfinite(xyz).all():
        raise ValueError("nonfinite_or_missing_coordinate")
    return xyz


def _norm_token(value: Any) -> str:
    if value is None or value is False:
        return ""
    text = str(value).strip()
    return "" if text in {"", ".", "?"} else text


def _atom_key(chain: Any, sequence: Any, insertion: Any, component: Any,
              atom_name: Any, element: Any) -> tuple[str, str, str, str, str, str]:
    return (
        _norm_token(chain),
        _norm_token(sequence),
        _norm_token(insertion),
        _norm_token(component).upper(),
        _norm_token(atom_name).upper(),
        _norm_token(element).title(),
    )


def _key_text(key: tuple[str, ...]) -> str:
    return "|".join(key)


def _pdb_heavy(path: Path, *, source_protein_only: bool = False) -> dict[
        tuple[str, str, str, str, str, str], np.ndarray]:
    atoms: dict[tuple[str, str, str, str, str, str], np.ndarray] = {}
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for number, line in enumerate(stream, 1):
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip().title() if len(line) >= 78 else ""
            if not element:
                letters = re.sub(r"[^A-Za-z]", "", atom_name)
                element = letters[:1].title()
            if element in {"H", "D"}:
                continue
            if source_protein_only and line.startswith("HETATM"):
                raise ValueError(f"source_pdb_nonprotein_heavy_record:{number}")
            try:
                xyz = _finite_xyz((line[30:38], line[38:46], line[46:54]))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid_pdb_atom_line:{number}") from exc
            key = _atom_key(
                line[21:22], line[22:26], line[26:27], line[17:20],
                atom_name, element,
            )
            if key in atoms:
                raise ValueError(f"duplicate_pdb_heavy_atom:{_key_text(key)}")
            atoms[key] = xyz
    if not atoms:
        raise ValueError("pdb_contains_no_heavy_atoms")
    return atoms


def _row_get(row: dict[str, Any], *names: str, default: Any = "") -> Any:
    lowered = {str(key).lower(): value for key, value in row.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return default


def _cif_heavy(path: Path) -> tuple[
        dict[tuple[str, str, str, str, str, str], np.ndarray],
        dict[tuple[str, str, str, str, str, str], np.ndarray],
        dict[str, Any]]:
    """Read a real mmCIF with gemmi, including semicolon-delimited text fields."""
    try:
        block = gemmi.cif.read_file(str(path)).sole_block()
        atom_rows = list(rows(block, "_atom_site."))
    except Exception as exc:
        raise ValueError("mmcif_parser_failure") from exc
    if not atom_rows:
        raise ValueError("cif_atom_site_missing")

    all_atoms: dict[tuple[str, str, str, str, str, str], np.ndarray] = {}
    protein_atoms: dict[tuple[str, str, str, str, str, str], np.ndarray] = {}
    skipped_model = 0
    skipped_alt = 0
    skipped_occupancy = 0
    skipped_hydrogen = 0

    for row_number, row in enumerate(atom_rows, 1):
        model = _norm_token(_row_get(row, "pdbx_PDB_model_num", default="1")) or "1"
        if model != "1":
            skipped_model += 1
            continue
        alt = _norm_token(_row_get(row, "label_alt_id", "auth_alt_id"))
        if alt not in {"", "A"}:
            skipped_alt += 1
            continue
        occupancy_text = _norm_token(_row_get(row, "occupancy", default="1")) or "1"
        try:
            occupancy = float(occupancy_text)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"invalid_cif_occupancy:{row_number}") from exc
        if not math.isfinite(occupancy):
            raise ValueError(f"invalid_cif_occupancy:{row_number}")
        if occupancy <= 0:
            skipped_occupancy += 1
            continue

        element = _norm_token(_row_get(row, "type_symbol")).title()
        if not element:
            raise ValueError(f"cif_atom_element_missing:{row_number}")
        if element in {"H", "D"}:
            skipped_hydrogen += 1
            continue
        try:
            xyz = _finite_xyz((
                _row_get(row, "Cartn_x"),
                _row_get(row, "Cartn_y"),
                _row_get(row, "Cartn_z"),
            ))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"invalid_cif_atom_coordinate:{row_number}") from exc

        key = _atom_key(
            _row_get(row, "auth_asym_id", "label_asym_id"),
            _row_get(row, "auth_seq_id", "label_seq_id"),
            _row_get(row, "pdbx_PDB_ins_code"),
            _row_get(row, "auth_comp_id", "label_comp_id"),
            _row_get(row, "auth_atom_id", "label_atom_id"),
            element,
        )
        if key in all_atoms:
            raise ValueError(f"duplicate_cif_source_identity:{_key_text(key)}")
        all_atoms[key] = xyz
        group = _norm_token(_row_get(row, "group_PDB", default="ATOM")).upper()
        if group == "ATOM":
            protein_atoms[key] = xyz

    if not all_atoms:
        raise ValueError("cif_contains_no_heavy_atoms")
    metadata = {
        "parser": "gemmi.cif.read_file",
        "data_block": block.name,
        "atom_site_rows": len(atom_rows),
        "eligible_model1_alt_blank_or_A_positive_occupancy_heavy_atoms": len(all_atoms),
        "eligible_protein_ATOM_records": len(protein_atoms),
        "skipped_non_model1_rows": skipped_model,
        "skipped_nonblank_nonA_alt_rows": skipped_alt,
        "skipped_nonpositive_occupancy_rows": skipped_occupancy,
        "skipped_hydrogen_rows": skipped_hydrogen,
    }
    return all_atoms, protein_atoms, metadata


def _subset_sha256(atoms: dict[tuple[str, str, str, str, str, str], np.ndarray]) -> str:
    payload = [
        {"identity": list(key), "xyz_A": [float(value) for value in atoms[key]]}
        for key in sorted(atoms)
    ]
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=False, allow_nan=False).encode("utf-8")
    return _sha256_bytes(encoded)


def _record_key(record: dict[str, Any]) -> tuple[str, str, str, str, str, str]:
    return _atom_key(
        record.get("label_asym_id", ""), record.get("auth_seq_id", ""),
        record.get("insertion_code", ""), record.get("label_comp_id", ""),
        record.get("label_atom_id", ""), record.get("type_symbol", ""),
    )


def _protein_id_from_key(key: tuple[str, str, str, str, str, str]) -> str:
    return ":".join((key[0], key[1], key[3], key[4]))


def _json_record(record: dict[str, Any]) -> dict[str, Any]:
    result = dict(record)
    result["xyz"] = [float(value) for value in _finite_xyz(result["xyz"])]
    return result


def _installed_pdb2pqr_version() -> str | None:
    for distribution in ("pdb2pqr", "PDB2PQR"):
        try:
            return importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            continue
    return None


def _actual_tool_versions(paths: list[Path]) -> dict[str, Any]:
    patterns = {
        "pdb2pqr": re.compile(r"\bPDB2PQR(?:\s+version)?\s*[:=]?\s*(v?[0-9][^\s,;]*)", re.I),
        "propka": re.compile(r"\bPROPKA(?:\s+version)?\s*[:=]?\s*(v?[0-9][^\s,;]*)", re.I),
    }
    result: dict[str, Any] = {
        "pdb2pqr_module_version_actual_environment": _installed_pdb2pqr_version(),
        "pdb2pqr_version_from_actual_output": None,
        "propka_version_from_actual_output": None,
    }
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="replace")
        for name, pattern in patterns.items():
            match = pattern.search(text)
            field = f"{name}_version_from_actual_output"
            if match and result[field] is None:
                result[field] = match.group(1)
    return result


def _copy_exclusive(source: Path, destination: Path) -> None:
    if not source.is_file() or source.is_symlink():
        raise ValueError("regular_non_symlink_source_required")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as incoming, destination.open("xb") as outgoing:
        shutil.copyfileobj(incoming, outgoing, length=1024 * 1024)


class ProteinHydrogenEvidence:
    """Validate PDB2PQR artifacts against an exact deposited source frame."""

    def __init__(self, source_pdb: os.PathLike[str] | str,
                 source_cif: os.PathLike[str] | str, receipt: dict[str, Any],
                 artifact_root: os.PathLike[str] | str,
                 *, coordinate_tolerance_A: float = _COORDINATE_TOLERANCE_A):
        self.source_pdb = Path(source_pdb).expanduser().resolve()
        self.source_cif = Path(source_cif).expanduser().resolve()
        self.receipt = receipt
        self.artifact_root = Path(artifact_root).expanduser().resolve()
        self.coordinate_tolerance_A = float(coordinate_tolerance_A)
        if not math.isfinite(self.coordinate_tolerance_A) or self.coordinate_tolerance_A < 0:
            raise ValueError("coordinate_tolerance_A must be finite and nonnegative")

    @staticmethod
    def _failure(failures: list[dict[str, Any]], code: str, **details: Any) -> None:
        failures.append({"code": code, **details})

    def evidence(self) -> dict[str, Any]:
        failures: list[dict[str, Any]] = []
        receipt = self.receipt if isinstance(self.receipt, dict) else {}
        if not isinstance(self.receipt, dict):
            self._failure(failures, "receipt_not_dictionary")

        for label, path in (("source_pdb", self.source_pdb), ("source_cif", self.source_cif)):
            if not _safe_regular(path, self.artifact_root):
                self._failure(failures, "unsafe_or_missing_file", artifact=label)

        source_pdb_atoms: dict[tuple[str, str, str, str, str, str], np.ndarray] = {}
        cif_all_atoms: dict[tuple[str, str, str, str, str, str], np.ndarray] = {}
        cif_protein_atoms: dict[tuple[str, str, str, str, str, str], np.ndarray] = {}
        selected_cif_atoms: dict[tuple[str, str, str, str, str, str], np.ndarray] = {}
        cif_scope: dict[str, Any] = {}
        if not failures:
            try:
                source_pdb_atoms = _pdb_heavy(self.source_pdb, source_protein_only=True)
                cif_all_atoms, cif_protein_atoms, cif_scope = _cif_heavy(self.source_cif)
            except (OSError, UnicodeError, ValueError, RuntimeError) as exc:
                self._failure(
                    failures, "source_parse_failure",
                    failure_type=type(exc).__name__, reason=str(exc),
                )

        source_frame_maximum: float | None = None
        missing_from_cif: list[tuple[str, str, str, str, str, str]] = []
        source_frame_changed: list[dict[str, Any]] = []
        if source_pdb_atoms and cif_protein_atoms:
            missing_from_cif = sorted(set(source_pdb_atoms) - set(cif_protein_atoms))
            if missing_from_cif:
                self._failure(
                    failures, "source_pdb_atom_missing_from_cif_protein_subset",
                    atoms=[_key_text(key) for key in missing_from_cif],
                )
            selected_cif_atoms = {
                key: cif_protein_atoms[key]
                for key in source_pdb_atoms if key in cif_protein_atoms
            }
            displacements = {
                key: float(np.linalg.norm(source_pdb_atoms[key] - selected_cif_atoms[key]))
                for key in selected_cif_atoms
            }
            source_frame_maximum = max(displacements.values(), default=0.0)
            source_frame_changed = [
                {"atom_id": _key_text(key), "displacement_A": distance}
                for key, distance in sorted(displacements.items())
                if distance > self.coordinate_tolerance_A
            ]
            if source_frame_changed:
                self._failure(
                    failures, "source_pdb_cif_coordinate_mismatch",
                    maximum_displacement_A=source_frame_maximum,
                    atoms=source_frame_changed,
                )

        expected_paths = {
            "input_path": receipt.get("input_sha256"),
            "output_pqr_path": receipt.get("output_sha256"),
            "output_pdb_path": receipt.get("output_pdb_sha256"),
            "stdout_path": receipt.get("stdout_sha256"),
            "stderr_path": receipt.get("stderr_sha256"),
        }
        verified_paths: dict[str, Path] = {}
        for field, expected_hash in expected_paths.items():
            raw = receipt.get(field)
            if not isinstance(raw, str) or not raw or not isinstance(expected_hash, str):
                self._failure(failures, "receipt_artifact_field_missing", field=field)
                continue
            path = Path(raw).expanduser().resolve()
            if not _safe_regular(path, self.artifact_root):
                self._failure(failures, "unsafe_or_missing_receipt_artifact", field=field)
                continue
            if _sha256(path) != expected_hash:
                self._failure(failures, "receipt_artifact_hash_mismatch", field=field)
                continue
            verified_paths[field] = path

        if verified_paths.get("input_path") != self.source_pdb:
            self._failure(failures, "receipt_input_is_not_registered_source_pdb")

        executable_raw = receipt.get("executable_path")
        executable_hash = receipt.get("executable_sha256")
        executable_path: Path | None = None
        if isinstance(executable_raw, str) and executable_raw and isinstance(executable_hash, str):
            executable_path = Path(executable_raw).expanduser().resolve()
            if not executable_path.is_file() or executable_path.is_symlink():
                self._failure(failures, "unsafe_or_missing_pdb2pqr_executable")
            elif _sha256(executable_path) != executable_hash:
                self._failure(failures, "pdb2pqr_executable_hash_mismatch")
        else:
            self._failure(failures, "pdb2pqr_executable_provenance_missing")

        fixed_heavy = receipt.get("preserve_heavy_requested") is True
        command = receipt.get("command")
        if not isinstance(command, list) or not all(isinstance(item, str) for item in command):
            self._failure(failures, "receipt_command_missing")
            command = []
        if fixed_heavy and not {"--noopt", "--nodebump"} <= set(command):
            self._failure(failures, "fixed_heavy_command_flags_missing")
        if not fixed_heavy and ({"--noopt", "--nodebump"} & set(command)):
            self._failure(failures, "optimized_command_contains_fixed_heavy_flags")
        if not any(item.startswith("--with-ph=") for item in command):
            self._failure(failures, "pH_command_context_missing")

        output_atoms: dict[tuple[str, str, str, str, str, str], np.ndarray] = {}
        output_pdb_path = verified_paths.get("output_pdb_path")
        if output_pdb_path is not None:
            try:
                output_atoms = _pdb_heavy(output_pdb_path)
            except (OSError, ValueError) as exc:
                self._failure(
                    failures, "prepared_pdb_parse_failure",
                    failure_type=type(exc).__name__, reason=str(exc),
                )

        missing: list[tuple[str, str, str, str, str, str]] = []
        extras: list[tuple[str, str, str, str, str, str]] = []
        changed: list[dict[str, Any]] = []
        maximum_displacement: float | None = None
        if source_pdb_atoms and output_atoms:
            source_keys, output_keys = set(source_pdb_atoms), set(output_atoms)
            missing = sorted(source_keys - output_keys)
            extras = sorted(output_keys - source_keys)
            displacements = {
                key: float(np.linalg.norm(output_atoms[key] - source_pdb_atoms[key]))
                for key in source_keys & output_keys
            }
            maximum_displacement = max(displacements.values(), default=0.0)
            changed = [
                {"atom_id": _key_text(key), "displacement_A": distance}
                for key, distance in sorted(displacements.items())
                if distance > self.coordinate_tolerance_A
            ]
            if missing:
                self._failure(
                    failures, "source_heavy_atom_missing",
                    atoms=[_key_text(key) for key in missing],
                )
            if changed:
                self._failure(failures, "source_heavy_atom_position_changed", atoms=changed)

        allowed_oxt = [key for key in extras if key[4] == "OXT" and key[5] == "O"]
        disallowed_extras = [key for key in extras if key not in allowed_oxt]
        if disallowed_extras:
            self._failure(
                failures, "unapproved_extra_heavy_atom",
                atoms=[_key_text(key) for key in disallowed_extras],
            )
        extra_annotations = [
            {
                "atom_id": _key_text(key),
                "annotation": "pending_terminal_OXT_scientific_review",
                "included_in_source_heavy_contact_groups": False,
            }
            for key in allowed_oxt
        ] + [
            {
                "atom_id": _key_text(key),
                "annotation": "unapproved_extra_heavy_atom",
                "included_in_source_heavy_contact_groups": False,
            }
            for key in disallowed_extras
        ]

        prepared: dict[str, Any] | None = None
        source_protein_atoms: list[dict[str, Any]] = []
        protein_hydrogens: list[dict[str, Any]] = []
        parent_uncertainties: list[dict[str, Any]] = []
        if output_pdb_path is not None and not any(
            item["code"] in {"receipt_artifact_hash_mismatch", "prepared_pdb_parse_failure"}
            for item in failures
        ):
            try:
                prepared = read_prepared_protein(receipt)
                source_protein_atoms = [
                    _json_record(record) for record in prepared["protein_atoms"]
                    if _record_key(record) in source_pdb_atoms
                ]
                heavy_for_validation = {
                    _protein_id_from_key(_record_key(record)): record
                    for record in source_protein_atoms
                }
                mutable_hydrogens = [dict(record) for record in prepared["protein_hydrogens"]]
                for hydrogen in mutable_hydrogens:
                    parent = hydrogen.get("parent_atom_id")
                    if parent not in heavy_for_validation:
                        hydrogen["parent_atom_id"] = None
                        hydrogen["parent_assignment_status"] = "unknown_excluded_from_directional_tests"

                valid_parent_map, actual_parent_uncertainties = (
                    chemical_states._protein_hydrogen_parents(
                        mutable_hydrogens, heavy_for_validation
                    )
                )
                valid_hydrogen_objects = {
                    id(hydrogen)
                    for attached in valid_parent_map.values() for hydrogen in attached
                }
                for hydrogen in mutable_hydrogens:
                    if id(hydrogen) not in valid_hydrogen_objects:
                        hydrogen["parent_atom_id"] = None
                        hydrogen["parent_assignment_status"] = "unknown_excluded_from_directional_tests"
                    else:
                        hydrogen["parent_assignment_status"] = "validated_source_parent"
                parent_uncertainties.extend(actual_parent_uncertainties)
                protein_hydrogens = [_json_record(record) for record in mutable_hydrogens]
                if len(source_protein_atoms) != len(source_pdb_atoms):
                    self._failure(failures, "prepared_source_heavy_mapping_incomplete")
            except (TypeError, ValueError, OSError) as exc:
                self._failure(
                    failures, "prepared_protein_read_failure",
                    failure_type=type(exc).__name__, reason=str(exc),
                )

        manifest: list[dict[str, Any]] = []
        generated_evidence_path = self.artifact_root / "protein_hydrogen_evidence.json"
        if self.artifact_root.is_dir() and not self.artifact_root.is_symlink():
            for path in sorted(self.artifact_root.rglob("*")):
                if path == generated_evidence_path:
                    if path.is_symlink() or not path.is_file():
                        self._failure(failures, "invalid_generated_evidence_json")
                    continue
                if path.is_symlink():
                    self._failure(failures, "symlink_in_artifact_tree")
                elif path.is_file():
                    manifest.append({
                        "path": path.relative_to(self.artifact_root).as_posix(),
                        "sha256": _sha256(path),
                        "size_bytes": path.stat().st_size,
                    })
                elif not path.is_dir():
                    self._failure(failures, "non_regular_artifact")
        else:
            self._failure(failures, "unsafe_artifact_root")

        pka_entries = receipt.get(
            "pKa_entries_from_actual_output",
            receipt.get("pka_entries_from_actual_output", []),
        )
        if not isinstance(pka_entries, list):
            self._failure(failures, "invalid_actual_pka_entry_format")
            pka_entries = []
        manifest_paths = {item["path"] for item in manifest}
        verified_pka_entries: list[dict[str, Any]] = []
        for entry in pka_entries:
            if not isinstance(entry, dict) or not isinstance(entry.get("source_path"), str):
                self._failure(failures, "unverifiable_actual_pka_entry")
                continue
            source_path = Path(entry["source_path"]).expanduser().resolve()
            if not _safe_regular(source_path, self.artifact_root):
                self._failure(failures, "pka_source_outside_verified_artifacts")
                continue
            relative = source_path.relative_to(self.artifact_root).as_posix()
            if relative not in manifest_paths:
                self._failure(failures, "pka_source_not_hash_registered")
                continue
            try:
                value = float(entry["pKa"])
            except (KeyError, TypeError, ValueError):
                self._failure(failures, "invalid_actual_pka_entry_format")
                continue
            if not math.isfinite(value):
                self._failure(failures, "invalid_actual_pka_entry_format")
                continue
            verified_pka_entries.append(dict(entry))

        version_paths = [
            path for field, path in verified_paths.items()
            if field in {"stdout_path", "stderr_path"}
        ]
        for item in verified_pka_entries:
            candidate = Path(item["source_path"]).resolve()
            if candidate not in version_paths:
                version_paths.append(candidate)
        tool_versions = _actual_tool_versions(version_paths)

        source_frame_accepted = bool(
            source_pdb_atoms
            and len(selected_cif_atoms) == len(source_pdb_atoms)
            and not source_frame_changed
        )
        prepared_frame_accepted = bool(source_pdb_atoms and not missing and not changed)
        orientation_review = bool(protein_hydrogens)
        state_flags = {
            "fixed_heavy_requested": fixed_heavy,
            "source_pdb_matches_selected_cif_protein_subset_exactly": source_frame_accepted,
            "source_heavy_coordinates_accepted": prepared_frame_accepted,
            "terminal_oxt_pending_review": bool(allowed_oxt),
            "added_hydrogen_orientation_optimized": not fixed_heavy,
            "added_hydrogen_orientation_requires_review": orientation_review,
            "unknown_hydrogen_parents_excluded_from_directional_tests": any(
                item.get("parent_atom_id") is None for item in protein_hydrogens
            ),
            "optimized_external_heavy_unchanged": (
                None if fixed_heavy else prepared_frame_accepted
            ),
            "computed_diagnostic": not bool(failures),
            "pending_human_review": True,
            "formal_scientific_approval": False,
            "pKa_prediction_performed_by_this_module": False,
            "population_prediction_performed_by_this_module": False,
        }
        status = "rejected" if failures else "computed_diagnostic_pending_human_review"
        prepared_uncertainties = list(prepared.get("uncertainties", [])) if prepared else []
        result = {
            "format": VERSION,
            "status": status,
            "protein_atoms": source_protein_atoms,
            "protein_hydrogens": protein_hydrogens,
            "receipt": receipt,
            "failures": failures,
            "state_flags": state_flags,
            "counts": {
                "source_pdb_protein_heavy_atoms": len(source_pdb_atoms),
                "source_pdb_heavy_atoms": len(source_pdb_atoms),
                "source_cif_eligible_heavy_atoms_full_deposition": len(cif_all_atoms),
                "source_cif_eligible_protein_heavy_atoms_full_deposition": len(cif_protein_atoms),
                "source_cif_selected_matching_protein_heavy_atoms": len(selected_cif_atoms),
                "source_cif_heavy_atoms": len(cif_all_atoms),
                "prepared_heavy_atoms_total": len(output_atoms),
                "source_heavy_atoms_for_contacts": len(source_protein_atoms),
                "prepared_hydrogens": len(protein_hydrogens),
                "hydrogens_with_validated_source_parent": sum(
                    item.get("parent_assignment_status") == "validated_source_parent"
                    for item in protein_hydrogens
                ),
                "hydrogens_with_known_source_parent": sum(
                    item.get("parent_atom_id") is not None for item in protein_hydrogens
                ),
                "extra_heavy_atoms": len(extras),
                "pending_terminal_oxt_atoms": len(allowed_oxt),
                "actual_pka_entries": len(verified_pka_entries),
            },
            "sha256": {
                "source_pdb_file": (
                    _sha256(self.source_pdb)
                    if _safe_regular(self.source_pdb, self.artifact_root) else None
                ),
                "source_cif_file": (
                    _sha256(self.source_cif)
                    if _safe_regular(self.source_cif, self.artifact_root) else None
                ),
                "source_pdb_selected_subset_canonical": (
                    _subset_sha256(source_pdb_atoms) if source_pdb_atoms else None
                ),
                "source_cif_selected_subset_canonical": (
                    _subset_sha256(selected_cif_atoms) if selected_cif_atoms else None
                ),
                "pdb2pqr_executable": (
                    executable_hash if executable_path is not None else None
                ),
            },
            "source_frame": {
                "identity": "auth_chain,auth_residue,insertion,component,atom_name,element",
                "coordinate_tolerance_A": self.coordinate_tolerance_A,
                "coordinate_policy": "exact at default tolerance 0.0 A",
                "source_pdb_scope": "heavy ATOM records only; heavy HETATM records are rejected",
                "source_cif_scope": {
                    **cif_scope,
                    "model": "1",
                    "alternate_location_policy": "blank or A",
                    "occupancy_policy": "finite and greater than zero",
                    "matching_record_group": "ATOM",
                    "selection": "exact identities present in source protein-only PDB",
                    "nonmatching_ligands_waters_and_additional_chains": "retained in full-file hash but excluded from subset equality",
                    "selected_auth_chains": sorted({key[0] for key in selected_cif_atoms}),
                    "unselected_eligible_heavy_atom_count": max(
                        0, len(cif_all_atoms) - len(selected_cif_atoms)
                    ),
                },
                "pdb_cif_maximum_displacement_A": source_frame_maximum,
                "prepared_source_maximum_displacement_A": maximum_displacement,
            },
            "extras": extra_annotations,
            "uncertainties": prepared_uncertainties + parent_uncertainties,
            "actual_pka_entries": verified_pka_entries,
            "pKa_population_policy": {
                "pKa_values_exposed": "only rows parsed from hash-registered actual tool output",
                "unreported_pKa_values": "no inference or replacement value is produced",
                "population_values": "not calculated or inferred by this module",
            },
            "method": {
                "worker": "packages.science.chemical_states.run_pdb2pqr",
                "reader": "packages.science.chemical_states.read_prepared_protein",
                "mmcif_parser": "gemmi via repository rows helper",
                "hydrogen_parent_validator": "packages.science.chemical_states._protein_hydrogen_parents",
                "fixed_heavy_flags": ["--noopt", "--nodebump"] if fixed_heavy else [],
                "orientation_interpretation": (
                    "Added hydrogen orientations are unoptimized and remain pending review."
                    if fixed_heavy else
                    "External optimization was requested; orientations still remain pending human review."
                ),
                "scientific_acceptance": "pending human review; no automatic approval",
            },
            "fingerprint": {
                "module_version": VERSION,
                "pdb2pqr_executable_sha256": executable_hash,
                **tool_versions,
            },
            "artifact_manifest": manifest,
            "api_cost": None,
        }
        json.dumps(result, ensure_ascii=False, allow_nan=False)
        return result

    def validate(self) -> dict[str, Any]:
        """Return validated diagnostic evidence or raise without repairing files."""
        result = self.evidence()
        if result["failures"]:
            raise EvidenceValidationError(result["failures"])
        return result


def prepare_evidence(source_pdb: os.PathLike[str] | str,
                     source_cif: os.PathLike[str] | str,
                     output_dir: os.PathLike[str] | str,
                     executable: os.PathLike[str] | str,
                     *, pH: float = 7.4, fixed_heavy: bool = True,
                     timeout: float = 300) -> dict[str, Any]:
    """Copy sources, execute the actual worker, and validate immutable evidence."""
    if not isinstance(fixed_heavy, bool):
        raise ValueError("fixed_heavy must be bool")
    destination = Path(output_dir).expanduser().resolve()
    if destination.exists():
        raise FileExistsError("output_directory_already_exists")
    source_pdb_path = Path(source_pdb).expanduser().resolve()
    source_cif_path = Path(source_cif).expanduser().resolve()
    executable_path = Path(executable).expanduser().resolve()
    for path in (source_pdb_path, source_cif_path, executable_path):
        if not path.is_file() or path.is_symlink():
            raise ValueError("regular_non_symlink_input_required")

    destination.mkdir(parents=True, exist_ok=False)
    copied_pdb = destination / "source" / "source.pdb"
    copied_cif = destination / "source" / "source.cif"
    _copy_exclusive(source_pdb_path, copied_pdb)
    _copy_exclusive(source_cif_path, copied_cif)
    preparation_dir = destination / "pdb2pqr"
    receipt = run_pdb2pqr(
        copied_pdb, preparation_dir, executable_path,
        pH=pH, timeout=timeout, preserve_heavy=fixed_heavy,
    )
    evidence = ProteinHydrogenEvidence(
        copied_pdb, copied_cif, receipt, destination
    ).validate()
    evidence_path = destination / "protein_hydrogen_evidence.json"
    with evidence_path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(evidence, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")
    return evidence


__all__ = [
    "VERSION", "EvidenceValidationError", "ProteinHydrogenEvidence",
    "prepare_evidence",
]
