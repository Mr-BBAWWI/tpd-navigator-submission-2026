#!/usr/bin/env python3
"""Export verified original docking poses and analog-specific state diagnostics."""
from __future__ import annotations

import argparse
import hashlib
import html
import io
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import Chem

from packages.science import expert_followup
from packages.science.analog_state_comparison import build_state_set
from scripts.export_expert_evidence import (
    ANALOG_IDS, DEFAULT_NEUTRAL, DEFAULT_PROTEIN, EXPECTED_SOURCE_HASHES,
    SOURCE_BINDING_PATHS, _checked_destination, _checked_file, _document,
    _field, _mol, _protein, _sha, _source_input_binding,
)


def _bytes_sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _read_bytes(port, reference) -> bytes:
    value = port.read(reference)
    if isinstance(value, str):
        value = value.encode("utf-8")
    if not isinstance(value, bytes):
        raise TypeError("ARTIFACT_BYTES_REQUIRED")
    return value


def _poses(value: bytes) -> list[Chem.Mol]:
    supplier = Chem.ForwardSDMolSupplier(
        io.BytesIO(value), removeHs=False, sanitize=True, strictParsing=True
    )
    records = list(supplier)
    if any(record is None for record in records) or len(records) != 5:
        raise ValueError("EXACTLY_FIVE_VALID_ORIGINAL_DOCKING_POSES_REQUIRED")
    for record in records:
        if record.GetNumConformers() != 1:
            raise ValueError("DOCKING_POSE_COORDINATE_CONFORMER_REQUIRED")
    return records


def _selected(document: dict) -> dict[str, dict]:
    rows = document.get("analogs")
    if not isinstance(rows, list):
        raise ValueError("ANALOG_LIST_REQUIRED")
    result = {}
    for analog_id in ANALOG_IDS:
        matches = [
            row for row in rows if isinstance(row, dict)
            and row.get("id") == analog_id and row.get("selected") is True
        ]
        if len(matches) != 1:
            raise ValueError(f"SELECTED_ANALOG_REQUIRED:{analog_id}")
        result[analog_id] = matches[0]
    return result


def _artifact_ref(analog: dict, key: str):
    docking = analog.get("docking")
    files = docking.get("files") if isinstance(docking, dict) else None
    reference = files.get(key) if isinstance(files, dict) else None
    if not isinstance(reference, dict):
        raise ValueError(f"ORIGINAL_DOCKING_ARTIFACT_REFERENCE_REQUIRED:{key}")
    artifact_id = reference.get("artifact_id")
    digest = reference.get("sha256")
    if not isinstance(artifact_id, str) or not artifact_id:
        raise ValueError(f"REGISTERED_DOCKING_ARTIFACT_ID_REQUIRED:{key}")
    if (not isinstance(digest, str) or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest.lower())):
        raise ValueError(f"REGISTERED_DOCKING_ARTIFACT_HASH_REQUIRED:{key}")
    media_type = reference.get("media_type")
    if key == "poses_sdf" and media_type != "chemical/x-mdl-sdfile":
        raise ValueError("ORIGINAL_DOCKING_POSES_SDF_MEDIA_TYPE_REQUIRED")
    if key == "poses_pdbqt" and media_type not in {
            "text/plain", "chemical/x-pdbqt"}:
        raise ValueError("ORIGINAL_DOCKING_POSES_PDBQT_MEDIA_TYPE_REQUIRED")
    return reference


def _verified_artifact_bytes(port, reference: dict, analog_id: str, key: str) -> bytes:
    raw = _read_bytes(port, reference)
    if _bytes_sha(raw) != reference["sha256"]:
        raise ValueError(f"DOCKING_ARTIFACT_HASH_MISMATCH:{analog_id}:{key}")
    return raw


def _required_json_dict(value, name: str) -> dict:
    if not isinstance(value, dict):
        raise ValueError(f"{name}_JSON_OBJECT_REQUIRED")
    return json.loads(json.dumps(value, ensure_ascii=False, allow_nan=False))


def _repair_pose_maps(service, raw_poses: list[Chem.Mol], restored_poses: list[Chem.Mol],
                      expected_mapped_smiles: str):
    """Restore each SDF record from only its same-index PDBQT model."""
    if len(raw_poses) != 5 or len(restored_poses) != 5:
        raise ValueError("EXACTLY_FIVE_SAME_INDEX_POSE_MODELS_REQUIRED")
    expected = Chem.MolFromSmiles(expected_mapped_smiles)
    if expected is None:
        raise ValueError("EXPECTED_MAPPED_SMILES_REQUIRED")
    expected = Chem.RemoveHs(expected)
    repaired = []
    receipts = []
    for index, (raw_pose, restored_pose) in enumerate(zip(raw_poses, restored_poses)):
        normalized, receipt = service._repair_actual_pose_maps(
            raw_pose, restored_pose, expected
        )
        receipt = dict(receipt)
        receipt.update({
            "pose_index_zero_based": index,
            "correspondence": "same_zero_based_sdf_record_and_pdbqt_model",
            "expected_mapped_smiles": expected_mapped_smiles,
            "nearest_pose_guess_used": False,
            "alignment_applied": False,
        })
        repaired.append(normalized)
        receipts.append(receipt)
    return repaired, receipts


def _prefix_state_artifacts(report: dict, analog_id: str, pose_index: int) -> None:
    prefix = f"states/{analog_id}/pose-{pose_index}/"
    for state in report.get("states", []):
        artifacts = state.get("artifacts")
        if isinstance(artifacts, dict):
            for key, value in list(artifacts.items()):
                if isinstance(value, str):
                    artifacts[key] = prefix + value


def _write_json(path: Path, value) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8", newline="\n",
    )


def _rejected_angles(profile: dict) -> list:
    result = []

    def visit(value):
        if isinstance(value, dict):
            text = " ".join(str(value.get(key, "")) for key in (
                "status", "reason", "decision", "kind"
            )).lower()
            angles = {
                str(key): item for key, item in value.items()
                if "angle" in str(key).lower()
            }
            if angles and ("reject" in text or value.get("accepted") is False):
                result.append({
                    "id": value.get("id"), "reason": value.get("reason"),
                    "angles": angles, "participants": value.get("participants"),
                })
            for key, item in value.items():
                if "rejected" in str(key).lower() and "angle" in str(key).lower():
                    result.append({str(key): item})
                visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)

    visit(profile)
    return result


def _hydrophobic_counts(state: dict) -> dict[str, int]:
    result = {"VAL1408": 0, "PHE1409": 0, "ILE1470": 0}
    for contact in state.get("hydrophobic_contacts_VAL1408_PHE1409_ILE1470", []):
        seen = set()
        for participant in contact.get("protein_atom_ids", []):
            fields = participant.split(":") if isinstance(participant, str) else []
            if len(fields) == 4:
                key = fields[2].upper() + fields[1]
                if key in result:
                    seen.add(key)
        for key in seen:
            result[key] += 1
    return result


def _index(comparison: dict) -> str:
    def escaped(value) -> str:
        if isinstance(value, (dict, list)):
            value = json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)
        return html.escape(str(value), quote=True)

    body = []
    for analog in comparison["analogs"]:
        state_rows = []
        for pose in analog["poses"]:
            report = pose["report"]
            for state in report["states"]:
                map17 = state.get("map17_ASN1464_OD1") or {}
                nitrogens = state.get("mapped_nitrogens") or {}
                links = [
                    f'<a href="{escaped(target)}">{escaped(label)}</a>'
                    for label, target in (state.get("artifacts") or {}).items()
                ]
                state_rows.append(
                    "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"
                    "<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"
                    "<td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
                        escaped(pose["pose_index"]), escaped(state["state_id"]),
                        escaped(state["status"]), escaped(state.get("formal_charge")),
                        escaped(report.get("shared_frame_core_rmsd_A")),
                        escaped(state.get("heavy_coordinate_max_displacement_A")),
                        escaped(map17.get("heavy_distance_A")),
                        escaped(map17.get("directional_hbond_observed")),
                        escaped(_rejected_angles(state.get("interaction_profile") or {})),
                        escaped(_hydrophobic_counts(state)),
                        escaped(nitrogens.get("19", "없음")),
                        escaped(nitrogens.get("5001", "없음")), " · ".join(links),
                    )
                )
        body.append(
            f"<h2>{escaped(analog['analog_id'])}</h2>"
            f"<p>원본 도킹 포즈 {escaped(analog['pose_count'])}개를 모두 유지했습니다. "
            "최적 포즈 또는 최종 상태를 주장하지 않습니다.</p>"
            "<table><thead><tr><th>포즈</th><th>상태</th><th>계산 상태</th>"
            "<th>총 전하</th><th>포즈 코어 RMSD (Å)</th><th>중원자 최대 이동 (Å)</th>"
            "<th>map17–ASN1464 OD1 거리 (Å)</th><th>map17 방향성 수소결합</th>"
            "<th>거부된 각도 근거</th><th>VAL1408/PHE1409/ILE1470 소수성 접촉 수</th>"
            "<th>N19 전하/H 상태</th><th>N5001 전하/H 상태</th><th>상태 SDF</th>"
            "</tr></thead><tbody>" + "".join(state_rows) + "</tbody></table>"
        )
    return """<!doctype html><html lang="ko"><head><meta charset="utf-8">
<title>아날로그 양성자화 상태 비교</title><style>
body{font-family:system-ui,sans-serif;max-width:1600px;margin:2rem auto;padding:0 1rem}
.pending{background:#fff4d6;border-left:5px solid #c77b00;padding:1rem}
table{border-collapse:collapse;width:100%;font-size:.9rem}th,td{border:1px solid #bbb;padding:.35rem;text-align:left;vertical-align:top;overflow-wrap:anywhere}
</style></head><body><h1>아날로그별 양성자화 상태 비교</h1>
<div class="pending"><strong>공식 과학 검토가 대기 중입니다.</strong> 상태와 포즈는 선택되지 않았고,
pKa 또는 개체군은 추정하지 않았으며, 이 보고서는 승인이 아닙니다.</div>
<p><a href="comparison.json">기계 판독 비교 JSON</a> ·
<a href="manifest.json">SHA-256 매니페스트</a></p>""" + "".join(body) + "</body></html>\n"


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--store", required=True)
    result.add_argument("--job-id", required=True)
    result.add_argument("--output", required=True)
    return result


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    destination = _checked_destination(Path(args.output))
    neutral_path = _checked_file(DEFAULT_NEUTRAL, "neutral parent SDF")
    protein_path = _checked_file(DEFAULT_PROTEIN, "protein CIF")
    if _sha(neutral_path) != EXPECTED_SOURCE_HASHES[neutral_path.name]:
        raise ValueError("ACTUAL_SOURCE_HASH_MISMATCH:SMARCA2-neutral-design.sdf")
    if _sha(protein_path) != EXPECTED_SOURCE_HASHES[protein_path.name]:
        raise ValueError("ACTUAL_SOURCE_HASH_MISMATCH:6HAZ.cif")

    from packages.platform.scientific_acceptance import ScientificAcceptanceService
    from packages.platform.store import Store

    service = ScientificAcceptanceService(Store(Path(args.store).resolve()), "local-research")
    with service.store.db() as db:
        job, result, binding, verified_source = service._verified_result(
            db, args.job_id, require_current=True
        )
        archive_reference = _field(_field(result, "files"), "json")
        archive_raw = _read_bytes(service.port, archive_reference)
        document = _document(json.loads(archive_raw))
        source_binding = _source_input_binding(document)
        selected = _selected(document)
        verified_binding = _required_json_dict(binding, "VERIFIED_BINDING")
        verified_source_json = _required_json_dict(verified_source, "VERIFIED_SOURCE")
        pose_archives = {}
        for analog_id, analog in selected.items():
            sdf_reference = _artifact_ref(analog, "poses_sdf")
            pdbqt_reference = _artifact_ref(analog, "poses_pdbqt")
            sdf_raw = _verified_artifact_bytes(
                service.port, sdf_reference, analog_id, "poses_sdf"
            )
            pdbqt_raw = _verified_artifact_bytes(
                service.port, pdbqt_reference, analog_id, "poses_pdbqt"
            )
            restored_pdbqt_poses = service._pdbqt_poses(pdbqt_raw)
            if len(restored_pdbqt_poses) != 5:
                raise ValueError(f"EXACTLY_FIVE_VALID_PDBQT_MODELS_REQUIRED:{analog_id}")
            pose_archives[analog_id] = (
                sdf_raw, sdf_reference, pdbqt_raw, pdbqt_reference,
                restored_pdbqt_poses,
            )

    verified_load_followup = _required_json_dict(
        expert_followup.load_followup(), "VERIFIED_LOAD_FOLLOWUP"
    )

    for name in (neutral_path.name, protein_path.name):
        record = source_binding.get(SOURCE_BINDING_PATHS[name])
        if not isinstance(record, dict) or record.get("status") != "configured_hash_verified":
            raise ValueError(f"SOURCE_BINDING_REQUIRED:{name}")
        if record.get("sha256") != EXPECTED_SOURCE_HASHES[name]:
            raise ValueError(f"SOURCE_BINDING_HASH_MISMATCH:{name}")

    parent = _mol(neutral_path)
    protein_atoms = _protein(protein_path, "A")
    destination.mkdir()
    sources = destination / "sources"
    sources.mkdir()
    archive_copy = sources / "verified-design-archive.json"
    archive_copy.write_bytes(archive_raw)
    shutil.copyfile(neutral_path, sources / neutral_path.name)
    shutil.copyfile(protein_path, sources / protein_path.name)
    followup_source = expert_followup.ROOT / expert_followup.SOURCE_RELATIVE_PATH
    followup_copy = sources / Path(expert_followup.SOURCE_RELATIVE_PATH).name
    shutil.copyfile(followup_source, followup_copy)
    if (_sha(followup_copy) != expert_followup.SOURCE_SHA256
            or followup_copy.stat().st_size != expert_followup.SOURCE_BYTES):
        raise ValueError("EXPERT_FOLLOWUP_SOURCE_COPY_MISMATCH")

    analog_reports = []
    flat_rows = []
    for analog_id in ANALOG_IDS:
        raw, reference, pdbqt_raw, pdbqt_reference, restored_pdbqt_poses = pose_archives[
            analog_id
        ]
        original_copy = sources / f"{analog_id}-original-docking-poses.sdf"
        original_pdbqt_copy = sources / f"{analog_id}-original-docking-poses.pdbqt"
        original_copy.write_bytes(raw)
        original_pdbqt_copy.write_bytes(pdbqt_raw)
        raw_poses = _poses(raw)
        poses, restoration_receipts = _repair_pose_maps(
            service, raw_poses, restored_pdbqt_poses,
            selected[analog_id].get("mapped_smiles", ""),
        )
        pose_reports = []
        for pose_index, pose in enumerate(poses):
            state_directory = destination / "states" / analog_id / f"pose-{pose_index}"
            report = build_state_set(
                analog_id, pose, parent, protein_atoms, output_dir=state_directory
            )
            _prefix_state_artifacts(report, analog_id, pose_index)
            pose_reports.append({
                "pose_index": pose_index,
                "pose_selection_policy": "all_original_poses_retained",
                "scientific_best_pose_claimed": False,
                "map_restoration_receipt": restoration_receipts[pose_index],
                "report": report,
            })
            for state in report["states"]:
                flat_rows.append({
                    "row_type": "analog_pose_state",
                    "analog_id": analog_id,
                    "pose_index": pose_index,
                    "state_index": state["state_index"],
                    "state_id": state["state_id"],
                    "status": state["status"],
                    "formal_charge": state.get("formal_charge"),
                    "shared_frame_core_rmsd_A": report["shared_frame_core_rmsd_A"],
                    "core_rmsd_within_1A": report["shared_frame_core_rmsd_within_limit"],
                    "heavy_coordinate_max_displacement_A": state.get(
                        "heavy_coordinate_max_displacement_A"
                    ),
                    "selected_state": None,
                    "microstate_gate": "pending",
                    "core_final": "pending",
                })
        analog_reports.append({
            "analog_id": analog_id,
            "pose_count": len(poses),
            "source_original_copy": original_copy.relative_to(destination).as_posix(),
            "source_original_pdbqt_copy": original_pdbqt_copy.relative_to(destination).as_posix(),
            "source_archive_sha256": _bytes_sha(raw),
            "source_pdbqt_sha256": _bytes_sha(pdbqt_raw),
            "source_artifact_id": reference.get("artifact_id"),
            "source_artifact_sha256": reference.get("sha256"),
            "source_pdbqt_artifact_id": pdbqt_reference.get("artifact_id"),
            "source_pdbqt_artifact_sha256": pdbqt_reference.get("sha256"),
            "map_restoration_policy": (
                "exact expected mapped SMILES and same-index PDBQT model; "
                "no nearest guess and no alignment"
            ),
            "map_restoration_receipts": restoration_receipts,
            "poses": pose_reports,
        })

    if len(flat_rows) != 40:
        raise ValueError("EXACTLY_FORTY_ACTUAL_ANALOG_POSE_STATES_REQUIRED")

    comparison = {
        "format": "analog-state-comparison-export/20261001.2",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "verified_job_id": args.job_id,
        "verified_job_state": _field(job, "state"),
        "verified_binding": verified_binding,
        "verified_source": verified_source_json,
        "verified_load_followup": verified_load_followup,
        "expert_followup_source": {
            "path": followup_copy.relative_to(destination).as_posix(),
            "original_path": expert_followup.SOURCE_RELATIVE_PATH,
            "sha256": _sha(followup_copy),
            "bytes": followup_copy.stat().st_size,
        },
        "source_archive_sha256": _bytes_sha(archive_raw),
        "source_original_copy": "sources/verified-design-archive.json",
        "parent_source": {
            "path": f"sources/{neutral_path.name}", "sha256": _sha(neutral_path)
        },
        "protein_source": {
            "path": f"sources/{protein_path.name}", "sha256": _sha(protein_path)
        },
        "pose_policy": "all five original docking poses per selected analog",
        "scientific_best_pose_claimed": False,
        "analogs": analog_reports,
        "rows": flat_rows,
        "selected_state": None,
        "microstate_gate": "pending",
        "core_final": "pending",
        "formal_status": "science_review_pending",
        "scientific_approved": False,
        "stored_job_mutated": False,
        "network_called": False,
    }
    _write_json(destination / "comparison.json", comparison)
    (destination / "index.html").write_text(
        _index(comparison), encoding="utf-8", newline="\n"
    )

    artifacts = []
    for path in sorted(destination.rglob("*")):
        if path.is_file() and not path.is_symlink() and path.name != "manifest.json":
            artifacts.append({
                "path": path.relative_to(destination).as_posix(),
                "sha256": _sha(path), "bytes": path.stat().st_size,
            })
    source_files = [
        {"path": f"sources/{neutral_path.name}", "sha256": _sha(neutral_path)},
        {"path": f"sources/{protein_path.name}", "sha256": _sha(protein_path)},
        {"path": "sources/verified-design-archive.json", "sha256": _bytes_sha(archive_raw)},
        {
            "path": followup_copy.relative_to(destination).as_posix(),
            "original_path": expert_followup.SOURCE_RELATIVE_PATH,
            "sha256": _sha(followup_copy),
            "bytes": followup_copy.stat().st_size,
        },
    ]
    for analog in analog_reports:
        source_files.extend([
            {
                "path": analog["source_original_copy"],
                "artifact_id": analog["source_artifact_id"],
                "sha256": analog["source_artifact_sha256"],
            },
            {
                "path": analog["source_original_pdbqt_copy"],
                "artifact_id": analog["source_pdbqt_artifact_id"],
                "sha256": analog["source_pdbqt_artifact_sha256"],
            },
        ])
    manifest = {
        "format": "analog-state-comparison-manifest/20261001.2",
        "source_archive_sha256": _bytes_sha(archive_raw),
        "verified_load_followup": verified_load_followup,
        "source_files": source_files,
        "artifact_hash_scope": "every output file except manifest.json itself",
        "artifacts": artifacts,
        "stored_job_mutated": False,
        "formal_status": "science_review_pending",
    }
    _write_json(destination / "manifest.json", manifest)
    print(json.dumps({
        "status": "exported_pending_review",
        "output": str(destination),
        "comparison": str(destination / "comparison.json"),
    }, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
