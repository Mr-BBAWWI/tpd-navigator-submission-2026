#!/usr/bin/env python3
"""Export a portable offline expert-evidence packet from verified local inputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import os
import re
import shutil
import sys
from pathlib import Path
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import Chem

from packages.science import novel_ternary
from packages.science.parent_fingerprint import build_parent_fingerprint
from packages.science.structures import atom_sites
from packages.science.synthesis_proposals import propose_synthesis, scheme_svg
from scripts.run_novel_ternary import _candidate

DEFAULT_NEUTRAL = ROOT / "cases" / "design_sources" / "SMARCA2-neutral-design.sdf"
DEFAULT_CHARGED = ROOT / "cases" / "design_sources" / "SMARCA2-parent.sdf"
DEFAULT_PROTEIN = ROOT / "cases" / "design_sources" / "6HAZ.cif"

PARENT_ID = "SMARCA2-FX5"
STRUCTURAL_FRAME = "6HAZ chain A"
PROTEIN_CHAIN = "A"
ANALOG_IDS = ("W-c2afc5e73c1a", "W-4c0a639c0a41", "W-80f8f4a11b5d")
E3_ORDER = ("VHL", "CRBN")
LINKER_ID = "alkyl_c6"
FIXED_ORIGINAL_DOCX_SHA256 = "42a78a5540f5549347bc8ee2af1fd9758cec43056c6e808f89ff9aedc6d521cb"
EXPECTED_SOURCE_HASHES = {
    "SMARCA2-neutral-design.sdf": "50e689f68d455083f2d7432b5928a93f6f7f971c4f2fc3cfda911bc0164ac6ee",
    "SMARCA2-parent.sdf": "7b39c8e3f36f1c7a94c02a4ddda6c1cda1d672d3746094bb8572b81cd298e9b9",
    "6HAZ.cif": "7752dda8c54b7e56ab45c84e03ca544e383d834fb30e4c0849fc3544278e2c9d",
}
SOURCE_BINDING_PATHS = {
    name: f"design_sources/{name}" for name in EXPECTED_SOURCE_HASHES
}
CONTACT_FIELDS = (
    "protein_atom_id",
    "chain",
    "auth_seq_id",
    "residue",
    "protein_atom",
    "ligand_atom_map",
    "interaction_kind",
    "heavy_distance_A",
    "directional_hydrogen",
)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bytes_sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _json_fingerprint(value) -> str:
    encoded = json.dumps(
        _jsonable(value),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return _bytes_sha(encoded)


def _field(value, name):
    if isinstance(value, dict):
        return value[name]
    return getattr(value, name)


def _jsonable(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if hasattr(value, "model_dump"):
        return _jsonable(value.model_dump())
    if hasattr(value, "__dict__"):
        return _jsonable(vars(value))
    return str(value)


def _document(value):
    if not isinstance(value, dict):
        raise ValueError("DESIGN_DOCUMENT_OBJECT_REQUIRED")
    if isinstance(value.get("protac_candidates"), list):
        return value
    for key in ("result", "job_result"):
        wrapped = value.get(key)
        if isinstance(wrapped, str):
            wrapped = json.loads(wrapped)
        if isinstance(wrapped, dict) and isinstance(wrapped.get("protac_candidates"), list):
            return wrapped
    raise ValueError("DESIGN_PROTAC_CANDIDATES_REQUIRED")


def _source_input_binding(document: dict) -> dict:
    input_binding = document.get("input_binding")
    if not isinstance(input_binding, dict):
        raise ValueError("OFFLINE_INPUT_BINDING_REQUIRED")
    source_files = input_binding.get("source_files")
    if not isinstance(source_files, dict):
        raise ValueError("OFFLINE_SOURCE_FILES_REQUIRED")

    selected = {}
    for name, expected in EXPECTED_SOURCE_HASHES.items():
        literal_path = SOURCE_BINDING_PATHS[name]
        record = source_files.get(literal_path)
        if not isinstance(record, dict):
            raise ValueError(f"SOURCE_BINDING_REQUIRED:{literal_path}")
        if record.get("status") != "configured_hash_verified":
            raise ValueError(f"SOURCE_BINDING_STATUS_MISMATCH:{literal_path}")
        if record.get("sha256") != expected:
            raise ValueError(f"SOURCE_BINDING_HASH_MISMATCH:{literal_path}")
        selected[literal_path] = {
            "status": record["status"],
            "sha256": record["sha256"],
        }
    return selected


def _verified_design(args) -> tuple[dict, bytes, dict, object]:
    if args.store:
        from packages.platform.scientific_acceptance import ScientificAcceptanceService
        from packages.platform.store import Store

        service = ScientificAcceptanceService(Store(Path(args.store).resolve()), "local-research")
        with service.store.db() as db:
            job, result, binding, source = service._verified_result(
                db, args.job_id, require_current=True
            )
            files = _field(result, "files")
            archive_raw = service.port.read(_field(files, "json"))
        if isinstance(archive_raw, str):
            archive_raw = archive_raw.encode("utf-8")
        if not isinstance(archive_raw, bytes):
            raise TypeError("VERIFIED_ARCHIVE_BYTES_REQUIRED")
        document = _document(json.loads(archive_raw))
        source_input_binding = _source_input_binding(document)
        metadata = {
            "binding_status": "server_verified_stored_job",
            "project": "local-research",
            "job_id": args.job_id,
            "state": _field(job, "state"),
            "binding": _jsonable(binding),
            "source_fingerprint": _json_fingerprint(source),
            "source_input_binding": source_input_binding,
        }
        return document, archive_raw, metadata, source_input_binding

    if not args.design_json:
        raise ValueError("--design-json or --store with --job-id is required")
    path = _checked_file(Path(args.design_json), "design JSON")
    archive_raw = path.read_bytes()
    document = _document(json.loads(archive_raw))
    source_input_binding = _source_input_binding(document)
    metadata = {
        "binding_status": "offline_source_verified",
        "project": None,
        "job_id": None,
        "design_source_path": str(path),
        "result_sha256": document.get("result_sha256", _bytes_sha(archive_raw)),
        "source_fingerprint": _json_fingerprint(source_input_binding),
        "source_input_binding": source_input_binding,
    }
    return document, archive_raw, metadata, source_input_binding


def _has_symlink_component(path: Path) -> bool:
    absolute = Path(os.path.abspath(os.path.expanduser(str(path))))
    current = Path(absolute.anchor)
    for part in absolute.parts[1:]:
        current = current / part
        if current.exists() and current.is_symlink():
            return True
    return False


def _checked_file(path: Path, role: str) -> Path:
    expanded = Path(os.path.abspath(os.path.expanduser(str(path))))
    if _has_symlink_component(expanded):
        raise ValueError(f"SYMLINK_PATH_DISALLOWED:{role}")
    if not expanded.is_file():
        raise FileNotFoundError(expanded)
    return expanded.resolve(strict=True)


def _checked_destination(path: Path) -> Path:
    expanded = Path(os.path.abspath(os.path.expanduser(str(path))))
    if expanded.exists():
        raise FileExistsError(expanded)
    if _has_symlink_component(expanded.parent):
        raise ValueError("SYMLINK_PATH_DISALLOWED:output")
    parent = expanded.parent.resolve(strict=True)
    return parent / expanded.name


def _validate_source_bindings(source, paths: dict[str, Path]) -> list[dict]:
    if not isinstance(source, dict):
        raise ValueError("SOURCE_INPUT_BINDING_OBJECT_REQUIRED")
    rows = []
    for name, expected in EXPECTED_SOURCE_HASHES.items():
        literal_path = SOURCE_BINDING_PATHS[name]
        record = source.get(literal_path)
        if not isinstance(record, dict):
            raise ValueError(f"SOURCE_BINDING_REQUIRED:{literal_path}")
        if record.get("status") != "configured_hash_verified":
            raise ValueError(f"SOURCE_BINDING_STATUS_MISMATCH:{literal_path}")
        if record.get("sha256") != expected:
            raise ValueError(f"SOURCE_BINDING_HASH_MISMATCH:{literal_path}")
        path = paths[name]
        actual = _sha(path)
        if actual != expected:
            raise ValueError(f"ACTUAL_SOURCE_HASH_MISMATCH:{name}")
        rows.append({
            "name": name,
            "source_path": str(path),
            "sha256": actual,
            "binding_status": record["status"],
        })
    return rows


def _mol(path: Path) -> Chem.Mol:
    supplier = Chem.ForwardSDMolSupplier(
        io.BytesIO(path.read_bytes()), removeHs=False, sanitize=True, strictParsing=True
    )
    records = list(supplier)
    if len(records) != 1 or records[0] is None:
        raise ValueError(f"exactly one valid SDF record required in {path}")
    return records[0]


def _protein(path: Path, chain: str) -> list[dict]:
    if chain != PROTEIN_CHAIN:
        raise ValueError("PROTEIN_CHAIN_MUST_BE_A")
    import gemmi

    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    if len(structure) == 0:
        raise ValueError("PROTEIN_MODEL_REQUIRED")
    matching_chains = [candidate for candidate in structure[0] if candidate.name == chain]
    if len(matching_chains) != 1:
        raise ValueError("PROTEIN_CHAIN_NOT_DECLARED")
    polymer = matching_chains[0].get_polymer()
    if not polymer or polymer.check_polymer_type() != gemmi.PolymerType.PeptideL:
        raise ValueError("PROTEIN_CHAIN_NOT_L_POLYPEPTIDE")
    selected = atom_sites(path, [chain])
    if not selected or {str(row["label_asym_id"]) for row in selected} != {chain}:
        raise ValueError("PROTEIN_CHAIN_ATOMS_INVALID")
    result = []
    for row in selected:
        result.append({
            "label_asym_id": str(row["label_asym_id"]),
            "auth_asym_id": str(row["auth_asym_id"]),
            "auth_seq_id": str(row["auth_seq_id"]),
            "label_seq_id": str(row["label_seq_id"]),
            "label_comp_id": str(row["label_comp_id"]),
            "label_atom_id": str(row["label_atom_id"]),
            "type_symbol": str(row["type_symbol"]),
            "xyz": [float(value) for value in row["xyz"]],
        })
    return result


def _linker(candidate: dict):
    value = candidate.get("linker_id", candidate.get("linker"))
    if isinstance(value, dict):
        return value.get("id", value.get("linker_id"))
    return value


def _representatives(document: dict) -> list[dict]:
    scope = document.get("parent_scope")
    if not isinstance(scope, dict):
        raise ValueError("PRIORITY_PARENT_SCOPE_REQUIRED")
    parent = scope.get("actual_design_parent_id")
    frame = scope.get("receptor_frame")
    if parent != PARENT_ID:
        raise ValueError("PRIORITY_PARENT_SCOPE_MISMATCH")
    if frame != STRUCTURAL_FRAME:
        raise ValueError("PRIORITY_STRUCTURAL_FRAME_MISMATCH")

    analogs = document.get("analogs")
    candidates = document.get("protac_candidates")
    if not isinstance(analogs, list) or not isinstance(candidates, list):
        raise ValueError("STRICT_DESIGN_LISTS_REQUIRED")
    selected = [row for row in analogs if isinstance(row, dict) and row.get("selected") is True]
    selected_ids = [row.get("id") for row in selected]
    if len(selected_ids) != 3 or set(selected_ids) != set(ANALOG_IDS):
        raise ValueError("EXPERT_SELECTED_ANALOG_SET_MISMATCH")
    if len(set(selected_ids)) != len(selected_ids):
        raise ValueError("EXPERT_SELECTED_ANALOG_DUPLICATE")
    for row in selected:
        if row.get("qualified") is False or row.get("pose_passed") is False:
            raise ValueError("EXPERT_SELECTED_ANALOG_NOT_QUALIFIED")

    result = []
    candidate_ids = set()
    graph_hashes = set()
    for analog_id in ANALOG_IDS:
        for e3 in E3_ORDER:
            matches = [
                row for row in candidates
                if isinstance(row, dict)
                and row.get("warhead_analog_id") == analog_id
                and row.get("e3_type") == e3
                and _linker(row) == LINKER_ID
            ]
            if len(matches) != 1:
                raise ValueError("PRIORITY_CANDIDATE_NOT_UNIQUE")
            identifier = matches[0].get("candidate_id")
            if not isinstance(identifier, str) or not identifier or identifier in candidate_ids:
                raise ValueError("PRIORITY_CANDIDATE_ID_INVALID_OR_DUPLICATE")
            checked = _candidate(document, identifier)
            if checked.get("warhead_analog_id") != analog_id or checked.get("e3_type") != e3:
                raise ValueError("PRIORITY_CANDIDATE_BINDING_MISMATCH")
            if _linker(checked) != LINKER_ID:
                raise ValueError("PRIORITY_LINKER_MISMATCH")
            graph = novel_ternary._mapped_graph(checked)
            graph_hash = graph.get("graph_sha256")
            if not isinstance(graph_hash, str) or not graph_hash or graph_hash in graph_hashes:
                raise ValueError("PRIORITY_CANDIDATE_GRAPH_INVALID_OR_DUPLICATE")
            candidate_ids.add(identifier)
            graph_hashes.add(graph_hash)
            result.append(checked)
    if len(result) != 6:
        raise ValueError("EXACTLY_SIX_REPRESENTATIVES_REQUIRED")
    return result


def _safe_name(value: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError("CANDIDATE_ID_REQUIRED")
    result = re.sub(r"[^A-Za-z0-9_-]", "_", value)[:64]
    if not result:
        raise ValueError("CANDIDATE_FILENAME_EMPTY")
    return result


def _contact_rows(parent: dict) -> list[dict]:
    """Flatten only contacts present in the parent-fingerprint schema."""
    rows = []

    def exact_contact(contact):
        if not isinstance(contact, dict):
            raise ValueError("PARENT_CONTACT_OBJECT_REQUIRED")
        missing = [field for field in CONTACT_FIELDS if field not in contact]
        if missing:
            raise ValueError("PARENT_CONTACT_FIELDS_REQUIRED:" + ",".join(missing))
        return {field: contact[field] for field in CONTACT_FIELDS}

    parents = parent.get("parents")
    if not isinstance(parents, list):
        raise ValueError("PARENT_FINGERPRINT_PARENTS_REQUIRED")
    for parent_record in parents:
        if not isinstance(parent_record, dict):
            raise ValueError("PARENT_FINGERPRINT_PARENT_OBJECT_REQUIRED")
        parent_role = parent_record.get("parent_role")
        if not isinstance(parent_role, str) or not parent_role:
            raise ValueError("PARENT_ROLE_REQUIRED")

        states = parent_record.get("states")
        if not isinstance(states, list):
            raise ValueError("PARENT_STATES_REQUIRED")
        for state in states:
            if not isinstance(state, dict) or "state_index" not in state:
                raise ValueError("PARENT_STATE_INDEX_REQUIRED")
            contacts = state.get("contacts")
            if not isinstance(contacts, list):
                raise ValueError("PARENT_STATE_CONTACTS_REQUIRED")
            for contact in contacts:
                rows.append({
                    "scope": "state",
                    "parent_role": parent_role,
                    "state_index": state["state_index"],
                    **exact_contact(contact),
                })

        proposal = parent_record.get("proposed_protected_interaction_set")
        if not isinstance(proposal, dict):
            raise ValueError("PARENT_PROPOSED_PROTECTED_INTERACTION_SET_REQUIRED")
        proposal_contacts = proposal.get("contacts")
        if not isinstance(proposal_contacts, list):
            raise ValueError("PARENT_PROPOSED_PROTECTED_CONTACTS_REQUIRED")
        for contact in proposal_contacts:
            rows.append({
                "scope": "proposed_protected_interaction_set",
                "parent_role": parent_role,
                "state_index": "common",
                **exact_contact(contact),
            })
    return rows


def _write_contacts(path: Path, parent: dict) -> None:
    columns = (
        "scope",
        "parent_role",
        "state_index",
        *CONTACT_FIELDS,
    )
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(_contact_rows(parent))


def _index(rows: list[dict]) -> str:
    links = []
    for row in rows:
        candidate = html.escape(row["candidate_id"], quote=True)
        report = html.escape("synthesis/" + row["report"], quote=True)
        scheme = html.escape("synthesis/" + row["scheme"], quote=True)
        links.append(
            f"<li><strong>{candidate}</strong> — "
            f"<a href=\"{report}\">합성 제안 JSON</a> · "
            f"<a href=\"{scheme}\">반응식 SVG</a></li>"
        )
    return """<!doctype html>
<html lang="ko"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>전문가 근거 패킷</title>
<style>body{font-family:system-ui,sans-serif;max-width:960px;margin:2rem auto;padding:0 1rem;line-height:1.6;color:#17212b}h1,h2{color:#173f5f}.pending{border-left:5px solid #d98e04;background:#fff8df;padding:1rem}li{margin:.55rem 0}code{background:#eef2f5;padding:.1rem .3rem}</style>
</head><body><h1>SMARCA2-FX5 전문가 검토용 계산 근거</h1>
<div class="pending"><strong>과학 검토 대기 / 승인 아님.</strong> 이 패킷은 고정된 6HAZ chain A 프레임의 기술적 근거와 합성 제안을 보존합니다. 효능, 분해능, 결합 우열, 합성 성공 또는 전문가 승인을 주장하지 않습니다.</div>
<h2>부모 구조 근거</h2><ul>
<li><a href="parent/index.html">부모 fingerprint 보고서</a></li>
<li><a href="parent-residue-contacts.csv">상태별 residue contact CSV</a></li>
</ul><h2>대표 합성 제안 6건</h2><p>선택 범위의 정확한 세 warhead와 각 VHL/CRBN 조합에 대해 개발자가 지정한 대표 linker <code>alkyl_c6</code>를 포함합니다. 제시된 합성 경로는 검증되지 않았습니다.</p><ul>""" + "".join(links) + """</ul>
<h2>추적성</h2><ul><li><a href="sources/design.json">원본 design archive</a></li><li><a href="manifest.json">SHA-256 manifest</a></li></ul>
</body></html>
"""


def _validate_local_links(destination: Path) -> None:
    from html.parser import HTMLParser

    class Links(HTMLParser):
        def __init__(self):
            super().__init__()
            self.values = []

        def handle_starttag(self, tag, attrs):
            for key, value in attrs:
                if key == "href" and value is not None:
                    self.values.append(value)

    for page in destination.rglob("*.html"):
        parser = Links()
        parser.feed(page.read_text(encoding="utf-8"))
        for href in parser.values:
            parsed = urlsplit(href)
            if parsed.scheme or parsed.netloc or href.startswith(("//", "/", "\\")):
                raise ValueError(f"EXTERNAL_OR_ABSOLUTE_LINK:{href}")
            relative = parsed.path
            if not relative:
                continue
            target = page.parent / relative
            try:
                resolved = target.resolve(strict=True)
                resolved.relative_to(destination.resolve(strict=True))
            except (OSError, ValueError) as exc:
                raise ValueError(f"BROKEN_OR_ESCAPING_LINK:{href}") from exc
            if target.is_symlink():
                raise ValueError(f"SYMLINK_LINK_TARGET:{href}")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--design-json")
    result.add_argument("--store", help="ScientificAcceptanceService Store path")
    result.add_argument("--job-id")
    result.add_argument("--output", required=True)
    result.add_argument("--neutral-sdf", default=str(DEFAULT_NEUTRAL))
    result.add_argument("--charged-sdf", default=str(DEFAULT_CHARGED))
    result.add_argument("--protein-cif", default=str(DEFAULT_PROTEIN))
    result.add_argument("--protein-chain", default=PROTEIN_CHAIN)
    result.add_argument("--original-docx")
    result.add_argument(
        "--expected-original-docx-sha256",
        default=FIXED_ORIGINAL_DOCX_SHA256,
        help=argparse.SUPPRESS,
    )
    result.add_argument("--max-states", type=int, default=8)
    return result


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    if bool(args.store) != bool(args.job_id):
        raise ValueError("--store and --job-id must be supplied together")
    if args.design_json and args.store:
        raise ValueError("DESIGN_SOURCE_MODES_ARE_MUTUALLY_EXCLUSIVE")
    if args.protein_chain != PROTEIN_CHAIN:
        raise ValueError("PROTEIN_CHAIN_MUST_BE_A")
    if args.max_states < 1 or args.max_states > 16:
        raise ValueError("MAX_STATES_OUT_OF_RANGE")
    if args.expected_original_docx_sha256 != FIXED_ORIGINAL_DOCX_SHA256:
        raise ValueError("ORIGINAL_DOCX_FIXED_HASH_CANNOT_BE_OVERRIDDEN")

    destination = _checked_destination(Path(args.output))
    neutral_path = _checked_file(Path(args.neutral_sdf), "neutral SDF")
    charged_path = _checked_file(Path(args.charged_sdf), "charged SDF")
    protein_path = _checked_file(Path(args.protein_cif), "protein CIF")
    document, design_raw, design_binding, source_binding = _verified_design(args)

    paths = {
        "SMARCA2-neutral-design.sdf": neutral_path,
        "SMARCA2-parent.sdf": charged_path,
        "6HAZ.cif": protein_path,
    }
    source_rows = _validate_source_bindings(source_binding, paths)
    representatives = _representatives(document)
    neutral_mol = _mol(neutral_path)
    charged_mol = _mol(charged_path)
    protein_atoms = _protein(protein_path, args.protein_chain)

    original_path = None
    if args.original_docx:
        original_path = _checked_file(Path(args.original_docx), "original DOCX")
        if _sha(original_path) != FIXED_ORIGINAL_DOCX_SHA256:
            raise ValueError("ORIGINAL_DOCX_FIXED_HASH_MISMATCH")

    stems = [_safe_name(str(candidate["candidate_id"])) for candidate in representatives]
    if len(stems) != len(set(stems)):
        raise ValueError("CANDIDATE_FILENAME_COLLISION")

    destination.mkdir()
    source_dir = destination / "sources"
    source_dir.mkdir()
    copied_sources = []
    for row in source_rows:
        source_path = Path(row["source_path"])
        target = source_dir / row["name"]
        shutil.copyfile(source_path, target)
        if _sha(target) != row["sha256"]:
            raise ValueError("SOURCE_COPY_HASH_MISMATCH")
        copied_sources.append({
            **row,
            "packet_path": target.relative_to(destination).as_posix(),
        })
    design_target = source_dir / "design.json"
    design_target.write_bytes(design_raw)
    if design_target.read_bytes() != design_raw:
        raise ValueError("DESIGN_ARCHIVE_COPY_MISMATCH")

    original_manifest = None
    if original_path is not None:
        original_target = source_dir / "original-opinion.docx"
        shutil.copyfile(original_path, original_target)
        original_manifest = {
            "source_path": str(original_path),
            "packet_path": original_target.relative_to(destination).as_posix(),
            "sha256": _sha(original_target),
            "role": "original_opinion_reference_only",
            "creates_or_changes_approval": False,
        }

    parent = build_parent_fingerprint(
        charged_mol,
        neutral_mol,
        protein_atoms,
        protein_hydrogens=None,
        output_dir=destination / "parent",
        pH_context=7.4,
        max_states=args.max_states,
    )
    _write_contacts(destination / "parent-residue-contacts.csv", parent)

    synthesis_dir = destination / "synthesis"
    synthesis_dir.mkdir()
    synth_rows = []
    for candidate, stem in zip(representatives, stems):
        identifier = str(candidate["candidate_id"])
        report_name = f"{stem}.json"
        scheme_name = f"scheme-{stem}.svg"
        report = propose_synthesis(candidate, exact_route_records=None, constraints=None)
        (synthesis_dir / report_name).write_text(
            json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        (synthesis_dir / scheme_name).write_text(
            scheme_svg(candidate), encoding="utf-8", newline="\n"
        )
        synth_rows.append({
            "candidate_id": identifier,
            "warhead_analog_id": candidate["warhead_analog_id"],
            "e3_type": candidate["e3_type"],
            "linker_id": LINKER_ID,
            "candidate_graph_sha256": novel_ternary._mapped_graph(candidate)["graph_sha256"],
            "report": report_name,
            "scheme": scheme_name,
        })

    (destination / "index.html").write_text(
        _index(synth_rows), encoding="utf-8", newline="\n"
    )

    producer_paths = [
        Path(__file__).resolve(),
        ROOT / "packages" / "science" / "parent_fingerprint.py",
        ROOT / "packages" / "science" / "structures.py",
        ROOT / "packages" / "science" / "synthesis_proposals.py",
        ROOT / "scripts" / "run_novel_ternary.py",
    ]
    producer_files = []
    for producer in producer_paths:
        checked = _checked_file(producer, "producer file")
        producer_files.append({"path": str(checked), "sha256": _sha(checked)})

    artifacts = []
    for path in sorted(destination.rglob("*")):
        if path.is_file() and not path.is_symlink() and path.name != "manifest.json":
            artifacts.append({
                "path": path.relative_to(destination).as_posix(),
                "sha256": _sha(path),
                "bytes": path.stat().st_size,
            })

    manifest = {
        "format": "expert-evidence-packet/20261002.2",
        "binding_status": design_binding["binding_status"],
        "design_binding": design_binding,
        "design_archive": {
            "packet_path": "sources/design.json",
            "sha256": _bytes_sha(design_raw),
            "bytes": len(design_raw),
            "exact_archive_bytes_preserved": True,
        },
        "source_files": copied_sources,
        "original_opinion": original_manifest,
        "scope": {
            "parent_id": PARENT_ID,
            "structural_frame": STRUCTURAL_FRAME,
            "protein_chain": PROTEIN_CHAIN,
            "selected_analog_ids": list(ANALOG_IDS),
            "representative_linker": LINKER_ID,
            "representative_count": 6,
        },
        "representatives": synth_rows,
        "parent_summary": {
            "protein_atom_count": len(protein_atoms),
            "protein_hydrogens_supplied": parent.get("protein_hydrogens_supplied"),
            "pH_context": parent.get("pH_context", 7.4),
            "fixed_receptor_alignment_hypothesis": True,
            "scientific_approval": False,
        },
        "claims": {
            "readiness_state": "science_review_pending",
            "expert_approved": False,
            "scientifically_approved": False,
            "efficacy_established": False,
            "degradation_established": False,
            "synthesis_success_established": False,
            "status": "science_review_pending",
        },
        "network_called": False,
        "stored_job_mutated": False,
        "producer_files": producer_files,
        "artifact_hash_scope": "every packet file except manifest.json, whose self-hash cannot be embedded immutably",
        "artifacts": artifacts,
    }
    (destination / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    _validate_local_links(destination)
    print(json.dumps({
        "status": "exported",
        "output": str(destination),
        "manifest": str(destination / "manifest.json"),
    }, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
