#!/usr/bin/env python3
"""Build a verified, portable, curated native-parent review pack."""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import re
import sys
from pathlib import Path, PurePosixPath
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

PARENT_ID = "SMARCA2-9D12-A1A1P"
DOI = "10.1021/acs.jmedchem.4c01903"
DEFAULT_PANEL = ROOT / ".localdata/final-sprint-20261001/parent/native-n3-panel-v1"
DEFAULT_BASELINE = ROOT / ".localdata/final-sprint-20261001/parent/native-explicit-remote-v3"
DEFAULT_AUDIT = ROOT / ".localdata/final-sprint-20261001/parent/source-pair-audit-v1"
DEFAULT_EXPORT = ROOT / ".localdata/expert-final-20261001/parent-exports" / PARENT_ID
DEFAULT_CIF = ROOT / ".localdata/acceptance-20260930/sources/9D12.cif"
DEFAULT_CSV = ROOT / ".localdata/acceptance-20260930/sources/jm4c01903_si_006.csv"
DEFAULT_DOC = ROOT / "docs/NATIVE_PARENT_SAR_REVIEW_20261001.md"
MANIFEST = "manifest.json"
BLOCKED_PARTS = {".git", ".venv", "venv", "__pycache__", ".cache", "cache", "node_modules"}
SECRET_NAMES = {".env", "credentials", "credentials.json", "keys", "keys.json"}


def _safe_resolved_path(value: Path, label: str) -> Path:
    """Reject symlinks in a lexical path before resolving it."""
    absolute = Path(os.path.abspath(os.fspath(value)))
    for candidate in (absolute, *absolute.parents):
        if candidate.is_symlink():
            raise ValueError(f"{label} path contains a symlink: {candidate}")
    return absolute.resolve()


def _reject_output_source_overlap(output: Path, sources: list[Path]) -> None:
    for source in sources:
        if output == source or output.is_relative_to(source) or source.is_relative_to(output):
            raise ValueError(f"Output and source paths must not contain one another: {output} / {source}")


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _read(path: Path) -> bytes:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Unsafe or missing file: {path}")
    return path.read_bytes()


def _json(path: Path) -> Any:
    return json.loads(_read(path).decode("utf-8"))


def _safe_relative(value: Any) -> PurePosixPath:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError("Manifest path must be a non-empty POSIX path")
    if re.match(r"^[A-Za-z]:($|/)", value):
        raise ValueError(f"Windows drive path is forbidden in manifest: {value!r}")
    rel = PurePosixPath(value)
    if rel.is_absolute() or any(part in ("", ".", "..") for part in rel.parts):
        raise ValueError(f"Unsafe manifest path: {value!r}")
    if rel.as_posix() == MANIFEST:
        raise ValueError("Manifest cannot declare itself")
    return rel


def _check_name(rel: PurePosixPath) -> None:
    lowered = [part.lower() for part in rel.parts]
    if any(part in BLOCKED_PARTS for part in lowered):
        raise ValueError(f"Forbidden cache/environment path: {rel}")
    if any(part in SECRET_NAMES or part.startswith(".env.") for part in lowered):
        raise ValueError(f"Possible secret path is forbidden: {rel}")


def _actual_files(root: Path, exclude_manifest: bool = True) -> dict[str, Path]:
    root = _safe_resolved_path(root, "Directory")
    result: dict[str, Path] = {}
    if not root.is_dir():
        raise ValueError(f"Unsafe or missing directory: {root}")
    for path in sorted(root.rglob("*")):
        rel = PurePosixPath(path.relative_to(root).as_posix())
        _check_name(rel)
        if path.is_symlink():
            raise ValueError(f"Symlink forbidden: {path}")
        if path.is_file() and not (exclude_manifest and rel.as_posix() == MANIFEST):
            result[rel.as_posix()] = path
    return result


def verify_manifest_tree(root: Path, required: set[str] | None = None) -> dict[str, Any]:
    """Verify byte/hash rows and exact all-file closure for a manifest directory."""
    root = _safe_resolved_path(Path(root), "Manifest root")
    manifest_path = root / MANIFEST
    if manifest_path.is_symlink() or not manifest_path.is_file():
        raise ValueError(f"Missing safe manifest: {manifest_path}")
    raw = _read(manifest_path)
    manifest = json.loads(raw.decode("utf-8"))
    rows = manifest.get("files")
    if not isinstance(rows, list):
        raise ValueError("Manifest files must be a list")
    declared: dict[str, dict[str, Any]] = {}
    for row in rows:
        if not isinstance(row, dict) or set(row) != {"path", "bytes", "sha256"}:
            raise ValueError("Malformed manifest row")
        rel = _safe_relative(row["path"])
        _check_name(rel)
        name = rel.as_posix()
        if name in declared:
            raise ValueError(f"Duplicate manifest path: {name}")
        path = root.joinpath(*rel.parts)
        cursor = root
        for part in rel.parts:
            cursor /= part
            if cursor.is_symlink():
                raise ValueError(f"Symlink forbidden: {name}")
        if not path.is_file() or not path.resolve().is_relative_to(root):
            raise ValueError(f"Manifest path missing or escaping: {name}")
        data = _read(path)
        if type(row["bytes"]) is not int or row["bytes"] != len(data) or row["sha256"] != _sha(data):
            raise ValueError(f"Manifest mismatch: {name}")
        declared[name] = row
    actual = _actual_files(root)
    if set(actual) != set(declared):
        raise ValueError("Manifest does not cover all and only files")
    if required and not required <= set(declared):
        raise ValueError(f"Manifest lacks required files: {sorted(required - set(declared))}")
    return {"root": root, "manifest": manifest, "manifest_sha256": _sha(raw), "files": actual}


def _source_pair_verifier(source_export: Path) -> Any:
    from packages.science.native_parent_docking import verify_source_pair
    return verify_source_pair(source_export, PARENT_ID)


def _native_probe_verifier(path: Path) -> Any:
    from packages.science.native_attachment_experiment import verify_native_probe
    return verify_native_probe(path)


def _copy_verified_tree(info: dict[str, Any], destination: Path) -> list[dict[str, Any]]:
    destination.mkdir(parents=True, exist_ok=False)
    records = []
    root: Path = info["root"]
    paths = dict(info["files"])
    paths[MANIFEST] = root / MANIFEST
    for rel, source in sorted(paths.items()):
        data = _read(source)
        before = _sha(data)
        target = destination.joinpath(*PurePosixPath(rel).parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        copied = _read(target)
        if copied != data or _sha(_read(source)) != before:
            raise ValueError(f"Source changed or copy mismatch: {source}")
        records.append({"path": target.as_posix(), "bytes": len(data), "sha256": before})
    return records


def _copy_unmanifested_verified(source: Path, destination: Path) -> list[dict[str, Any]]:
    root = _safe_resolved_path(Path(source), "Unmanifested source")
    files = _actual_files(root, exclude_manifest=False)
    destination.mkdir(parents=True, exist_ok=False)
    records = []
    for rel, path in sorted(files.items()):
        data = _read(path)
        digest = _sha(data)
        target = destination.joinpath(*PurePosixPath(rel).parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        if _read(target) != data or _sha(_read(path)) != digest:
            raise ValueError(f"Source changed or copy mismatch: {path}")
        records.append({"path": target.as_posix(), "bytes": len(data), "sha256": digest})
    return records


def _copy_file(source: Path, target: Path) -> dict[str, Any]:
    source = _safe_resolved_path(Path(source), "File source")
    data = _read(source)
    digest = _sha(data)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    if _read(target) != data or _sha(_read(source)) != digest:
        raise ValueError(f"Source changed or copy mismatch: {source}")
    return {"path": target.as_posix(), "bytes": len(data), "sha256": digest}


def _case_path(case: dict[str, Any], case_id: str, value: Any) -> str | None:
    if not isinstance(value, str) or not value:
        return None
    rel = _safe_relative(value)
    path = PurePosixPath("n3-panel/cases") / case_id / rel
    return path.as_posix()


def _link(path: str | None, label: str) -> str:
    return f'<a href="{html.escape(path, quote=True)}">{html.escape(label)}</a>' if path else "—"


def _fmt(value: Any, digits: int = 3) -> str:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return f"{value:.{digits}f}"
    return "—"


def _molecule_svg(mapped_smiles: Any) -> str:
    from rdkit import Chem
    from rdkit.Chem import Draw
    if not isinstance(mapped_smiles, str):
        raise ValueError("Case mapped_smiles is missing")
    mol = Chem.MolFromSmiles(mapped_smiles)
    if mol is None:
        raise ValueError("Case mapped_smiles is invalid")
    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    drawer = Draw.MolDraw2DSVG(360, 220)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = re.sub(r"<\?xml[^>]*>\s*", "", svg)
    svg = re.sub(r"<!DOCTYPE[^>]*>\s*", "", svg)
    return svg


def _panel_data(panel_root: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    panel = _json(panel_root / "panel.json")
    protocol = _json(panel_root / "protocol.json")
    compact = _json(panel_root / "compact-summary.json")
    cases = panel.get("cases")
    requested = panel.get("requested_count")
    if panel.get("status") != "completed_with_limits" or requested != 12 or panel.get("retained_count") != 12:
        raise ValueError("Panel must be completed_with_limits with requested/retained count 12")
    if not isinstance(cases, list) or len(cases) != requested or protocol.get("requested_candidate_count") != 12:
        raise ValueError("Panel/protocol must contain exactly twelve cases")
    ids = [row.get("id") for row in cases if isinstance(row, dict)]
    protocol_ids = [row.get("id") for row in protocol.get("cases", []) if isinstance(row, dict)]
    if len(ids) != 12 or len(set(ids)) != 12 or ids != protocol_ids:
        raise ValueError("Panel and protocol case identities differ")
    return panel, protocol, compact


def _calibration_counts(compact: dict[str, Any]) -> dict[str, tuple[int, int]]:
    counts: dict[str, list[int]] = {}
    for run in compact.get("native_parent_and_Br_runs", []):
        if not isinstance(run, dict):
            continue
        ligand = str(run.get("ligand", "unknown"))
        passed, total = counts.setdefault(ligand, [0, 0])
        counts[ligand] = [passed + int((run.get("pose_pass_count", 0) or 0) > 0), total + 1]
    return {ligand: (values[0], values[1]) for ligand, values in counts.items()}


def _calibration_label(compact: dict[str, Any]) -> str:
    counts = _calibration_counts(compact)
    parent = counts.get("native_parent_Cl", (0, 0))
    br = counts.get("measured_probe_Br", (0, 0))
    return f"Parent {parent[0]}/{parent[1]}, Br {br[0]}/{br[1]} seed-level pass."


def _render_html(panel: dict[str, Any], compact: dict[str, Any]) -> str:
    cases = panel["cases"]
    pass_cases = sum(int((c.get("pose_metrics", {}).get("pose_preservation_pass_count", 0) or 0) > 0) for c in cases)
    passes = sum(int(c.get("pose_metrics", {}).get("pose_preservation_pass_count", 0) or 0) for c in cases)
    assemblies = sum(len(c.get("assemblies", [])) for c in cases)
    cards = []
    for case in cases:
        cid = str(case["id"])
        metrics = case.get("pose_metrics", {})
        branches = {b.get("e3_type"): b for b in case.get("assemblies", []) if isinstance(b, dict)}
        links = [_link(_case_path(case, cid, case.get("candidate_sdf")), "candidate SDF")]
        if case.get("docking") is not None:
            links.append(_link(f"n3-panel/cases/{html.escape(cid, quote=True)}/dock-receipt.json", "dock receipt"))
        for e3 in ("CRBN", "VHL"):
            branch = branches.get(e3)
            links.append(_link(_case_path(case, cid, branch.get("sdf") if branch else None), f"{e3} SDF"))
        cards.append(f'''<article class="case"><h3>{html.escape(cid)}</h3>
<div class="mol">{_molecule_svg(case.get("mapped_smiles"))}</div>
<p class="status">{html.escape(str(case.get("status", "unknown")))}</p>
<dl><dt>candidate</dt><dd>{html.escape(str(case.get("candidate_id", "—")))}</dd>
<dt>best geometry RMSD Å</dt><dd>{_fmt(metrics.get("best_geometry_core_RMSD_A"))}</dd>
<dt>same-pose contact</dt><dd>{_fmt(metrics.get("best_geometry_contact_retention"))}</dd>
<dt>pass poses</dt><dd>{html.escape(str(metrics.get("pose_preservation_pass_count", 0)))}</dd>
<dt>assemblies</dt><dd>{html.escape(str(len(case.get("assemblies", []))))}</dd></dl>
<p class="links">{' · '.join(links)}</p></article>''')
    runs = compact.get("native_parent_and_Br_runs", [])
    run_rows = "".join(f"<tr><td>{html.escape(str(r.get('ligand')))}</td><td>{html.escape(str(r.get('seed')))}</td><td>{html.escape(str(r.get('status')))}</td><td>{html.escape(str(r.get('pose_pass_count')))}</td><td>{_fmt(r.get('minimum_core_RMSD_A'), 4)}</td></tr>" for r in runs)
    return f'''<!doctype html><html lang="ko"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Native parent N3 12-case review</title><style>
:root{{--paper:#f5f0e6;--ink:#1d2925;--green:#1f5c4a;--line:#c9c0ae}}*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.5 system-ui,sans-serif}}main{{max-width:1240px;margin:auto;padding:32px}}h1,h2,h3{{font-family:Georgia,serif}}.hero,.case,.note{{background:#fffdf8;border:1px solid var(--line);border-radius:8px;padding:18px}}.metrics{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}}.metric b{{display:block;font-size:28px;color:var(--green)}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:16px}}.mol svg{{width:100%;height:220px}}dl{{display:grid;grid-template-columns:1fr 1fr}}dt{{font-weight:700}}dd{{margin:0;text-align:right}}table{{width:100%;border-collapse:collapse;background:#fffdf8}}th,td{{padding:8px;border-bottom:1px solid var(--line);text-align:left}}a{{color:var(--green)}}code{{overflow-wrap:anywhere}}@media(max-width:700px){{.metrics{{grid-template-columns:1fr 1fr}}main{{padding:14px}}}}
</style></head><body><main><section class="hero"><p>검토용 · exploratory source-extrapolation · strict count 0</p><h1>Native parent N3 attachment panel</h1><div class="metrics"><div class="metric"><b>{len(cases)}</b>completed/retained</div><div class="metric"><b>{pass_cases}/{len(cases)}</b>pose-pass cases</div><div class="metric"><b>{passes}</b>total pass poses</div><div class="metric"><b>{assemblies}</b>E3 branches</div></div></section>
<section><h2>Source evidence</h2><div class="note">Raw source map 19 → selected-parent map 3 (N3), state <b>UNKNOWN</b>. 이는 source extrapolation이며 measured potency transfer가 아니다. DOI <a href="https://doi.org/{DOI}">{DOI}</a> · <a href="sources/9D12.cif">9D12 CIF</a> · <a href="sources/jm4c01903_si_006.csv">primary CSV</a> · <a href="source-pair-audit/evidence-audit.json">source-pair audit</a></div></section>
<section><h2>12 cases</h2><div class="grid">{''.join(cards)}</div></section>
<section><h2>Native parent / Br calibration</h2><p>{html.escape(_calibration_label(compact))} 6HAZ 개선 주장은 하지 않는다: 비교 protocol과 atom mask가 다르다.</p><table><thead><tr><th>ligand</th><th>seed</th><th>status</th><th>pass poses</th><th>min RMSD Å</th></tr></thead><tbody>{run_rows}</tbody></table></section>
<section><h2>Limits</h2><div class="note">실패·filtered·all pose 자료는 원본 패널에 그대로 보존된다. Geometry와 docking heuristic은 결합력·효능·승인·합성경로 또는 사람 승인이 아니다. Scientific review remains pending.</div></section></main></body></html>'''


def _readme(panel: dict[str, Any], compact: dict[str, Any]) -> str:
    cases = panel["cases"]
    pass_cases = sum((c.get("pose_metrics", {}).get("pose_preservation_pass_count", 0) or 0) > 0 for c in cases)
    pass_poses = sum(c.get("pose_metrics", {}).get("pose_preservation_pass_count", 0) or 0 for c in cases)
    branches = sum(len(c.get("assemblies", [])) for c in cases)
    filtered = sum(c.get("status") == "filtered_not_replaced" for c in cases)
    docking_failed = sum(c.get("status") == "docking_failed_not_replaced" for c in cases)
    return f"""# Native parent 검토 패키지 — 최신 12-case 결과

`index.html`을 브라우저에서 직접 열면 검토할 수 있다. 이 README와 `n3-panel/panel.json`이 현재 12-case 상태의 최종 기준이다. 별도 과거 문서는 배경 snapshot일 뿐이며 N3가 진행 중이라고 적힌 문구는 현재 결과를 대체하지 않는다.

## 실제 집계
- requested / completed-retained: **12 / {len(cases)}**; 상태 `completed_with_limits`
- pose-pass case / total pass poses: **{pass_cases} / {pass_poses}**
- E3 assembly: **{branches} branches** (CRBN/VHL 실제 생성분만 집계)
- filtered / docking failed: **{filtered} / {docking_failed}**; 실패와 limits는 교체 없이 보존
- strict-qualified: **0**. 기존 9-parent aggregate는 4 broad family였고, 새 N3 패널은 별도 exploratory 1 family이다. 이를 합산하지 않으며 parent별 six-family gate는 변하지 않았다.
- native calibration: **{_calibration_label(compact)}** 현재 mask가 달라 6HAZ 수치 개선을 주장하지 않는다.
- raw source map 19 → selected-parent map 3 N3, 상태 **UNKNOWN**. source extrapolation이며 measured potency transfer가 아니다.
- scientific review pending; geometry/docking/assembly는 효능·결합·승인·합성경로 또는 사람 승인이 아니다.
- native context는 원격 incomplete residue `A:1385`를 명시적으로 제외했다. 개발자 guard는 **20 Å**이고, 기록된 실제 거리는 별도로 **21.55 Å**이다.

## 1차 근거와 구성
- DOI: `{DOI}`; 구조: `9D12`
- `sources/9D12.cif`, `sources/jm4c01903_si_006.csv`
- `native-baseline/`: 6 baseline runs와 30 poses의 검증된 전체 산출물
- `n3-panel/`: 12 cases, candidate SDF, docking receipts/all poses, 성공한 CRBN/VHL SDF, compact summary
- `source-pair-audit/`: raw map 19 / selected map 3 근거 감사
- `source-export/{PARENT_ID}/`: `verify_source_pair`로 검증한 실제 export 전체; fake official export가 아님
- `manifest.json`: README/index를 포함한 패키지 파일의 bytes/SHA-256. 자체 manifest는 제외

## 재실행 (active repository root에서)
이 패키지는 Python installer가 아니다. core modules, cases, RDKit/Vina와 가상환경 의존성은 live v6 active app/repository에 있다. HTML 검토에는 실행이 필요 없고, 아래 실행에는 active 환경이 필요하다.

```powershell
$pack = '실제 검토 패키지 경로'
.venv/Scripts/python.exe -X utf8 scripts/run_native_parent_probe.py --native-cif "$pack/sources/9D12.cif" --source-export "$pack/source-export/{PARENT_ID}" --output .localdata/native-parent-new-run --exclude-remote-incomplete-residue A:1385
.venv/Scripts/python.exe -X utf8 scripts/run_native_attachment_panel.py --native-probe "$pack/native-baseline" --native-cif "$pack/sources/9D12.cif" --source-export "$pack/source-export/{PARENT_ID}" --output .localdata/native-n3-new-run
```

## 해석 경계
이 자료는 source SAR, parent/receptor calibration, docking diagnostics, panel graph와 dual-E3 graph hypothesis에 대한 사실적 기여를 제공한다. **14 acceptance가 모두 pass했다는 뜻이 아니다.** API 요청·worker 실행·secret·human-approval override는 패키징 과정이나 패키지에 포함되지 않는다.
"""


def build_pack(
    output: Path,
    panel: Path = DEFAULT_PANEL,
    baseline: Path = DEFAULT_BASELINE,
    source_pair_audit: Path = DEFAULT_AUDIT,
    source_export: Path = DEFAULT_EXPORT,
    cif: Path = DEFAULT_CIF,
    csv: Path = DEFAULT_CSV,
    doc: Path | None = DEFAULT_DOC,
    source_pair_verify: Callable[[Path], Any] = _source_pair_verifier,
    native_probe_verify: Callable[[Path], Any] = _native_probe_verifier,
) -> dict[str, Any]:
    """Verify sources and create a new review pack without running workers."""
    output = _safe_resolved_path(Path(output), "Output")
    panel = _safe_resolved_path(Path(panel), "Panel source")
    baseline = _safe_resolved_path(Path(baseline), "Baseline source")
    source_pair_audit = _safe_resolved_path(Path(source_pair_audit), "Audit source")
    source_export = _safe_resolved_path(Path(source_export), "Export source")
    cif = _safe_resolved_path(Path(cif), "CIF source")
    csv = _safe_resolved_path(Path(csv), "CSV source")
    doc = _safe_resolved_path(Path(doc), "Document source") if doc is not None else None
    sources = [panel, baseline, source_pair_audit, source_export, cif, csv]
    if doc is not None:
        sources.append(doc)
    _reject_output_source_overlap(output, sources)
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    panel_info = verify_manifest_tree(panel, {"panel.json", "protocol.json", "compact-summary.json"})
    baseline_info = verify_manifest_tree(baseline, {"native-probe-result.json", "protocol.json"})
    audit_info = verify_manifest_tree(source_pair_audit)
    native_receipt = native_probe_verify(baseline)
    if native_receipt.get("manifest_sha256") != baseline_info["manifest_sha256"]:
        raise ValueError("Native helper and generic manifest verification disagree")
    source_receipt = source_pair_verify(source_export)
    panel_data, _protocol, compact = _panel_data(panel_info["root"])

    output.mkdir(parents=True, exist_ok=False)
    copied: list[dict[str, Any]] = []
    try:
        copied += _copy_verified_tree(baseline_info, output / "native-baseline")
        copied += _copy_verified_tree(panel_info, output / "n3-panel")
        copied += _copy_verified_tree(audit_info, output / "source-pair-audit")
        copied += _copy_unmanifested_verified(source_export, output / "source-export" / PARENT_ID)
        copied.append(_copy_file(cif, output / "sources/9D12.cif"))
        copied.append(_copy_file(csv, output / "sources/jm4c01903_si_006.csv"))
        if doc is not None and doc.is_file():
            copied.append(_copy_file(doc, output / "background/NATIVE_PARENT_SAR_REVIEW_20261001.md"))
        (output / "README.md").write_text(_readme(panel_data, compact), encoding="utf-8", newline="\n")
        (output / "index.html").write_text(_render_html(panel_data, compact), encoding="utf-8", newline="\n")
        rows = []
        for rel, path in sorted(_actual_files(output).items()):
            data = _read(path)
            rows.append({"path": rel, "bytes": len(data), "sha256": _sha(data)})
        manifest = {
            "format_version": "native-parent-review-pack-v1",
            "hash_algorithm": "sha256",
            "files": rows,
            "manifest_self_hash_policy": "manifest.json is excluded",
            "original_manifest_hashes": {
                "native-baseline": baseline_info["manifest_sha256"],
                "n3-panel": panel_info["manifest_sha256"],
                "source-pair-audit": audit_info["manifest_sha256"],
            },
            "verification_receipts": {
                "native_probe": {k: native_receipt.get(k) for k in ("manifest_sha256", "result_sha256", "protocol_sha256")},
                "source_pair_verified": True,
                "source_pair_receipt_sha256": _sha(json.dumps(source_receipt, sort_keys=True, default=str).encode()),
                "verified_copy_count": len(copied),
            },
        }
        (output / MANIFEST).write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        verify_manifest_tree(output)
        return manifest
    except Exception as exc:
        failure = {
            "exception_type": type(exc).__name__,
            "exception_text": str(exc),
        }
        try:
            (output / "packaging-failure.json").write_text(
                json.dumps(failure, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
                newline="\n",
            )
        except Exception:
            pass
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--source-pair-audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--source-export", type=Path, default=DEFAULT_EXPORT)
    parser.add_argument("--cif", type=Path, default=DEFAULT_CIF)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--doc", type=Path, default=DEFAULT_DOC)
    args = parser.parse_args(argv)
    build_pack(args.output, args.panel, args.baseline, args.source_pair_audit, args.source_export, args.cif, args.csv, args.doc)
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
