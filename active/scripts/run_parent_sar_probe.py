from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import Chem

from packages.science.parent_sar_probe import (
    FORMAT,
    PARENT_ID,
    build_exact_halogen_probe,
    sha256_file,
    summarise_exports,
    unresolved_potential_stereo,
    validate_predeclared_seeds,
)


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _find_9d12_result(audit: dict) -> tuple[dict, dict]:
    rows = [row for row in audit["parents"] if row["parent_id"] == PARENT_ID]
    if len(rows) != 1:
        raise ValueError("EXACTLY_ONE_9D12_EXPORT_REQUIRED")
    row = rows[0]
    receipt_path = Path(row["receipt_path"])
    if sha256_file(receipt_path) != row["raw_receipt_sha256"]:
        raise ValueError("9D12_JOB_RECEIPT_CHANGED_AFTER_AUDIT")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    result = receipt.get("result")
    if not isinstance(result, dict):
        raise ValueError("9D12_RECEIPT_RESULT_REQUIRED")
    return row, result


def _linkage_signature(linkage: dict) -> tuple:
    matched = linkage.get("matched_record")
    matched_compound = matched.get("compound_id") if isinstance(matched, dict) else None
    pairs = []
    for row in linkage.get("source_sar", []):
        if isinstance(row, dict):
            pairs.append((
                row.get("pair_id"),
                row.get("requested_change_label"),
                row.get("other_compound_id"),
                tuple(row.get("selected_parent_affected_atom_maps", [])),
            ))
    return (
        linkage.get("status"),
        linkage.get("selected_parent_id"),
        matched_compound,
        tuple(sorted(pairs)),
    )


def _run_docking(output: Path, proposal: dict, sites: dict, seeds: list[int], exhaustiveness: int) -> list[dict]:
    from packages.science.design_docking import dock
    from packages.science.dual_e3 import REFERENCE_SOURCE
    from packages.science.reference_parents import load_reference_parent
    from packages.science.structures import atom_sites

    module_root = Path(__file__).resolve().parents[1]
    parent, _ = load_reference_parent(PARENT_ID, REFERENCE_SOURCE)
    protein = atom_sites(module_root / "cases" / "design_sources" / "6HAZ.cif", ["A"])
    parent_ligand = Chem.Mol(parent)
    parent_ligand.RemoveAllConformers()
    candidate = Chem.MolFromSmiles(proposal["mapped_smiles"])
    if candidate is None:
        raise ValueError("PROPOSAL_MOLECULE_UNREADABLE")
    parent_mask = sorted(int(row["atom_map"]) for row in sites.get("atoms", []) if row.get("state") == "PROTECTED")
    candidate_mask = [value for value in parent_mask if value != 8]
    proposal["parent_diagnostic_mask"] = parent_mask
    proposal["candidate_diagnostic_mask"] = candidate_mask
    proposal["candidate_mask_exclusion"] = {"atom_map": 8, "reason": "only actually changed atom"}
    proposal["strict_gate_comparable"] = False

    stereo_reason = unresolved_potential_stereo(parent_ligand)
    rows: list[dict] = []
    for label, molecule, mask in (("parent", parent_ligand, parent_mask), ("observed_br", candidate, candidate_mask)):
        for seed in seeds:
            run_dir = output / "docking" / label / f"seed-{seed}"
            run_dir.mkdir(parents=True, exist_ok=False)
            if stereo_reason:
                result = {
                    "status": "not_run",
                    "reason": "SOURCE_PARENT_UNRESOLVED_POTENTIAL_STEREOCHEMISTRY",
                    "source_parent_unresolved_potential_stereo": stereo_reason,
                    "seed": seed,
                    "ligand": label,
                }
            else:
                try:
                    result = dock(
                        molecule, parent, protein, mask, run_dir,
                        exhaustiveness=exhaustiveness, seed=seed,
                        selected_parent_id=PARENT_ID,
                    )
                except Exception as exc:
                    result = {"status": "failed_exception", "error": str(exc), "seed": seed}
                result["ligand"] = label
                result["seed"] = seed
            result["diagnostic_mask"] = mask
            result["strict_gate_comparable"] = False
            result["creates_acceptance_pass"] = False
            _write_json(run_dir / "result.json", result)
            rows.append({"ligand": label, "seed": seed, "directory": str(run_dir.relative_to(output)), "result": result})
    validate_predeclared_seeds(seeds, rows)
    return rows


def _report(audit: dict, proposal: dict, probe: dict) -> str:
    lines = [
        "# 부모 SAR 소스-결속 프로브",
        "",
        f"- 검증된 부모 export: {audit['parent_count']}개",
        "- pooled family union은 부모별 ≥6개 distinct strict broad-family 잠재 coverage 게이트를 충족시키지 않습니다.",
        "- 과학 승인: 아니오; 설계 정책 변경: 아니오",
        "",
        "## 부모별 감사",
    ]
    for row in audit["parents"]:
        lines.append(
            f"- {row['parent_id']}: cheap-pass {row['cheap_passed_count']}, "
            f"strict candidates {len(row['strict_qualified_candidate_ids'])}, "
            f"strict families {row['strict_qualified_family_count']}, "
            f"families={', '.join(row['strict_qualified_families']) or '-'}, "
            f"docking errors={row['retained_docking_error_count']}, "
            f"attachment exclusions={row['attachment_exclusion_count']}"
        )
    lines += [
        "",
        "## 정확 관측 소스 쌍",
        f"- {proposal['source_pair_id']}: map 9 방향족 탄소에 결합한 terminal map 8 Cl을 Br로 정확히 1회 치환",
        f"- 보호 map 충돌: {proposal['protected_map_conflict']}",
        "- 후보 비교 mask는 실제 변경 map 8만 제외하며 strict gate와 비교 불가입니다.",
        "- CSV 입체화학은 미지정이고 parent 입체화학은 reference/model 할당입니다. 측정 potency를 열거 입체형 또는 docking 후보로 이전하지 않습니다.",
        "",
        "## 도킹",
        f"- 모드: {'실행 요청' if probe['dock_requested'] else '순수 보고서(미실행)'}",
        f"- 보존된 실행 레코드: {len(probe['docking_runs'])}",
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Offline source-bound parent SAR probe")
    parser.add_argument("--exports-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--dock", action="store_true")
    parser.add_argument("--seeds", default="23,41,61")
    parser.add_argument("--exhaustiveness", type=int, default=16)
    args = parser.parse_args(argv)
    if not 1 <= args.exhaustiveness <= 64:
        parser.error("--exhaustiveness must be between 1 and 64")
    try:
        seeds = [int(value.strip()) for value in args.seeds.split(",") if value.strip()]
    except ValueError:
        parser.error("--seeds must be comma-separated integers")
    if not seeds or len(seeds) != len(set(seeds)) or any(seed < 0 or seed > 0xFFFFFFFF for seed in seeds):
        parser.error("--seeds must contain unique nonnegative 32-bit integers")
    if args.output.exists():
        parser.error("output already exists; refusing overwrite")

    audit = summarise_exports(args.exports_root)
    parent_row, result = _find_9d12_result(audit)
    archived_linkage = result.get("selected_parent_measured_evidence")
    source_catalog = result.get("source_catalog")
    sources = source_catalog.get("sources") if isinstance(source_catalog, dict) else None
    medchem = sources.get("medchem_evidence.json") if isinstance(sources, dict) else None
    if not isinstance(archived_linkage, dict) or not isinstance(medchem, dict):
        raise ValueError("9D12 archived linkage and source_catalog medchem evidence required")
    sites = result.get("sites") if isinstance(result.get("sites"), dict) else {"atoms": []}

    from packages.science import medchem_sar
    from packages.science.dual_e3 import REFERENCE_SOURCE
    from packages.science.reference_parents import load_reference_parent
    parent, _ = load_reference_parent(PARENT_ID, REFERENCE_SOURCE)
    recomputed_linkage = medchem_sar.link_selected_parent(parent, PARENT_ID, medchem)
    if not isinstance(recomputed_linkage, dict) or _linkage_signature(recomputed_linkage) != _linkage_signature(archived_linkage):
        raise ValueError("ARCHIVED_SOURCE_LINKAGE_PAIR_OR_ATOM_MAP_MISMATCH")
    linkage = dict(recomputed_linkage)
    linkage["sites"] = archived_linkage.get("sites", sites.get("atoms", []))
    proposal = build_exact_halogen_probe(parent, linkage, medchem)

    args.output.mkdir(parents=True, exist_ok=False)
    _write_json(args.output / "evidence-audit.json", audit)
    _write_json(args.output / "source-pair-proposals.json", {"format": FORMAT, "proposals": [proposal]})
    docking_runs = _run_docking(args.output, proposal, sites, seeds, args.exhaustiveness) if args.dock else []
    probe = {
        "format": FORMAT,
        "receipt": {
            "operation": "source_pair_probe",
            "official_M2_job": False,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "parameters": {"dock": args.dock, "seeds": seeds, "exhaustiveness": args.exhaustiveness},
            "source_parent_job_id": parent_row["job_id"],
            "source_parent_id": PARENT_ID,
        },
        "status": "completed_with_limits",
        "dock_requested": args.dock,
        "docking_runs": docking_runs,
        "strict_acceptance_created": False,
        "qualified_for_counts": False,
        "scientific_approval": False,
        "changes_design_policy": False,
    }
    _write_json(args.output / "probe-result.json", probe)
    (args.output / "report.md").write_text(_report(audit, proposal, probe), encoding="utf-8")

    records = []
    for path in sorted(p for p in args.output.rglob("*") if p.is_file() and p.name != "manifest.json"):
        records.append({
            "path": str(path.relative_to(args.output)).replace("\\", "/"),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    manifest = {
        "format": FORMAT,
        "job_metadata": probe["receipt"],
        "files": records,
        "manifest_excludes_itself": True,
        "scientific_approval": False,
        "changes_design_policy": False,
    }
    _write_json(args.output / "manifest.json", manifest)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"parent-sar-probe refused: {exc}", file=sys.stderr)
        raise
