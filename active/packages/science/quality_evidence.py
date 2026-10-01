"""Join quality diagnostics to immutable candidate evidence by exact sample hashes."""
from __future__ import annotations

import html
from pathlib import Path

from .evidence_cli import read_export
from .evidence_common import file_index, safe_relative, schema_check, seal, verify_seal, write_new_directory
from .handoff import check, encoded, parse, sha
from .structure_quality import read_quality, strict_json

FORMAT = "tpd-quality-evidence/0.1.0-draft"


def tree_files(directory):
    directory = Path(directory)
    files = {}
    for path in sorted(directory.rglob("*")):
        check(not path.is_symlink(), "QUALITY_EXPORT_SYMLINK")
        if path.is_file():
            name = path.relative_to(directory).as_posix()
            files[name] = safe_relative(directory, name).read_bytes()
    return files


def sample_binding(run, model, result):
    return {"compound_id": run["compound_id"], "molecule_id": run["molecule_id"],
            "run_id": run["run_id"], "model_rank": model["rank"], "input_sha256": run["input_sha256"],
            "structure_sha256": result["source_manifest"]["files"]["structure"]["sha256"]}


def project(evidence, base_files, reports):
    collection = parse(base_files["collection.json"])
    by_sample = {}
    for path, report in reports:
        key = (report["binding"]["run_id"], report["binding"]["model_rank"])
        check(key not in by_sample, "QUALITY_DUPLICATE_SAMPLE_REPORT")
        by_sample[key] = (path, report)
    used, candidates = set(), []
    for candidate in evidence["candidates"]:
        cid = candidate["compound"]["compound_id"]
        samples = []
        for run in collection["runs"]:
            if run["compound_id"] != cid:
                continue
            for model in run["models"]:
                key = (run["run_id"], model["rank"])
                row = {"run_id": run["run_id"], "rank": model["rank"], "structure_status": model["status"],
                       "source_warnings": [] if model["summary"] is None else model["summary"]["warnings"],
                       "input_atom_identity": "checked_by_existing_reader" if model["status"] == "parsed" else "not_assessed",
                       "ligand_stereochemistry": "not_assessed", "quality_status": "not_assessed",
                       "quality_report": None, "quality_summary": None, "analyses": []}
                if model["status"] == "parsed":
                    row["ligand_stereochemistry"] = "requires_review" if "STEREOCHEMISTRY_REQUIRES_REVIEW" in model["summary"]["warnings"] else "checked_by_existing_reader"
                if key in by_sample:
                    check(model["status"] == "parsed", "QUALITY_CANNOT_ATTACH_TO_UNPARSED_MODEL")
                    path, report = by_sample[key]
                    result = parse(base_files[model["comparison_file"]["path"]])
                    check(report["binding"] == sample_binding(run, model, result), "QUALITY_SAMPLE_BINDING_MISMATCH")
                    expected_mode = "synthetic_test" if run["origin"] == "synthetic_test" else "real"
                    check(report["data_mode"] == expected_mode, "QUALITY_SAMPLE_ORIGIN_MISMATCH")
                    used.add(key)
                    row.update(quality_status=report["summary"]["interpretation_status"],
                               quality_report={"path": path + "/report.json", "sha256": sha(encoded(report))},
                               quality_summary=report["summary"], analyses=report["analyses"])
                samples.append(row)
        candidates.append({"compound_id": cid, "molecule_id": candidate["compound"]["molecule_id"],
                           "prediction_status": candidate["prediction"]["status"], "samples": samples,
                           "no_prediction_reason": "NO_PREDICTION_RECORD_SUPPLIED" if not samples else None})
    check(used == set(by_sample), "QUALITY_ORPHAN_SAMPLE_REPORT")
    bundle = seal({"format": FORMAT, "data_mode": evidence["data_mode"], "case_id": evidence["case_id"],
                 "base_evidence_digest": evidence["digest"], "base_manifest_sha256": sha(base_files["manifest.json"]),
                 "candidates": candidates,
                 "authority": {"human_review": "pending", "approval_record_created": False, "dispatch_authorized": False,
                               "public_release_ready": False, "efficacy_established": False},
                 "meaning": "Quality observations and missing coverage; no pooled score, automatic repair, candidate selection or scientific acceptance."})
    schema_check(bundle, "b_quality_evidence.schema.json")
    return bundle


def markdown(bundle):
    def text(value):
        return html.escape(str(value), quote=False).replace("|", "\\|").replace("\n", " ").replace("\r", " ")
    states = {"not_assessed": "미검증", "review_required": "검토 필요", "no_findings_in_declared_scope": "선택한 검사 범위 내 발견 없음"}
    run_states = {"completed": "실행 완료", "partial": "일부 실행", "not_run": "미실행", "failed": "검증 실패"}
    warnings = {"HYDROGENS_UNPLACED": "일부 수소 추가 미완료", "ATOM_TYPING_UNVERIFIED": "원자형 적합성 미확인",
                "NONSTANDARD_RESIDUE_TYPING_UNVERIFIED": "비표준 잔기의 원자형 미확인",
                "HYDROGEN_PREPARATION_UNVERIFIED": "수소 준비 결과 확인 불가",
                "BIDIRECTIONAL_COVERAGE_INCOMPLETE": "양방향 접촉 검사 미완료"}
    hydrogen_names = {"xray": "X선 기준", "nuclear": "핵 위치 기준"}
    lines = ["# 구조 검증을 연결한 후보 근거표", "", "**내부 검토 자료 · 전문가 승인 전**", "",
             "구조 파일 판독, 입체화학 검사, 접촉 검증의 실행 상태와 한계를 구분합니다. 충돌 개수는 효능 점수나 전체 MolProbity clashscore가 아닙니다.", "",
             "[기존 후보별 문헌·분자·예측 근거표](evidence/README.md)는 원본 그대로 보존했습니다.", "",
             f'자료 모드: `{bundle["data_mode"]}`. 연결 버전: `{bundle["digest"]}`.', "",
             "| 후보 | 실행 / rank | 구조 판독 | 원자·입체화학 검사 | 접촉 검증 |", "|---|---|---|---|---|"]
    for candidate in bundle["candidates"]:
        if not candidate["samples"]:
            lines.append(f'| {text(candidate["compound_id"])} | 제공된 예측 없음 | 미실행 | 미검증 | 미검증 |')
        for sample in candidate["samples"]:
            identity = "기존 판독기 검사 기록 있음" if sample["input_atom_identity"] != "not_assessed" else "미검증"
            if sample["ligand_stereochemistry"] == "requires_review":
                identity += " / 입체화학 검토 필요"
            lines.append(f'| {text(candidate["compound_id"])} | {text(sample["run_id"])} / {sample["rank"]} | {text(sample["structure_status"])} | {identity} | {states[sample["quality_status"]]} |')
    lines += ["", "검증 자료가 없거나 도구가 실패한 결과를 ‘충돌 없음’으로 표시하지 않습니다. 검사가 끝나도 전체 구조의 과학적 타당성이 승인된 것은 아닙니다.", ""]
    for candidate in bundle["candidates"]:
        for sample in candidate["samples"]:
            lines += [f'## {text(candidate["compound_id"])} · {text(sample["run_id"])} / {sample["rank"]}', ""]
            if sample["source_warnings"]:
                lines.append("기존 판독 경고: " + text(", ".join(sample["source_warnings"])))
            if sample["quality_report"] is None:
                lines += ["접촉 검증 자료 미제공. 후속 기하 해석에 필요한 검증이 남아 있습니다.", ""]
                continue
            lines += [f'[검증 원자료와 상세 판독]({sample["quality_report"]["path"]})', "",
                      "| 검사 범위 | 수소 거리 설정 | 수소 준비 / 접촉 검사 | 충돌 신호 원자쌍 | 검사 범위 상태 |", "|---|---|---|---|---|"]
            for a in sample["analyses"]:
                count = str(a["bad_overlap_pair_count"]) if a["finding_status"] != "not_assessed" else "판정 불가"
                coverage = "불완전" if a["coverage"] == "incomplete" else "명시한 범위만 검사"
                lines.append(f'| {text(a["scope"]["description"])} | {hydrogen_names[a["hydrogen_convention"]]} | {run_states[a["preparation_status"]]} / {run_states[a["execution_status"]]} | {count} | {coverage} |')
            lines.append("")
            for a in sample["analyses"]:
                if a["limitations"]:
                    limitations = ", ".join(warnings.get(code, code) for code in a["limitations"])
                    lines.append(f'- {text(a["scope"]["description"])} ({hydrogen_names[a["hydrogen_convention"]]}) 제한: {text(limitations)}')
                if a["unplaced_hydrogens"]:
                    lines.append("  수소 추가 미완료: " + text(", ".join(a["unplaced_hydrogens"])))
                if a["preparation_error"]:
                    lines.append("  수소 준비 오류: " + text(a["preparation_error"]))
                for direction, status in a["directions"].items():
                    if status["status"] == "failed":
                        lines.append(f'  {text(direction)} 검증 실패: {text(status["reason"])}')
            lines.append("")
    lines += ["원 구조와 진단 사본을 보존합니다. 이 묶음은 모델 가중치·좌표·공식 승인 상태를 변경하지 않습니다.", ""]
    return "\n".join(lines).encode("utf-8")


def export_quality_evidence(evidence_directory, quality_directories, destination, *, allow_synthetic=False):
    check(not Path(destination).exists(), "OUTPUT_EXISTS")
    base = Path(evidence_directory)
    evidence = read_export(base, allow_synthetic=allow_synthetic)
    base_files = tree_files(base)
    files = {"evidence/" + name: data for name, data in base_files.items()}
    reports, paths = [], []
    for index, directory in enumerate(quality_directories):
        report = read_quality(directory, allow_synthetic=allow_synthetic)
        path = f"quality/q{index:04d}"
        paths.append(path)
        reports.append((path, report))
        files.update({path + "/" + name: data for name, data in tree_files(directory).items()})
    bundle = project(evidence, base_files, reports)
    files.update({"quality-evidence.json": encoded(bundle), "README.md": markdown(bundle)})
    files["manifest.json"] = encoded(seal({"format": FORMAT, "quality_directories": paths,
                                          "files": file_index(files), "bundle_digest": bundle["digest"]}))
    write_new_directory(destination, files)
    return bundle


def read_quality_evidence(directory, *, allow_synthetic=False):
    directory = Path(directory)
    manifest = strict_json((directory / "manifest.json").read_bytes())
    verify_seal(manifest)
    check(manifest["format"] == FORMAT, "QUALITY_EXPORT_FORMAT")
    files = tree_files(directory)
    check(set(files) == set(manifest["files"]) | {"manifest.json"}, "QUALITY_EXPORT_UNINDEXED_FILES")
    for name, entry in manifest["files"].items():
        check(sha(files[name]) == entry["sha256"] and len(files[name]) == entry["bytes"], "QUALITY_EXPORT_FILE_HASH")
    paths = manifest["quality_directories"]
    check(len(paths) == len(set(paths)), "QUALITY_DUPLICATE_DIRECTORY")
    reports = [(p, read_quality(safe_relative(directory, p), allow_synthetic=allow_synthetic)) for p in paths]
    evidence = read_export(directory / "evidence", allow_synthetic=allow_synthetic)
    base_files = {n[9:]: d for n, d in files.items() if n.startswith("evidence/")}
    expected = project(evidence, base_files, reports)
    check(encoded(expected) == files["quality-evidence.json"] and expected["digest"] == manifest["bundle_digest"], "QUALITY_EXPORT_PROJECTION")
    check(markdown(expected) == files["README.md"], "QUALITY_EXPORT_REPORT_PROJECTION")
    return expected
