"""Core-only, version-bound view of B candidate facts and successive diagnostics."""
from __future__ import annotations

import copy
import html
from pathlib import Path
import re

from .evidence_cli import read_export
from .evidence_common import file_index, safe_relative, schema_check, seal, verify_seal, write_new_directory
from .handoff import check, encoded, sha
from .ligand_contacts import read_contacts
from .ligand_preparation import read_preparation
from .quality_evidence import read_quality_evidence, sample_binding, tree_files
from .structure_quality import strict_json

FORMAT = "tpd-b-review-evidence/0.1.0-draft"
BASE = "source/quality"
ROOTS = {"candidate_evidence": BASE + "/evidence", "prior_quality": BASE}


def _key(report):
    return report["binding"]["run_id"], report["binding"]["model_rank"]


def _index(reports, stage):
    result = {}
    for path, report, context in reports:
        key = _key(report)
        check(key not in result, "REVIEW_DUPLICATE_" + stage.upper())
        result[key] = (path, report, context)
    return result


def _stage(item):
    if item is None:
        return {"status": "not_provided", "report_ref": None, "report": None}
    path, report, _ = item
    return {"status": report["status"], "report_ref": {"path": path + "/report.json",
            "sha256": sha(encoded(report)), "digest": report["digest"]}, "report": copy.deepcopy(report)}


def _issues(row):
    issues = []
    def add(code, stage, detail, mode=None, direction=None):
        issues.append({"code": code, "stage": stage, "detail": detail, "mode": mode, "direction": direction})
    quality = row["prior_quality"]
    if quality["structure_status"] != "parsed":
        add("STRUCTURE_" + quality["structure_status"].upper(), "prediction", "예측 구조 판독 미완료; 원 예측 상태와 실패/결측 기록을 확인하세요.")
    for warning in quality["source_warnings"]:
        add("PREDICTION_WARNING", "prediction", warning)
    if quality["quality_report"] is None:
        add("PRIOR_QUALITY_NOT_PROVIDED", "prior_quality", "기존 구조 품질 자료 미제공.")
    elif quality["quality_status"] == "review_required":
        add("PRIOR_QUALITY_REVIEW_REQUIRED", "prior_quality", "기존 방법의 충돌·실패·검사 범위 제한은 해당 원 보고서에 보존됩니다.")
    for analysis in quality["analyses"]:
        if analysis["preparation_error"]:
            add("PRIOR_PREPARATION_ERROR", "prior_quality", analysis["analysis_id"] + ": " + analysis["preparation_error"], analysis["hydrogen_convention"])
        for direction, result in analysis["directions"].items():
            if result["status"] == "failed":
                add("PRIOR_CONTACT_DIRECTION_FAILED", "prior_quality", analysis["analysis_id"] + ": " + result["reason"], analysis["hydrogen_convention"], direction)
    prep = row["ligand_preparation"]
    if prep["status"] == "not_provided":
        add("LIGAND_PREPARATION_NOT_PROVIDED", "ligand_preparation", "리간드 준비 자료 미제공.")
    elif prep["status"] == "not_prepared":
        add("LIGAND_NOT_PREPARED", "ligand_preparation", prep["report"]["error"])
    contact = row["ligand_contacts"]
    if contact["status"] == "not_provided":
        add("LIGAND_CONTACTS_NOT_PROVIDED", "ligand_contacts", "후속 복합체 접촉 검사 자료 미제공; 충돌 0으로 해석할 수 없습니다.")
    else:
        report = contact["report"]
        for limit in report["limitations"]:
            add(limit, "ligand_contacts", limit)
        for model in report["models"]:
            mode = model["hydrogen_convention"]
            add("CONTACT_COVERAGE_INCOMPLETE", "ligand_contacts", "선택한 리간드–단백질 범위의 진단이며 전체 구조 검증이 아닙니다.", mode)
            if model["bad_overlap_pair_count"]:
                add("CONTACT_BAD_OVERLAPS", "ligand_contacts", str(model["bad_overlap_pair_count"]) + " bo pairs; 전체 clashscore나 효능 점수가 아닙니다.", mode)
            # Deduplicate shared atom findings across directions without discarding failures.
            roles, missing = set(), set()
            for direction, result in model["directions"].items():
                if result["execution_status"] != "completed":
                    add("CONTACT_DIRECTION_FAILED", "ligand_contacts", result["error"], mode, direction)
                if result["typing"]:
                    roles.update(r["atom_map"] for r in result["typing"]["role_disagreements"])
                    missing.update(tuple(a) for a in result["typing"]["missing_atom_types"])
            if roles:
                add("CONTACT_ROLE_DISAGREEMENT", "ligand_contacts", "Probe/Lipinski 역할 차이 atom maps: " + ", ".join(map(str, sorted(roles))), mode)
            if missing:
                add("CONTACT_ATOM_TYPES_MISSING", "ligand_contacts", str(len(missing)) + " atoms; 원자형 정보 미완료.", mode)
    add("HUMAN_REVIEW_PENDING", "review", "전문가 검수·후속 판단은 별도입니다.")
    return issues


def project_review(evidence, quality, base_files, preparations, contacts):
    """Project already validated sources; archive readers validate the supplied bytes."""
    check(evidence["digest"] == quality["base_evidence_digest"] and evidence["data_mode"] == quality["data_mode"], "REVIEW_BASE_MISMATCH")
    collection = strict_json(base_files["evidence/collection.json"])
    runs = {r["run_id"]: r for r in collection["runs"]}
    check(len(runs) == len(collection["runs"]), "REVIEW_DUPLICATE_RUN")
    prep_index = _index(preparations, "preparation"); contact_index = _index(contacts, "contacts")
    used_p, used_c, seen = set(), set(), set()
    candidates = []
    qualities = {c["compound_id"]: c for c in quality["candidates"]}
    check(len(qualities) == len(quality["candidates"]) == len(evidence["candidates"]), "REVIEW_CANDIDATE_SET")
    for facts in evidence["candidates"]:
        cid = facts["compound"]["compound_id"]
        check(cid in qualities and qualities[cid]["molecule_id"] == facts["compound"]["molecule_id"], "REVIEW_CANDIDATE_BINDING")
        samples = []
        for old in qualities[cid]["samples"]:
            key = (old["run_id"], old["rank"])
            check(key not in seen, "REVIEW_DUPLICATE_SAMPLE"); seen.add(key)
            run = runs[old["run_id"]]
            model = next(m for m in run["models"] if m["rank"] == old["rank"])
            binding = None
            if model["status"] == "parsed":
                result = strict_json(base_files["evidence/" + model["comparison_file"]["path"]])
                binding = sample_binding(run, model, result)
            prep, contact = prep_index.get(key), contact_index.get(key)
            for item in (prep, contact):
                if item:
                    check(binding is not None, "REVIEW_UNPARSED_SAMPLE")
                    check(item[1]["binding"] == binding, "REVIEW_SAMPLE_BINDING")
                    check(item[1]["data_mode"] == ("synthetic_test" if run["origin"] == "synthetic_test" else "real"), "REVIEW_SAMPLE_ORIGIN")
            if prep:
                used_p.add(key)
            if contact:
                used_c.add(key)
                check(prep is not None and prep[1]["status"] == "prepared_for_review", "REVIEW_CONTACT_REQUIRES_PREPARATION")
                check(old["quality_report"] is not None, "REVIEW_CONTACT_REQUIRES_PRIOR_QUALITY")
                prior = strict_json(base_files[old["quality_report"]["path"]])
                report, context = contact[1], contact[2]
                check(context["ligand_preparation_digest"] == prep[1]["digest"], "REVIEW_STALE_PREPARATION")
                check(report["prior_quality_digest"] == context["prior_quality_digest"] == prior["digest"], "REVIEW_STALE_QUALITY")
                check(report["input_digest"] == context["digest"] and context["binding"] == binding, "REVIEW_CONTACT_CONTEXT")
            row = {"run_id": old["run_id"], "rank": old["rank"], "binding": binding,
                   "prior_quality": copy.deepcopy(old), "ligand_preparation": _stage(prep), "ligand_contacts": _stage(contact),
                   "review_status": "review_required"}
            row["review_issues"] = _issues(row)
            samples.append(row)
        candidates.append({"compound_id": cid, "molecule_id": facts["compound"]["molecule_id"],
                           "evidence": copy.deepcopy(facts), "samples": samples,
                           "synthesis_assessment": {"status": "not_assessed", "reason": "NO_DEDICATED_SYNTHESIS_ASSESSMENT_IN_SOURCE_CONTRACT"},
                           "review_status": "review_required"})
    check(used_p == set(prep_index) and used_c == set(contact_index), "REVIEW_ORPHAN_REPORT")
    rows = [s for c in candidates for s in c["samples"]]
    summary = {"candidate_count": len(candidates), "sample_count": len(rows),
               "structures": {k: sum(s["prior_quality"]["structure_status"] == k for s in rows) for k in ("parsed", "missing", "invalid")},
               "preparations": {k: sum(s["ligand_preparation"]["status"] == k for s in rows) for k in ("not_provided", "prepared_for_review", "not_prepared")},
               "contacts": {k: sum(s["ligand_contacts"]["status"] == k for s in rows) for k in ("not_provided", "completed_with_limits", "failed_or_partial")}}
    bundle = seal({"format": FORMAT, "data_mode": evidence["data_mode"], "case_id": evidence["case_id"],
                   "base_quality_digest": quality["digest"], "base_manifest_sha256": sha(base_files["manifest.json"]),
                   "source_roots": copy.deepcopy(ROOTS), "starting_ligand": copy.deepcopy(evidence["starting_ligand"]),
                   "shared_hypothesis": copy.deepcopy(evidence["shared_hypothesis"]),
                   "literature_reconciliation": copy.deepcopy(evidence["literature_reconciliation"]),
                   "candidates": candidates, "summary": summary, "authority": copy.deepcopy(evidence["authority"]),
                   "meaning": "Stage-scoped source records and review issues; no replacement of observations, automatic acceptance, ranking or synthesis/efficacy claim."})
    schema_check(bundle, "b_review_evidence.schema.json")
    return bundle


def markdown(bundle):
    def text(value):
        return html.escape(str(value), quote=False).replace("|", "\\|").replace("\n", " ").replace("\r", " ")
    labels = {"not_provided": "자료 미제공", "prepared_for_review": "준비됨·검토 전", "not_prepared": "준비 실패/거부",
              "completed_with_limits": "계산 완료·제한 있음", "failed_or_partial": "계산 실패/부분 결과"}
    lines = ["# B 후보별 최종 검토 근거 묶음", "", "**내부 검토 자료 · 전문가 승인 전 · 공개 준비 완료 아님**", "",
             "문헌 관측·기존 계산 물성·모델 예측·검사 결과를 함께 읽습니다. 원자형 역할 차이와 수소 방향 미최적화를 유지하며 자동 합격·효능 순위는 제공하지 않습니다.", "",
             f'[문헌·분자·예측 원 근거표]({BASE}/evidence/README.md) · [이전 구조 품질 근거]({BASE}/README.md)', "",
             "원 단계의 미검증 표시는 그 당시 기록입니다. 이후 접촉 계산이 있어도 준비 보고서나 과거 충돌 결과를 덮어쓰지 않습니다.", "",
             "| 후보 | 실행 / rank | 구조 판독 | 리간드 준비 | 후속 접촉 검사 |", "|---|---|---|---|---|"]
    for candidate in bundle["candidates"]:
        if not candidate["samples"]:
            lines.append(f'| {text(candidate["compound_id"])} | 제공된 예측 없음 | 미실행 | 자료 미제공 | 자료 미제공 |')
        for row in candidate["samples"]:
            lines.append(f'| {text(candidate["compound_id"])} | {text(row["run_id"])} / {row["rank"]} | {text(row["prior_quality"]["structure_status"])} | {labels[row["ligand_preparation"]["status"]]} | {labels[row["ligand_contacts"]["status"]]} |')
    for candidate in bundle["candidates"]:
        facts = candidate["evidence"]
        lines += ["", f'## {text(candidate["compound_id"])}', "", "계산 물성: " + text(facts["computed_properties"]["values"]),
                  "", "문헌 합성 근거: 별도 검토 자료 없음·미평가. 조립 성공을 합성 가능성으로 해석하지 않습니다.", ""]
        for item in facts["review_items"]:
            lines.append("- " + text(item))
        for row in candidate["samples"]:
            lines += ["", f'### {text(row["run_id"])} / {row["rank"]}', ""]
            for field, title in (("ligand_preparation", "리간드 준비 원 보고서"), ("ligand_contacts", "후속 접촉 원 보고서")):
                stage = row[field]
                if stage["report_ref"]:
                    lines.append(f'[{title}]({stage["report_ref"]["path"]}) · {labels[stage["status"]]}')
            contact = row["ligand_contacts"]["report"]
            if contact:
                lines += ["", "| H 거리 기준 | 이전 bo 쌍 | 이번 bo 쌍 | 양방향 완료 | 검사 범위 |", "|---|---|---|---|---|"]
                for model in contact["models"]:
                    mode = model["hydrogen_convention"]
                    # Prior analysis IDs are producer-defined. Use no guessed pairing;
                    # the original full analyses remain available in JSON and linked report.
                    count = "판정 불가" if model["bad_overlap_pair_count"] is None else str(model["bad_overlap_pair_count"])
                    lines.append(f'| {text(mode)} | 원 보고서 참조 | {count} | {model["both_directions_completed"]} | 불완전 |')
                lines += ["", "전체 리간드 H 방향과 거리 규약이 달라질 수 있어, 이전 대비 수치 차이를 구조 개선이나 누락 H만의 효과로 단정하지 않습니다."]
            lines.append("")
            for issue in row["review_issues"]:
                scope = "/".join(str(issue[k]) for k in ("stage", "mode", "direction") if issue[k])
                lines.append(f'- [{text(scope)}] {text(issue["code"])}: {text(issue["detail"])}')
    lines += ["", "이 reader는 파일 결합·버전·투영을 검사합니다. 과학 계산 재실행·실행자 인증·사람 승인·공식 M2 상태 갱신은 수행하지 않습니다.", ""]
    return "\n".join(lines).encode()


def _load(base, preparation_paths, contact_paths, *, allow_synthetic):
    quality = read_quality_evidence(base / BASE, allow_synthetic=allow_synthetic)
    evidence = read_export(base / BASE / "evidence", allow_synthetic=allow_synthetic)
    base_files = tree_files(base / BASE)
    preparations = [(p, read_preparation(safe_relative(base, p), allow_synthetic=allow_synthetic), None) for p in preparation_paths]
    contacts = []
    for path in contact_paths:
        directory = safe_relative(base, path)
        report = read_contacts(directory, allow_synthetic=allow_synthetic)
        # read_contacts verified this context and all its source bytes.
        context = strict_json((directory / "input/context.json").read_bytes())
        contacts.append((path, report, context))
    return project_review(evidence, quality, base_files, preparations, contacts)


def export_review_evidence(quality_directory, preparation_directories, contact_directories, destination, *, allow_synthetic=False):
    """Write a new archive only. No science runtime, network or subprocess is used."""
    destination = Path(destination)
    check(not destination.exists(), "OUTPUT_EXISTS")
    quality = read_quality_evidence(quality_directory, allow_synthetic=allow_synthetic)
    evidence = read_export(Path(quality_directory) / "evidence", allow_synthetic=allow_synthetic)
    base_files = tree_files(quality_directory)
    files = {BASE + "/" + n: d for n, d in base_files.items()}
    preparations, contacts = [], []
    for stage, directories, reader, entries in (("preparation", preparation_directories, read_preparation, preparations),
                                                ("contacts", contact_directories, read_contacts, contacts)):
        loaded = [(Path(d), reader(d, allow_synthetic=allow_synthetic)) for d in directories]
        loaded.sort(key=lambda pair: (*_key(pair[1]), pair[1]["digest"]))
        for i, (directory, report) in enumerate(loaded):
            path = f"{stage}/{i:04d}"
            context = strict_json((directory / "input/context.json").read_bytes()) if stage == "contacts" else None
            entries.append((path, report, context))
            files.update({path + "/" + n: d for n, d in tree_files(directory).items()})
    bundle = project_review(evidence, quality, base_files, preparations, contacts)
    files.update({"review-evidence.json": encoded(bundle), "README.md": markdown(bundle)})
    files["manifest.json"] = encoded(seal({"format": FORMAT, "preparation_directories": [p for p, _, _ in preparations],
                                          "contact_directories": [p for p, _, _ in contacts], "files": file_index(files), "bundle_digest": bundle["digest"]}))
    write_new_directory(destination, files)
    return bundle


def read_review_evidence(directory, *, allow_synthetic=False):
    directory = Path(directory)
    manifest = strict_json((directory / "manifest.json").read_bytes()); verify_seal(manifest)
    check(set(manifest) == {"format", "preparation_directories", "contact_directories", "files", "bundle_digest", "digest"}
          and manifest["format"] == FORMAT, "REVIEW_ARCHIVE_FORMAT")
    files = tree_files(directory)
    check(set(files) == set(manifest["files"]) | {"manifest.json"}, "REVIEW_UNINDEXED_FILES")
    for name, ref in manifest["files"].items():
        check(name != "manifest.json" and sha(files[name]) == ref["sha256"] and len(files[name]) == ref["bytes"], "REVIEW_FILE_HASH")
    for name, stage in (("preparation_directories", "preparation"), ("contact_directories", "contacts")):
        paths = manifest[name]
        check(isinstance(paths, list) and len(paths) == len(set(paths))
              and all(isinstance(p, str) and re.fullmatch(stage + r"/[0-9]{4,}", p) for p in paths), "REVIEW_SOURCE_DIRECTORIES")
    prefixes = [BASE + "/"] + [p + "/" for p in manifest["preparation_directories"] + manifest["contact_directories"]]
    check(all(n in ("manifest.json", "review-evidence.json", "README.md") or any(n.startswith(p) for p in prefixes) for n in files), "REVIEW_UNEXPECTED_FILES")
    expected = _load(directory, manifest["preparation_directories"], manifest["contact_directories"], allow_synthetic=allow_synthetic)
    check(encoded(expected) == files["review-evidence.json"] and expected["digest"] == manifest["bundle_digest"], "REVIEW_PROJECTION")
    check(markdown(expected) == files["README.md"], "REVIEW_MARKDOWN_PROJECTION")
    return expected
