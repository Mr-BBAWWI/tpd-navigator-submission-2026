#!/usr/bin/env python3
"""입력 산출물만 읽어 휴대 가능한 오프라인 최종 보고서를 만든다."""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

NOW = lambda: datetime.now(timezone.utc).isoformat()
URL_RE = re.compile(r"https?://[^\s<>'\"]+", re.I)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def enc(value: Any, limit: int = 3000) -> str:
    if isinstance(value, str):
        text = value
    else:
        try:
            text = json.dumps(value, ensure_ascii=False, sort_keys=True)
        except TypeError:
            text = str(value)
    text = URL_RE.sub("[외부 URL 생략]", text)
    if len(text) > limit:
        text = text[:limit] + "… (원문은 raw/assessment.json 참조)"
    return html.escape(text, quote=True)


def plain(value: Any, limit: int = 1000) -> str:
    if isinstance(value, str):
        text = value
    else:
        text = json.dumps(value, ensure_ascii=False, sort_keys=True)
    text = URL_RE.sub("[외부 URL 생략]", text).replace("\r", " ").replace("\n", " ")
    return text[:limit] + ("…" if len(text) > limit else "")


def load_json(path: Path) -> tuple[bytes, Any]:
    raw = path.read_bytes()
    return raw, json.loads(raw.decode("utf-8"))


def absolute_unresolved(path: Path) -> Path:
    return Path(os.path.abspath(os.fspath(path)))


def is_link_entry(path: Path) -> bool:
    junction_check = getattr(path, "is_junction", None)
    return path.is_symlink() or bool(junction_check and junction_check())


def reject_symlink_chain(path: Path, allow_missing_leaf: bool = False) -> None:
    p = absolute_unresolved(path)
    parts = p.parts
    cur = Path(parts[0])
    for i, part in enumerate(parts[1:], 1):
        cur /= part
        if not os.path.lexists(cur):
            if allow_missing_leaf:
                return
            raise FileNotFoundError(str(cur))
        if is_link_entry(cur):
            raise ValueError(f"link path is forbidden: {cur.name}")


def validate_tree(root: Path) -> Path:
    root = absolute_unresolved(root)
    reject_symlink_chain(root)
    if not root.is_dir():
        raise NotADirectoryError(str(root))
    for base, dirs, files in os.walk(root, followlinks=False):
        for name in dirs + files:
            p = Path(base) / name
            if is_link_entry(p):
                raise ValueError(f"link in source tree: {name}")
        for name in files:
            if not (Path(base) / name).is_file():
                raise ValueError(f"non-regular source entry: {name}")
    return root


def validate_file(path: Path) -> Path:
    path = absolute_unresolved(path)
    reject_symlink_chain(path)
    if not path.is_file():
        raise FileNotFoundError(str(path))
    return path


def is_within(child: Path, parent: Path) -> bool:
    try:
        child.relative_to(parent)
        return True
    except ValueError:
        return False


def file_record(path: Path, logical: str) -> dict[str, Any]:
    raw = path.read_bytes()
    stat = path.stat()
    return {"logical_path": logical, "sha256": digest(raw), "bytes": len(raw),
            "modified_at": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()}


def tree_records(root: Path, label: str) -> list[dict[str, Any]]:
    out = []
    for p in sorted(x for x in root.rglob("*") if x.is_file()):
        rel = p.relative_to(root).as_posix()
        if ".." in Path(rel).parts or Path(rel).is_absolute():
            raise ValueError("unsafe relative source path")
        out.append(file_record(p, f"{label}/{rel}"))
    return out


def copy_tree_verified(source: Path, target: Path) -> None:
    target.mkdir()
    for p in sorted(source.rglob("*")):
        rel = p.relative_to(source)
        if rel.is_absolute() or ".." in rel.parts or is_link_entry(p):
            raise ValueError("unsafe tree entry")
        dest = target / rel
        if p.is_dir():
            dest.mkdir(exist_ok=True)
        elif p.is_file():
            dest.parent.mkdir(parents=True, exist_ok=True)
            data = p.read_bytes()
            dest.write_bytes(data)
            if digest(dest.read_bytes()) != digest(data):
                raise OSError(f"copy verification failed: {rel}")
        else:
            raise ValueError(f"unsupported tree entry: {rel}")


def find_first(value: Any, keys: set[str]) -> Any:
    if isinstance(value, dict):
        for key, item in value.items():
            if key.lower() in keys and isinstance(item, (str, int, float, bool)):
                return item
        for item in value.values():
            found = find_first(item, keys)
            if found is not None:
                return found
    elif isinstance(value, list):
        for item in value:
            found = find_first(item, keys)
            if found is not None:
                return found
    return None


def number(value: Any) -> int | float | None:
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def test_result(text: str) -> dict[str, Any]:
    matches = list(re.finditer(r"(?m)^Ran\s+(\d+)\s+tests?\b", text))
    if not matches:
        return {"status": "unverified", "tests_run": None, "skipped": None,
                "passed": None, "failed": None,
                "reason": "terminal 'Ran N tests' 기록 없음"}
    match = matches[-1]
    tests_run = int(match.group(1))
    tail = text[match.end():]
    failed_match = re.search(r"(?m)^FAILED(?:\s|\()([^\n]*)", tail)
    ok_match = re.search(r"(?m)^OK(?:\s|$)([^\n]*)", tail)
    terminal = failed_match or ok_match
    suffix = terminal.group(0) if terminal else ""
    skipped_match = re.search(r"skipped=(\d+)", suffix)
    skipped = int(skipped_match.group(1)) if skipped_match else 0
    if failed_match and (not ok_match or failed_match.start() < ok_match.start()):
        failures = sum(int(value) for value in re.findall(r"(?:failures|errors)=(\d+)", suffix))
        return {"status": "failed", "tests_run": tests_run, "skipped": skipped,
                "passed": None, "failed": failures,
                "reason": "terminal FAILED 기록"}
    if ok_match:
        return {"status": "pass", "tests_run": tests_run, "skipped": skipped,
                "passed": tests_run - skipped, "failed": 0,
                "reason": "terminal Ran N tests 및 OK 기록"}
    return {"status": "unverified", "tests_run": tests_run, "skipped": None,
            "passed": None, "failed": None,
            "reason": "terminal OK/FAILED 기록 없음"}


def parent_rows(evidence: Any) -> list[dict[str, Any]]:
    if not isinstance(evidence, dict) or not isinstance(evidence.get("parents"), list):
        return []
    fields = {
        "attempt": "docking_attempts",
        "completed": "docking_completed",
        "failure": "failed_docking",
        "selected": "selected_analogs",
        "assembled": "assemblies",
    }
    rows = []
    for i, row in enumerate(evidence["parents"]):
        if not isinstance(row, dict):
            continue
        ident = next((row.get(k) for k in ("parent_id", "id", "parent_key", "name")
                      if isinstance(row.get(k), (str, int))), f"parent-{i + 1}")
        counts = row.get("counts_recomputed_from_records")
        counts = counts if isinstance(counts, dict) else {}
        out: dict[str, Any] = {"parent": str(ident)}
        for label, source_key in fields.items():
            out[label] = number(counts.get(source_key))
        rows.append(out)
    return rows


def authoritative_human_approval(assessment: dict[str, Any]) -> bool | str:
    decision = assessment.get("decision")
    if not isinstance(decision, dict) or decision.get("formal_confirmation") is not True:
        return "unverified"
    approval = decision.get("human_approval")
    return approval if isinstance(approval, bool) else "unverified"


def explicit_usage_summary(usage: Any) -> dict[str, Any]:
    keys = ("calls", "known_usage_calls", "unknown_usage_calls",
            "input_tokens", "output_tokens", "total_tokens")
    result: dict[str, Any] = {}
    for ledger_name in ("developer_calls", "product_calls"):
        ledger = usage.get(ledger_name) if isinstance(usage, dict) else None
        result[ledger_name] = {
            key: number(ledger.get(key)) if isinstance(ledger, dict) else None
            for key in keys
        }
    result["unknown_usage_failures"] = (
        usage.get("unknown_usage_failures") if isinstance(usage, dict) else None
    )
    return result


def campaign_call_inventory(campaign_root: Path) -> dict[str, Any]:
    source_roles = ("proposer", "adversarial_critic", "experiment_judge")
    successful = 0
    failed = 0
    rounds_root = campaign_root / "rounds"
    for round_dir in sorted(rounds_root.glob("round-*")):
        if not round_dir.is_dir() or is_link_entry(round_dir):
            continue
        for role in source_roles:
            metadata = round_dir / f"{role}.response.metadata.json"
            if metadata.is_file() and not is_link_entry(metadata):
                try:
                    record = json.loads(metadata.read_text(encoding="utf-8"))
                except (OSError, UnicodeDecodeError, json.JSONDecodeError):
                    record = None
                if record:
                    successful += 1
                    continue
            failure_files = [p for p in round_dir.glob(f"{role}*.failure*.json")
                             if p.is_file() and not is_link_entry(p) and p.stat().st_size > 0]
            if failure_files:
                failed += 1
    return {"successful_call_records": successful,
            "failed_call_records": failed,
            "actual_call_count": successful + failed}


def diagnostic_inventory(root: Path, output_name: str) -> dict[str, Any]:
    files = [p.relative_to(root).as_posix() for p in sorted(root.rglob("*")) if p.is_file()]
    indexes = [p for p in files if p.lower().endswith("index.html")]
    manifests = [p for p in files if "manifest" in Path(p).name.lower() and p.lower().endswith(".json")]
    fallbacks = [p for p in files if p.lower().endswith((".json", ".html"))]
    links = indexes or (manifests + [p for p in fallbacks if p not in manifests])[:20]
    return {"name": output_name, "links": links, "file_count": len(files)}


def criteria_from(assessment: Any) -> list[dict[str, Any]]:
    rows = assessment.get("criteria", []) if isinstance(assessment, dict) else []
    categories = {
        "pass": "해당 기준 충족",
        "pending": "근거·검토 대기",
        "failed": "계산·근거 기준 미충족",
        "blocked": "선행조건 차단",
        "unverified": "근거·검토 대기",
    }
    result = []
    for i, row in enumerate(rows):
        if not isinstance(row, dict):
            continue
        status = str(row.get("status", "unverified")).lower()
        if status not in categories:
            status = "unverified"
        required_key = next((key for key in ("required", "requirement", "threshold") if key in row), None)
        observed_key = next((key for key in ("observed", "actual", "value") if key in row), None)
        required = row.get(required_key) if required_key else "입력에 별도 required 필드 없음"
        observed = row.get(observed_key) if observed_key else "입력에 별도 observed 필드 없음"
        base_locator = f"raw/assessment.json#/criteria/{i}"
        result.append({
            "id": str(row.get("id", f"criterion-{i + 1}")),
            "title": str(row.get("title", row.get("id", f"기준 {i + 1}"))),
            "status": status,
            "observed": plain(observed, 3000),
            "observed_locator": f"{base_locator}/{observed_key}" if observed_key else base_locator,
            "required": plain(required, 1500),
            "required_locator": f"{base_locator}/{required_key}" if required_key else base_locator,
            "reason": plain(row.get("reason", "입력에 판정 사유 없음"), 1500),
            "category": categories[status],
        })
    return result


def write_redirect(path: Path, target: str, label: str) -> None:
    text = ("<!doctype html><meta charset='utf-8'><title>historical v5</title>"
            f"<p>historical v5(현재 판정 아님): <a href='{target}'>{html.escape(label)}</a></p>"
            f"<meta http-equiv='refresh' content='0;url={target}'>")
    path.parent.mkdir()
    path.write_text(text, encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ("assessment", "campaign-dir", "dashboard-dir", "parent-pack", "science-pack",
                 "prior-report-dir", "demo-dir", "usage-json", "test-log", "output-dir"):
        ap.add_argument("--" + name, required=True, type=Path)
    args = ap.parse_args()

    files = {"assessment": validate_file(args.assessment), "usage": validate_file(args.usage_json),
             "test_log": validate_file(args.test_log)}
    trees = {"campaign": validate_tree(args.campaign_dir), "dashboard": validate_tree(args.dashboard_dir),
             "parent": validate_tree(args.parent_pack), "science": validate_tree(args.science_pack),
             "prior": validate_tree(args.prior_report_dir), "demo": validate_tree(args.demo_dir)}
    output = absolute_unresolved(args.output_dir)
    reject_symlink_chain(output, allow_missing_leaf=True)
    if os.path.lexists(output):
        raise FileExistsError(str(output))
    for source in list(files.values()) + list(trees.values()):
        if is_within(output, source) or is_within(source, output):
            raise ValueError("output and source paths must not overlap")
    demo_video = trees["demo"] / "TPD_Navigator_시연.mp4"
    if not demo_video.is_file() or demo_video.is_symlink():
        raise FileNotFoundError("TPD_Navigator_시연.mp4")
    campaign_result = trees["campaign"] / "rounds" / "result.json"
    evidence_path = trees["campaign"] / "evidence.json"
    validate_file(campaign_result)
    validate_file(evidence_path)

    assessment_raw, assessment = load_json(files["assessment"])
    usage_raw, usage = load_json(files["usage"])
    campaign_raw, campaign = load_json(campaign_result)
    evidence_raw, evidence = load_json(evidence_path)
    token_report_path = trees["campaign"] / "token-report.json"
    if token_report_path.exists():
        validate_file(token_report_path)
        token_report_raw, token_report = load_json(token_report_path)
    else:
        token_report_raw, token_report = None, None
    log_raw = files["test_log"].read_bytes()
    log_text = log_raw.decode("utf-8", errors="replace")
    if not isinstance(assessment, dict):
        raise ValueError("assessment root must be an object")

    inputs = [file_record(files["assessment"], "assessment.json"),
              file_record(files["usage"], "api-usage.json"),
              file_record(files["test_log"], "tests.log")]
    for label, root in trees.items():
        inputs.extend(tree_records(root, label))
    script_data = Path(__file__).read_bytes()
    source_info = {"name": Path(__file__).name, "sha256": digest(script_data), "bytes": len(script_data)}

    output.mkdir(parents=True)
    mapping = {"dashboard": "research-campaign", "parent": "parent-research",
               "science": "science-research", "prior": "previous-v5", "demo": "demo"}
    for key, name in mapping.items():
        copy_tree_verified(trees[key], output / name)
    raw_dir = output / "raw"
    raw_dir.mkdir()
    raw_copies = {"assessment.json": assessment_raw, "api-usage.json": usage_raw, "tests.log": log_raw,
                  "campaign-result.json": campaign_raw, "evidence.json": evidence_raw}
    if token_report_raw is not None:
        raw_copies["token-report.json"] = token_report_raw
    for name, data in raw_copies.items():
        (raw_dir / name).write_bytes(data)

    write_redirect(output / "expert-packet" / "index.html", "../previous-v5/expert-packet/index.html", "expert packet")
    write_redirect(output / "analog-states" / "index.html", "../previous-v5/analog-states/index.html", "analog states")

    criteria = criteria_from(assessment)
    tests = test_result(log_text)
    all_gate_pass = len(criteria) == 14 and all(x["status"] == "pass" for x in criteria)
    scientific_all_pass = bool(assessment.get("scientific_accepted") is True and all_gate_pass)
    approval = authoritative_human_approval(assessment)
    usage_summary = explicit_usage_summary(usage)
    call_inventory = campaign_call_inventory(trees["campaign"])
    observed_usage = token_report.get("observed_usage") if isinstance(token_report, dict) else None
    unknown_usage_calls = token_report.get("unknown_usage_calls") if isinstance(token_report, dict) else None
    ai_summary = {
        "rounds_completed": number(campaign.get("rounds_completed")) if isinstance(campaign, dict) else None,
        "roles": ["proposer", "critic", "judge"],
        **call_inventory,
        "observed_usage": observed_usage,
        "unknown_usage_calls": unknown_usage_calls,
        "ready_for_experiment": campaign.get("ready_for_experiment") if isinstance(campaign, dict) else None,
    }
    parents = parent_rows(evidence)
    counts = evidence.get("campaign_counts_recomputed_from_records", {}) if isinstance(evidence, dict) else {}
    diagnostics = [diagnostic_inventory(trees["parent"], "parent-research"),
                   diagnostic_inventory(trees["science"], "science-research")]
    provided = {"public_url": "unverified/미제공",
                "youtube": "unverified/미제공",
                "slides_pdf": "unverified/미제공"}
    facts = {
        "format": "offline-final-report-facts/1", "generated_at": NOW(),
        "assessment_id": assessment.get("id"), "assessment_created_at": assessment.get("created_at"),
        "criteria": criteria, "criterion_count": len(criteria),
        "artifact_status": {"software_tests_pass": tests["status"] == "pass",
                            "scientific_all_pass": scientific_all_pass,
                            "human_approval": approval if approval is not None else "unverified"},
        "test_result": tests, "campaign_counts": counts, "parents": parents,
        "ai_campaign": ai_summary, "api_usage": usage_summary,
        "diagnostics": diagnostics, "submission_items_provided": provided,
        "source": source_info, "inputs": inputs,
    }
    (output / "facts.json").write_text(json.dumps(facts, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    gate_rows = "".join(
        f"<tr><td>{enc(x['id'],300)}</td><td>{enc(x['title'],500)}</td><td>{enc(x['status'])}</td>"
        f"<td>{enc(x['observed'],3000)}<br><a href='{html.escape(x['observed_locator'], quote=True)}'>raw locator</a></td>"
        f"<td>{enc(x['required'],1500)}<br><a href='{html.escape(x['required_locator'], quote=True)}'>raw locator</a></td>"
        f"<td>{enc(x['reason'],1500)}</td><td>{enc(x['category'])}</td></tr>" for x in criteria)
    parent_table = "".join("<tr>" + "".join(f"<td>{enc(row.get(k),300)}</td>" for k in
                                           ("parent", "attempt", "completed", "failure", "selected", "assembled")) + "</tr>"
                           for row in parents) or "<tr><td colspan='6'>파싱 가능한 부모별 count 없음</td></tr>"
    diag_html = ""
    for diag in diagnostics:
        links = "".join(f"<li><a href='{diag['name']}/{html.escape(p, quote=True)}'>{html.escape(Path(p).name)}</a></li>"
                        for p in diag["links"])
        diag_html += f"<h3>{enc(diag['name'])}</h3><p>파일 {diag['file_count']}개</p><ul>{links or '<li>JSON/HTML 링크 없음</li>'}</ul>"
    requests = "".join(f"<li><b>{enc(x['title'],500)}</b>: {enc(x['required'],1000)} — {enc(x['reason'],1200)}</li>"
                       for x in criteria if x["status"] != "pass") or "<li>미통과 항목이 0개로 파싱되었으나, 관측된 테스트·검토 없이 문제 없음으로 간주하지 않음</li>"
    status = facts["artifact_status"]
    compact_diag = {"campaign_counts_recomputed_from_records": counts,
                    "parent_count": len(parents), "diagnostic_packs": diagnostics}
    report_html = f"""<!doctype html><html lang='ko'><head><meta charset='utf-8'><title>TPD 최종 오프라인 보고서</title>
<style>body{{font-family:system-ui,sans-serif;max-width:1200px;margin:auto;padding:24px;line-height:1.55}}.table-wrap{{max-width:100%;overflow-x:auto}}table{{border-collapse:collapse;width:100%;table-layout:fixed;font-size:.9rem}}th,td{{border:1px solid #bbb;padding:7px;vertical-align:top;overflow-wrap:anywhere;word-break:break-word}}th{{background:#eef}}code,pre{{white-space:pre-wrap;overflow-wrap:anywhere}}.warn{{background:#fff4d6;padding:12px}}video{{max-width:100%}}</style></head><body>
<h1>최종 오프라인 보고서</h1><p class='warn'>이 보고서는 제공된 입력의 현재 판정을 보존한다. 진단 자료나 사람 승인은 다른 실패·대기 기준을 덮어쓰지 않는다.</p>
<h2>개요와 실행 순서</h2><ol><li>04: 파일 해시 확인</li><li>01: 프로그램 시작</li><li>10: 연구 비교</li><li>02: CPU 실험</li><li>07: 이전 C01/C02 확인</li><li>03: 종료</li></ol>
<p>오프라인 질의 3종:</p><ul><li>CRBN·VHL 후보와 linker 비교</li><li>도킹 실패와 재실험 항목 확인</li><li>14개 기준과 AI개선제안 확인</li></ul><p><a href='research-campaign/index.html'>연구 비교실에서 3개 조회 실행</a></p>
<h2>현재 산출물 상태</h2><ul><li>software_tests_pass: <b>{enc(status['software_tests_pass'])}</b> ({enc(tests)})</li><li>scientific_all_pass: <b>{enc(status['scientific_all_pass'])}</b></li><li>human_approval: <b>{enc(status['human_approval'])}</b></li></ul>
<h2>14 gate 현재 표</h2><p>입력에서 읽힌 기준 수: {len(criteria)}. 정확히 14개가 모두 pass여야 전체 기준 통과로 계산하며 누락 기준을 만들어내지 않는다.</p><div class='table-wrap'><table><thead><tr><th>ID</th><th>기준</th><th>현재 상태</th><th>관측</th><th>요구</th><th>정확한 사유</th><th>분류</th></tr></thead><tbody>{gate_rows}</tbody></table></div>
<h2>현재 변경 및 campaign 집계</h2><p>엄격 기준과 탐색 조립은 별개다. 다음 집계는 evidence의 실제 레코드에서 읽은 값이며 탐색 조립은 strict 통과를 뜻하지 않는다.</p><pre>{enc(counts,3000)}</pre>
<div class='table-wrap'><table><tr><th>부모</th><th>attempt</th><th>completed</th><th>failure</th><th>selected</th><th>assembled</th></tr>{parent_table}</table></div>
<h2>진단 팩</h2>{diag_html}<details><summary>추가 계산·진단 (원 판정과 별도)</summary><pre>{enc(compact_diag,3000)}</pre></details>
<h2>AI campaign 3역할 기록</h2><p>역할은 proposer/critic/judge이며, call 수는 round 디렉터리의 역할별 metadata 또는 failure 기록을 역할당 한 번만 세었다. token-report.json의 직접 관측값만 표시한다.</p><pre>{enc(ai_summary,2000)}</pre><p>AI 역할 기록이나 ready_for_experiment=true는 측정된 과학 기준 통과 또는 실험 실행의 근거가 아니다.</p>
<h2>과학 장벽 및 전문가 요청</h2><ul>{requests}</ul><p>전문가 답변은 이 보고서가 생성하지 않는다. 공식 accept도 다른 failed/blocked/pending 기준을 넓게 면제하지 않는다.</p>
<h2>API 사용</h2><pre>{enc(usage_summary,1200)}</pre><p>호출·토큰 값이 없으면 unknown이다. quota 스냅샷은 비단조적일 수 있어 사용량의 완전한 누적 증거로 간주하지 않는다.</p>
<h2>테스트 로그 판정</h2><p>{enc(tests)}</p><p>terminal 'Ran N tests'와 OK/FAILED가 없으면 unverified이며 테스트 코드 파일 수를 통과 수로 세지 않는다.</p>
<h2>시연 영상</h2><video controls src='demo/TPD_Navigator_시연.mp4'></video><p>실제 앱 캡처를 편집하고 내레이션을 결합한 자료이며 연속 화면 녹화라고 주장하지 않는다. demo 폴더의 SRT·narration·manifest도 원본 그대로 포함한다.</p>
<h2>이전 v5</h2><p><a href='expert-packet/index.html'>historical v5 expert packet</a> · <a href='analog-states/index.html'>historical v5 analog states</a>. 현재 판정이 아니다.</p>
<h2>제출 잔여 항목</h2><p>이 CLI에서는 artifact나 URL을 실제 검증하지 않는다.</p><ul><li>필수 public URL: {enc(provided['public_url'])}</li><li>YouTube 링크: {enc(provided['youtube'])}</li><li>slides PDF: {enc(provided['slides_pdf'])}</li></ul>
<p><a href='research-campaign/index.html'>현재 campaign dashboard</a> · <a href='facts.json'>facts.json</a> · <a href='report.md'>report.md</a></p></body></html>"""
    (output / "report.html").write_text(report_html, encoding="utf-8")

    failed_lines = [f"- {plain(x['title'])}: {x['status']} — {plain(x['reason'])}" for x in criteria if x["status"] != "pass"]
    md = (f"# 최종 오프라인 보고서\n\n생성: {facts['generated_at']}\n\n"
          f"- software_tests_pass: {status['software_tests_pass']}\n- scientific_all_pass: {status['scientific_all_pass']}\n"
          f"- human_approval: {plain(status['human_approval'])}\n- 입력 기준 수: {len(criteria)}\n\n"
          "## 현재 미통과 과학 기준\n" + ("\n".join(failed_lines) or "- 미통과 항목이 0개로 파싱됨; 관측된 검토 없이 문제 없음으로 간주하지 않음") +
          "\n\n## 주의\n탐색 조립과 추가 진단은 strict 통과 또는 실험 실행을 뜻하지 않는다. 이전 v5는 현재 판정이 아니다.\n")
    (output / "report.md").write_text(md, encoding="utf-8")

    outputs = []
    manifest_path = output / "report-manifest.json"
    for p in sorted(x for x in output.rglob("*") if x.is_file() and x != manifest_path):
        rel = p.relative_to(output).as_posix()
        raw = p.read_bytes()
        outputs.append({"path": rel, "sha256": digest(raw), "bytes": len(raw)})
    manifest = {"format": "offline-final-report-manifest/1", "created_at": NOW(),
                "compiler_source": source_info, "inputs": inputs, "outputs": outputs,
                "manifest_excludes_itself": True}
    (output / "report-manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"ok": True, "output_leaf": output.name, "files": len(outputs) + 1}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({"ok": False, "error": type(exc).__name__, "message": str(exc)}, ensure_ascii=False), file=sys.stderr)
        raise SystemExit(1)
