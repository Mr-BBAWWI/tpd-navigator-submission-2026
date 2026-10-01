#!/usr/bin/env python3
import argparse
import hashlib
import html
import json
import math
import os
import shutil
from pathlib import Path, PurePosixPath

IDS = [
    "parent_funnel", "expert_parent_selection", "modifiable_sites",
    "distinct_constitutional_graphs", "actual_broad_families", "qualified_panel",
    "both_e3_assembly", "core_interaction_preservation", "microstates_h_direction",
    "known_crbn_calibration", "novel_ternary_repeats", "novel_ternary_geometry",
    "exact_synthesis_review", "formal_expert_decision",
]
STATUSES = {"pass", "failed", "pending"}
SOURCE_PATH = "analog-states/sources/reply.docx"
RAW_FIXED = {
    "raw/assessment.json", "raw/policy.json", "raw/test-summary.json",
    "raw/api-usage.json", "raw/runtime-strict.json", "raw/tests.log",
}
MANIFEST_NAMES = {"manifest.json", "artifact-manifest.json"}
def die(message):
    raise ValueError(message)


def load_json(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        die(f"invalid JSON {path}: {exc}")


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(131072), b""):
            digest.update(block)
    return digest.hexdigest()


def is_link_or_junction(path):
    try:
        if path.is_symlink():
            return True
        checker = getattr(path, "is_junction", None)
        return bool(checker and checker())
    except OSError as exc:
        die(f"cannot inspect path {path}: {exc}")


def safe_rel(value):
    if not isinstance(value, str) or not value:
        die("path must be a non-empty string")
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or ".." in path.parts
        or "\\" in value
        or ":" in value
        or any(part in ("", ".") for part in path.parts)
    ):
        die(f"unsafe path: {value}")
    return path.as_posix()


def finite(value, where="facts"):
    if isinstance(value, float) and not math.isfinite(value):
        die(f"non-finite number at {where}")
    if isinstance(value, dict):
        for key, item in value.items():
            finite(item, f"{where}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            finite(item, f"{where}[{index}]")


def original_notes(assessment):
    notes = []
    seen = set()
    for field, text_key in (
        ("expert_questions", "question"),
        ("pending_actions", "action"),
    ):
        entries = assessment.get(field)
        if not isinstance(entries, list):
            die(f"assessment {field} must be a list")
        for entry in entries:
            if not isinstance(entry, dict):
                die(f"assessment {field} entries must be objects")
            text = entry.get(text_key)
            if not isinstance(text, str) or not text:
                die(f"assessment {field}.{text_key} must be non-empty text")
            if text not in seen:
                seen.add(text)
                notes.append(text)
    return notes


def validate_source_record(source):
    if not isinstance(source, dict):
        die("missing source")
    if source.get("path") != SOURCE_PATH:
        die(f"source.path must be {SOURCE_PATH}")
    if not isinstance(source.get("bytes"), int) or isinstance(
        source.get("bytes"), bool
    ) or source["bytes"] < 0:
        die("source.bytes is invalid")
    digest = source.get("sha256")
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(char not in "0123456789abcdef" for char in digest)
    ):
        die("source.sha256 is invalid")


def validate_facts(facts, evidence):
    finite(facts)
    if facts.get("format") != "tpd-final-delivery-facts/1":
        die("unsupported facts format")

    criteria = facts.get("criteria")
    if not isinstance(criteria, list) or len(criteria) != 14:
        die("criteria must contain exactly 14 entries")
    if [item.get("id") for item in criteria if isinstance(item, dict)] != IDS:
        die("criteria IDs or ordering are invalid")
    for item in criteria:
        if not isinstance(item, dict) or {
            "id", "title", "status", "reason", "observed_text", "required_text",
        } - set(item):
            die(f"incomplete criterion: {item.get('id') if isinstance(item, dict) else None}")
        if item["status"] not in STATUSES:
            die(f"invalid criterion status: {item['status']}")
        for key in ("title", "reason", "observed_text", "required_text"):
            if not isinstance(item[key], str):
                die(f"{item['id']}.{key} must be text")

    if facts.get("scientific_accepted") is not False:
        die("scientific_accepted must be false")
    if all(item["status"] == "pass" for item in criteria):
        die("all-pass criteria conflict with scientific_accepted=false")
    if facts.get("delivery_status") != "DEMO/INCOMPLETE":
        die("delivery_status must be DEMO/INCOMPLETE")

    tests = facts.get("tests")
    required_tests = {
        "command", "tests_run", "failures", "errors", "skipped",
        "exit_code", "log_sha256",
    }
    if not isinstance(tests, dict) or set(tests) != required_tests:
        die("invalid tests object")
    if not isinstance(tests["command"], str):
        die("tests.command must be text")
    for key in ("tests_run", "failures", "errors", "skipped", "exit_code"):
        if not isinstance(tests[key], int) or isinstance(tests[key], bool):
            die(f"tests.{key} must be an integer")
    if min(
        tests["tests_run"], tests["failures"], tests["errors"], tests["skipped"]
    ) < 0:
        die("test counts cannot be negative")
    if tests["failures"] + tests["errors"] + tests["skipped"] > tests["tests_run"]:
        die("test result counts exceed tests_run")
    if (
        not isinstance(tests["log_sha256"], str)
        or len(tests["log_sha256"]) != 64
        or any(char not in "0123456789abcdef" for char in tests["log_sha256"])
    ):
        die("tests.log_sha256 is invalid")

    tests_log = evidence / "raw/tests.log"
    if not tests_log.is_file() or sha(tests_log) != tests["log_sha256"]:
        die("tests.log checksum mismatch")
    if load_json(evidence / "raw/test-summary.json") != tests:
        die("test-summary does not exactly match facts.tests")

    features = facts.get("software_features")
    if not isinstance(features, list) or not all(isinstance(x, str) for x in features):
        die("software_features must be strings")

    measurements = facts.get("measurements")
    if not isinstance(measurements, list):
        die("measurements must be a list")
    for table in measurements:
        if (
            not isinstance(table, dict)
            or not isinstance(table.get("title"), str)
            or not isinstance(table.get("headers"), list)
            or not table["headers"]
            or not isinstance(table.get("rows"), list)
        ):
            die("invalid measurement table")
        if not all(isinstance(x, str) for x in table["headers"]):
            die("measurement headers must be strings")
        for row in table["rows"]:
            if not isinstance(row, list) or len(row) != len(table["headers"]):
                die("measurement row width mismatch")
            if any(isinstance(x, (dict, list)) for x in row):
                die("measurement cells must be scalar")

    api = facts.get("api_usage")
    required_api = {
        "development_known_tokens", "product_known_tokens",
        "total_known_tokens", "unknown_usage_calls",
        "latest_quota_estimate", "quota_note",
    }
    if not isinstance(api, dict) or set(api) != required_api:
        die("invalid api_usage object")
    for key in (
        "development_known_tokens", "product_known_tokens",
        "total_known_tokens", "unknown_usage_calls",
    ):
        if not isinstance(api[key], int) or isinstance(api[key], bool) or api[key] < 0:
            die(f"api_usage.{key} must be a non-negative integer")
    if api["total_known_tokens"] != (
        api["development_known_tokens"] + api["product_known_tokens"]
    ):
        die("api_usage.total_known_tokens must equal the known-token sum")
    if (
        not isinstance(api["latest_quota_estimate"], (int, float))
        or isinstance(api["latest_quota_estimate"], bool)
        or api["latest_quota_estimate"] < 0
    ):
        die("api_usage.latest_quota_estimate must be non-negative")
    if not isinstance(api["quota_note"], str):
        die("api_usage.quota_note must be text")

    freshness = facts.get("freshness")
    if not isinstance(freshness, dict):
        die("missing freshness")
    for key in (
        "current_policy", "current_implementation", "source_runtime_current",
    ):
        if not isinstance(freshness.get(key), bool):
            die(f"freshness.{key} must be boolean")

    source = facts.get("source")
    validate_source_record(source)
    source_path = evidence / safe_rel(source["path"])
    if not source_path.is_file():
        die("source file is missing")
    if source_path.stat().st_size != source["bytes"] or sha(source_path) != source["sha256"]:
        die("source factual-input binding mismatch")

    docfacts = facts.get("docfacts")
    if docfacts is not None and (
        not isinstance(docfacts, dict) or docfacts.get("source") != source
    ):
        die("docfacts.source does not bind to source")

    assessment = load_json(evidence / "raw/assessment.json")
    assessment_format = assessment.get("format") if isinstance(assessment, dict) else None
    if (
        not isinstance(assessment_format, str)
        or not assessment_format.startswith("scientific-acceptance/")
    ):
        die("invalid assessment format")
    if assessment.get("scientific_accepted") is not facts["scientific_accepted"]:
        die("assessment scientific_accepted does not exactly match facts")

    assessed_criteria = assessment.get("criteria")
    if not isinstance(assessed_criteria, list):
        die("assessment criteria must be a list")
    if not all(isinstance(item, dict) for item in assessed_criteria):
        die("assessment criteria entries must be objects")
    assessed = [
        {"id": item.get("id"), "status": item.get("status")}
        for item in assessed_criteria
    ]
    expected = [{"id": item["id"], "status": item["status"]} for item in criteria]
    if assessed != expected:
        die("assessment criteria do not exactly match facts")

    followup = assessment.get("expert_followup")
    assessment_source = followup.get("source") if isinstance(followup, dict) else None
    if not isinstance(assessment_source, dict):
        die("assessment expert_followup.source is missing")
    if (
        assessment_source.get("bytes") != source["bytes"]
        or assessment_source.get("sha256") != source["sha256"]
    ):
        die("assessment source hash or bytes do not exactly match facts.source")

    original_notes(assessment)
    return assessment


def manifest_entries(data):
    entries = data.get("files", data.get("artifacts")) if isinstance(data, dict) else None
    if not isinstance(entries, list):
        die("artifact manifest requires files/artifacts list")
    result = {}
    for entry in entries:
        if not isinstance(entry, dict):
            die("invalid artifact manifest entry")
        path = safe_rel(entry.get("path"))
        size, digest = entry.get("bytes"), entry.get("sha256")
        if (
            path in result
            or not isinstance(size, int)
            or isinstance(size, bool)
            or size < 0
        ):
            die("invalid or duplicate manifest entry")
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(char not in "0123456789abcdef" for char in digest)
        ):
            die("invalid manifest hash")
        result[path] = (size, digest)
    return result


def validate_evidence(root):
    if is_link_or_junction(root):
        die(f"symlink or junction root refused: {root}")
    if not root.is_dir():
        die("evidence-dir is not a directory")

    for current, directories, files in os.walk(root, topdown=True, followlinks=False):
        current_path = Path(current)
        for name in directories + files:
            path = current_path / name
            if is_link_or_junction(path):
                die(f"symlink or junction refused: {path}")
            lowered = name.lower()
            if any(
                word in lowered
                for word in (".env", "secret", "password", ".pem", ".key")
            ):
                die(f"secret-like artifact refused: {path}")

    allowed = set(RAW_FIXED)
    receipts = root / "raw/parent-receipts"
    if not receipts.is_dir():
        die("raw/parent-receipts is required")
    for path in receipts.iterdir():
        if not path.is_file() or path.suffix != ".json":
            die("parent-receipts may contain only JSON files")
        allowed.add(path.relative_to(root).as_posix())

    for subtree in ("expert-packet", "analog-states"):
        base = root / subtree
        manifests = [
            path for path in base.iterdir()
            if path.is_file() and path.name in MANIFEST_NAMES
        ] if base.is_dir() else []
        if len(manifests) != 1 or not (base / "index.html").is_file():
            die(f"{subtree} requires index.html and exactly one manifest")
        entries = manifest_entries(load_json(manifests[0]))
        actual = {
            path.relative_to(base).as_posix(): path
            for path in base.rglob("*")
            if path.is_file() and path != manifests[0]
        }
        if set(entries) != set(actual):
            die(f"{subtree} manifest has missing or extra paths")
        for rel, path in actual.items():
            size, digest = entries[rel]
            if path.stat().st_size != size or sha(path) != digest:
                die(f"{subtree} artifact mismatch: {rel}")
            allowed.add(f"{subtree}/{rel}")
        allowed.add(f"{subtree}/{manifests[0].name}")

    actual = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file()
    }
    if actual != allowed:
        die(f"unexpected or missing evidence files: {sorted(actual ^ allowed)}")


def esc(value):
    return html.escape(str(value), quote=True)


def md(value):
    return (
        html.escape(str(value), quote=True)
        .replace("|", "\\|")
        .replace("\r", " ")
        .replace("\n", " ")
    )


def status_counts(criteria):
    return {
        name: sum(item["status"] == name for item in criteria)
        for name in ("pass", "failed", "pending")
    }


def freshness_text(freshness):
    labels = {
        "current_policy": "정책",
        "current_implementation": "구현",
        "source_runtime_current": "source runtime",
    }
    return [
        f"{labels[key]}: {'현재 기준과 일치' if freshness[key] else '현재 기준과 불일치'}"
        for key in labels
    ]


def render(facts, assessment):
    counts = status_counts(facts["criteria"])
    tests = facts["tests"]
    cards = "".join(
        f"<div class=card><b>{esc(key.upper())}</b><strong>{value}</strong></div>"
        for key, value in counts.items()
    )
    cards += (
        f"<div class=card><b>TESTS</b><strong>{tests['tests_run']}</strong></div>"
        f"<div class=card><b>TEST SKIPPED</b><strong>{tests['skipped']}</strong></div>"
    )
    rows = "".join(
        "<tr>" + "".join(
            f"<td>{esc(item[key])}</td>"
            for key in (
                "id", "title", "status", "reason", "observed_text", "required_text",
            )
        ) + "</tr>"
        for item in facts["criteria"]
    )
    tables = ""
    for table in facts["measurements"]:
        head = "".join(f"<th>{esc(value)}</th>" for value in table["headers"])
        body = "".join(
            "<tr>" + "".join(f"<td>{esc(value)}</td>" for value in row) + "</tr>"
            for row in table["rows"]
        )
        tables += (
            f"<h3>{esc(table['title'])}</h3><div class=scroll><table>"
            f"<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"
        )
    outcome = (
        f"실패 {tests['failures']}건, 오류 {tests['errors']}건"
        if tests["failures"] or tests["errors"]
        else "실패 및 오류 0건"
    )
    feature_items = "".join(
        f"<li>{esc(feature)}</li>" for feature in facts["software_features"]
    )
    note_items = "".join(
        f"<li>{esc(note)}</li>" for note in original_notes(assessment)
    )
    freshness_items = "".join(
        f"<li>{esc(item)}</li>" for item in freshness_text(facts["freshness"])
    )
    api = facts["api_usage"]
    api_rows = [
        ("개발 known tokens", api["development_known_tokens"]),
        ("제품 known tokens", api["product_known_tokens"]),
        ("known tokens 합계", api["total_known_tokens"]),
        ("사용량을 알 수 없는 호출", api["unknown_usage_calls"]),
        ("최근 quota 추정값", api["latest_quota_estimate"]),
        ("quota 설명", api["quota_note"]),
    ]
    api_body = "".join(
        f"<tr><th>{esc(label)}</th><td>{esc(value)}</td></tr>"
        for label, value in api_rows
    )
    css = """body{margin:0;font:15px/1.55 system-ui;color:#16322b;background:#f4f8f6}header{background:#102f46;color:white;padding:28px}main{max-width:1180px;margin:auto;padding:22px}.cards{display:flex;gap:10px;flex-wrap:wrap}.card{background:white;border-top:4px solid #23805f;padding:12px 20px;min-width:110px;box-shadow:0 2px 9px #ccd}.card strong{display:block;font-size:28px}section{background:white;margin:18px 0;padding:20px;border-radius:8px}.scroll{overflow:auto}table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccd8d3;padding:8px;vertical-align:top}th{background:#173f55;color:white}a{color:#176b50}.response{height:42px;border-bottom:1px solid #777;margin:8px 0 16px}.warn{border-left:5px solid #d19a22;padding:10px}@media(max-width:650px){header{padding:18px}main{padding:10px}th,td{min-width:120px}}"""
    return f"""<!doctype html><html lang=ko><meta charset=utf-8><meta name=viewport content="width=device-width,initial-scale=1"><title>최종 전달 보고서</title><style>{css}</style><header><h1>최종 전달 보고서</h1><p>{esc(facts['delivery_status'])} · 과학적 승인 아님</p></header><main><div class=cards>{cards}</div>
<section><h2>판정 요약</h2><p class=warn>알려진 calibration 실패만으로 새 후보의 비효능을 판정할 수 없으며 E3 winner는 선정되지 않았습니다. 테스트 통과 여부와 기능 목록은 과학적 완결성을 의미하지 않습니다.</p><p>테스트 명령: {esc(tests['command'])}; {outcome}; skipped {tests['skipped']}건; exit {tests['exit_code']}.</p></section>
<section><h2>소프트웨어 기능</h2><ul>{feature_items}</ul></section>
<section><h2>API 사용량</h2><p>known 사용량 합계와 unknown 호출 수를 분리했습니다. quota 값은 추정치이며 검증된 잔액이 아닙니다.</p><table><tbody>{api_body}</tbody></table></section>
<section><h2>14개 기준</h2><div class=scroll><table><thead><tr><th>ID</th><th>제목</th><th>상태</th><th>이유</th><th>관찰</th><th>요구</th></tr></thead><tbody>{rows}</tbody></table></div></section>
<section><h2>측정값</h2>{tables}</section>
<section><h2>전문가 후속 질문 및 조치</h2><ul>{note_items}</ul></section>
<section><h2>현재성</h2><ul>{freshness_items}</ul><p>archive 자료는 열람 근거이며 현재 구현 또는 formal approval 여부는 위 현재성 값과 별도로 판단합니다.</p></section>
<section><h2>실행 순서</h2><ol><li>04 파일 검증</li><li>01 UI 실행</li><li>실행된 UI에서 기존 examples 중 하나 선택</li><li>02 no-dock exploratory CPU test</li><li>08 state</li><li>09 report</li><li>03 stop</li></ol><p>실제 dock 및 검증 사실 수치는 측정표에 기록된 값만 사용합니다.</p></section>
<section><h2>UI 예시</h2><p>먼저 01로 애플리케이션을 실행한 뒤 UI의 기존 <code>/design</code> examples에서 예시를 선택합니다. 이 오프라인 보고서에는 실행 중인 애플리케이션 경로로 연결되는 링크가 없습니다.</p></section>
<section><h2>입력 소스</h2><p><a href="{esc(facts['source']['path'])}">reply.docx</a> · {facts['source']['bytes']} bytes · SHA-256 {esc(facts['source']['sha256'])}</p></section>
<section><h2>근거</h2><a href="expert-packet/index.html">expert packet</a> · <a href="analog-states/index.html">states</a> · <a href="raw/assessment.json">assessment</a></section></main></html>"""


def markdown(facts, assessment):
    counts = status_counts(facts["criteria"])
    tests = facts["tests"]
    api = facts["api_usage"]
    lines = [
        "# 최종 전달 보고서",
        "",
        "**DEMO/INCOMPLETE — 과학적 승인 아님**",
        "",
        " · ".join(f"{key}: {value}" for key, value in counts.items()),
        "",
        "테스트 통과 여부와 기능 목록은 과학적 완결성을 의미하지 않습니다.",
        "",
        "## 테스트",
        "",
        "|항목|값|",
        "|---|---|",
        f"|명령|{md(tests['command'])}|",
        f"|실행 수|{tests['tests_run']}|",
        f"|실패|{tests['failures']}|",
        f"|오류|{tests['errors']}|",
        f"|skipped|{tests['skipped']}|",
        f"|exit code|{tests['exit_code']}|",
        f"|tests.log SHA-256|{md(tests['log_sha256'])}|",
        "",
        "## 소프트웨어 기능",
        "",
    ]
    lines.extend(f"- {md(feature)}" for feature in facts["software_features"])
    lines += [
        "",
        "## API 사용량",
        "",
        "known 사용량 합계와 unknown 호출 수를 분리했습니다. quota 값은 추정치이며 검증된 잔액이 아닙니다.",
        "",
        "|항목|값|",
        "|---|---|",
        f"|개발 known tokens|{api['development_known_tokens']}|",
        f"|제품 known tokens|{api['product_known_tokens']}|",
        f"|known tokens 합계|{api['total_known_tokens']}|",
        f"|사용량을 알 수 없는 호출|{api['unknown_usage_calls']}|",
        f"|최근 quota 추정값|{md(api['latest_quota_estimate'])}|",
        f"|quota 설명|{md(api['quota_note'])}|",
        "",
        "## 14개 기준",
        "",
        "|ID|제목|상태|이유|관찰|요구|",
        "|---|---|---|---|---|---|",
    ]
    lines.extend(
        "|" + "|".join(
            md(item[key])
            for key in (
                "id", "title", "status", "reason", "observed_text", "required_text",
            )
        ) + "|"
        for item in facts["criteria"]
    )
    lines += ["", "## 측정값", ""]
    for table in facts["measurements"]:
        lines += [
            f"### {md(table['title'])}",
            "",
            "|" + "|".join(md(value) for value in table["headers"]) + "|",
            "|" + "|".join("---" for _ in table["headers"]) + "|",
        ]
        lines.extend(
            "|" + "|".join(md(value) for value in row) + "|"
            for row in table["rows"]
        )
        lines.append("")
    lines += ["## 전문가 후속 질문 및 조치", ""]
    lines.extend(f"- {md(note)}" for note in original_notes(assessment))
    lines += ["", "## 현재성", ""]
    lines.extend(f"- {md(item)}" for item in freshness_text(facts["freshness"]))
    lines += [
        "",
        "archive 자료는 열람 근거이며 현재 구현 또는 formal approval 여부는 위 현재성 값과 별도로 판단합니다.",
        "",
        "## 실행 순서",
        "",
        "1. 04 파일 검증",
        "2. 01 UI 실행",
        "3. 실행된 UI에서 기존 examples 중 하나 선택",
        "4. 02 no-dock exploratory CPU test",
        "5. 08 state",
        "6. 09 report",
        "7. 03 stop",
        "",
        "실제 dock 및 검증 사실 수치는 측정표에 기록된 값만 사용합니다.",
        "",
        "## UI 예시",
        "",
        "먼저 01로 애플리케이션을 실행한 뒤 UI의 기존 `/design` examples에서 예시를 선택합니다. 이 오프라인 보고서에는 실행 중인 애플리케이션 경로로 연결되는 링크가 없습니다.",
        "",
        "## 입력 소스",
        "",
        f"[reply.docx]({md(facts['source']['path'])}) · {facts['source']['bytes']} bytes · SHA-256 {md(facts['source']['sha256'])}",
        "",
        "## 근거",
        "",
        "[expert packet](expert-packet/index.html) · [states](analog-states/index.html) · [assessment](raw/assessment.json)",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--facts", required=True)
    parser.add_argument("--evidence-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    facts_path = Path(args.facts)
    evidence = Path(args.evidence_dir)
    output = Path(args.output)
    if output.exists():
        die("output must be new and nonexistent")

    facts = load_json(facts_path)
    validate_evidence(evidence)
    assessment = validate_facts(facts, evidence)
    report_html = render(facts, assessment)
    report_md = markdown(facts, assessment)

    output.mkdir(parents=True)
    for child in evidence.iterdir():
        if child.is_dir():
            shutil.copytree(child, output / child.name)
        else:
            shutil.copy2(child, output / child.name)
    shutil.copyfile(facts_path, output / "facts.json")
    (output / "report.html").write_text(report_html, encoding="utf-8")
    (output / "report.md").write_text(report_md, encoding="utf-8")
    follow = report_html.replace(
        "<title>최종 전달 보고서</title>",
        "<title>전문가 후속 응답</title>",
    )
    (output / "expert-followup.html").write_text(follow, encoding="utf-8")
    files = {
        path.relative_to(output).as_posix(): {
            "bytes": path.stat().st_size,
            "sha256": sha(path),
        }
        for path in sorted(output.rglob("*"))
        if path.is_file() and path.name != "report-manifest.json"
    }
    (output / "report-manifest.json").write_text(
        json.dumps({"files": files}, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
