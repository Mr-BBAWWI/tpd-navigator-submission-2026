"""Summarize explicitly supplied API audit files and read-only workbench stores."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sqlite3

QUOTA_HEADERS = (
    "x-team-remaining-quota-tokens", "x-team-tokens-consumed",
    "x-team-remaining-tokens", "x-team-remaining-requests",
)
SAFE_ID = re.compile(r"^[A-Za-z0-9._:/-]{1,200}$")
SAFE_ARTIFACT = re.compile(r"^[A-Za-z0-9._-]{1,200}$")
_API_INDICATORS = {"response_id", "requested_model", "returned_model", "usage", "quota_headers", "request_sha256"}
_FAILURE_INDICATORS = {"error_type", "exception_code", "request_sha256", "failure_marker"}


def _dict(value):
    if isinstance(value, dict):
        return value
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
            return parsed if isinstance(parsed, dict) else None
        except (ValueError, TypeError):
            pass
    return None


def _iso(value, fallback):
    if not isinstance(value, str) or len(value) > 40:
        return fallback
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
        if parsed.utcoffset() != timezone.utc.utcoffset(parsed):
            return fallback
        return parsed.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
    except ValueError:
        return fallback


def _stamp(value):
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp()
    except (ValueError, TypeError):
        return 0.0


def _usage(raw):
    result = {}
    if isinstance(raw, dict):
        for name in ("input_tokens", "output_tokens", "total_tokens"):
            value = raw.get(name)
            if type(value) is int and value >= 0:
                result[name] = value
    if all(name in result for name in ("input_tokens", "output_tokens", "total_tokens")):
        if result["total_tokens"] < result["input_tokens"] + result["output_tokens"]:
            result.pop("total_tokens")
    return result


def _extract(value, fallback):
    root = _dict(value) or {}
    candidates = []

    def visit(item):
        item = _dict(item)
        if not item:
            return
        candidates.append(item)
        for name in (
            "api_metadata", "metadata", "provider_metadata", "response", "api_response",
            "before_failure", "after_failure", "before", "after",
        ):
            visit(item.get(name))

    visit(root)
    merged = {}
    for candidate in candidates:
        for name in ("response_id", "usage", "quota_headers", "observed_at"):
            if name in candidate:
                merged[name] = candidate[name]
        if "response_id" not in merged and isinstance(candidate.get("id"), str) and (
                isinstance(candidate.get("usage"), dict) or "object" in candidate):
            merged["response_id"] = candidate["id"]
    response_id = merged.get("response_id")
    if not isinstance(response_id, str) or not SAFE_ID.fullmatch(response_id):
        response_id = None
    quota = {}
    raw_quota = merged.get("quota_headers")
    if isinstance(raw_quota, dict):
        quota = {
            name: value for name, value in raw_quota.items()
            if name in QUOTA_HEADERS and type(value) is int and value >= 0
        }
    return {
        "response_id": response_id,
        "usage": _usage(merged.get("usage")),
        "quota_headers": quota,
        "observed_at": _iso(merged.get("observed_at"), fallback),
        "conflict": False,
        "source_error": None,
    }


def _merge(records, identity, item):
    previous = records.get(identity)
    if previous is None:
        records[identity] = item
        return
    for name, value in item["usage"].items():
        if name in previous["usage"] and previous["usage"][name] != value:
            previous["conflict"] = True
        else:
            previous["usage"][name] = value
    usage = previous["usage"]
    if all(name in usage for name in ("input_tokens", "output_tokens", "total_tokens")):
        if usage["total_tokens"] < usage["input_tokens"] + usage["output_tokens"]:
            previous["conflict"] = True
    if item["conflict"]:
        previous["conflict"] = True
    if item.get("source_error"):
        previous["source_error"] = item["source_error"]
    if _stamp(item["observed_at"]) >= _stamp(previous["observed_at"]):
        previous["observed_at"] = item["observed_at"]
        if item["quota_headers"]:
            previous["quota_headers"] = item["quota_headers"]


def _totals(records):
    result = {
        "calls": len(records), "known_usage_calls": 0, "unknown_usage_calls": 0,
        "partial_usage_calls": 0, "conflict_calls": 0,
        "input_known_calls": 0, "output_known_calls": 0, "total_known_calls": 0,
        "input_tokens": 0, "output_tokens": 0, "total_tokens": 0,
        "source_error_calls": 0,
    }
    for item in records.values():
        if item.get("source_error"):
            result["source_error_calls"] += 1
        if item["conflict"]:
            result["conflict_calls"] += 1
            continue
        usage = item["usage"]
        if not usage:
            result["unknown_usage_calls"] += 1
            continue
        result["known_usage_calls"] += 1
        if len(usage) < 3:
            result["partial_usage_calls"] += 1
        for name, count_name in (
            ("input_tokens", "input_known_calls"),
            ("output_tokens", "output_known_calls"),
            ("total_tokens", "total_known_calls"),
        ):
            if name in usage:
                result[count_name] += 1
                result[name] += usage[name]
    return result


def _store_paths(explicit):
    seen = set()
    for supplied in explicit:
        supplied = Path(supplied)
        candidates = [supplied] if supplied.is_file() else (
            list(supplied.rglob("*.sqlite")) + list(supplied.rglob("*.sqlite3")) + list(supplied.rglob("*.db"))
            if supplied.is_dir() else []
        )
        for candidate in candidates:
            try:
                resolved = candidate.resolve()
            except OSError:
                continue
            if resolved not in seen and resolved.is_file():
                seen.add(resolved)
                yield resolved


def _columns(connection, table):
    return {row[1] for row in connection.execute("PRAGMA table_info(" + table + ")")}


def _read_store(path):
    connection = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        required = {"workbench_calls", "workbench_jobs", "artifacts"}
        if not required <= tables:
            return [], {}, None
        if not {"job_id", "ordinal", "body"} <= _columns(connection, "workbench_calls"):
            return [], {}, None
        jobs = {}
        for row in connection.execute("SELECT id,body FROM workbench_jobs"):
            body = _dict(row["body"])
            if body:
                jobs[row["id"]] = body
        calls = [dict(row) for row in connection.execute("SELECT job_id,ordinal,body FROM workbench_calls")]
        artifacts = {}
        for row in connection.execute("SELECT id,project,metadata FROM artifacts"):
            artifacts[row["id"]] = (row["project"], _dict(row["metadata"]))
        return calls, jobs, artifacts
    finally:
        connection.close()


def _artifact_metadata(root, reference, project, artifacts, add_source):
    if not isinstance(reference, dict) or not isinstance(project, str) or artifacts is None:
        return None, "corruptedref"
    artifact_id = reference.get("artifact_id")
    if not isinstance(artifact_id, str) or not SAFE_ARTIFACT.fullmatch(artifact_id) or artifact_id in {".", ".."}:
        return None, "corruptedref"
    registered = artifacts.get(artifact_id)
    if not registered or registered[0] != project or registered[1] != reference:
        return None, "corruptedref"
    digest = reference.get("sha256")
    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
        return None, "corruptedref"
    blobs = (root / "blobs").resolve()
    path = (blobs / artifact_id).resolve()
    if path.parent != blobs or path.name != artifact_id or not path.is_file():
        return None, "corruptedref"
    try:
        raw = path.read_bytes()
    except OSError:
        return None, "corruptedref"
    add_source(path, raw)
    if hashlib.sha256(raw).hexdigest() != digest:
        return None, "corruptedref"
    try:
        blob = json.loads(raw)
    except (ValueError, UnicodeDecodeError):
        return None, "corruptedref"
    if not isinstance(blob, dict) or not isinstance(blob.get("metadata"), dict):
        return None, "corruptedref"
    return blob["metadata"], None


def summarize(audit_dirs, store_dirs):
    developer, product = {}, {}
    failure_ids = set()
    observations, observation_keys = [], set()
    source_hashes, source_paths = [], set()
    source_errors = {}
    parse_errors = []

    def add_source(path, raw):
        resolved = Path(path).resolve()
        if resolved in source_paths:
            return
        source_paths.add(resolved)
        source_hashes.append({"path": str(resolved), "sha256": hashlib.sha256(raw).hexdigest()})

    def observe(identity, item):
        if not item["quota_headers"]:
            return
        key = (identity, item["observed_at"], json.dumps(item["quota_headers"], sort_keys=True))
        if key not in observation_keys:
            observation_keys.add(key)
            observations.append(item)

    audit_paths = set()
    for directory in map(Path, audit_dirs):
        if directory.is_file() and directory.name.endswith((".metadata.json", ".failure.json")):
            audit_paths.add(directory.resolve())
        elif directory.is_dir():
            audit_paths.update(path.resolve() for path in directory.rglob("*.metadata.json"))
            audit_paths.update(path.resolve() for path in directory.rglob("*.failure.json"))
    for path in sorted(audit_paths):
        try:
            raw = path.read_bytes()
        except OSError:
            continue
        try:
            document = json.loads(raw)
        except (ValueError, UnicodeDecodeError) as error:
            add_source(path, raw)
            parse_errors.append({"path": str(path), "error": type(error).__name__})
            continue
        if not isinstance(document, dict):
            continue
        is_failure = path.name.endswith(".failure.json")
        nested_api = document.get("api_metadata") if isinstance(document.get("api_metadata"), dict) else {}
        recognized = bool(_API_INDICATORS & (set(document) | set(nested_api)))
        if is_failure:
            recognized = recognized and bool(_FAILURE_INDICATORS & set(document))
        if not recognized:
            continue
        add_source(path, raw)
        fallback = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat().replace("+00:00", "Z")
        item = _extract(document, fallback)
        identity = item["response_id"] or "file:" + str(path)
        _merge(developer, identity, item)
        observe("developer:" + identity, item)
        if is_failure:
            failure_ids.add(identity)

    for path in _store_paths(store_dirs):
        try:
            raw = path.read_bytes()
            add_source(path, raw)
            calls, jobs, artifacts = _read_store(path)
        except (OSError, sqlite3.Error):
            source_errors["database"] = source_errors.get("database", 0) + 1
            continue
        if artifacts is None:
            source_errors["schema"] = source_errors.get("schema", 0) + 1
            continue
        fallback = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat().replace("+00:00", "Z")
        for row in calls:
            body = _dict(row.get("body")) or {}
            job = jobs.get(row.get("job_id"), {})
            project = job.get("project_id")
            observed = body.get("completed_at", body.get("requested_at"))
            row_fallback = _iso(observed, fallback)
            body_item = _extract(body, row_fallback)
            item = body_item
            error = None
            if body.get("answer_ref") is not None:
                metadata, error = _artifact_metadata(path.parent, body.get("answer_ref"), project, artifacts, add_source)
                if metadata is not None:
                    authoritative = _extract(metadata, row_fallback)
                    if authoritative["usage"]:
                        item = authoritative
                    else:
                        body_item["response_id"] = authoritative["response_id"] or body_item["response_id"]
                        body_item["quota_headers"] = authoritative["quota_headers"] or body_item["quota_headers"]
                        item = body_item
            if error:
                item["source_error"] = error
                source_errors[error] = source_errors.get(error, 0) + 1
            identity = item["response_id"] or "call:%s:%s:%s" % (
                project or "?", row.get("job_id"), row.get("ordinal")
            )
            _merge(product, identity, item)
            observe("product:" + identity, item)

    observations.sort(key=lambda item: _stamp(item["observed_at"]))
    nonmonotonic = False
    for before, after in zip(observations, observations[1:]):
        for name in ("x-team-remaining-quota-tokens", "x-team-remaining-tokens"):
            if name in before["quota_headers"] and name in after["quota_headers"]:
                if after["quota_headers"][name] > before["quota_headers"][name]:
                    nonmonotonic = True
    latest = observations[-1] if observations else None
    unresolved_failures = sum(
        1 for identity in failure_ids
        if identity not in developer or developer[identity]["conflict"] or not developer[identity]["usage"]
    )
    warnings = ["quota observations are not monotonic"] if nonmonotonic else []
    return {
        "developer_calls": _totals(developer),
        "product_calls": _totals(product),
        "unknown_usage_failures": unresolved_failures,
        "latest_quota": ({
            "values": latest["quota_headers"],
            "observed_at": latest["observed_at"],
            "label": "estimated",
            "nonmonotonic_observations": nonmonotonic,
        } if latest else None),
        "warnings": warnings,
        "quota_caveat": "Estimated observations only; counters may be nonmonotonic and no balance is inferred without a server baseline.",
        "source_errors": source_errors,
        "parse_errors": parse_errors,
        "source_hashes": sorted(source_hashes, key=lambda item: item["path"]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-dir", action="append", default=[])
    parser.add_argument("--store-dir", action="append", default=[])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = summarize(args.audit_dir, args.store_dir)
    rendered = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
        print(json.dumps({
            "developer_calls": result["developer_calls"]["calls"],
            "product_calls": result["product_calls"]["calls"],
            "unknown_usage_failures": result["unknown_usage_failures"],
        }, sort_keys=True))
    else:
        print(rendered)


if __name__ == "__main__":
    main()
