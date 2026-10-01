"""Stable import bridge to the independently versioned I1 contracts."""
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
path = ROOT / "contracts/payloads/v0.1.0/validation.py"
spec = importlib.util.spec_from_file_location("tpd_i1_contracts", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
ArtifactReader = module.ArtifactReader
ContractError = module.ContractError
validate_payload = module.validate_payload
validate_exchange = module.validate_exchange
parse_json = module.parse_json
DOMAIN_URI = module.DOMAIN_URI


def encoded(value):
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def now():
    return datetime.now(timezone.utc).isoformat()


def known(value):
    return {"state": "known", "value": str(value), "reason": None}


def missing(reason="원문에서 확인하지 못했습니다.", state="not_reported"):
    return {"state": state, "value": None, "reason": reason}


def issue(code, message, category="invalid_output", ids=None, retry="do_not_retry"):
    return {"code": code, "category": category, "message": message,
            "affected_artifact_ids": ids or [], "retry_hint": retry}


def payload(kind, **fields):
    value = {"payload_type": kind, "payload_version": "0.1.0", **fields}
    validate_payload(kind, value)
    return value
