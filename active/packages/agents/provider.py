"""Small DACON Responses adapter; no SDK assumptions or automatic retries."""
from __future__ import annotations

import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
import re
import time
from dataclasses import dataclass

import httpx

BASE_URL = "https://dacon-apim-hackathon-0903.azure-api.net/hackathon/openai/v1"
MODELS = ("gpt-5.6-sol", "gpt-5.6-terra", "gpt-5.6-luna")
DEFAULT_KEY_FILE = Path(__file__).resolve().parents[3] / ".localdata/secrets/dacon-api-key.local.md"
QUOTA_HEADERS = (
    "x-team-remaining-quota-tokens", "x-team-tokens-consumed",
    "x-team-remaining-tokens", "x-team-remaining-requests",
)
_SAFE_TEXT = re.compile(r"^[A-Za-z0-9._:/\[\]-]{1,200}$")


def _utc_now():
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _redact(value, key):
    return value.replace(key, "[REDACTED]") if isinstance(value, str) and key else value


def _safe_text(value, key):
    value = _redact(value, key)
    return value if isinstance(value, str) and _SAFE_TEXT.fullmatch(value) else None


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


def _quota(headers):
    result = {}
    for name in QUOTA_HEADERS:
        value = headers.get(name)
        if value is not None and re.fullmatch(r"\d+", value):
            result[name] = int(value)
    return result


def _metadata(raw, key=""):
    """Return only redacted, allowlisted audit fields."""
    if not isinstance(raw, dict):
        return {}
    clean = {}
    for name in ("requested_model", "returned_model", "response_id", "status"):
        value = _safe_text(raw.get(name), key)
        if value is not None:
            clean[name] = value
    status = raw.get("http_status")
    if type(status) is int and 0 <= status <= 999:
        clean["http_status"] = status
    observed = raw.get("observed_at")
    if isinstance(observed, str):
        try:
            parsed = datetime.fromisoformat(observed.replace("Z", "+00:00"))
            if parsed.utcoffset() == timezone.utc.utcoffset(parsed):
                clean["observed_at"] = parsed.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
        except ValueError:
            pass
    clean["usage"] = _usage(raw.get("usage"))
    clean["usage_known"] = "total_tokens" in clean["usage"]
    quota = raw.get("quota_headers")
    clean["quota_headers"] = {
        name: value for name, value in (quota.items() if isinstance(quota, dict) else ())
        if name in QUOTA_HEADERS and type(value) is int and value >= 0
    }
    clean["quota_is_estimate"] = True
    elapsed = raw.get("elapsed_seconds")
    if type(elapsed) in (int, float) and not isinstance(elapsed, bool) and elapsed >= 0 and math.isfinite(elapsed):
        clean["elapsed_seconds"] = round(elapsed, 3)
    return clean


class AgentError(RuntimeError):
    """A fixed public code plus strictly allowlisted API observations."""

    def __init__(self, code, metadata=None, key=""):
        super().__init__(code)
        self.code = code
        self.metadata = _metadata(metadata, key)


def load_key(path: Path | None = None) -> str:
    key = os.environ.get("DACON_API_KEY", "").strip()
    if not key:
        try:
            content = (path or DEFAULT_KEY_FILE).read_text(encoding="utf-8-sig")
        except OSError:
            raise AgentError("DACON_KEY_NOT_CONFIGURED") from None
        matches = re.findall(r"^```dacon-api-key\s*\n([^\r\n]+)\r?\n```", content, re.M)
        if len(matches) != 1:
            raise AgentError("DACON_KEY_FILE_FORMAT")
        key = matches[0].strip()
    if key == "PASTE_YOUR_DACON_API_KEY_HERE" or not key or not key.isascii() or any(c.isspace() for c in key):
        raise AgentError("DACON_KEY_NOT_CONFIGURED")
    return key


@dataclass
class Completion:
    text: str
    metadata: dict


class DaconProvider:
    mode = "live_dacon_api"

    def __init__(self, key: str, transport=None):
        self._key = key
        self._client = httpx.Client(
            timeout=httpx.Timeout(90, connect=15), follow_redirects=False,
            trust_env=False, transport=transport,
        )

    def close(self):
        self._client.close()

    def _observed(self, response, model, started):
        metadata = {
            "requested_model": model,
            "http_status": response.status_code,
            "observed_at": _utc_now(),
            "usage": {},
            "quota_headers": _quota(response.headers),
            "elapsed_seconds": time.monotonic() - started,
        }
        try:
            data = response.json()
        except (ValueError, TypeError):
            return None, _metadata(metadata, self._key)
        if not isinstance(data, dict):
            return None, _metadata(metadata, self._key)
        for source, target in (("id", "response_id"), ("model", "returned_model"), ("status", "status")):
            value = _safe_text(data.get(source), self._key)
            if value is not None:
                metadata[target] = value
        metadata["usage"] = _usage(data.get("usage"))
        return data, _metadata(metadata, self._key)

    def complete(self, *, model: str, instructions: str, context: dict, max_output_tokens: int) -> Completion:
        if model not in MODELS:
            raise AgentError("MODEL_NOT_ALLOWED")
        payload = {
            "model": model, "instructions": instructions,
            "input": json.dumps(context, ensure_ascii=False),
            "max_output_tokens": max_output_tokens, "store": False,
        }
        started = time.monotonic()
        try:
            response = self._client.post(
                BASE_URL + "/responses", json=payload,
                headers={"api-key": self._key, "Authorization": "Bearer " + self._key},
            )
        except httpx.HTTPError:
            raise AgentError("API_TRANSPORT_ERROR_USAGE_UNKNOWN") from None

        data, metadata = self._observed(response, model, started)
        if response.status_code != 200:
            raise AgentError(
                f"API_HTTP_{response.status_code}_NO_RETRY_USAGE_UNRECONCILED",
                metadata, self._key,
            )
        if data is None:
            raise AgentError("API_RESPONSE_FORMAT", metadata, self._key)
        if data.get("status") != "completed":
            raise AgentError("API_INCOMPLETE_USAGE_UNRECONCILED", metadata, self._key)

        parts = []
        output_items = data.get("output", [])
        if not isinstance(output_items, list):
            raise AgentError("API_UNEXPECTED_OUTPUT", metadata, self._key)
        for item in output_items:
            if not isinstance(item, dict):
                raise AgentError("API_UNEXPECTED_OUTPUT", metadata, self._key)
            if item.get("type") == "reasoning":
                continue
            if item.get("type") != "message" or item.get("role") != "assistant":
                raise AgentError("API_UNEXPECTED_OUTPUT", metadata, self._key)
            content = item.get("content", [])
            if not isinstance(content, list):
                raise AgentError("API_UNEXPECTED_CONTENT", metadata, self._key)
            for part in content:
                if not isinstance(part, dict):
                    raise AgentError("API_UNEXPECTED_CONTENT", metadata, self._key)
                if part.get("type") == "refusal":
                    raise AgentError("API_REFUSAL", metadata, self._key)
                if part.get("type") != "output_text" or not isinstance(part.get("text"), str):
                    raise AgentError("API_UNEXPECTED_CONTENT", metadata, self._key)
                parts.append(_redact(part["text"], self._key))
        output = "".join(parts)
        if not output:
            raise AgentError("API_EMPTY_OUTPUT", metadata, self._key)
        if "total_tokens" not in metadata["usage"]:
            raise AgentError("API_USAGE_UNKNOWN", metadata, self._key)
        return Completion(output, metadata)
