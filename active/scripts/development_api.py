"""Explicit, audited API assistance for this development bundle; never executes model text."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import httpx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from packages.agents.provider import AgentError, DaconProvider, load_key


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("request", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    request_bytes = args.request.read_bytes()
    request_hash = hashlib.sha256(request_bytes).hexdigest()
    request = json.loads(request_bytes.decode("utf-8"))
    provider = DaconProvider(load_key())
    provider._client.timeout = httpx.Timeout(300, connect=15)
    try:
        result = provider.complete(
            model="gpt-5.6-sol", instructions=request["instructions"],
            context=request["context"],
            max_output_tokens=request.get("max_output_tokens", 16000),
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result.text, encoding="utf-8")
        metadata = dict(
            result.metadata,
            request_sha256=request_hash,
            output_sha256=hashlib.sha256(result.text.encode()).hexdigest(),
            applied_automatically=False,
        )
        args.output.with_suffix(".metadata.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(metadata, ensure_ascii=True, sort_keys=True))
    except Exception as error:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        api_metadata = error.metadata if isinstance(error, AgentError) else {}
        usage = api_metadata.get("usage", {}) if isinstance(api_metadata.get("usage", {}), dict) else {}
        usage_known = api_metadata.get("usage_known") is True or type(usage.get("total_tokens")) is int
        failure = {
            "exception_code": error.code if isinstance(error, AgentError) else "DEVELOPMENT_REQUEST_ERROR",
            "usage_known": usage_known,
            "usage": usage,
            "automatic_retry": False,
            "request_sha256": request_hash,
            "api_metadata": api_metadata,
            "failure_marker": "after_request",
        }
        args.output.with_suffix(".failure.json").write_text(
            json.dumps(failure, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        raise
    finally:
        provider.close()


if __name__ == "__main__":
    main()
