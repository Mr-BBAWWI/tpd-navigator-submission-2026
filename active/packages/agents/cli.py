"""Run from active/: python -m packages.agents.cli --help."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .bundle import build_bundle, digest
from .provider import AgentError, DaconProvider, MODELS, load_key
from .runtime import ReviewRuntime


def main():
    parser = argparse.ArgumentParser(description="Internal agent review experiment; no official approvals.")
    parser.add_argument("command", choices=["prepare", "run", "compare"])
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--architecture", choices=["hierarchical", "single", "no_critic"], default="hierarchical")
    parser.add_argument("--model", choices=MODELS, default="gpt-5.6-sol")
    parser.add_argument("--key-file", type=Path)
    args = parser.parse_args()
    provider = None
    try:
        if args.output.exists():
            raise AgentError("OUTPUT_ALREADY_EXISTS")
        bundle = build_bundle(args.case_dir, args.raw_dir)
        if args.command == "prepare":
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x", encoding="utf-8") as stream:
                json.dump(bundle, stream, ensure_ascii=False, indent=2)
            print(json.dumps({"status": "prepared", "input_digest": bundle["input_digest"], "facts": len(bundle["facts"])}))
            return 0
        provider = DaconProvider(load_key(args.key_file))
        architectures = ["single", "no_critic", "hierarchical"] if args.command == "compare" else [args.architecture]
        if args.command == "compare":
            args.output.mkdir(parents=True, exist_ok=False)
        results = []
        for architecture in architectures:
            output = args.output / architecture if args.command == "compare" else args.output
            runtime = ReviewRuntime(provider, output, model=args.model)
            packet = runtime.run(bundle, architecture=architecture,
                                 revalidate=lambda: build_bundle(args.case_dir, args.raw_dir)["input_digest"])
            summary = {k: packet.get(k) for k in ("architecture", "execution_mode", "status", "input_digest", "calls", "total_tokens", "elapsed_seconds", "revision_rounds", "error")}
            summary["claims"] = len(packet["claims"])
            results.append(summary)
            print(json.dumps(summary, ensure_ascii=False), flush=True)
            # Do not turn a shared authentication/quota/transport problem into more requests.
            if packet.get("error", "").startswith("API_"):
                break
        implementation = {p.name: digest(p.read_bytes()) for p in Path(__file__).parent.glob("*.py")}
        manifest = {"model": args.model, "runs": results, "implementation_sha256": implementation,
                    "contract_sha256": digest((Path(__file__).resolve().parents[2] / "contracts/drafts/agent_review.schema.json").read_bytes()),
                    "interpretation": "Feasibility and cost observations only; no expert-scored quality or architecture superiority established."}
        (args.output / "experiment.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
        return 0 if len(results) == len(architectures) and all(r["status"] != "held" for r in results) else 2
    except (AgentError, OSError, ValueError) as error:
        code = str(error) if isinstance(error, AgentError) else "LOCAL_INPUT_OR_IO_ERROR"
        print(json.dumps({"status": "held", "error": code}))
        return 2
    finally:
        if provider:
            provider.close()


if __name__ == "__main__":
    raise SystemExit(main())
