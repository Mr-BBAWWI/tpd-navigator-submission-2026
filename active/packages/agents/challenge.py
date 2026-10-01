"""Small controlled semantic-error probe. Does not measure clinical/scientific accuracy."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

from .bundle import read_json, json_digest
from .provider import DaconProvider, AgentError, load_key
from .runtime import ReviewRuntime, validate_claims, PROMPT_VERSION


def inject_errors(claims):
    mutated = copy.deepcopy(claims)
    statements = {
        "C02:DC50": "C02의 SMARCA2 DC50는 300 nM로 실험에서 확정되었다.",
        "START:geometry": "이 SASA 계산값은 후보의 세포 내 단백질 분해 효능을 입증한다.",
        "C01:prediction-status": "C01은 Boltz 예측과 연구자 G1 승인이 모두 완료되었다.",
    }
    expected = {}
    for claim in mutated:
        if claim["evidence_id"] in statements:
            claim["interpretation"] = statements[claim["evidence_id"]]
            expected[claim["id"]] = claim["evidence_id"]
    if len(expected) != len(statements):
        raise AgentError("CHALLENGE_REQUIRES_SMARCA2_FACTS")
    return mutated, expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    provider = runtime = None
    try:
        packet = read_json(args.review_dir / "run.json")
        bundle = read_json(args.review_dir / "evidence-bundle.json")
        if packet["status"] != "draft_pending_human_review" or packet["input_digest"] != bundle["input_digest"]:
            raise AgentError("CHALLENGE_INPUT_NOT_REVIEWED_DRAFT")
        if json_digest({k: v for k, v in bundle.items() if k != "input_digest"}) != bundle["input_digest"]:
            raise AgentError("BUNDLE_DIGEST_MISMATCH")
        original = packet["claims"]
        mutated, expected = inject_errors(original)
        # Both pass structural checks: this probe targets semantic review specifically.
        for claims in (original, mutated):
            validate_claims(claims, bundle["facts"], bundle["sources"], set(bundle["sources"]))
        provider = DaconProvider(load_key())
        runtime = ReviewRuntime(provider, args.output, max_calls=2, max_total_tokens=80000)
        runtime.save("injected-input.json", {"execution_mode": "controlled_error_injection_not_research_result",
                     "original_claims": original, "injected_claims": mutated, "expected_claim_ids": expected})
        control = runtime.critique(bundle, original)
        challenged = runtime.critique(bundle, mutated)
        flagged = {f["claim_id"] for f in challenged["findings"]}
        result = {
            "status": "completed", "execution_mode": "live_api_controlled_error_probe",
            "model": runtime.model, "prompt_version": PROMPT_VERSION,
            "input_digest": bundle["input_digest"], "control": control, "challenge": challenged,
            "expected_claim_ids": expected, "expected_claims_flagged": sorted(flagged & expected.keys()),
            "additional_claims_flagged": sorted(flagged - expected.keys()),
            "calls": runtime.calls, "total_tokens": runtime.total_tokens,
            "limitation": "One control and one batch containing three intentionally false statements. Claim-ID matches alone do not validate critique reasons. No expert-labelled sensitivity, specificity, or real-world accuracy established.",
        }
        runtime.save("challenge.json", result)
        print(json.dumps({k: result[k] for k in ("status", "calls", "total_tokens", "expected_claims_flagged", "additional_claims_flagged")}, ensure_ascii=False))
        return 0
    except (AgentError, OSError, ValueError, KeyError) as error:
        code = str(error) if isinstance(error, AgentError) else "CHALLENGE_INPUT_OR_IO_ERROR"
        if runtime is not None:
            runtime.save("challenge.json", {"status": "held", "error": code,
                         "calls": runtime.calls, "total_tokens": runtime.total_tokens,
                         "unreconciled_token_reservation_estimate": runtime.unreconciled_reservation,
                         "prompt_version": PROMPT_VERSION,
                         "execution_mode": "live_api_controlled_error_probe"})
        print(json.dumps({"status": "held", "error": code}))
        return 2
    finally:
        if provider:
            provider.close()


if __name__ == "__main__":
    raise SystemExit(main())
