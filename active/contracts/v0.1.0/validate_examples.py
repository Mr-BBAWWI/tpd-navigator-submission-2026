"""Validate boundary examples. No network, API, GPU or scientific calculation."""

import copy
import json
from pathlib import Path

from jsonschema import Draft202012Validator


ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main():
    schema = read("module_boundary.schema.json")
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)
    request = read("examples/literature.request.json")
    partial = read("examples/literature.partial-result.json")
    success = copy.deepcopy(partial)
    success.update(status="succeeded", issues=[])
    failed = copy.deepcopy(partial)
    failed.update(status="failed", output_artifacts=[])
    positive = [request, partial, success, failed]
    for example in positive:
        validator.validate(example)

    negative = []

    def case(name, base, mutate):
        item = copy.deepcopy(base)
        mutate(item)
        negative.append((name, item))

    case("unknown field", request, lambda x: x.update(api_key="fixture-only"))
    case("unregistered operation", request, lambda x: x.update(operation_id="run_shell"))
    case("negative budget", request, lambda x: x["limits"].update(max_gpu_seconds=-1))
    case("public GPU budget", request, lambda x: (x.update(execution_mode="public_explore"), x["limits"].update(max_gpu_seconds=1)))
    case("public GPU operation", request, lambda x: x.update(execution_mode="public_explore", operation_id="predict_ternary"))
    case("replay dispatch", request, lambda x: x.update(execution_mode="public_replay"))
    case("GPU without approval", request, lambda x: (x.update(operation_id="predict_ternary"), x["limits"].update(max_gpu_seconds=1)))
    case("partial without issue", partial, lambda x: x.update(issues=[]))
    case("success without artifact", success, lambda x: x.update(output_artifacts=[]))
    case("invalid hash", partial, lambda x: x["output_artifacts"][0].update(sha256="unknown"))
    case("run status used as result", partial, lambda x: x.update(status="held"))
    for name, example in negative:
        if validator.is_valid(example):
            raise AssertionError(f"Invalid example accepted: {name}")

    operations = read("operations.json")["operations"]
    op_ids = [item["operation_id"] for item in operations]
    assert len(op_ids) == len(set(op_ids)), "Duplicate operation IDs"
    for kind in ("JobSpec", "ModuleResult"):
        assert schema["$defs"][kind]["properties"]["operation_id"]["enum"] == op_ids
    ready = [item for item in operations if item["live_execution_ready"]]
    assert {item["operation_id"] for item in ready} == {"resolve_target", "collect_literature"}
    for item in ready:
        assert item["payload_schema_status"] == "schema_validated"
        assert item["execution_profiles"] == ["local-single-user-i2"]
        assert item["public_dispatch_ready"] is False
        assert (ROOT / item["verification_file"]).is_file()
    for key in ("contract_version", "project_id", "run_id", "job_id", "attempt", "input_revision", "input_digest", "operation_id", "data_mode"):
        assert request[key] == partial[key], f"Example correlation mismatch: {key}"
    assert request["input_digest"] == request["input_manifest"]["sha256"]
    print(f"PASS: schema valid; {len(positive)} positive and {len(negative)} negative cases; {len(op_ids)} operation IDs consistent.")
    print("Envelope validation only. Payload schemas, artifact access, runtime authorization and scientific correctness remain untested.")


if __name__ == "__main__":
    main()
