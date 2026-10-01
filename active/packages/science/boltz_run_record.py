"""Record a completed attempt's declared settings and file hashes; never launch inference."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
from pathlib import Path

from .evidence_common import schema_check
from .handoff import check, encoded, parse


def file_hash(path):
    if path is None:
        return None
    check(path.is_file() and not path.is_symlink(), "RUN_RECORD_REGULAR_FILE_REQUIRED")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def environment():
    # A device/version query only. No model or checkpoint is loaded.
    import torch
    return {"boltz": importlib.metadata.version("boltz"), "rdkit": importlib.metadata.version("rdkit"),
            "torch": torch.__version__, "cuda": torch.version.cuda,
            "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "unavailable_at_recording"}


def record(*, run_id, compound_id, requested_samples, exit_code, settings, runtime, checkpoint=None, log=None, msa_dir=None):
    hashes = {}
    if msa_dir is not None:
        check(msa_dir.is_dir() and not msa_dir.is_symlink(), "RUN_RECORD_MSA_DIRECTORY_REQUIRED")
        for path in sorted(msa_dir.rglob("*")):
            check(not path.is_symlink(), "RUN_RECORD_MSA_SYMLINK")
            if path.is_file():
                hashes[path.relative_to(msa_dir).as_posix()] = file_hash(path)
    value = {"format": "tpd-boltz-run/0.1.0-draft", "run_id": run_id, "compound_id": compound_id,
             "origin": "boltz_inference_reported", "requested_samples": requested_samples, "exit_code": exit_code,
             "settings": settings, "environment": runtime, "weights_sha256": file_hash(checkpoint),
             "msa_files_sha256": hashes, "log_sha256": file_hash(log)}
    schema_check(value, "b_boltz_run.schema.json")
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--compound-id", choices=("C01", "C02"), required=True)
    parser.add_argument("--samples", type=int, required=True)
    parser.add_argument("--exit-code", type=int, required=True)
    parser.add_argument("--settings", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--log", type=Path)
    parser.add_argument("--msa-dir", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        check(not args.out.exists(), "OUTPUT_EXISTS")
        value = record(run_id=args.run_id, compound_id=args.compound_id, requested_samples=args.samples,
                       exit_code=args.exit_code, settings=parse(args.settings.read_bytes()), runtime=environment(),
                       checkpoint=args.checkpoint, log=args.log, msa_dir=args.msa_dir)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("xb") as stream:
            stream.write(encoded(value))
        print("Recorded caller-reported attempt; no inference/approval authenticated: " + str(args.out))
    except (OSError, ValueError, ImportError) as error:
        parser.exit(2, f"{type(error).__name__}: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
