"""6BOY BRD4–CRBN Boltz2 독립 구조 재현 실행기."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.boltz_worker import BoltzWorkerError, run_seed


def _root_path(value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--reference", default="cases/design_sources/6BOY.cif",
                        help="저장소 루트 기준 6BOY mmCIF 경로")
    result.add_argument("--metadata", default="cases/design_sources/crbn_benchmark.json",
                        help="저장소 루트 기준 참조 메타데이터")
    result.add_argument("--output", required=True, help="새 seed 디렉터리와 evidence JSON을 쓸 로컬 경로")
    result.add_argument("--boltz-executable", required=True, help="격리 환경의 boltz 실행 파일 경로")
    result.add_argument("--checkpoint", required=True, help="사용자 제공 Boltz checkpoint 경로")
    result.add_argument("--cache", required=True, help="Boltz 로컬 cache 경로")
    result.add_argument("--seed", type=int, action="append", dest="seeds",
                        help="반복 가능. 생략하면 seed 0 한 번 실행")
    result.add_argument("--msa-mode", choices=("server", "single_sequence"), required=True,
                        help="server는 명시적 MSA 서버 사용, single_sequence는 msa: empty")
    result.add_argument("--accelerator", choices=("gpu", "cpu"), default="gpu")
    result.add_argument("--max-msa-seqs", type=int, default=8192)
    result.add_argument("--timeout-seconds", type=float, default=21600.0)
    result.add_argument("--ligand-groups", help="선택 사항: group 이름 -> RN6 CCD atom 이름 배열 JSON")
    result.add_argument("--frozen-seed-directory",
                        help="검증된 이전 seed receipt 및 boltz_output에서 msa/processed만 재사용")
    return result


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    seeds = [0] if args.seeds is None else args.seeds
    if any(type(seed) is not int or seed < 0 for seed in seeds) or len(seeds) != len(set(seeds)):
        parser().error("seed는 중복 없는 0 이상의 정수여야 합니다")
    paths = {name: _root_path(getattr(args, name)) for name in
             ("reference", "metadata", "output", "boltz_executable", "checkpoint", "cache")}
    groups = _root_path(args.ligand_groups) if args.ligand_groups else None
    frozen = _root_path(args.frozen_seed_directory) if args.frozen_seed_directory else None
    output = paths["output"]
    output.mkdir(parents=True, exist_ok=True)
    evidence_path = output / "crbn_calibration_evidence.json"
    if evidence_path.exists():
        parser().error(f"기존 evidence를 덮어쓰지 않습니다: {evidence_path}")

    receipts = []
    try:
        for seed in seeds:
            receipts.append(run_seed(reference=paths["reference"], metadata=paths["metadata"],
                                     output_root=output, executable=paths["boltz_executable"],
                                     checkpoint=paths["checkpoint"], cache=paths["cache"], seed=seed,
                                     msa_mode=args.msa_mode, accelerator=args.accelerator,
                                     max_msa_seqs=args.max_msa_seqs, timeout=args.timeout_seconds,
                                     ligand_groups=groups, frozen_seed_directory=frozen))
    except (OSError, ValueError, BoltzWorkerError) as error:
        print(f"실행 준비 실패: {type(error).__name__}: {error}", file=sys.stderr)
        return 2

    success_count = sum(receipt.get("status") == "success" for receipt in receipts)
    evidence = {
        "format": "tpd-crbn-calibration-evidence/1.0",
        "case": "CRBN-dBET6-6BOY",
        "reference": {"path": str(paths["reference"]), "sha256": _hash(paths["reference"]),
                      "known_structure_training_overlap_possible": True},
        "msa_mode": args.msa_mode,
        "msa_limitation_ko": ("단일 서열 입력이며 공진화 MSA 정보가 없다. VHL 실행과 같은 MSA라고 주장하지 않는다."
                              if args.msa_mode == "single_sequence" else
                              "명시적으로 서버 MSA를 요청했다. 생성된 MSA hash는 각 receipt를 검토해야 하며 VHL과 동일하다고 간주하지 않는다."),
        "requested_seeds": seeds, "successful_runs": success_count,
        "failed_runs": len(receipts) - success_count,
        "receipts": [{"seed": item["seed"], "status": item["status"],
                      "path": str(output / f"seed-{item['seed']}" / "receipt.json"),
                      "sha256": _hash(output / f"seed-{item['seed']}" / "receipt.json")}
                     for item in receipts],
        "calibration_complete": False,
        "M2_integration_status": "검토 대기",
        "해석": "알려진 한 구조에 대한 모델 예측 비교이며 효능 또는 분해 성공의 측정값이 아니다.",
        "limitations": [
            "6BOY가 모델 학습 자료와 겹쳤을 가능성을 배제할 수 없다.",
            "단일 seed는 광범위한 보정이나 seed 분포를 확립하지 않는다." if len(seeds) == 1 else
            "여러 seed도 단일 알려진 구조에 국한되며 일반화를 확립하지 않는다.",
            "모델 confidence는 측정된 구조 재현도와 구분된다.",
            "graph assembly 또는 docking 결과는 효능 승인 근거가 아니다.",
        ],
    }
    with evidence_path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(evidence, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")
    print(f"증거 JSON 기록: {evidence_path}")
    return 0 if success_count == len(receipts) else 1


if __name__ == "__main__":
    raise SystemExit(main())
