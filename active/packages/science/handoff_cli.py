"""Export an offline H1/H2 review package and optionally inspect Boltz drafts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .handoff import Snapshot, artifact_ref, build_material, encoded, match_compound, compare_observation


def render(material, depictions=False):
    lines = ["# B H1/H2 인계 자료", "", "사람 검수 전 자료입니다. G1 승인·GPU 예측 결과가 아닙니다.", "",
             f'사례: `{material["case_id"]}`', f'자료 버전: `{material["material_digest"]}`', "",
             "| 관리 ID | 논문 이름 | 화합물 번호 | PDB / CCD | 부착 원자 |", "|---|---|---|---|---|"]
    for c in material["h1"]["compounds"]:
        lines.append(f'| {c["compound_id"]} | {c["paper"]["name"]} | {c["paper"]["compound_number"]} | {c["structure"]["pdb"]} / {c["structure"]["ccd"]} | {c["attachment_atom"]["ccd_atom_id"]} |')
    lines += ["", "`material.json`에 원문 위치·화합물별 실험값/결측·조건·분자 버전·원자 대응을 담았습니다.",
              "`consumer-example.json`은 B가 준비한 소비 예제이며 동료 A의 실제 추출 결과가 아닙니다.",
              "파일 참조는 입력 CPU snapshot의 basename입니다. 함께 제공된 기존 snapshot과 SHA-256을 대조합니다.", "",
              "G1 검토 대상은 출발 리간드의 부착 가설입니다. C01/C02는 알려진 기준 PROTAC의 재구성 자료이며, 이를 G1 전 Phase1 후보로 자동 등록하지 않습니다.",
              "부착점의 노출·거리 수치는 설명 자료이며 효능·부착 허용 cutoff가 아닙니다.", "",
              "| 기준 후보 | 대칭 대응 수 | 매칭에만 적용한 전하 처리 |", "|---|---:|---|"]
    for r in material["h2"]["reference_reconstructions"]:
        m = r["warhead_mapping"]
        n = m["normalization_for_matching_only"]
        lines.append(f'| {r["compound_id"]} | {m["equivalent_mapping_count"]} | {n["atom"]}: {n["before"]["formal_charge"]} → {n["after_formal_charge"]}; 원본 보존 |')
    lines += ["", "실험값과 조건:", ""]
    for r in material["h1"]["observations"]:
        o = r["original"]
        lines.append(f'- {o["candidate_id"]} {o["endpoint"]}: {o["value"] if o["value"] is not None else "미추출"} {o["unit"]}; {o.get("missing_reason", "")}' )
    lines += ["", "미확인·미연결:", ""] + ["- " + x for x in material["h2"]["unresolved"]]
    if depictions:
        lines += ["", "부착 원자 확인 그림 (2D):", ""]
        for c in material["h1"]["compounds"]:
            lines.append(f'![{c["compound_id"]} {c["attachment_atom"]["ccd_atom_id"]}]({c["compound_id"]}-attachment.svg)')
    return "\n".join(lines) + "\n"


def export(directory, raw_directory, output, include_preflight=False):
    output = Path(output)
    if output.exists():
        raise ValueError("OUTPUT_EXISTS: choose a new output directory")
    snapshot = Snapshot(directory, raw_directory)
    material = build_material(snapshot)
    candidate = next(c for c in material["h1"]["compounds"] if c["role"] == "known_reference_protac")
    query = {"doi": candidate["paper"]["doi"], "paper_name": candidate["paper"]["name"],
             "compound_number": candidate["paper"]["compound_number"], "expected_digest": material["material_digest"]}
    example = {"producer": "B_prepared_consumer_example_not_A_extraction", "query": query,
               "result": match_compound(material, **query)}
    original = next((r["original"] for r in material["h1"]["observations"] if r["original"]["candidate_id"] == candidate["compound_id"]), None)
    if original is not None:
        example["observation_comparison"] = compare_observation(material, query, original)
    files = {"material.json": encoded(material), "consumer-example.json": encoded(example),
             "README.md": render(material, include_preflight).encode("utf-8")}
    if include_preflight:
        from .boltz_preflight import preflight
        from .depictions import attachment_depictions
        files["boltz-preflight.json"] = encoded(preflight(snapshot))
        files.update(attachment_depictions(snapshot, material))
    # All validation precedes output creation. A concurrent existing export is refused.
    output.mkdir(parents=True, exist_ok=False)
    for name, data in files.items():
        with (output / name).open("xb") as stream:
            stream.write(data)
    manifest = {"format": "tpd-b-review-export/0.1.0-draft", "source_snapshot_sha256": material["source_snapshot"]["handoff_sha256"],
                "material_ref": artifact_ref(files["material.json"], schema="urn:tpd-navigator:b-review-material:0.1.0-draft"),
                "files": [{"path": n, "bytes": len(d), "ref": artifact_ref(d, {".md": "text/markdown", ".svg": "image/svg+xml"}.get(Path(n).suffix, "application/json"))}
                          for n,d in files.items()], "registered_in_M2": False}
    with (output / "manifest.json").open("xb") as stream:
        stream.write(encoded(manifest))
    return material


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", required=True, type=Path)
    parser.add_argument("--raw", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args(argv)
    try:
        material = export(args.snapshot, args.raw, args.out, args.preflight)
    except (OSError, ValueError, KeyError, TypeError, ImportError) as error:
        parser.exit(2, f"Handoff not exported: {error}\n")
    print(json.dumps({"case_id": material["case_id"], "material_digest": material["material_digest"],
                      "human_review": "pending", "dispatch_authorized": False}, ensure_ascii=False))


if __name__ == "__main__":
    main()
