"""Candidate evidence projections and Korean review report; no scientific execution."""
from __future__ import annotations

import copy
import html

from .boltz_collection import empty_collection, validate_collection
from .evidence_common import schema_check, seal, verify_seal
from .handoff import artifact_ref, check, compare_observation, encoded, validate_material


def reference(path, value, schema_id):
    return {"path": path, "ref": artifact_ref(encoded(value), schema=schema_id)}


def reconcile(material, exchange, *, allow_synthetic=False):
    if exchange is None:
        return {"status": "not_provided", "producer": None, "records": [],
                "meaning": "No A-produced example supplied; B example tests do not establish A integration."}
    schema_check(exchange, "b_literature_exchange.schema.json")
    check(exchange["material_digest"] == material["material_digest"], "LITERATURE_STALE_MATERIAL")
    check(allow_synthetic or exchange["data_mode"] != "synthetic_test", "SYNTHETIC_LITERATURE_NOT_ALLOWED")
    check(len({r["record_id"] for r in exchange["records"]}) == len(exchange["records"]), "DUPLICATE_LITERATURE_RECORD")
    records = []
    for row in exchange["records"]:
        check(row["source"]["doi"] == row["query"]["doi"], "LITERATURE_SOURCE_DOI_MISMATCH")
        comparison = compare_observation(material, {**row["query"], "expected_digest": material["material_digest"]}, row["observation"])
        flags = []
        if row["observation"].get("locator") != row["source"]["locator"]:
            flags.append("INCOMING_LOCATOR_CONFLICT")
        known_hashes = {o["citation"]["source_sha256"] for o in material["h1"]["observations"] if o["citation"]}
        if row["source"]["sha256"] not in known_hashes:
            flags.append("INCOMING_SOURCE_VERSION_NOT_VERIFIED")
        records.append({"record_id": row["record_id"], "query": copy.deepcopy(row["query"]),
                        "source_reported": copy.deepcopy(row["source"]), "comparison": comparison,
                        "source_flags": flags, "human_review": "pending"})
    return {"status": "compared_pending_review", "producer": exchange["producer"], "records": records,
            "meaning": "Comparison only; incoming observations and source bytes are not independently extracted or approved here."}


def prediction_view(cid, collection):
    runs = [r for r in collection["runs"] if r["compound_id"] == cid]
    groups = {}
    samples = []
    for run in runs:
        groups.setdefault(run["condition_digest"], []).append(run["run_id"])
        samples += [{"run_id": run["run_id"], "condition_digest": run["condition_digest"],
                     "origin": run["origin"], **copy.deepcopy(m)} for m in run["models"]]
    return {"status": "not_run" if not runs else "available_unreviewed" if all(r["readiness"] == "available_unreviewed" for r in runs) else "needs_review",
            "meaning": "State of supplied records only; no official M2 execution/approval status.",
            "run_ids": [r["run_id"] for r in runs],
            "condition_groups": [{"condition_digest": k, "run_ids": v} for k, v in sorted(groups.items())],
            "counts": {"requested": sum(r["reported_execution"]["requested_samples"] for r in runs),
                       **{s: sum(m["status"] == s for m in samples) for s in ("parsed", "missing", "invalid")},
                       "unexpected": sum(not m["expected"] for m in samples)},
            "samples": samples, "automatic_selection": None,
            "run_warnings": {r["run_id"]: r["warnings"] for r in runs}}


def build_evidence(material, collection=None, exchange=None, *, allow_synthetic=False):
    validate_material(material)
    collection = empty_collection(material) if collection is None else collection
    validate_collection(collection, material, allow_synthetic=allow_synthetic)
    reconciliation = reconcile(material, exchange, allow_synthetic=allow_synthetic)
    start = next(c for c in material["h1"]["compounds"] if c["role"] == "starting_ligand")
    candidates = []
    for c in material["h1"]["compounds"]:
        if c["role"] != "known_reference_protac":
            continue
        cid = c["compound_id"]
        observations = [o for o in material["h1"]["observations"] if o["original"]["candidate_id"] == cid]
        reconstruction = next(r for r in material["h2"]["reference_reconstructions"] if r["compound_id"] == cid)
        prediction = prediction_view(cid, collection)
        review = ["약학 검수와 실제 G1/G2/G3 판단은 별도 기록 필요", "보충자료의 해당 화합물 합성 근거 미검수"]
        for obs in observations:
            o = obs["original"]
            if o["kind"] == "missing" or o.get("missing_reason"):
                review.append(f'{o["endpoint"]}: {o.get("missing_reason", "관측 미확인")}')
        if prediction["status"] == "not_run":
            review.append("Boltz 실제 예측 미실행: 예측 기반 구조 판단 대기")
        elif prediction["status"] == "needs_review":
            review.append("계산 누락·실패·출처·입체화학 또는 신뢰도 결측 확인 필요")
        if reconciliation["status"] == "not_provided":
            review.append("동료 A의 실제 문헌 출력 대조 대기")
        candidates.append({"compound": copy.deepcopy(c), "literature_observations": copy.deepcopy(observations),
                           "reference_reconstruction": copy.deepcopy(reconstruction),
                           "computed_properties": {"kind": "computed", "values": copy.deepcopy(c["identity"]["properties"]),
                                                   "source": "Existing verified CPU snapshot; not newly measured ADME or synthesis evidence."},
                           "prediction": prediction, "review_items": review})
    synthetic = collection["data_mode"] == "synthetic_test" or exchange is not None and exchange["data_mode"] == "synthetic_test"
    bundle = seal({"format": "tpd-candidate-evidence/0.1.0-draft", "case_id": material["case_id"],
                   "data_mode": "synthetic_test" if synthetic else "real", "material_digest": material["material_digest"],
                   "collection_digest": collection["digest"], "source_files": {
                       "material": reference("material.json", material, "urn:tpd-navigator:b-review-material:0.1.0-draft"),
                       "collection": reference("collection.json", collection, "urn:tpd-navigator:b-boltz-collection:0.1.0-draft"),
                       "literature_exchange": reference("literature-exchange.json", exchange, "urn:tpd-navigator:b-literature-exchange:0.1.0-draft") if exchange else None},
                   "starting_ligand": copy.deepcopy(start), "shared_hypothesis": copy.deepcopy(material["h2"]),
                   "candidates": candidates, "literature_reconciliation": reconciliation,
                   "authority": {"human_review": "pending", "approval_record_created": False, "dispatch_authorized": False,
                                 "public_release_ready": False, "efficacy_claim": "not_established"}})
    schema_check(bundle, "b_candidate_evidence.schema.json")
    return bundle


def validate_evidence(bundle, material, collection, exchange=None, *, allow_synthetic=False):
    schema_check(bundle, "b_candidate_evidence.schema.json")
    verify_seal(bundle)
    expected = build_evidence(material, collection, exchange, allow_synthetic=allow_synthetic)
    check(encoded(bundle) == encoded(expected), "EVIDENCE_PROJECTION_MISMATCH")
    return bundle


def markdown(bundle):
    """Render deterministic review text from already validated evidence; no new claims."""
    def text(value):
        if value is None:
            return "미확인"
        return html.escape(str(value), quote=False).replace("|", "\\|").replace("\n", " ")

    status = {"not_run": "미실행", "available_unreviewed": "자료 있음·검토 전", "needs_review": "누락/오류/출처 검토 필요"}
    lines = ["# SMARCA2 후보별 통합 근거표", "", "**내부 검토 자료 · 사람 승인 전 · 공개 결과 확정 전**", "",
             "실측 문헌 관측, 기존 CPU 계산, 보고된 모델 예측, 미확인 사항을 구분합니다. 합성 가능성·분해 효능을 자동 판정하지 않습니다.", "",
             f'자료 모드: `{text(bundle["data_mode"])}`. 자료 버전: `{bundle["digest"]}`.', "",
             "## 전체 상태", "", "| 후보 | 논문 화합물 | 분자/부착 원자 | 예측 | 요청 / 판독 / 누락 / 오류 |", "|---|---|---|---|---|"]
    for candidate in bundle["candidates"]:
        c, p = candidate["compound"], candidate["prediction"]
        n = p["counts"]
        lines.append(f'| {text(c["compound_id"])} | {text(c["paper"]["name"])} · 번호 {text(c["paper"]["compound_number"])} | {text(c["structure"]["ccd"])} / {text(c["attachment_atom"]["ccd_atom_id"])} | {status[p["status"]]} | {n["requested"]} / {n["parsed"]} / {n["missing"]} / {n["invalid"]} |')
    lines += ["", "요청/판독 수는 제공된 실행 기록의 개수입니다. 재시도·여러 샘플을 독립 실험으로 세거나 후보 총점으로 합산하지 않습니다.", ""]
    for candidate in bundle["candidates"]:
        c, p = candidate["compound"], candidate["prediction"]
        lines += [f'## {text(c["compound_id"])} — {text(c["paper"]["name"])}', "",
                  f'알려진 문헌 PROTAC의 재구성 자료. 분자 버전: `{text(c["molecule_id"])}`.', "",
                  f'출처 DOI: {text(c["paper"]["doi"])} · 실험 구조: {text(c["structure"]["pdb"])} / {text(c["structure"]["ccd"])}.', "",
                  f'정확한 입체화학 SMILES: `{text(c["identity"]["canonical_isomeric_smiles"])}`', "",
                  f'[구조 SDF](cpu/{c["files"]["sdf_2d"]["path"]}) · [2D 그림](cpu/{c["files"]["depiction_2d"]["path"]}) · [원자 대응](cpu/{c["files"]["atom_map"]["path"]})', "",
                  "### 문헌 관측", "", "| 항목 | 원래 값 | 조건·결측 | 원문 위치 |", "|---|---|---|---|"]
        for observation in candidate["literature_observations"]:
            o = observation["original"]
            value = "미추출" if o["kind"] == "missing" else f'{text(o["relation"])} {text(o["value"])} {text(o["unit"])}'
            conditions = f'세포 {text(o.get("cell_line"))}; 노출 시간 {text(o.get("exposure_time_h"))}; 분석 {text(o.get("assay"))}; 반복 {text(o.get("biological_replicates"))}. {text(o.get("missing_reason", ""))}'
            lines.append(f'| {text(o["endpoint"])} | {value} | {conditions} | {text(o["locator"])} |')
        lines += ["", "문헌값은 기존 수동 전사이며 약학 검수 전입니다. 원래 단위·부등호·조건과 결측을 유지합니다.", "",
                  "### 기존 CPU 물성 계산", "", "| 항목 | 값 |", "|---|---|"]
        lines += [f'| {text(k)} | {text(v)} |' for k, v in sorted(candidate["computed_properties"]["values"].items())]
        lines += ["", "실측 ADME나 합성 검증 결과가 아닙니다.", "", "### 구조 예측 검토", "", f'상태: **{status[p["status"]]}**. 자동 후보 선택 없음.', ""]
        if p["samples"]:
            lines += ["| 실행 / rank | 판독 | confidence | 표적 RMSD Å | 전체 리간드 RMSD 범위 Å | VHL 상대 RMSD Å | 비고 |", "|---|---|---|---|---|---|---|"]
            for m in p["samples"]:
                s = m["summary"] or {}
                lines.append(f'| {text(m["run_id"])} / {m["rank"]} | {text(m["status"])} | {text(s.get("confidence_score"))} | {text(s.get("target_rmsd_A"))} | {text(s.get("ligand_rmsd_range_A"))} | {text(s.get("vhl_rmsd_after_target_alignment_A"))} | {text(m["reason"] or ", ".join(s.get("warnings", [])))} |')
            lines.append("")
            for group in p["condition_groups"]:
                lines.append(f'- 조건 그룹 `{group["condition_digest"]}`: {text(", ".join(group["run_ids"]))}')
            for run_id, warnings in sorted(p["run_warnings"].items()):
                if warnings:
                    lines.append(f'- {text(run_id)}: {text(", ".join(warnings))}')
        else:
            lines += ["GPU 계산을 보류한 상태입니다. 예측 좌표·신뢰도·구조 비교값은 아직 없습니다."]
        lines += ["", "### 확인할 항목", ""] + ["- " + text(v) for v in candidate["review_items"]] + [""]
    h2 = bundle["shared_hypothesis"]
    lines += ["## 공통 출발 리간드의 부착점 근거", "",
              "아래 기하는 START 실험 구조에서 계산한 값입니다. 완성 PROTAC이나 Boltz 예측값으로 옮겨 적지 않습니다.", "",
              f'부착 원자: {text(h2["descriptive_geometry"]["attachment_ccd_atom"])}. SASA: {text(h2["descriptive_geometry"]["SASA_A2"])} Å².', "",
              text(h2["attachment_literature"]["summary"]), "", "## 동료 문헌 결과 대조", ""]
    reconciliation = bundle["literature_reconciliation"]
    if reconciliation["status"] == "not_provided":
        lines.append("동료 A의 실제 출력은 아직 제공되지 않았습니다. B 측 예제 검사를 실제 A/B 연결 완료로 표시하지 않습니다.")
    else:
        lines += [f'제공자: {text(reconciliation["producer"])}. 모든 대조는 사람 검토 전입니다.', ""]
        for row in reconciliation["records"]:
            result = row["comparison"]
            lines.append(f'- {text(row["record_id"])}: {text(result["status"])}; 차이: {text(", ".join(result["differences"]))}; 출처: {text(", ".join(row["source_flags"]))}')
    lines += ["", "원자료와 상세 결과는 [evidence.json](evidence.json), [material.json](material.json), [collection.json](collection.json)에 연결됩니다. 이 자료는 공식 G1/G2/G3 승인이나 최종 공개 보고서를 대신하지 않습니다.", ""]
    return "\n".join(lines).encode("utf-8")
