import copy
import re

from packages.contracts import known, missing, parse_json, payload
from packages.transport import FetchError
from packages.ports import Artifacts, SourceTransport

URL = "https://rest.uniprot.org/uniprotkb/search"


def safe_term(text):
    return re.sub(r'[^\w .\-]', ' ', text, flags=re.UNICODE).strip()


def resolve_target(query_ref: dict, artifacts: Artifacts, transport: SourceTransport) -> dict:
    query = artifacts.json(query_ref)
    term = safe_term(query["query"])
    if not term:
        raise ValueError("표적명 또는 식별자를 입력해 주세요.")
    field = {"uniprot_accession": "accession", "gene_symbol": "gene_exact", "name": "protein_name"}[query["query_kind"]]
    expression = f'({field}:"{term}") AND (organism_id:{query["taxon_id"]})'
    response = transport.get(URL, {"query": expression, "format": "json", "size": 11,
        "fields": "accession,protein_name,gene_names,organism_id,ft_transmem,xref_pdb"}, max_bytes=2_000_000)
    source = artifacts.put_raw(response.body)
    data = parse_json(response.body)
    if not isinstance(data.get("results"), list):
        raise FetchError("INVALID_PROVIDER_RESPONSE", "UniProt 응답 형식을 확인할 수 없습니다.")
    total = int(response.headers.get("x-total-results", len(data["results"])))
    candidates = []
    for row in data["results"]:
        accession = row["primaryAccession"]
        description = row.get("proteinDescription", {})
        protein = description.get("recommendedName", {}).get("fullName", {}).get("value")
        if not protein:
            protein = next((n.get("fullName", {}).get("value") for n in description.get("submissionNames", []) if n.get("fullName")), accession)
        gene = next((g.get("geneName", {}).get("value") for g in row.get("genes", []) if g.get("geneName")), None)
        identity = {"uniprot_accession": accession, "taxon_id": row["organism"]["taxonId"], "label": protein,
                    "gene_symbol": known(gene) if gene else missing(), "isoform": missing("isoform은 별도로 확인하지 않았습니다.", "unknown"),
                    "mutation": missing("변이 상태는 별도로 확인하지 않았습니다.", "unknown"), "construct": missing("사용할 construct는 정하지 않았습니다.", "unknown")}
        candidates.append({"candidate_id": accession, "identity": identity, "source_refs": [source], "match_reason": "UniProt의 요청 필드·종 검색 결과"})
    count = len(candidates)
    # Truncated broad queries still have at least two candidates; never take the first.
    state = "not_found" if count == 0 else "resolved" if count == 1 and total == 1 else "ambiguous"
    if state == "ambiguous" and count < 2:
        raise FetchError("INCOMPLETE_TARGET_RESPONSE", "전체 후보 수와 받은 후보가 달라 선택할 수 없습니다. 입력을 구체화해 주세요.")
    reasons = []
    if total > count:
        reasons.append(f"전체 {total}건 중 {count}건만 표시합니다. 더 정확한 표적명/접근번호로 다시 검색할 수 있습니다.")
    scope = "undetermined"
    selection = None
    chosen = None
    if state == "resolved":
        chosen = candidates[0]["candidate_id"]
        scope, scope_reasons = assess_scope(data["results"][0])
        reasons += scope_reasons
        selection = {"kind": "unique_match", "annotation_ref": None, "reason": "요청 범위에서 후보가 정확히 한 건입니다."}
    elif not reasons:
        reasons.append("표적을 식별한 뒤 제품 범위를 확인합니다.")
    return payload("TargetResolution", query_ref=query_ref, resolution_status=state, candidates=candidates,
                   selected_candidate_id=chosen, selection=selection, scope_status=scope, scope_reasons=reasons, issues=[])


def assess_scope(row):
    if row["organism"]["taxonId"] != 9606:
        return "outside_current_scope", ["현재 제품은 사람 표적을 대상으로 합니다."]
    if any(f.get("type", "").lower() == "transmembrane" for f in row.get("features", [])):
        return "outside_current_scope", ["UniProt에 막관통 구간이 주석되어 있어 현재 제품 범위에서 제외합니다."]
    return "within_current_scope", ["사람 표적이며 조회 주석에서 막관통 구간을 찾지 못했습니다. 실제 구조·모델 적용성 검증은 이후 단계입니다."]


def select_target(previous: dict, candidate_id: str, annotation_ref: dict, artifacts: Artifacts) -> dict:
    value = copy.deepcopy(previous)
    candidate = next((c for c in value["candidates"] if c["candidate_id"] == candidate_id), None)
    if candidate is None:
        raise ValueError("제시된 표적 후보 중에서 선택해 주세요.")
    data = artifacts.json(candidate["source_refs"][0])
    row = next(r for r in data["results"] if r["primaryAccession"] == candidate_id)
    scope, reasons = assess_scope(row)
    value.update(resolution_status="resolved", selected_candidate_id=candidate_id,
        selection={"kind": "human_selection", "annotation_ref": annotation_ref, "reason": "표시된 후보에서 사용자가 선택했습니다."},
        scope_status=scope, scope_reasons=reasons)
    return value
