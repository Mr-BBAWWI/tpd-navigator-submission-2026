"""Europe PMC + optional PDF -> I1 EvidenceBundle; no scientific claim extraction."""
import uuid
import re
from urllib.parse import quote

from packages.contracts import issue, known, missing, parse_json, payload, now, encoded
from packages.transport import FetchError
from packages.ports import Artifacts, SourceTransport
from packages.literature.extract import VERSION, extract_pdf, extract_xml, strip_markup

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
WARNING_TEXT = {
    "SEGMENT_LIMIT": "이번 실행의 추출 구간 상한에 도달했습니다.",
    "SEGMENT_TEXT_TRUNCATED": "긴 구간은 일부만 추출했습니다. 원문을 함께 확인해 주세요.",
    "FIGURES_NOT_INTERPRETED": "그림 설명은 추출했지만 그림 속 구조·수치는 해석하지 않았습니다.",
    "TABLE_LAYOUT_NOT_INTERPRETED": "표 텍스트는 추출했지만 행·열 대응은 검증하지 않았습니다.",
    "PDF_LAYOUT_AND_IMAGES_NOT_INTERPRETED": "PDF 텍스트만 추출했습니다. 도표·분자 구조·열 배치는 해석하지 않았습니다.",
    "PDF_PAGE_OR_TEXT_LIMIT": "PDF 페이지/구간 길이 한도로 일부만 추출했습니다.",
    "PDF_NO_EXTRACTABLE_TEXT": "텍스트를 추출하지 못했습니다. 스캔 PDF는 OCR이 필요할 수 있습니다.",
}


def questions_for(identity):
    labels = [identity["uniprot_accession"], identity["label"]]
    if identity["gene_symbol"]["state"] == "known":
        labels.insert(0, identity["gene_symbol"]["value"])
    labels = list(dict.fromkeys(label.replace('"', ' ').strip() for label in labels))
    terms = " OR ".join('"' + label + '"' for label in labels)
    return [{"question_id": "q-target-tpd", "question": "이 표적의 TPD·PROTAC 및 부착점 관련 문헌은 무엇인가?",
             "query": "(" + terms + ') AND (PROTAC OR "targeted protein degradation")', "purpose": "attachment_sar"}]


class Collector:
    def __init__(self, artifacts: Artifacts, transport: SourceTransport):
        self.artifacts, self.transport = artifacts, transport

    def collect(self, request_ref, progress=lambda message: None):
        request = self.artifacts.json(request_ref)
        policy = self.artifacts.json(request["collection_policy_ref"])
        self.remaining_bytes = policy["max_source_bytes"]
        self.remaining_segments = policy["max_segments"]
        documents, segments, notices = [], [], []
        for index, reference in enumerate(request["attached_source_refs"]):
            progress(f"첨부 PDF {index + 1}/{len(request['attached_source_refs'])} 확인 중")
            raw = self.artifacts.read(reference)
            self.remaining_bytes -= len(raw)
            document, extracted = self.document(raw, "pdf", known(f"첨부 PDF {index + 1}"), [],
                [{"kind": "upload", "provider": "researcher", "acquired_at": now(), "external_uri": None}], reference)
            documents.append(document)
            segments += extracted
        searches = []
        for question in request["questions"]:
            progress("Europe PMC 문헌 검색 중")
            slots = policy["max_documents"] - len(documents)
            if slots < 1:
                notice = issue("DOCUMENT_LIMIT", "첨부 자료가 문서 한도를 사용하여 검색을 실행하지 못했습니다.")
                searches.append({"question_id": question["question_id"], "query": question["query"], "provider": "europe_pmc", "status": "not_run",
                    "total_hits": {"state": "unknown", "value": None, "reason": "검색 미실행"}, "search_result_ref": None,
                    "retrieved_document_ids": [], "termination": "not_run", "issues": [notice]})
                notices.append(notice)
                continue
            try:
                response = self.transport.get(BASE + "/search", {"query": question["query"], "format": "json", "resultType": "core", "pageSize": min(slots, 12)}, max_bytes=2_000_000)
                search_ref = self.artifacts.put_raw(response.body)
                data = parse_json(response.body)
                records = data.get("resultList", {}).get("result")
                if not isinstance(records, list) or not isinstance(data.get("hitCount"), int):
                    raise FetchError("INVALID_PROVIDER_RESPONSE", "문헌 검색 응답을 해석하지 못했습니다.")
                found = []
                search_issues = []
                termination = "record_limit" if data["hitCount"] > len(records[:slots]) else "source_exhausted"
                seen = set()
                for index, row in enumerate(records[:slots]):
                    dedup = (row.get("source"), row.get("id"))
                    if dedup in seen:
                        continue
                    seen.add(dedup)
                    progress(f"검색 문헌 {index + 1}/{len(records[:slots])} 원문 확인 중")
                    try:
                        document, extracted = self.from_record(row)
                    except FetchError as exc:
                        notice = issue(exc.code, str(exc))
                        notices.append(notice)
                        search_issues.append(notice)
                        termination = "byte_limit" if exc.code == "BYTE_LIMIT" else "time_limit" if exc.code == "TIME_LIMIT" else "record_limit"
                        break
                    documents.append(document)
                    segments += extracted
                    found.append(document["document_id"])
                searches.append({"question_id": question["question_id"], "query": question["query"], "provider": "europe_pmc", "status": "completed",
                    "total_hits": {"state": "known", "value": data["hitCount"], "reason": None}, "search_result_ref": search_ref,
                    "retrieved_document_ids": found, "termination": termination, "issues": search_issues})
            except FetchError as exc:
                notice = issue(exc.code, str(exc), "transient", retry="same_input_may_retry")
                searches.append({"question_id": question["question_id"], "query": question["query"], "provider": "europe_pmc", "status": "failed",
                    "total_hits": {"state": "unknown", "value": None, "reason": "검색 응답을 확보하지 못했습니다."}, "search_result_ref": None,
                    "retrieved_document_ids": [], "termination": "time_limit" if exc.code == "TIME_LIMIT" else "network_error" if exc.code == "NETWORK_ERROR" else "provider_error", "issues": [notice]})
                notices.append(notice)
        notices += [item for d in documents for item in d["issues"] + d["extraction"]["issues"]]
        coverage = {"stored_document_count": len(documents), "fulltext_acquired_count": sum(d["access_level"] == "fulltext" for d in documents),
            "extracted_segment_count": len(segments), "llm_read_segment_count": 0, "evidence_record_count": 0,
            "reviewed_record_count": 0, "unavailable_document_count": sum(d["acquisition_status"] != "retrieved" for d in documents),
            "independent_experiment_count": {"state": "known", "value": 0, "reason": None}}
        return payload("EvidenceBundle", project_id=self.artifacts.project, request_ref=request_ref,
            target_resolution_ref=request["target_resolution_ref"], documents=documents, segments=segments, compounds=[],
            evidence_records=[], searches=searches, coverage=coverage, issues=notices)

    def from_record(self, row):
        identifiers = [{"namespace": scheme, "value": str(row[key])} for key, scheme in [("doi", "doi"), ("pmid", "pmid"), ("pmcid", "pmcid")] if row.get(key)]
        if row.get("source") == "MED" and not row.get("pmid") and row.get("id"):
            identifiers.append({"namespace": "pmid", "value": str(row["id"])})
        title = known(strip_markup(row["title"])) if row.get("title") else missing()
        uri = "https://europepmc.org/article/" + quote(str(row.get("source", "MED")), safe="") + "/" + quote(str(row.get("id", "")), safe="")
        paths = [{"kind": "search", "provider": "europe_pmc", "acquired_at": now(), "external_uri": uri}]
        notices = []
        if row.get("isOpenAccess") == "Y" and re.fullmatch(r"PMC\d+", str(row.get("pmcid", ""))) and self.remaining_segments:
            try:
                response = self.transport.get(BASE + "/" + row["pmcid"] + "/fullTextXML", max_bytes=max(1, self.remaining_bytes - 20000))
                document, segments = self.document(response.body, "jats_xml", title, identifiers, paths)
                if document["extraction"]["status"] in {"complete", "partial"}:
                    self.remaining_bytes -= len(response.body)
                    return document, segments
                notices += document["extraction"]["issues"]
            except FetchError as exc:
                notices.append(issue(exc.code, str(exc), "transient"))
        if row.get("abstractText"):
            raw = row["abstractText"].encode("utf-8")
            kind = "abstract_text"
        else:
            raw = encoded(row)
            kind = "metadata"
        if len(raw) > self.remaining_bytes:
            # A metadata stub still preserves the acquisition failure without inventing text.
            raw = encoded({"id": row.get("id"), "source": row.get("source"), "title": row.get("title", "")[:300]})
            kind = "metadata"
            notices.append(issue("SOURCE_BYTE_LIMIT", "원문 바이트 한도로 메타데이터 일부만 보존했습니다."))
        if len(raw) > self.remaining_bytes:
            raise FetchError("BYTE_LIMIT", "메타데이터를 저장할 원문 용량도 남지 않았습니다.")
        document, segments = self.document(raw, kind, title, identifiers, paths)
        document["issues"] += notices
        self.remaining_bytes -= len(raw)
        return document, segments

    def document(self, raw, kind, title, identifiers, paths, reference=None):
        media = {"pdf": "application/pdf", "jats_xml": "application/xml", "abstract_text": "text/plain", "metadata": "application/json"}[kind]
        source = reference or self.artifacts.put_raw(raw, media)
        document_id = "d-" + uuid.uuid4().hex
        issues, chunks, warnings = [], [], []
        try:
            if self.remaining_segments <= 0:
                chunks, warnings = [], ["SEGMENT_LIMIT"]
            elif kind == "jats_xml":
                chunks, warnings = extract_xml(raw, self.remaining_segments)
            elif kind == "pdf":
                chunks, warnings = extract_pdf(raw, self.remaining_segments)
            elif kind == "abstract_text":
                chunks = [{"text": strip_markup(raw.decode("utf-8")), "content_kind": "abstract"}]
            else:
                metadata = parse_json(raw)
                chunks = [{"text": "제목: " + str(metadata.get("title", "미제공")), "content_kind": "metadata"}]
            extraction_state = "partial" if warnings else "complete"
            for code in warnings:
                issues.append(issue(code, WARNING_TEXT.get(code, code), ids=[source["artifact_id"]]))
        except Exception:
            extraction_state = "failed"
            issues.append(issue("SOURCE_EXTRACTION_FAILED", "내용을 추출하지 못했습니다. 암호화·손상·형식 또는 처리 시간 한도를 확인해 주세요.", ids=[source["artifact_id"]], retry="requires_changed_input"))
        raw_segments = []
        for chunk in chunks[:self.remaining_segments]:
            if not chunk["text"].strip():
                continue
            identifier = "s-" + uuid.uuid4().hex
            text_ref = self.artifacts.put_raw(chunk["text"].encode("utf-8"), "text/plain", "computed")
            raw_segments.append({"segment_id": identifier, "text_ref": text_ref, "hint": {k: v for k, v in chunk.items() if k != "text"}})
        snapshot = None
        if extraction_state != "failed":
            snapshot = self.artifacts.put_json({"extractor_version": VERSION, "source_ref": source, "segments": raw_segments, "warnings": warnings})
        self.remaining_segments -= len(raw_segments)
        segments = []
        for segment in raw_segments:
            hint = segment["hint"]
            locator = {"document_id": document_id, "document_revision": 1, "source_artifact_ref": source,
                       "extraction_snapshot_ref": snapshot, "extractor_version": VERSION, "kind": kind}
            if kind == "jats_xml":
                locator.update(block_id=segment["segment_id"], xpath=hint["xpath"], table_id=hint["table_id"], figure_id=hint["figure_id"])
            elif kind == "pdf":
                locator.update(page=hint["page"], segment_id=segment["segment_id"], bbox=None, table_id=None, figure_id=None)
            elif kind == "abstract_text":
                locator.update(record_id=document_id, field="abstract", segment_id=segment["segment_id"])
            else:
                locator.update(record_id=document_id, field="title")
            segments.append({"segment_id": segment["segment_id"], "locator": locator, "text_artifact_ref": segment["text_ref"],
                             "content_kind": hint["content_kind"], "read_events": []})
        level = "unverified" if kind in {"pdf", "jats_xml"} and (extraction_state == "failed" or not chunks) else "fulltext" if kind in {"pdf", "jats_xml"} else "abstract" if kind == "abstract_text" else "metadata"
        document = {"document_id": document_id, "revision": 1, "project_id": self.artifacts.project, "title": title,
            "identifiers": identifiers, "acquisition_paths": paths, "acquisition_status": "retrieved", "access_level": level,
            "format": kind, "source_artifact_ref": source,
            "extraction": {"status": extraction_state, "extractor_version": VERSION, "snapshot_ref": snapshot, "issues": issues},
            "relations": [], "replay_rights": {"status": "undetermined", "basis_ref": None, "reason": "공개 재생을 위한 원문 권리는 별도 검토 전입니다."}, "issues": []}
        return document, segments
