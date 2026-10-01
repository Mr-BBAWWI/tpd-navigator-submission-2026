"""Source-bound scientific assessment, trusted computation, and human decisions."""
from __future__ import annotations

import copy
import hashlib
import io
import json
import math
import uuid
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator
from rdkit import Chem

from packages.contracts import ContractError, encoded, now, parse_json, validate_payload
from packages.platform.design_panel import input_binding
from packages.platform.dossiers import require
from packages.platform.review_identity import ReviewIdentityService
from packages.platform.saved_results import SavedResultsService
from packages.platform import protocol_reassessment
from packages.science import calibration_import, expert_followup, expert_opinion, interaction_review, novel_ternary, scientific_assessment, synthesis_review
from packages.science.protein_hydrogen_evidence import ProteinHydrogenEvidence

VERSION = "scientific-acceptance-service/20261001.3"
ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = ROOT / "contracts" / "drafts" / "scientific_acceptance.schema.json"
MAX_RESULT_BYTES = 32 * 1024 * 1024
MAX_STRUCTURED_STATEMENT_BYTES = 50000
MAX_MICROSTATE_SOURCE_BYTES = 1024 * 1024
MAX_MICROSTATE_SOURCE_ROWS = 256
POLICY_KINDS = {
    "numerical_criteria",
    "parent_funnel",
    "site_policy",
    "interaction_requirements",
    "calibration_criterion",
    "microstate_decision",
    "synthesis_policy",
    "ternary_criterion",
    "protein_hydrogen_review",
    "followup_geometry_review",
}
ENGINE_FILES = (
    "packages/platform/scientific_acceptance.py",
    "packages/platform/computed_ligand_pka.py",
    "packages/platform/review_identity.py",
    "packages/platform/protocol_reassessment.py",
    "packages/science/interaction_review.py",
    "packages/science/scientific_assessment.py",
    "packages/science/benchmark_distribution.py",
    "packages/science/boltz_worker.py",
    "packages/science/calibration_import.py",
    "packages/science/synthesis_review.py",
    "packages/science/novel_ternary.py",
    "packages/science/protein_hydrogen_evidence.py",
    "packages/science/expert_opinion.py",
    "packages/science/expert_followup.py",
    "cases/expert_opinions/20260930/opinion.json",
    "cases/expert_opinions/20260930/reply.docx",
    "cases/expert_opinions/20261001/reply.docx",
    "contracts/drafts/scientific_acceptance.schema.json",
)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _digest(value: Any) -> str:
    return _sha(encoded(value))


def _copy(value: Any) -> Any:
    return json.loads(json.dumps(value, ensure_ascii=False, allow_nan=False))


def _ref(value: Any) -> bool:
    return isinstance(value, dict) and {
        "artifact_id", "version", "sha256", "media_type", "schema_id", "provenance"
    } <= set(value)


def engine_fingerprint() -> str:
    digest = hashlib.sha256()
    for relative in ENGINE_FILES:
        path = ROOT / relative
        if not path.is_file() or path.is_symlink():
            raise ContractError("SCIENTIFIC_IMPLEMENTATION_FILE_MISSING")
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


class ScientificAcceptanceService:
    def __init__(self, store, project):
        self.store, self.project, self.port = store, project, store.scope(project)
        self.auth = ReviewIdentityService(store, project)
        self.saved = SavedResultsService(store, project)
        schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
        self.validator = Draft202012Validator(schema)
        with store.db() as db:
            db.executescript("""
                CREATE TABLE IF NOT EXISTS scientific_assessments(
                    id TEXT PRIMARY KEY,
                    project TEXT NOT NULL,
                    job_id TEXT NOT NULL,
                    revision INTEGER NOT NULL,
                    policy_revision INTEGER NOT NULL,
                    policy_digest TEXT NOT NULL,
                    fingerprint TEXT NOT NULL,
                    ref TEXT NOT NULL,
                    UNIQUE(project,job_id,revision),
                    UNIQUE(project,job_id,fingerprint));
                CREATE TABLE IF NOT EXISTS scientific_decisions(
                    id TEXT PRIMARY KEY,
                    project TEXT NOT NULL,
                    job_id TEXT NOT NULL,
                    assessment_id TEXT NOT NULL,
                    policy_revision INTEGER NOT NULL,
                    kind TEXT NOT NULL,
                    actor_id TEXT NOT NULL,
                    withdrawn INTEGER NOT NULL,
                    fingerprint TEXT NOT NULL,
                    body TEXT NOT NULL,
                    UNIQUE(project,job_id,fingerprint));
                CREATE TABLE IF NOT EXISTS scientific_policy_events(
                    project TEXT NOT NULL,
                    job_id TEXT NOT NULL,
                    revision INTEGER NOT NULL,
                    event_kind TEXT NOT NULL,
                    decision_id TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    PRIMARY KEY(project,job_id,revision));
                CREATE TABLE IF NOT EXISTS scientific_evidence(
                    id TEXT PRIMARY KEY,
                    project TEXT NOT NULL,
                    job_id TEXT NOT NULL,
                    kind TEXT NOT NULL,
                    identity_key TEXT NOT NULL,
                    fingerprint TEXT NOT NULL,
                    binding TEXT NOT NULL,
                    ref TEXT NOT NULL,
                    superseded INTEGER NOT NULL,
                    created_at TEXT NOT NULL,
                    UNIQUE(project,job_id,fingerprint));
                CREATE TABLE IF NOT EXISTS scientific_statements(
                    id TEXT PRIMARY KEY,
                    project TEXT NOT NULL,
                    actor_id TEXT NOT NULL,
                    ref TEXT NOT NULL,
                    created_at TEXT NOT NULL);
            """)

    def _validate(self, value: Any, definition: str) -> None:
        wrapper = {"kind": definition, "payload": value}
        errors = sorted(self.validator.iter_errors(wrapper), key=lambda error: list(error.path))
        if errors:
            path = ".".join(str(item) for item in errors[0].absolute_path)
            raise ContractError("SCIENTIFIC_SCHEMA_INVALID" + (":" + path if path else ""))

    def _authenticate(self, token, db):
        return self.auth.authenticate(token, db)

    def _job(self, db, job_id: str) -> dict:
        row = db.execute(
            "SELECT project,state,body FROM workbench_jobs WHERE id=? AND project=?",
            (job_id, self.project),
        ).fetchone()
        if row is None:
            raise KeyError(job_id)
        body = json.loads(row["body"])
        require(row["project"] == self.project and body.get("project_id") == self.project,
                "SCIENTIFIC_JOB_PROJECT")
        require(row["state"] == "completed" and body.get("state") == "completed",
                "SCIENTIFIC_JOB_NOT_COMPLETED")
        require(body.get("operation") == "design_panel", "SCIENTIFIC_JOB_KIND")
        require(isinstance(body.get("result"), dict), "SCIENTIFIC_JOB_RESULT")
        return body

    def _source_state(self, job: dict) -> dict:
        current_binding = input_binding()
        inputs_current = current_binding.get("digest") == job.get("input_digest")
        try:
            from packages.platform.workbench import identity
            runtime_current = job.get("runtime") == identity()
        except Exception:
            runtime_current = False
        return {
            "source_inputs_current": inputs_current,
            "source_inputs_reason": None if inputs_current else "DESIGN_SOURCE_INPUT_DIGEST_CHANGED",
            "source_runtime_current": runtime_current,
            "source_runtime_reason": None if runtime_current else "DESIGN_JOB_RUNTIME_CHANGED",
            "current_input_binding": current_binding,
        }

    def _verified_result(self, db, job_id: str, *, require_current: bool = False):
        job = self._job(db, job_id)
        result = job["result"]
        require(len(encoded(result)) <= MAX_RESULT_BYTES, "SCIENTIFIC_RESULT_TOO_LARGE")
        require(result.get("format") == "design-panel/20260930.4", "SCIENTIFIC_RESULT_FORMAT")
        result_binding = result.get("input_binding")
        require(isinstance(result_binding, dict), "SCIENTIFIC_RESULT_BINDING")
        require(result_binding.get("digest") == job.get("input_digest"),
                "SCIENTIFIC_RESULT_JOB_INPUT_MISMATCH")
        require(job.get("binding", {}).get("digest") == job.get("input_digest"),
                "SCIENTIFIC_JOB_BINDING_MISMATCH")

        files = result.get("files")
        require(isinstance(files, dict) and _ref(files.get("json")) and _ref(files.get("report")),
                "SCIENTIFIC_RESULT_FILES")
        outputs = job.get("outputs")
        require(isinstance(outputs, list), "SCIENTIFIC_JOB_OUTPUTS")
        output_refs = [item.get("ref") for item in outputs if isinstance(item, dict)]
        require(sum(ref == files["json"] for ref in output_refs) == 1,
                "SCIENTIFIC_JSON_OUTPUT_NOT_EXACT")
        require(sum(ref == files["report"] for ref in output_refs) == 1,
                "SCIENTIFIC_REPORT_OUTPUT_NOT_EXACT")

        json_raw = self._read_registered_ref(db, files["json"])
        self._read_registered_ref(db, files["report"])
        require(len(json_raw) <= MAX_RESULT_BYTES, "SCIENTIFIC_RESULT_TOO_LARGE")
        archived = parse_json(json_raw)
        expected = copy.deepcopy(result)
        expected["files"].pop("json", None)
        expected["files"].pop("report", None)
        require(encoded(archived) == encoded(expected), "SCIENTIFIC_DESIGN_JSON_MISMATCH")
        require(files["json"]["sha256"] == _sha(json_raw), "SCIENTIFIC_DESIGN_JSON_HASH")
        parent = result.get("parent_scope", {}).get("actual_design_parent_id")
        require(isinstance(parent, str) and parent, "SCIENTIFIC_PARENT_SCOPE")

        source = self._source_state(job)
        if require_current:
            require(source["source_inputs_current"], "SCIENTIFIC_DESIGN_SOURCE_CHANGED")
            require(source["source_runtime_current"], "SCIENTIFIC_DESIGN_RUNTIME_CHANGED")
        binding = {
            "project": self.project,
            "job_id": job_id,
            "input_sha256": job["input_digest"],
            "result_sha256": files["json"]["sha256"],
            "result_binding_kind": "exact_archived_design_json_bytes",
            "parent_id": parent,
        }
        return job, result, binding, source

    def _verified_reinspection_snapshot(
            self, db, job_id: str, *,
            allow_archived_runtime: bool = False):
        job, result, binding, source = self._verified_result(
            db, job_id, require_current=False
        )
        require(source["source_inputs_current"],
                "SCIENTIFIC_DESIGN_SOURCE_CHANGED")
        if not allow_archived_runtime:
            require(source["source_runtime_current"],
                    "SCIENTIFIC_DESIGN_RUNTIME_CHANGED")
        return job, result, binding, source

    def _reinspection_report_provenance(self, source: dict,
                                        actor: dict) -> dict:
        runtime_current = source["source_runtime_current"]
        return {
            "kind": "authenticated_immutable_archived_result_reinspection",
            "operator_id": actor["id"],
            "source_inputs_current": source["source_inputs_current"],
            "source_runtime_current": runtime_current,
            "historical_runtime_reinspection": not runtime_current,
            "note": (
                "reinspection of immutable archived result; "
                "source_runtime_current=false"
                if not runtime_current else
                "reinspection of immutable archived result under current runtime"
            ),
            "scientific_approval": False,
            "formal_acceptance": False,
        }

    def _write(self, db, paths, value, media="application/json", provenance="computed"):
        if media == "application/json" and isinstance(value, bytes):
            parse_json(value)
            raw = value
        else:
            if media == "application/json":
                self._verify_refs(db, value)
            raw = encoded(value) if media == "application/json" else value
        require(isinstance(raw, bytes), "SCIENTIFIC_ARTIFACT_BYTES")
        return self.saved._write(db, paths, raw, media, provenance)

    def _refs(self, value: Any):
        if _ref(value):
            yield value
        elif isinstance(value, dict):
            for item in value.values():
                yield from self._refs(item)
        elif isinstance(value, list):
            for item in value:
                yield from self._refs(item)

    def _read_registered_ref(self, db, ref: dict) -> bytes:
        validate_payload("ArtifactRef", ref)
        identifier = ref["artifact_id"]
        row = db.execute(
            "SELECT id,project,metadata FROM artifacts WHERE id=?",
            (identifier,),
        ).fetchone()
        require(row is not None, "SCIENTIFIC_ARTIFACT_NOT_REGISTERED")
        require(row["id"] == identifier, "SCIENTIFIC_ARTIFACT_ID")
        require(row["project"] == self.project, "SCIENTIFIC_ARTIFACT_PROJECT")
        try:
            registered = json.loads(row["metadata"])
        except (TypeError, ValueError) as error:
            raise ContractError("SCIENTIFIC_ARTIFACT_METADATA") from error
        validate_payload("ArtifactRef", registered)
        require(registered == ref, "SCIENTIFIC_ARTIFACT_METADATA")

        relative = Path(identifier)
        require(
            not relative.is_absolute()
            and len(relative.parts) == 1
            and relative.name == identifier
            and identifier not in {".", ".."}
            and "/" not in identifier
            and "\\" not in identifier,
            "SCIENTIFIC_ARTIFACT_PATH",
        )
        blob_directory = self.store.root / "blobs"
        require(blob_directory.is_dir() and not blob_directory.is_symlink(),
                "SCIENTIFIC_ARTIFACT_PATH")
        try:
            strict_root = blob_directory.resolve(strict=True)
            candidate = blob_directory / identifier
            require(not candidate.is_symlink(), "SCIENTIFIC_ARTIFACT_PATH")
            strict_path = candidate.resolve(strict=True)
            require(strict_path.parent == strict_root and strict_path.is_file()
                    and not strict_path.is_symlink(),
                    "SCIENTIFIC_ARTIFACT_PATH")
            raw = strict_path.read_bytes()
        except ContractError:
            raise
        except OSError as error:
            raise ContractError("SCIENTIFIC_ARTIFACT_PATH") from error
        require(_sha(raw) == registered["sha256"], "SCIENTIFIC_ARTIFACT_HASH")
        return raw

    def _verify_refs(self, db, value: Any, actor_id: str | None = None) -> None:
        for ref in self._refs(value):
            raw = self._read_registered_ref(db, ref)
            try:
                document = parse_json(raw)
            except Exception:
                document = None
            if isinstance(document, dict) and document.get("type") == "human_statement":
                row = db.execute(
                    "SELECT actor_id FROM scientific_statements WHERE project=? AND ref=?",
                    (self.project, json.dumps(ref, sort_keys=True)),
                ).fetchone()
                require(row is not None, "SCIENTIFIC_UNREGISTERED_HUMAN_STATEMENT")
                if actor_id is not None:
                    require(row["actor_id"] == actor_id, "SCIENTIFIC_STATEMENT_AUTHOR")

    def _require_human_statement_ref(self, db, ref: Any,
                                     actor_id: str | None = None) -> dict:
        require(_ref(ref), "SCIENTIFIC_SITE_SOURCE_REF")
        self._verify_refs(db, ref, actor_id)
        document = parse_json(self._read_registered_ref(db, ref))
        require(isinstance(document, dict) and document.get("type") == "human_statement",
                "SCIENTIFIC_SITE_HUMAN_STATEMENT_REQUIRED")
        row = db.execute(
            "SELECT actor_id FROM scientific_statements WHERE project=? AND ref=?",
            (self.project, json.dumps(ref, sort_keys=True)),
        ).fetchone()
        require(row is not None, "SCIENTIFIC_UNREGISTERED_HUMAN_STATEMENT")
        if actor_id is not None:
            require(row["actor_id"] == actor_id, "SCIENTIFIC_STATEMENT_AUTHOR")
            active = db.execute(
                "SELECT active FROM review_actors WHERE project=? AND id=?",
                (self.project, actor_id),
            ).fetchone()
            require(active is not None and active["active"] == 1,
                    "SCIENTIFIC_REVIEWER_INACTIVE")
        return document

    def _strict_quantitative_document(self, raw: bytes) -> dict:
        def pairs(items):
            value = {}
            for key, item in items:
                require(key not in value,
                        "SCIENTIFIC_MICROSTATE_SOURCE_DUPLICATE_KEY")
                value[key] = item
            return value

        def nonfinite(_value):
            raise ContractError("SCIENTIFIC_MICROSTATE_SOURCE_NONFINITE")

        def finite_float(value):
            try:
                number = float(value)
            except (ValueError, OverflowError) as error:
                raise ContractError("SCIENTIFIC_MICROSTATE_SOURCE_NONFINITE") from error
            require(math.isfinite(number),
                    "SCIENTIFIC_MICROSTATE_SOURCE_NONFINITE")
            return number

        def finite_int(value):
            try:
                number = int(value)
                converted = float(number)
            except (ValueError, OverflowError) as error:
                raise ContractError("SCIENTIFIC_MICROSTATE_SOURCE_NONFINITE") from error
            require(math.isfinite(converted),
                    "SCIENTIFIC_MICROSTATE_SOURCE_NONFINITE")
            return number

        try:
            value = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs,
                               parse_constant=nonfinite, parse_float=finite_float,
                               parse_int=finite_int)
        except ContractError:
            raise
        except (UnicodeError, TypeError, ValueError, json.JSONDecodeError) as error:
            raise ContractError("SCIENTIFIC_MICROSTATE_SOURCE_JSON") from error
        require(isinstance(value, dict), "SCIENTIFIC_MICROSTATE_SOURCE_OBJECT")
        return value

    def _validate_quantitative_source_record(
            self, db, job_id: str, result: dict, record: Any) -> None:
        require(isinstance(record, dict),
                "SCIENTIFIC_MICROSTATE_SOURCE_RECORD")
        result_sha256 = result.get("files", {}).get("json", {}).get("sha256")
        require(isinstance(result_sha256, str) and result_sha256,
                "SCIENTIFIC_MICROSTATE_RESULT_BINDING")
        require(record.get("project_id") == self.project
                and record.get("job_id") == job_id
                and record.get("result_sha256") == result_sha256,
                "SCIENTIFIC_MICROSTATE_SOURCE_SCOPE")

        analog_id = record.get("analog_id")
        analog = self._analog_map(result).get(analog_id)
        require(analog is not None, "SCIENTIFIC_MICROSTATE_SCOPE")
        state_index = record.get("selected_state_index")
        require(type(state_index) is int and state_index >= 0,
                "SCIENTIFIC_MICROSTATE_STATE_INDEX")
        analog_graph = scientific_assessment._isomeric_record(analog)
        require(analog_graph is not None, "SCIENTIFIC_MICROSTATE_ANALOG_GRAPH")
        require(scientific_assessment._canonical_isomeric(
                    record.get("analog_canonical_isomeric_graph")) == analog_graph,
                "SCIENTIFIC_MICROSTATE_SOURCE_ANALOG_GRAPH")
        state_graph = scientific_assessment._canonical_isomeric(
            record.get("state_canonical_isomeric_graph"))
        require(state_graph is not None,
                "SCIENTIFIC_MICROSTATE_SOURCE_STATE_GRAPH")
        source_pH = record.get("pH")
        require(type(source_pH) in (int, float)
                and math.isfinite(float(source_pH))
                and 0.0 <= float(source_pH) <= 14.0,
                "SCIENTIFIC_MICROSTATE_SOURCE_PH")

        interaction_ref = record.get("interaction_ref")
        require(_ref(interaction_ref),
                "SCIENTIFIC_MICROSTATE_SOURCE_INTERACTION_REF")
        source_match = next((item for item in self._interaction_evidence_rows(
            db, job_id, result_sha256, current_only=False)
            if item[1] == interaction_ref), None)
        require(source_match is not None,
                "SCIENTIFIC_MICROSTATE_SOURCE_INTERACTION_REF")
        report = source_match[2]
        require(report.get("analog_id") == analog_id,
                "SCIENTIFIC_MICROSTATE_REPORT")
        report_id = report.get("report_id") or report.get("job_id")
        require(isinstance(report_id, str) and report_id,
                "SCIENTIFIC_MICROSTATE_REPORT")
        declared_report_id = record.get("report_id")
        if declared_report_id is not None:
            require(declared_report_id == report_id,
                    "SCIENTIFIC_MICROSTATE_REPORT")
        self._microstate_report_state(
            report, analog_id, report_id, state_index, state_graph,
            float(source_pH))

        current_match = next((item for item in self._interaction_evidence_rows(
            db, job_id, result_sha256, current_only=True)
            if item[2].get("analog_id") == analog_id), None)
        require(current_match is not None,
                "SCIENTIFIC_MICROSTATE_REPORT")
        current_report = current_match[2]
        current_report_id = (current_report.get("report_id")
                             or current_report.get("job_id"))
        require(current_report_id == report_id,
                "SCIENTIFIC_MICROSTATE_REPORT")
        self._microstate_report_state(
            current_report, analog_id, report_id, state_index, state_graph,
            float(source_pH))

        required_text = ("method", "protocol", "conditions",
                         "uncertainty_description", "interpretation")
        require(all(isinstance(record.get(key), str) and record[key].strip()
                    for key in required_text),
                "SCIENTIFIC_MICROSTATE_SOURCE_METADATA")
        source_id = record.get("source_identifier")
        source_url = record.get("source_url")
        require((isinstance(source_id, str) and source_id.strip()) or
                (isinstance(source_url, str) and source_url.strip()),
                "SCIENTIFIC_MICROSTATE_SOURCE_IDENTIFIER")
        require(record.get("source_kind") in {"measured", "computed"},
                "SCIENTIFIC_MICROSTATE_SOURCE_KIND")
        has_pka = "pKa" in record
        has_population = "state_population" in record
        require(has_pka != has_population,
                "SCIENTIFIC_MICROSTATE_SOURCE_QUANTITY")
        quantity = record["pKa" if has_pka else "state_population"]
        require(type(quantity) in (int, float)
                and math.isfinite(float(quantity)),
                "SCIENTIFIC_MICROSTATE_SOURCE_NONFINITE")
        if has_population:
            require(0.0 <= float(quantity) <= 1.0,
                    "SCIENTIFIC_MICROSTATE_POPULATION_RANGE")

    def _validate_quantitative_source_document(
            self, db, job_id: str, result: dict, document: dict) -> None:
        require(document.get("type") != "human_statement",
                "SCIENTIFIC_MICROSTATE_QUANTITATIVE_SOURCE_REQUIRED")

        def approved(value: Any) -> bool:
            if isinstance(value, dict):
                return value.get("scientific_approved") is True or any(
                    approved(item) for item in value.values())
            if isinstance(value, list):
                return any(approved(item) for item in value)
            return False

        require(not approved(document),
                "SCIENTIFIC_MICROSTATE_SOURCE_APPROVAL_FORBIDDEN")
        rows = document.get("records")
        require(isinstance(rows, list)
                and 1 <= len(rows) <= MAX_MICROSTATE_SOURCE_ROWS,
                "SCIENTIFIC_MICROSTATE_SOURCE_ROWS")
        for record in rows:
            self._validate_quantitative_source_record(
                db, job_id, result, record)

    def _registered_population_source(self, db, job_id: str,
                                      result: dict, ref: Any) -> dict:
        require(_ref(ref), "SCIENTIFIC_MICROSTATE_SOURCE_REF")
        result_sha256 = result.get("files", {}).get("json", {}).get("sha256")
        rows = db.execute("""SELECT id,binding,ref,created_at
            FROM scientific_evidence WHERE project=? AND job_id=?
            AND kind='microstate_population_source' AND superseded=0
            ORDER BY rowid DESC""", (self.project, job_id)).fetchall()
        matches = []
        for row in rows:
            try:
                stored_ref = json.loads(row["ref"])
                binding = json.loads(row["binding"])
            except (TypeError, ValueError, json.JSONDecodeError):
                continue
            if stored_ref != ref or binding.get("result_sha256") != result_sha256:
                continue
            try:
                raw = self._read_registered_ref(db, stored_ref)
                document = self._strict_quantitative_document(raw)
                digest = _digest(document)
                require(digest == binding.get("evidence_sha256")
                        and digest == binding.get("source_document_sha256"),
                        "SCIENTIFIC_MICROSTATE_SOURCE_DIGEST")
                self._validate_quantitative_source_document(
                    db, job_id, result, document)
            except (ContractError, TypeError, ValueError, UnicodeError,
                    json.JSONDecodeError):
                continue
            matches.append({"id": row["id"], "binding": binding,
                            "ref": stored_ref, "created_at": row["created_at"],
                            "document": document})
        require(len(matches) == 1,
                "SCIENTIFIC_MICROSTATE_SOURCE_NOT_REGISTERED")
        return matches[0]

    def register_quantitative_source(self, job_id: str, text: str,
                                     session_token: str) -> dict:
        require(isinstance(text, str), "SCIENTIFIC_MICROSTATE_SOURCE_TEXT")
        try:
            raw = text.encode("utf-8", errors="strict")
        except UnicodeError as error:
            raise ContractError("SCIENTIFIC_MICROSTATE_SOURCE_JSON") from error
        require(len(raw) <= MAX_MICROSTATE_SOURCE_BYTES,
                "SCIENTIFIC_MICROSTATE_SOURCE_SIZE")
        document = self._strict_quantitative_document(raw)
        paths = []
        try:
            with self.store.db() as db:
                db.execute("BEGIN IMMEDIATE")
                actor = self._authenticate(session_token, db)
                _, result, binding, _ = self._verified_result(
                    db, job_id, require_current=True)
                self._validate_quantitative_source_document(
                    db, job_id, result, document)
                canonical_sha256 = _digest(document)
                rows = db.execute("""SELECT ref FROM scientific_evidence
                    WHERE project=? AND job_id=?
                    AND kind='microstate_population_source'
                    AND identity_key=? AND superseded=0
                    ORDER BY rowid DESC""",
                    (self.project, job_id, canonical_sha256)).fetchall()
                for row in rows:
                    try:
                        existing_ref = json.loads(row["ref"])
                        existing = self._registered_population_source(
                            db, job_id, result, existing_ref)
                    except (ContractError, TypeError, ValueError,
                            json.JSONDecodeError):
                        continue
                    existing_binding = existing["binding"]
                    return {
                        "id": existing["id"],
                        "ref": existing["ref"],
                        "created": False,
                        "source_upload_sha256": existing_binding[
                            "source_upload_sha256"],
                        "source_document_sha256": existing_binding[
                            "source_document_sha256"],
                        "state": "unverified_quantity_pending_expert_review",
                        "scientific_approved": False,
                    }
                evidence_binding = {
                    **binding,
                    "actor_id": actor["id"],
                    "source_upload_sha256": _sha(raw),
                    "source_document_sha256": canonical_sha256,
                    "source_state": "unverified_quantity_pending_expert_review",
                    "scientific_approval": False,
                    "formal_acceptance": False,
                }
                made = self._register_evidence(
                    db, paths, job_id, "microstate_population_source",
                    canonical_sha256, evidence_binding, document,
                    provenance="source")
                return {
                    "id": made["id"],
                    "ref": made["ref"],
                    "created": made["created"],
                    "source_upload_sha256": made["binding"][
                        "source_upload_sha256"],
                    "source_document_sha256": made["binding"][
                        "source_document_sha256"],
                    "state": "unverified_quantity_pending_expert_review",
                    "scientific_approved": False,
                }
        except BaseException:
            for path in paths:
                path.unlink(missing_ok=True)
            raise

    def list_quantitative_sources(self, job_id: str,
                                  session_token: str) -> dict:
        with self.store.db() as db:
            actor = self._authenticate(session_token, db)
            _, result, _, source_state = self._verified_result(
                db, job_id, require_current=False)
            result_sha256 = result.get("files", {}).get("json", {}).get("sha256")
            rows = db.execute("""SELECT ref FROM scientific_evidence
                WHERE project=? AND job_id=?
                AND kind='microstate_population_source' AND superseded=0
                ORDER BY rowid DESC""", (self.project, job_id)).fetchall()
            sources = []
            seen = set()
            for row in rows:
                try:
                    ref = json.loads(row["ref"])
                    key = json.dumps(ref, sort_keys=True, separators=(",", ":"))
                    if key in seen:
                        continue
                    verified = self._registered_population_source(
                        db, job_id, result, ref)
                except (ContractError, TypeError, ValueError,
                        json.JSONDecodeError):
                    continue
                seen.add(key)
                binding = verified["binding"]
                sources.append({
                    "id": verified["id"],
                    "ref": verified["ref"],
                    "created_at": verified["created_at"],
                    "source_document_sha256": binding[
                        "source_document_sha256"],
                    "state": "unverified_quantity_pending_expert_review",
                })

            statement_rows = db.execute("""SELECT id,ref,created_at
                FROM scientific_statements
                WHERE project=? AND actor_id=? ORDER BY rowid DESC""",
                (self.project, actor["id"])).fetchall()
            statements = []
            for row in statement_rows:
                try:
                    ref = json.loads(row["ref"])
                    document = parse_json(self._read_registered_ref(db, ref))
                except (ContractError, TypeError, ValueError,
                        json.JSONDecodeError):
                    continue
                if not (isinstance(document, dict)
                        and document.get("type") == "human_statement"):
                    continue
                text = document.get("text", "")
                statements.append({
                    "id": row["id"], "ref": ref,
                    "locator": document.get("locator"),
                    "text": text[:500] if isinstance(text, str) else "",
                    "created_at": row["created_at"],
                })

            interactions = []
            for _, ref, report in self._interaction_evidence_rows(
                    db, job_id, result_sha256, current_only=True):
                interactions.append({
                    "analog_id": report.get("analog_id"),
                    "report_id": report.get("report_id") or report.get("job_id"),
                    "ref": ref,
                })
            return {
                "project_id": self.project,
                "job_id": job_id,
                "result_sha256": result_sha256,
                "source_inputs_current": source_state["source_inputs_current"],
                "source_runtime_current": source_state["source_runtime_current"],
                "sources": sources,
                "review_statements": statements,
                "interaction_reports": interactions,
            }

    def _interaction_evidence_rows(self, db, job_id: str,
                                   result_sha256: str,
                                   *, current_only: bool = False) -> list[tuple[dict, dict, dict]]:
        suffix = " AND superseded=0" if current_only else ""
        rows = db.execute("""SELECT binding,ref FROM scientific_evidence
            WHERE project=? AND job_id=? AND kind='interaction'""" + suffix +
            " ORDER BY rowid DESC", (self.project, job_id)).fetchall()
        values = []
        for row in rows:
            binding = json.loads(row["binding"])
            if binding.get("result_sha256") != result_sha256:
                continue
            ref = json.loads(row["ref"])
            report = parse_json(self._read_registered_ref(db, ref))
            require(_digest(report) == binding.get("evidence_sha256"),
                    "SCIENTIFIC_EVIDENCE_BINDING")
            values.append((binding, ref, report))
        return values

    def _microstate_report_state(self, report: dict, analog_id: str,
                                 report_id: str, state_index: int,
                                 expected_graph: str, expected_pH: float) -> dict:
        require(report.get("analog_id") == analog_id,
                "SCIENTIFIC_MICROSTATE_INTERACTION_ANALOG")
        require((report.get("report_id") or report.get("job_id")) == report_id,
                "SCIENTIFIC_MICROSTATE_REPORT")
        states = report.get("state_alternatives")
        require(isinstance(states, list) and state_index < len(states),
                "SCIENTIFIC_MICROSTATE_STATE_INDEX")
        state = states[state_index]
        require(isinstance(state, dict) and state.get("status") == "computed",
                "SCIENTIFIC_MICROSTATE_STATE_NOT_COMPUTED")
        graph = scientific_assessment._canonical_isomeric(
            state.get("mapped_smiles", state.get("canonical_smiles")))
        require(graph == expected_graph, "SCIENTIFIC_MICROSTATE_STATE_GRAPH")
        report_pH = scientific_assessment._state_pH(state, report)
        require(report_pH is not None and report_pH == expected_pH,
                "SCIENTIFIC_MICROSTATE_PH_BINDING")
        return state

    def _microstate_population_binding(self, db, job_id: str, result: dict,
                                       row: dict, actor_id: str) -> dict:
        evidence_ref, review_ref = row["evidence_ref"], row["review_ref"]
        require(evidence_ref != review_ref,
                "SCIENTIFIC_MICROSTATE_SOURCE_REVIEW_DISTINCT")
        self._registered_population_source(db, job_id, result, evidence_ref)
        source = self._strict_quantitative_document(
            self._read_registered_ref(db, evidence_ref))
        require(source.get("type") != "human_statement",
                "SCIENTIFIC_MICROSTATE_QUANTITATIVE_SOURCE_REQUIRED")
        record = self._resolve_locator(source, row["locator"])
        require(isinstance(record, dict), "SCIENTIFIC_MICROSTATE_SOURCE_RECORD")
        self._require_human_statement_ref(db, review_ref, actor_id)
        require(isinstance(row.get("interpretation"), str)
                and row["interpretation"].strip(),
                "SCIENTIFIC_MICROSTATE_INTERPRETATION")

        state_index = row["state_index"]
        require(type(state_index) is int and state_index >= 0,
                "SCIENTIFIC_MICROSTATE_STATE_INDEX")
        analog_id = row["analog_id"]
        analog = self._analog_map(result).get(analog_id)
        require(analog is not None, "SCIENTIFIC_MICROSTATE_SCOPE")
        analog_graph = scientific_assessment._isomeric_record(analog)
        require(analog_graph is not None, "SCIENTIFIC_MICROSTATE_ANALOG_GRAPH")
        result_sha256 = result.get("files", {}).get("json", {}).get("sha256")
        require(isinstance(result_sha256, str) and result_sha256,
                "SCIENTIFIC_MICROSTATE_RESULT_BINDING")
        require(record.get("project_id") == self.project
                and record.get("job_id") == job_id
                and record.get("result_sha256") == result_sha256,
                "SCIENTIFIC_MICROSTATE_SOURCE_SCOPE")
        require(record.get("analog_id") == analog_id
                and type(record.get("selected_state_index")) is int
                and record["selected_state_index"] == state_index,
                "SCIENTIFIC_MICROSTATE_SOURCE_SELECTION")

        report_id = row["report_id"]
        declared_pH = row["pH"]
        require(type(declared_pH) in (int, float)
                and math.isfinite(float(declared_pH)),
                "SCIENTIFIC_MICROSTATE_PH_BINDING")
        declared_pH = float(declared_pH)
        state_graph = scientific_assessment._canonical_isomeric(
            record.get("state_canonical_isomeric_graph"))
        require(state_graph is not None, "SCIENTIFIC_MICROSTATE_SOURCE_STATE_GRAPH")
        require(scientific_assessment._canonical_isomeric(
                    record.get("analog_canonical_isomeric_graph")) == analog_graph,
                "SCIENTIFIC_MICROSTATE_SOURCE_ANALOG_GRAPH")
        source_pH = record.get("pH")
        require(type(source_pH) in (int, float)
                and math.isfinite(float(source_pH))
                and float(source_pH) == declared_pH,
                "SCIENTIFIC_MICROSTATE_SOURCE_PH")

        source_interaction_ref = record.get("interaction_ref")
        require(_ref(source_interaction_ref),
                "SCIENTIFIC_MICROSTATE_SOURCE_INTERACTION_REF")
        source_match = next((item for item in self._interaction_evidence_rows(
            db, job_id, result_sha256) if item[1] == source_interaction_ref), None)
        require(source_match is not None,
                "SCIENTIFIC_MICROSTATE_SOURCE_INTERACTION_REF")
        source_binding, _, source_report = source_match
        self._microstate_report_state(
            source_report, analog_id, report_id, state_index,
            state_graph, declared_pH)

        latest = next((item for item in self._interaction_evidence_rows(
            db, job_id, result_sha256, current_only=True)
            if item[2].get("analog_id") == analog_id), None)
        require(latest is not None,
                "SCIENTIFIC_MICROSTATE_CURRENT_INTERACTION_REQUIRED")
        self._microstate_report_state(
            latest[2], analog_id, report_id, state_index,
            state_graph, declared_pH)

        required_text = ("method", "protocol", "conditions",
                         "uncertainty_description", "interpretation")
        require(all(isinstance(record.get(key), str) and record[key].strip()
                    for key in required_text),
                "SCIENTIFIC_MICROSTATE_SOURCE_METADATA")
        source_id, source_url = record.get("source_identifier"), record.get("source_url")
        require((isinstance(source_id, str) and source_id.strip()) or
                (isinstance(source_url, str) and source_url.strip()),
                "SCIENTIFIC_MICROSTATE_SOURCE_IDENTIFIER")
        require(record.get("source_kind") in {"measured", "computed"},
                "SCIENTIFIC_MICROSTATE_SOURCE_KIND")
        has_pka, has_population = "pKa" in record, "state_population" in record
        require(has_pka != has_population,
                "SCIENTIFIC_MICROSTATE_SOURCE_QUANTITY")
        quantity_kind = "pKa" if has_pka else "state_population"
        quantity = record[quantity_kind]
        require(type(quantity) in (int, float) and math.isfinite(float(quantity)),
                "SCIENTIFIC_MICROSTATE_SOURCE_NONFINITE")
        if quantity_kind == "state_population":
            require(0.0 <= float(quantity) <= 1.0,
                    "SCIENTIFIC_MICROSTATE_POPULATION_RANGE")
        return {
            "_validated_by_platform": True,
            "project_id": self.project,
            "job_id": job_id,
            "analog_id": analog_id,
            "report_id": report_id,
            "state_index": state_index,
            "pH": declared_pH,
            "analog_canonical_isomeric_graph": analog_graph,
            "state_canonical_isomeric_graph": state_graph,
            "bound_result_sha256": result_sha256,
            "source_interaction_ref": _copy(source_interaction_ref),
            "source_interaction_sha256": source_binding["evidence_sha256"],
            "quantity_kind": quantity_kind,
            "value": float(quantity),
            "source": (source_id.strip() if isinstance(source_id, str)
                       and source_id.strip() else source_url.strip()),
            "source_kind": record["source_kind"],
            "method": record["method"].strip(),
            "protocol": record["protocol"].strip(),
            "conditions": record["conditions"].strip(),
            "uncertainty_description": record["uncertainty_description"].strip(),
            "scientist_interpretation": row["interpretation"].strip(),
            "source_interpretation": record["interpretation"].strip(),
            "evidence_ref": _copy(evidence_ref),
            "source_locator": row["locator"],
            "review_ref": _copy(review_ref),
        }

    def _validated_population_evidence(self, db, job_id: str, result: dict,
                                       rows: Any, actor_id: str) -> list[dict]:
        if rows is None:
            return []
        require(isinstance(rows, list), "SCIENTIFIC_MICROSTATE_POPULATION_ROWS")
        keys, trusted = set(), []
        for row in rows:
            require(isinstance(row, dict), "SCIENTIFIC_MICROSTATE_POPULATION_ROWS")
            key = (
                json.dumps(row.get("evidence_ref"), sort_keys=True),
                row.get("locator"), row.get("analog_id"), row.get("report_id"),
                row.get("state_index"), row.get("pH"),
            )
            require(key not in keys, "SCIENTIFIC_MICROSTATE_POPULATION_DUPLICATE")
            keys.add(key)
            trusted.append(self._microstate_population_binding(
                db, job_id, result, row, actor_id))
        return trusted

    def _active_decisions(self, db, job_id: str) -> list[dict]:
        rows = db.execute("""SELECT d.body FROM scientific_decisions d
            JOIN review_actors a ON a.project=d.project AND a.id=d.actor_id
            WHERE d.project=? AND d.job_id=? AND d.withdrawn=0 AND a.active=1
            ORDER BY d.policy_revision,d.rowid""", (self.project, job_id)).fetchall()
        values = [json.loads(row[0]) for row in rows]
        for value in values:
            self._verify_refs(db, value)
        return values

    def _latest_policy_revision(self, db, job_id: str) -> int:
        row = db.execute(
            "SELECT COALESCE(MAX(revision),0) FROM scientific_policy_events WHERE project=? AND job_id=?",
            (self.project, job_id),
        ).fetchone()
        return int(row[0])

    def _next_policy_revision(self, db, job_id: str, event_kind: str, decision_id: str) -> int:
        revision = self._latest_policy_revision(db, job_id) + 1
        db.execute("""INSERT INTO scientific_policy_events(
            project,job_id,revision,event_kind,decision_id,created_at)
            VALUES(?,?,?,?,?,?)""",
                   (self.project, job_id, revision, event_kind, decision_id, now()))
        return revision

    def _resolve_locator(self, document: Any, locator: str) -> Any:
        require(isinstance(locator, str) and locator.startswith("/"),
                "SCIENTIFIC_SOURCE_LOCATOR")
        value = document
        for token in locator.split("/")[1:]:
            token = token.replace("~1", "/").replace("~0", "~")
            if isinstance(value, list):
                require(token.isdigit() and int(token) < len(value), "SCIENTIFIC_SOURCE_LOCATOR")
                value = value[int(token)]
            elif isinstance(value, dict):
                require(token in value, "SCIENTIFIC_SOURCE_LOCATOR")
                value = value[token]
            else:
                require(False, "SCIENTIFIC_SOURCE_LOCATOR")
        return _copy(value)

    def _policy(self, db, job_id: str, result: dict) -> dict:
        numerical = dict(scientific_assessment._DEFAULT_NUMERICAL)
        state = {
            "_validated_by_platform": True,
            "numerical_criteria": numerical,
            "scope": {"project_id": self.project, "job_id": job_id},
            "decisions": {},
        }
        sites, interactions = {}, {}
        routes, constraints = {}, {}
        active = self._active_decisions(db, job_id)
        for decision in active:
            if decision.get("action") != "accept" or decision.get("kind") not in POLICY_KINDS:
                continue
            kind, data = decision["kind"], decision["data"]
            if kind == "numerical_criteria":
                numerical.update(data["values"])
                if data.get("minimum_reduction_reviewed") is True:
                    state["decisions"]["repeat_minimum_reduction_reviewed"] = True
            elif kind == "parent_funnel":
                state["decisions"]["parent_funnel"] = data["parents"]
                state["decisions"]["parent_selection_confirmed"] = True
            elif kind == "site_policy":
                for site in data["sites"]:
                    sites[(data["parent_id"], site["atom_map"])] = site
            elif kind == "interaction_requirements":
                interactions[(data["parent_id"], data["analog_id"])] = data
            elif kind == "calibration_criterion":
                state["decisions"]["calibration_criterion"] = {
                    **data, "expert_accept_failed_distribution": False,
                }
            elif kind == "microstate_decision":
                actor_id = decision.get("actor", {}).get("id")
                require(isinstance(actor_id, str) and actor_id,
                        "SCIENTIFIC_MICROSTATE_REVIEWER")
                state["decisions"]["microstate_expert_decision"] = {
                    "accepted": True,
                    "state_choices": data["state_choices"],
                    "pH_conditions": data["pH_conditions"],
                    "population_evidence": self._validated_population_evidence(
                        db, job_id, result, data.get("population_evidence"), actor_id),
                }
            elif kind == "synthesis_policy":
                state["decisions"]["synthesis_criterion"] = {
                    "chosen_candidate_ids": data["chosen_candidate_ids"]
                }
                for item in data["candidate_policies"]:
                    constraints[item["candidate_id"]] = {
                        "blocked_smarts": item["blocked_smarts"]
                    }
                for record in data["exact_route_records"]:
                    source = self.port.json(record["source_ref"])
                    resolved = self._resolve_locator(source, record["locator"])
                    require(isinstance(resolved, dict), "SCIENTIFIC_ROUTE_SOURCE_RECORD")
                    resolved.setdefault("candidate_id", record["candidate_id"])
                    require(resolved.get("candidate_id") == record["candidate_id"],
                            "SCIENTIFIC_ROUTE_CANDIDATE")
                    routes.setdefault(record["candidate_id"], []).append(resolved)
            elif kind == "ternary_criterion":
                state["decisions"]["ternary_criterion"] = {
                    "chosen_candidates": data["chosen_candidates"],
                    "expected_seeds": data["expected_seeds"],
                }
                geometry = _copy(data["geometry"])
                geometry["disposition"] = (
                    "accept" if geometry.pop("accepted") is True else "reject"
                )
                state["decisions"]["ternary_geometry_criterion"] = geometry
            elif kind == "followup_geometry_review":
                trusted = _copy(data)
                for branch in trusted["branches"]:
                    for review in branch["seed_reviews"]:
                        review["evidence_binding"] = (
                            self._followup_geometry_evidence_binding(
                                db, job_id, result, review["evidence_ref"],
                                branch["candidate_id"], branch["e3_type"],
                                review["seed"], require_current=False,
                            )
                        )
                trusted.update({
                    "_validated_by_platform": True,
                    "decision_id": decision["id"],
                    "actor": _copy(decision["actor"]),
                    "scope": {"project_id": self.project, "job_id": job_id},
                    "source_binding": {
                        "project_id": self.project,
                        "job_id": job_id,
                        "result_sha256": result.get("files", {}).get(
                            "json", {}).get("sha256"),
                        "source_ref": _copy(data["source_ref"]),
                        "source_sha256": data["source_ref"]["sha256"],
                        "followup_source_policy_digest":
                            expert_followup.SOURCE_POLICY_DIGEST,
                    },
                })
                state["decisions"]["followup_geometry_review"] = trusted
            elif kind == "protein_hydrogen_review":
                state["decisions"]["protein_hydrogen_review"] = {
                    "evidence_ref": _copy(data["evidence_ref"]),
                    "accepted": data["accepted"],
                    "rationale": data["rationale"],
                    "source_ref": _copy(data["source_ref"]),
                    "actor": _copy(decision["actor"]),
                    "decision_id": decision["id"],
                }
        state["decisions"]["site_policy"] = [sites[key] for key in sorted(sites)]
        state["decisions"]["interaction_requirements"] = [
            interactions[key] for key in sorted(interactions)
        ]
        expert_opinion.apply_bundled_opinion(
            state,
            result,
            active,
            archived_result_verified=True,
        )
        if state.get("opinion", {}).get("scope", {}).get("fx5_bound_scope") is True:
            state["followup"] = expert_followup.load_followup()
        state["_resolved_synthesis"] = {
            "constraints": constraints,
            "exact_route_records": routes,
        }
        state["revision"] = self._latest_policy_revision(db, job_id)
        state["active_policy_decision_ids"] = [
            item["id"] for item in active
            if item.get("kind") in POLICY_KINDS and item.get("action") == "accept"
        ]
        return state

    def _policy_digest(self, policy: dict) -> str:
        return scientific_assessment._digest(policy)

    def _public_policy(self, policy: dict) -> dict:
        value = _copy(policy)
        value.pop("_validated_by_platform", None)
        value.pop("_resolved_synthesis", None)
        return value

    def policy(self, job_id: str) -> dict:
        with self.store.db() as db:
            _, result, _, source = self._verified_result(db, job_id)
            policy = self._policy(db, job_id, result)
            value = self._public_policy(policy)
            value.update(source_inputs_current=source["source_inputs_current"],
                         source_runtime_current=source["source_runtime_current"])
            return value

    def policy_for_design(self, assessment_id: str) -> dict:
        with self.store.db() as db:
            assessment = self._assessment(db, assessment_id)
            _, result, binding, source = self._verified_result(
                db, assessment["job_id"], require_current=True)
            policy = self._policy(db, assessment["job_id"], result)
            digest = self._policy_digest(policy)
            require(policy["revision"] >= 1, "SCIENTIFIC_ACTIVE_POLICY_REQUIRED")
            require(assessment.get("policy_revision") == policy["revision"] and
                    assessment.get("policy_digest") == digest,
                    "SCIENTIFIC_STALE_ASSESSMENT")
            require(source["source_inputs_current"] and source["source_runtime_current"],
                    "SCIENTIFIC_DESIGN_SOURCE_CHANGED")
            parent_id = binding["parent_id"]
            active = self._active_decisions(db, assessment["job_id"])
            site_rows = {}
            parent_funnel = []
            for decision in active:
                if decision.get("action") != "accept":
                    continue
                if decision.get("kind") == "parent_funnel":
                    parent_funnel = _copy(decision["data"]["parents"])
                elif (decision.get("kind") == "site_policy" and
                      decision.get("data", {}).get("parent_id") == parent_id):
                    for site in decision["data"]["sites"]:
                        release = site.get("allow_release_protected") is True
                        source_ref = (site.get("release_source_ref") if release
                                      else site.get("source_ref"))
                        rationale = (site.get("release_justification") if release
                                     else site.get("rationale"))
                        require(_ref(source_ref), "SCIENTIFIC_SITE_SOURCE_REF")
                        require(isinstance(rationale, str) and rationale.strip(),
                                "SCIENTIFIC_SITE_RATIONALE")
                        resolved = {
                            "atom_map": site["atom_map"],
                            "state": site["state"],
                            "allowed_transforms": list(site["allowed_rule_ids"]),
                            "rationale": rationale,
                            "source_ref": _copy(source_ref),
                            "allow_release_protected": release,
                        }
                        site_rows[site["atom_map"]] = resolved
            require(bool(site_rows), "SCIENTIFIC_ACCEPTED_SITE_POLICY_REQUIRED")
            return {
                "scope": {
                    "project_id": self.project,
                    "job_id": assessment["job_id"],
                    "parent_id": parent_id,
                    "assessment_id": assessment_id,
                    "policy_revision": policy["revision"],
                    "policy_digest": digest,
                },
                "site_policy": [site_rows[key] for key in sorted(site_rows)],
                "parent_funnel": parent_funnel,
            }

    def _supplements(self, db, job_id: str, binding: dict, policy: dict) -> dict:
        policy_digest = self._policy_digest(policy)
        result = {"interactions": [], "synthesis": [], "ternary": [], "calibration": []}
        rows = db.execute("""SELECT kind,binding,ref FROM scientific_evidence
            WHERE project=? AND job_id=? AND superseded=0 ORDER BY rowid""",
            (self.project, job_id)).fetchall()
        for row in rows:
            stored = json.loads(row["binding"])
            if stored.get("result_sha256") != binding["result_sha256"]:
                continue
            if stored.get("evaluation_policy_digest") != policy_digest:
                continue
            ref = json.loads(row["ref"])
            value = self.port.json(ref)
            require(_digest(value) == stored["evidence_sha256"],
                    "SCIENTIFIC_EVIDENCE_BINDING")
            if row["kind"] == "calibration_distribution":
                measurement = stored.get("measurement_binding")
                if not isinstance(measurement, dict):
                    continue
                if measurement.get("source_id") != "CRBN_6BOY":
                    continue
                if measurement.get("module_version") != calibration_import.VERSION:
                    continue
                if measurement.get("module_fingerprint") != calibration_import.measurement_fingerprint():
                    continue
            wrapped = {
                **value,
                "_validated_by_platform": True,
                "job_id": job_id,
                "current_job_id": job_id,
                "policy_digest": policy_digest,
                "measurement_binding": stored.get("measurement_binding"),
                "evaluation_policy": stored.get("evaluation_policy"),
                "bound_result_sha256": stored.get("result_sha256"),
                "interaction_evidence_sha256": stored.get("evidence_sha256"),
            }
            if row["kind"] == "interaction":
                result["interactions"].append(wrapped)
            elif row["kind"] == "synthesis":
                result["synthesis"].append(wrapped)
            elif row["kind"] == "ternary":
                result["ternary"].append(wrapped)
            elif row["kind"] == "calibration_distribution":
                result["calibration"].append(wrapped)
        return result

    def _registered_evidence(self, db, job_id: str, result: dict) -> list[dict]:
        result_sha256 = result.get("files", {}).get("json", {}).get("sha256")
        rows = db.execute("""SELECT id,kind,binding,ref FROM scientific_evidence
            WHERE project=? AND job_id=? AND superseded=0 ORDER BY rowid""",
            (self.project, job_id)).fetchall()
        registered = []
        for row in rows:
            binding = json.loads(row["binding"])
            if binding.get("result_sha256") != result_sha256:
                continue
            ref = json.loads(row["ref"])
            raw = self.port.read(ref)
            require(_sha(raw) == ref["sha256"], "SCIENTIFIC_ARTIFACT_HASH")
            value = parse_json(raw)
            require(_digest(value) == binding.get("evidence_sha256"),
                    "SCIENTIFIC_EVIDENCE_BINDING")
            registered.append({
                "evidence_id": row["id"],
                "registered_type": row["kind"],
                "ref": ref,
            })
        return registered

    def create(self, job_id, expected_result_sha256=None, *, reviewer_id=None,
               evidence_refs=None, protocol_ids=None):
        paths = []
        evidence_refs = [] if evidence_refs is None else evidence_refs
        if protocol_ids is not None:
            require(isinstance(protocol_ids, list)
                    and 1 <= len(protocol_ids) <= protocol_reassessment.MAX_PROTOCOL_IDS
                    and all(isinstance(item, str) and item for item in protocol_ids)
                    and len(set(protocol_ids)) == len(protocol_ids),
                    "SCIENTIFIC_CREATE_PROTOCOL_IDS")
        require(isinstance(evidence_refs, list) and all(_ref(item) for item in evidence_refs),
                "SCIENTIFIC_CREATE_EVIDENCE_REFS")
        try:
            with self.store.db() as db:
                db.execute("BEGIN IMMEDIATE")
                job, result, binding, source = self._verified_result(db, job_id)
                if expected_result_sha256 is not None:
                    require(expected_result_sha256 == binding["result_sha256"],
                            "SCIENTIFIC_RESULT_CHANGED")
                self._verify_refs(db, evidence_refs)
                if reviewer_id is not None:
                    actor = self.auth.actor(db, reviewer_id)
                    require(actor is not None and actor["active"], "SCIENTIFIC_REVIEWER_INACTIVE")
                policy = self._policy(db, job_id, result)
                policy_digest = self._policy_digest(policy)
                protocol_closure = None
                if protocol_ids is not None:
                    protocol_closure = protocol_reassessment.load(
                        db, self.store, self.project, job_id, protocol_ids,
                        binding, result, policy, policy_digest,
                    )
                supplements = self._supplements(db, job_id, binding, policy)
                assessment = scientific_assessment.evaluate_design(
                    result, supplements=supplements, policy=policy)
                if protocol_closure is not None:
                    protocol_reassessment.apply(assessment, protocol_closure)
                implementation = engine_fingerprint()
                assessment_engine = scientific_assessment.compute_engine_fingerprint()
                assessment.update({
                    "id": "assessment-" + uuid.uuid4().hex,
                    "project_id": self.project,
                    "job_id": job_id,
                    "created_at": now(),
                    "source_binding": binding,
                    "policy_revision": policy["revision"],
                    "policy_digest": policy_digest,
                    "assigned_reviewer_id": reviewer_id,
                    "evidence_refs": evidence_refs,
                    "service_version": VERSION,
                    "service_engine_fingerprint": implementation,
                    "assessment_engine_fingerprint": assessment_engine,
                    "archived_source_inputs_current": source["source_inputs_current"],
                    "archived_source_runtime_current": source["source_runtime_current"],
                })
                fingerprint = _digest({
                    "binding": binding,
                    "policy_digest": policy_digest,
                    "supplements": supplements,
                    "implementation": implementation,
                    "assessment_engine": assessment_engine,
                    "reviewer_id": reviewer_id,
                    "evidence_refs": evidence_refs,
                    "protocol_evaluation": protocol_closure,
                })
                prior = db.execute("""SELECT ref FROM scientific_assessments
                    WHERE project=? AND job_id=? AND fingerprint=?""",
                    (self.project, job_id, fingerprint)).fetchone()
                if prior:
                    archived = self.port.json(json.loads(prior["ref"]))
                    return self._effective(db, archived), False
                revision = db.execute("""SELECT COALESCE(MAX(revision),0)+1
                    FROM scientific_assessments WHERE project=? AND job_id=?""",
                    (self.project, job_id)).fetchone()[0]
                assessment["revision"] = revision
                ref = self._write(db, paths, assessment)
                db.execute("""INSERT INTO scientific_assessments(
                    id,project,job_id,revision,policy_revision,policy_digest,
                    fingerprint,ref) VALUES(?,?,?,?,?,?,?,?)""",
                    (assessment["id"], self.project, job_id, revision,
                     policy["revision"], policy_digest, fingerprint, json.dumps(ref)))
                return self._effective(db, assessment), True
        except BaseException:
            for path in paths:
                path.unlink(missing_ok=True)
            raise

    def _assessment(self, db, identifier: str) -> dict:
        row = db.execute(
            "SELECT ref FROM scientific_assessments WHERE id=? AND project=?",
            (identifier, self.project),
        ).fetchone()
        if row is None:
            raise KeyError(identifier)
        value = self.port.json(json.loads(row["ref"]))
        self._verified_result(db, value["job_id"], require_current=False)
        return value

    def _effective(self, db, assessment: dict) -> dict:
        value = _copy(assessment)
        job, result, _, source = self._verified_result(db, value["job_id"])
        policy = self._policy(db, value["job_id"], result)
        current_policy_digest = self._policy_digest(policy)
        value["current_policy_revision"] = policy["revision"]
        value["current_policy_digest"] = current_policy_digest
        value["current_policy"] = (
            value.get("policy_revision") == policy["revision"]
            and value.get("policy_digest") == current_policy_digest
        )
        value["source_inputs_current"] = source["source_inputs_current"]
        value["source_inputs_reason"] = source["source_inputs_reason"]
        value["source_runtime_current"] = source["source_runtime_current"]
        value["source_runtime_reason"] = source["source_runtime_reason"]
        value["current_implementation"] = (
            value.get("service_engine_fingerprint") == engine_fingerprint()
            and value.get("assessment_engine_fingerprint") ==
                scientific_assessment.compute_engine_fingerprint()
        )
        value["current_job_runtime"] = source["source_runtime_current"]
        value["archived"] = True
        protocol_current, protocol_reason = protocol_reassessment.freshness(
            db, self, value, policy, current_policy_digest
        )
        value["protocol_evidence_current"] = protocol_current
        value["protocol_evidence_freshness_reason"] = protocol_reason
        if not protocol_current:
            protocol_reassessment.mark_stale(value, protocol_reason)
        value["registered_evidence"] = self._registered_evidence(
            db, value["job_id"], result
        )

        decisions = self._active_decisions(db, value["job_id"])
        formal = [item for item in decisions if item.get("kind") == "formal_decision"
                  and item.get("assessment_id") == value["id"]]
        latest = formal[-1] if formal else None
        all_current = all((
            value["current_policy"], value["source_inputs_current"],
            value["source_runtime_current"], value["current_implementation"],
            value["protocol_evidence_current"],
        ))
        for criterion in value.get("criteria", []):
            if criterion.get("id") == "formal_expert_decision":
                action = latest.get("action") if latest and all_current else None
                criterion["observed"] = action
                criterion["status"] = (
                    "pass" if action == "accept" else
                    "failed" if action == "reject" else "pending"
                )
        failures = [item["id"] for item in value.get("criteria", [])
                    if item.get("status") != "pass"]
        value["scientific_accepted"] = bool(
            all_current and latest and latest.get("action") == "accept" and not failures)
        value["confirmation_status"] = (
            "accepted" if value["scientific_accepted"] else
            "rejected" if latest and latest.get("action") == "reject" else "unconfirmed"
        )
        value["active_formal_decision"] = latest
        value["criterion_failures"] = failures
        return value

    def view(self, identifier: str) -> dict:
        with self.store.db() as db:
            return self._effective(db, self._assessment(db, identifier))

    def list(self, job_id: str, *, compact: bool = False) -> list[dict]:
        with self.store.db() as db:
            if compact:
                self._job(db, job_id)
                rows = db.execute("""SELECT id,revision FROM scientific_assessments
                    WHERE project=? AND job_id=? ORDER BY revision DESC""",
                    (self.project, job_id)).fetchall()
                return [{
                    "id": row["id"],
                    "revision": row["revision"],
                    "summary_only": True,
                } for row in rows]
            self._verified_result(db, job_id, require_current=False)
            ids = [row[0] for row in db.execute("""SELECT id FROM scientific_assessments
                WHERE project=? AND job_id=? ORDER BY revision DESC""",
                (self.project, job_id))]
            return [self._effective(db, self._assessment(db, identifier)) for identifier in ids]

    def _candidate_map(self, result: dict) -> dict[str, dict]:
        return {row["candidate_id"]: row for row in result.get("protac_candidates", [])
                if isinstance(row, dict) and isinstance(row.get("candidate_id"), str)}

    def _analog_map(self, result: dict) -> dict[str, dict]:
        return {row["id"]: row for row in result.get("analogs", [])
                if isinstance(row, dict) and isinstance(row.get("id"), str)}

    def _parent_mol(self, parent_id: str):
        from packages.science.dual_e3 import REFERENCE_SOURCE, SOURCE
        if parent_id == "SMARCA2-FX5":
            path = SOURCE / "SMARCA2-neutral-design.sdf"
            require(path.is_file() and not path.is_symlink(),
                    "SCIENTIFIC_PARENT_STRUCTURE_MISSING")
            mol = Chem.MolFromMolBlock(path.read_text(encoding="utf-8"), removeHs=False)
        else:
            from packages.science.reference_parents import load_reference_parent
            mol, _ = load_reference_parent(parent_id, REFERENCE_SOURCE)
        require(mol is not None and mol.GetNumConformers() == 1,
                "SCIENTIFIC_PARENT_STRUCTURE_INVALID")
        return mol

    def _protein_atoms(self):
        from packages.science.dual_e3 import SOURCE
        from packages.science.structures import atom_sites
        return atom_sites(SOURCE / "6HAZ.cif", ["A"])

    def _followup_geometry_evidence_binding(
            self, db, job_id: str, result: dict, evidence_ref: dict,
            candidate_id: str, e3_type: str, seed: int,
            *, require_current: bool) -> dict:
        require(_ref(evidence_ref),
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_EVIDENCE_REF")
        rows = db.execute("""SELECT id,binding,ref,superseded
            FROM scientific_evidence
            WHERE project=? AND job_id=? AND kind='ternary'
            ORDER BY rowid""", (self.project, job_id)).fetchall()
        matches = [row for row in rows
                   if json.loads(row["ref"]) == evidence_ref]
        require(len(matches) == 1,
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_EVIDENCE")
        row = matches[0]
        if require_current:
            require(row["superseded"] == 0,
                    "SCIENTIFIC_FOLLOWUP_GEOMETRY_EVIDENCE_NOT_CURRENT")

        raw = self._read_registered_ref(db, evidence_ref)
        require(_sha(raw) == evidence_ref["sha256"],
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_ARTIFACT_HASH")
        report = parse_json(raw)
        require(isinstance(report, dict),
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_EVIDENCE_JSON")
        stored = json.loads(row["binding"])
        require(_digest(report) == stored.get("evidence_sha256"),
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_EVIDENCE_HASH")
        self._verify_refs(db, report)
        require(report.get("derived_interpretation") is not True and
                report.get("candidate_id") == candidate_id and
                report.get("e3_type") == e3_type and
                report.get("seed") == seed and
                report.get("actual_computation") is True and
                report.get("execution_success") is True and
                report.get("reference_free") is True and
                report.get("failure") is None,
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_REPORT_IDENTITY")

        raw_artifacts = report.get("registered_artifacts")
        require(isinstance(raw_artifacts, dict) and raw_artifacts and
                all(_ref(ref) for ref in raw_artifacts.values()),
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_RAW_ARTIFACTS")
        raw_hashes = {}
        for name, ref in sorted(raw_artifacts.items()):
            artifact_raw = self._read_registered_ref(db, ref)
            require(_sha(artifact_raw) == ref["sha256"],
                    "SCIENTIFIC_FOLLOWUP_GEOMETRY_RAW_ARTIFACT_HASH")
            raw_hashes[name] = ref["sha256"]
        for name in ("plan_ref", "aggregate_receipt_ref"):
            ref = report.get(name)
            require(_ref(ref),
                    "SCIENTIFIC_FOLLOWUP_GEOMETRY_RAW_ARTIFACTS")
            artifact_raw = self._read_registered_ref(db, ref)
            require(_sha(artifact_raw) == ref["sha256"],
                    "SCIENTIFIC_FOLLOWUP_GEOMETRY_RAW_ARTIFACT_HASH")
            raw_hashes[name] = ref["sha256"]

        result_ref = result.get("files", {}).get("json", {})
        parent_id = result.get("parent_scope", {}).get(
            "actual_design_parent_id")
        require(stored.get("project") == self.project and
                stored.get("job_id") == job_id and
                stored.get("input_sha256") and
                stored.get("result_sha256") == result_ref.get("sha256") and
                stored.get("parent_id") == parent_id,
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_SOURCE_BINDING")
        measurement = stored.get("measurement_binding")
        protocol = report.get("measurement_protocol")
        producer = report.get("producer_policy_binding")
        require(isinstance(measurement, dict) and
                isinstance(protocol, dict) and
                isinstance(producer, dict) and
                measurement.get("producer_policy") == producer and
                measurement.get("candidate_graph_sha256") ==
                    protocol.get("candidate_graph_sha256") and
                measurement.get("plan_digest") == protocol.get("plan_digest") and
                protocol.get("seed") == seed,
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_MEASUREMENT_BINDING")

        plan = self.port.json(report["plan_ref"])
        plan_job = plan.get("bindings", {}).get("job", {})
        plan_graph = plan.get("candidate_graph", {})
        require(plan.get("plan_digest") == measurement.get("plan_digest") and
                plan_job.get("project") == self.project and
                plan_job.get("job_id") == job_id and
                plan_job.get("input_sha256") == stored.get("input_sha256") and
                plan_job.get("result_sha256") == stored.get("result_sha256") and
                plan_graph.get("candidate_id") == candidate_id and
                plan_graph.get("e3_type") == e3_type,
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_PLAN_BINDING")

        candidates = self._candidate_map(result)
        require(candidate_id in candidates and
                candidates[candidate_id].get("e3_type") == e3_type,
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_CANDIDATE")
        from scripts.run_novel_ternary import _candidate as qualified_candidate
        graph = novel_ternary._mapped_graph(
            qualified_candidate(result, candidate_id))
        require(graph.get("graph_sha256") ==
                    measurement.get("candidate_graph_sha256") ==
                    plan_graph.get("graph_sha256") ==
                    plan_job.get("candidate_graph_sha256"),
                "SCIENTIFIC_FOLLOWUP_GEOMETRY_GRAPH_BINDING")
        return {
            "evidence_id": row["id"],
            "evidence_ref": _copy(evidence_ref),
            "evidence_sha256": stored["evidence_sha256"],
            "artifact_sha256": evidence_ref["sha256"],
            "raw_artifact_sha256": raw_hashes,
            "candidate_graph_sha256": graph["graph_sha256"],
            "plan_digest": measurement["plan_digest"],
            "project_id": self.project,
            "job_id": job_id,
            "input_sha256": stored["input_sha256"],
            "result_sha256": stored["result_sha256"],
            "parent_id": parent_id,
            "candidate_id": candidate_id,
            "e3_type": e3_type,
            "seed": seed,
            "actual_computation": True,
            "execution_success": True,
            "reference_free": True,
        }

    def _validate_policy_data(self, db, job_id: str, kind: str,
                              data: dict, result: dict,
                              actor_id: str | None = None) -> None:
        if kind == "numerical_criteria":
            values = data["values"]
            defaults = scientific_assessment._DEFAULT_NUMERICAL
            require(set(values) <= set(defaults), "SCIENTIFIC_NUMERICAL_KEY")
            require(all(type(value) is int and value >= 1 for value in values.values()),
                    "SCIENTIFIC_NUMERICAL_VALUE")
            lowered = sorted(key for key, value in values.items() if value < defaults[key])
            if lowered:
                require(data.get("minimum_reduction_reviewed") is True
                        and isinstance(data.get("reduction_justification"), str)
                        and data["reduction_justification"].strip()
                        and _ref(data.get("reduction_source_ref")),
                        "SCIENTIFIC_NUMERICAL_REDUCTION_REVIEW")
            merged = {**defaults, **values}
            require(merged["parent_funnel_min"] <= merged["parent_funnel_max"] <= 10,
                    "SCIENTIFIC_PARENT_RANGE")
            require(merged["panel_min"] <= merged["panel_max"] <= 20,
                    "SCIENTIFIC_PANEL_RANGE")
            require(merged["novel_ternary_repeats_min"] >= 2,
                    "SCIENTIFIC_TERNARY_MINIMUM")
        elif kind == "parent_funnel":
            from packages.science.dual_e3 import DEFAULT_PARENT_ID, catalog
            available = {DEFAULT_PARENT_ID: True}
            for row in catalog()["reference_parents"]["available"]:
                if isinstance(row, dict) and isinstance(row.get("id"), str):
                    available[row["id"]] = bool(row.get("ready") or row.get("technical_ready"))
            identifiers = [row["parent_id"] for row in data["parents"]]
            require(len(identifiers) == len(set(identifiers)), "SCIENTIFIC_PARENT_SELECTION")
            require(all(identifier in available and available[identifier]
                        for identifier in identifiers), "SCIENTIFIC_PARENT_NOT_READY")
        elif kind == "site_policy":
            actual_parent = result["parent_scope"]["actual_design_parent_id"]
            require(data["parent_id"] == actual_parent, "SCIENTIFIC_SITE_PARENT")
            atoms = result.get("sites", {}).get("atoms", [])
            require(isinstance(atoms, list), "SCIENTIFIC_SITE_MAP")
            actual = {row.get("atom_map"): row for row in atoms if isinstance(row, dict)}
            require(len(actual) == len(atoms), "SCIENTIFIC_SITE_MAP")
            submitted_maps = [row["atom_map"] for row in data["sites"]]
            require(len(submitted_maps) == len(set(submitted_maps)),
                    "SCIENTIFIC_SITE_MAP_DUPLICATE")
            from packages.science.dual_e3 import catalog
            rules = catalog()["rule_catalog"]
            rule_ids = {row["rule_id"] for row in rules
                        if isinstance(row, dict) and isinstance(row.get("rule_id"), str)}
            allowed_changes = {
                ("UNKNOWN", "MODIFIABLE"),
                ("UNKNOWN", "PROTECTED"),
                ("MODIFIABLE", "PROTECTED"),
                ("PROTECTED", "MODIFIABLE"),
            }
            for row in data["sites"]:
                require(row["atom_map"] in actual, "SCIENTIFIC_SITE_MAP")
                original = actual[row["atom_map"]].get("state")
                desired = row["state"]
                require(original in {"MODIFIABLE", "PROTECTED", "UNKNOWN"}
                        and desired in {"MODIFIABLE", "PROTECTED", "UNKNOWN"},
                        "SCIENTIFIC_SITE_STATE")
                if "original_state" in row:
                    require(row["original_state"] == original,
                            "SCIENTIFIC_SITE_ORIGINAL_STATE_MISMATCH")
                require(set(row["allowed_rule_ids"]) <= rule_ids,
                        "SCIENTIFIC_SITE_RULE")
                changed = desired != original
                if changed:
                    require((original, desired) in allowed_changes,
                            "SCIENTIFIC_SITE_STATE_TRANSITION")
                    require(row.get("original_state") == original,
                            "SCIENTIFIC_SITE_ORIGINAL_STATE_REQUIRED")
                    require(row.get("expert_classification") is True,
                            "SCIENTIFIC_SITE_EXPERT_CLASSIFICATION_REQUIRED")
                    require(isinstance(row.get("rationale"), str)
                            and row["rationale"].strip(),
                            "SCIENTIFIC_SITE_RATIONALE")
                    self._require_human_statement_ref(
                        db, row.get("source_ref"), actor_id
                    )
                else:
                    require(row.get("expert_classification") is not True,
                            "SCIENTIFIC_SITE_EXPERT_CLASSIFICATION_NOT_A_CHANGE")
                protected_release = original == "PROTECTED" and desired == "MODIFIABLE"
                if protected_release:
                    require(row.get("allow_release_protected") is True
                            and isinstance(row.get("release_justification"), str)
                            and row["release_justification"].strip()
                            and _ref(row.get("release_source_ref")),
                            "SCIENTIFIC_PROTECTED_SITE_RELEASE")
                    self._require_human_statement_ref(
                        db, row["release_source_ref"], actor_id
                    )
                elif row.get("allow_release_protected") is True:
                    # A legacy no-op PROTECTED declaration may retain release evidence,
                    # but release can never classify UNKNOWN or alter another state.
                    require(original == desired == "PROTECTED"
                            and isinstance(row.get("release_justification"), str)
                            and row["release_justification"].strip()
                            and _ref(row.get("release_source_ref")),
                            "SCIENTIFIC_SITE_RELEASE_FLAG")
                    self._require_human_statement_ref(
                        db, row["release_source_ref"], actor_id
                    )
                else:
                    require(row.get("allow_release_protected") is False,
                            "SCIENTIFIC_SITE_RELEASE_FLAG")
        elif kind == "interaction_requirements":
            parent_id = result["parent_scope"]["actual_design_parent_id"]
            require(data["parent_id"] == parent_id, "SCIENTIFIC_INTERACTION_PARENT")
            analog = self._analog_map(result).get(data["analog_id"])
            require(analog is not None, "SCIENTIFIC_INTERACTION_ANALOG")
            parent = self._parent_mol(parent_id)
            analog_mol = Chem.MolFromSmiles(analog.get("mapped_smiles", ""))
            require(analog_mol is not None, "SCIENTIFIC_INTERACTION_ANALOG_GRAPH")
            parent_maps = {atom.GetAtomMapNum() for atom in parent.GetAtoms()
                           if atom.GetAtomicNum() > 1 and atom.GetAtomMapNum() > 0}
            analog_maps = {atom.GetAtomMapNum() for atom in analog_mol.GetAtoms()
                           if atom.GetAtomicNum() > 1 and atom.GetAtomMapNum() > 0}
            shared = parent_maps & analog_maps
            protein_ids = {
                ":".join((str(atom["label_asym_id"]), str(atom["auth_seq_id"]),
                          str(atom["label_comp_id"]).upper(),
                          str(atom["label_atom_id"]).upper()))
                for atom in self._protein_atoms()
            }
            for requirement in data["requirements"]:
                require(set(requirement["ligand_maps"]) <= shared,
                        "SCIENTIFIC_INTERACTION_LIGAND_MAP")
                require(set(requirement["protein_atom_ids"]) <= protein_ids,
                        "SCIENTIFIC_INTERACTION_PROTEIN_ATOM")
        elif kind == "calibration_criterion":
            receipts = result.get("calibration", {}).get("CRBN", {}).get("seed_receipts", [])
            metrics = set()
            for receipt in receipts if isinstance(receipts, list) else []:
                if not isinstance(receipt, dict):
                    continue
                raw = receipt.get("comparison", {}).get("metrics", {})
                if not raw and isinstance(receipt.get("inspection"), dict):
                    raw = receipt["inspection"].get("comparison", {}).get("metrics", {})
                if isinstance(raw, dict):
                    metrics.update(raw)
            require(data["metric"] in metrics, "SCIENTIFIC_CALIBRATION_METRIC")
        elif kind == "microstate_decision":
            analogs = self._analog_map(result)
            require(set(data["state_choices"]) <= set(analogs),
                    "SCIENTIFIC_MICROSTATE_SCOPE")
            require(all(type(index) is int and index >= 0
                        for index in data["state_choices"].values()),
                    "SCIENTIFIC_MICROSTATE_STATE_INDEX")
            pH = data["pH_conditions"]["pH"]
            require(type(pH) in (int, float) and math.isfinite(float(pH)),
                    "SCIENTIFIC_MICROSTATE_PH")
            trusted = self._validated_population_evidence(
                db, job_id, result, data.get("population_evidence"), actor_id or "")
            require(all(row["pH"] == float(pH) for row in trusted),
                    "SCIENTIFIC_MICROSTATE_PH_BINDING")
            require(all(data["state_choices"].get(row["analog_id"]) ==
                        row["state_index"] for row in trusted),
                    "SCIENTIFIC_MICROSTATE_STATE_CHOICE_BINDING")
        elif kind == "synthesis_policy":
            candidates = self._candidate_map(result)
            chosen = data["chosen_candidate_ids"]
            require(set(chosen) <= set(candidates), "SCIENTIFIC_SYNTHESIS_CANDIDATE")
            require({item["candidate_id"] for item in data["candidate_policies"]} == set(chosen),
                    "SCIENTIFIC_SYNTHESIS_POLICY_COVERAGE")
            for record in data["exact_route_records"]:
                require(record["candidate_id"] in chosen, "SCIENTIFIC_ROUTE_CANDIDATE")
                source = self.port.json(record["source_ref"])
                resolved = self._resolve_locator(source, record["locator"])
                require(isinstance(resolved, dict), "SCIENTIFIC_ROUTE_SOURCE_RECORD")
                require(resolved.get("candidate_id") == record["candidate_id"],
                        "SCIENTIFIC_ROUTE_CANDIDATE")
                expected = candidates[record["candidate_id"]].get("canonical_smiles")
                route_mol = Chem.MolFromSmiles(resolved.get("canonical_smiles", ""))
                candidate_mol = Chem.MolFromSmiles(expected or "")
                require(route_mol is not None and candidate_mol is not None and
                        Chem.MolToSmiles(route_mol, True) == Chem.MolToSmiles(candidate_mol, True),
                        "SCIENTIFIC_ROUTE_GRAPH")
                require(isinstance(resolved.get("steps"), list) and resolved["steps"],
                        "SCIENTIFIC_ROUTE_STEPS")
        elif kind == "ternary_criterion":
            candidates = self._candidate_map(result)
            chosen = data["chosen_candidates"]
            require(set(chosen) == {"CRBN", "VHL"}, "SCIENTIFIC_TERNARY_E3_SCOPE")
            for e3, candidate_ids in chosen.items():
                require(isinstance(candidate_ids, list) and candidate_ids,
                        "SCIENTIFIC_TERNARY_CANDIDATE")
                require(all(candidate_id in candidates and
                            candidates[candidate_id].get("e3_type") == e3
                            for candidate_id in candidate_ids),
                        "SCIENTIFIC_TERNARY_CANDIDATE")
            require(len(data["expected_seeds"]) >= 2,
                    "SCIENTIFIC_TERNARY_SEEDS")
            require(data["geometry"]["reference_free"] is True,
                    "SCIENTIFIC_TERNARY_REFERENCE_FREE")
        elif kind == "followup_geometry_review":
            require(result.get("parent_scope", {}).get(
                        "actual_design_parent_id") == "SMARCA2-FX5",
                    "SCIENTIFIC_FOLLOWUP_GEOMETRY_SCOPE")
            followup = expert_followup.load_followup()
            require(followup.get("source_policy_digest") ==
                        expert_followup.SOURCE_POLICY_DIGEST,
                    "SCIENTIFIC_FOLLOWUP_GEOMETRY_SOURCE")
            self._require_human_statement_ref(
                db, data["source_ref"], actor_id)
            expected = {
                "CRBN": "D-99b12e64986a",
                "VHL": "D-b39273b7a53b",
            }
            branches = data["branches"]
            require(len(branches) == 2 and
                    [branch["e3_type"] for branch in branches].count("CRBN") == 1 and
                    [branch["e3_type"] for branch in branches].count("VHL") == 1 and
                    {branch["candidate_id"] for branch in branches} ==
                        set(expected.values()),
                    "SCIENTIFIC_FOLLOWUP_GEOMETRY_BRANCHES")
            for branch in branches:
                e3_type = branch["e3_type"]
                candidate_id = branch["candidate_id"]
                require(candidate_id == expected[e3_type] and
                        isinstance(branch["clash_definition"], str) and
                        branch["clash_definition"].strip() and
                        isinstance(branch["rationale"], str) and
                        branch["rationale"].strip(),
                        "SCIENTIFIC_FOLLOWUP_GEOMETRY_BRANCH")
                reviews = branch["seed_reviews"]
                seeds = [review["seed"] for review in reviews]
                require(len(reviews) == 5 and len(set(seeds)) == 5 and
                        set(seeds) == set(expert_followup.EXPECTED_SEEDS),
                        "SCIENTIFIC_FOLLOWUP_GEOMETRY_SEEDS")
                for review in reviews:
                    require(isinstance(review["topology_label"], str) and
                            review["topology_label"].strip() and
                            type(review["severe_clash_acceptable"]) is bool and
                            isinstance(review["rationale"], str) and
                            review["rationale"].strip(),
                            "SCIENTIFIC_FOLLOWUP_GEOMETRY_SEED_REVIEW")
                    self._followup_geometry_evidence_binding(
                        db, job_id, result, review["evidence_ref"],
                        candidate_id, e3_type, review["seed"],
                        require_current=True,
                    )
        elif kind == "protein_hydrogen_review":
            require(_ref(data["evidence_ref"]) and _ref(data["source_ref"]),
                    "SCIENTIFIC_PROTEIN_H_REVIEW_REFS")
            require(isinstance(data["rationale"], str) and data["rationale"].strip(),
                    "SCIENTIFIC_PROTEIN_H_REVIEW_RATIONALE")
            row = db.execute("""SELECT id,ref,binding FROM scientific_evidence
                WHERE project=? AND job_id=? AND kind='protein_h' AND superseded=0
                ORDER BY rowid DESC LIMIT 1""",
                (self.project, job_id)).fetchone()
            require(row is not None and json.loads(row["ref"]) == data["evidence_ref"],
                    "SCIENTIFIC_PROTEIN_H_REVIEW_EVIDENCE")
            stored = json.loads(row["binding"])
            result_ref = result.get("files", {}).get("json", {})
            require(stored.get("project") == self.project and
                    stored.get("job_id") == job_id and
                    stored.get("result_sha256") == result_ref.get("sha256"),
                    "SCIENTIFIC_PROTEIN_H_REVIEW_BINDING")
            evidence = self.port.json(data["evidence_ref"])
            require(_digest(evidence) == stored.get("evidence_sha256"),
                    "SCIENTIFIC_PROTEIN_H_EVIDENCE_CHANGED")

    def record_decision(self, assessment_id, submission, session_token):
        self._validate(submission, "decision")
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            actor = self._authenticate(session_token, db)
            assessment = self._assessment(db, assessment_id)
            assigned = assessment.get("assigned_reviewer_id")
            require(assigned is None or assigned == actor["id"],
                    "SCIENTIFIC_ASSIGNED_REVIEWER")
            job_id = assessment["job_id"]
            _, result, _, _ = self._verified_result(db, job_id, require_current=True)
            fingerprint = _digest({
                "assessment_id": assessment_id,
                "actor": actor["id"],
                "submission": submission,
            })
            prior = db.execute("""SELECT body FROM scientific_decisions
                WHERE project=? AND job_id=? AND fingerprint=?""",
                (self.project, job_id, fingerprint)).fetchone()
            if prior:
                return json.loads(prior["body"]), False
            latest = self._latest_policy_revision(db, job_id)
            require(submission["expected_policy_revision"] == latest,
                    "SCIENTIFIC_STALE_POLICY")
            kind = submission["kind"]
            if kind in POLICY_KINDS:
                require(assessment["policy_revision"] == latest,
                        "SCIENTIFIC_STALE_ASSESSMENT")
                self._verify_refs(db, submission, actor["id"])
                self._validate_policy_data(
                    db, job_id, kind, submission["data"], result, actor["id"]
                )
                decision_id = "decision-" + uuid.uuid4().hex
                revision = self._next_policy_revision(db, job_id, "decision", decision_id)
            else:
                decision_id = "decision-" + uuid.uuid4().hex
                revision = latest
            if kind == "formal_decision" and submission["action"] == "accept":
                effective = self._effective(db, assessment)
                require(effective["current_policy"] and effective["source_inputs_current"]
                        and effective["source_runtime_current"]
                        and effective["current_implementation"]
                        and effective.get("protocol_evidence_current", True),
                        "SCIENTIFIC_STALE_ASSESSMENT")
                policy = self._policy(db, job_id, result)
                require(bool(policy["active_policy_decision_ids"]),
                        "SCIENTIFIC_ACTIVE_POLICY_REQUIRED")
                failures = [item["id"] for item in effective["criteria"]
                            if item["id"] != "formal_expert_decision"
                            and item["status"] != "pass"]
                require(not failures, "SCIENTIFIC_ACCEPTANCE_GATES")
            body = {
                "id": decision_id,
                "project_id": self.project,
                "job_id": job_id,
                "assessment_id": assessment_id,
                "policy_revision": revision,
                "kind": kind,
                "action": submission["action"],
                "intent": submission["intent"],
                "reason": submission["reason"],
                "data": submission.get("data", {}),
                "actor": actor,
                "created_at": now(),
                "withdrawn": False,
            }
            db.execute("""INSERT INTO scientific_decisions(
                id,project,job_id,assessment_id,policy_revision,kind,actor_id,
                withdrawn,fingerprint,body) VALUES(?,?,?,?,?,?,?,?,?,?)""",
                       (body["id"], self.project, job_id, assessment_id, revision,
                        kind, actor["id"], 0, fingerprint, json.dumps(body)))
            return body, True

    def withdraw(self, decision_id, session_token, expected_revision=None):
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            actor = self._authenticate(session_token, db)
            row = db.execute("""SELECT actor_id,body,withdrawn,kind,job_id
                FROM scientific_decisions WHERE id=? AND project=?""",
                (decision_id, self.project)).fetchone()
            if row is None:
                raise KeyError(decision_id)
            require(row["actor_id"] == actor["id"], "SCIENTIFIC_WITHDRAW_AUTHOR")
            body = json.loads(row["body"])
            if row["withdrawn"]:
                return body
            latest = self._latest_policy_revision(db, row["job_id"])
            if expected_revision is not None:
                require(expected_revision == latest, "SCIENTIFIC_STALE_POLICY")
            if row["kind"] in POLICY_KINDS:
                body["withdrawal_policy_revision"] = self._next_policy_revision(
                    db, row["job_id"], "withdrawal", decision_id)
            body.update(withdrawn=True, withdrawn_at=now())
            db.execute("UPDATE scientific_decisions SET withdrawn=1,body=? WHERE id=? AND project=?",
                       (json.dumps(body), decision_id, self.project))
            return body

    def withdraw_decision(self, decision_id, session_token, expected_revision=None):
        return self.withdraw(decision_id, session_token, expected_revision)

    def _structured_statement_document(self, text: str) -> dict:
        require(len(text.encode("utf-8")) <= MAX_STRUCTURED_STATEMENT_BYTES,
                "SCIENTIFIC_STRUCTURED_STATEMENT_SIZE")

        def pairs(items):
            value = {}
            for key, item in items:
                require(key not in value, "SCIENTIFIC_STRUCTURED_STATEMENT_DUPLICATE_KEY")
                value[key] = item
            return value

        def nonfinite(_value):
            raise ContractError("SCIENTIFIC_STRUCTURED_STATEMENT_NONFINITE")

        def finite_float(value):
            number = float(value)
            require(math.isfinite(number), "SCIENTIFIC_STRUCTURED_STATEMENT_NONFINITE")
            return number

        try:
            document = json.loads(text, object_pairs_hook=pairs,
                                  parse_constant=nonfinite, parse_float=finite_float)
        except ContractError:
            raise
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ContractError("SCIENTIFIC_STRUCTURED_STATEMENT_JSON") from exc
        require(isinstance(document, dict), "SCIENTIFIC_STRUCTURED_STATEMENT_OBJECT")
        return document

    def create_statement(self, text, locator, session_token):
        submission = {"text": text, "locator": locator}
        if (isinstance(text, str) and isinstance(locator, str)
                and locator.startswith("/document/")):
            require(len(text.encode("utf-8")) <= MAX_STRUCTURED_STATEMENT_BYTES,
                    "SCIENTIFIC_STRUCTURED_STATEMENT_SIZE")
        self._validate(submission, "statement")
        paths = []
        try:
            with self.store.db() as db:
                db.execute("BEGIN IMMEDIATE")
                actor = self._authenticate(session_token, db)
                value = {
                    "type": "human_statement",
                    "text": text,
                    "locator": locator,
                    "author": actor,
                    "created_at": now(),
                    "authority": "authenticated local reviewer statement; not signed PKI",
                }
                if locator.startswith("/document/"):
                    document = self._structured_statement_document(text)
                    self._verify_refs(db, document, actor["id"])
                    selected = self._resolve_locator({"document": document}, locator)
                    require(isinstance(selected, dict),
                            "SCIENTIFIC_STRUCTURED_STATEMENT_LOCATOR_TARGET")
                    value.update({
                        "document": document,
                        "structured_evidence": True,
                        "evidence_provenance": "authenticated human-supplied",
                        "source_character": "author-provided route proposal or documented experiment; not independently verified",
                        "institutional_or_laboratory_proof": False,
                        "scientific_approved": False,
                        "route_approved": False,
                        "synthesized": False,
                    })
                ref = self._write(db, paths, value, provenance="source")
                identifier = "statement-" + uuid.uuid4().hex
                db.execute("""INSERT INTO scientific_statements(
                    id,project,actor_id,ref,created_at) VALUES(?,?,?,?,?)""",
                           (identifier, self.project, actor["id"],
                            json.dumps(ref, sort_keys=True), value["created_at"]))
                return {"id": identifier, "ref": ref, **value}
        except BaseException:
            for path in paths:
                path.unlink(missing_ok=True)
            raise

    def _register_evidence(self, db, paths, job_id, kind, identity_key,
                           binding, value, *, provenance="computed"):
        fingerprint = _digest({
            "kind": kind,
            "identity": identity_key,
            "binding": binding,
            "value": value,
            "engine": engine_fingerprint(),
        })
        prior = db.execute("""SELECT id,ref,binding FROM scientific_evidence
            WHERE project=? AND job_id=? AND fingerprint=?""",
            (self.project, job_id, fingerprint)).fetchone()
        if prior:
            return {
                "id": prior["id"], "kind": kind,
                "ref": json.loads(prior["ref"]),
                "binding": json.loads(prior["binding"]), "created": False,
            }
        ref = self._write(db, paths, value, provenance=provenance)
        stored_binding = {**binding, "evidence_sha256": _digest(value)}
        identifier = "evidence-" + uuid.uuid4().hex
        db.execute("""UPDATE scientific_evidence SET superseded=1
            WHERE project=? AND job_id=? AND kind=? AND identity_key=? AND superseded=0""",
            (self.project, job_id, kind, identity_key))
        db.execute("""INSERT INTO scientific_evidence(
            id,project,job_id,kind,identity_key,fingerprint,binding,ref,
            superseded,created_at) VALUES(?,?,?,?,?,?,?,?,?,?)""",
                   (identifier, self.project, job_id, kind, identity_key, fingerprint,
                    json.dumps(stored_binding), json.dumps(ref), 0, now()))
        return {"id": identifier, "kind": kind, "ref": ref,
                "binding": stored_binding, "created": True}

    def _protein_source_paths(self) -> tuple[Path, Path]:
        from packages.science.dual_e3 import SOURCE
        pdb = SOURCE / "SMARCA2-receptor.pdb"
        cif = SOURCE / "6HAZ.cif"
        require(pdb.is_file() and cif.is_file() and not pdb.is_symlink() and not cif.is_symlink(),
                "SCIENTIFIC_PROTEIN_H_REGISTERED_SOURCES")
        return pdb.resolve(), cif.resolve()

    def _protein_h_for_compute(self, db, job_id: str, binding: dict,
                               policy: dict) -> dict | None:
        row = db.execute("""SELECT id,binding,ref FROM scientific_evidence
            WHERE project=? AND job_id=? AND kind='protein_h' AND superseded=0
            ORDER BY rowid DESC LIMIT 1""", (self.project, job_id)).fetchone()
        if row is None:
            return None
        stored = json.loads(row["binding"])
        require(stored.get("result_sha256") == binding["result_sha256"] and
                stored.get("input_sha256") == binding["input_sha256"],
                "SCIENTIFIC_PROTEIN_H_JOB_BINDING")
        ref = json.loads(row["ref"])
        value = self.port.json(ref)
        require(_digest(value) == stored.get("evidence_sha256"),
                "SCIENTIFIC_PROTEIN_H_EVIDENCE_CHANGED")
        self._verify_refs(db, value)
        pdb, cif = self._protein_source_paths()
        source_hashes = stored.get("source_sha256", {})
        require(_sha(pdb.read_bytes()) == source_hashes.get("pdb") and
                _sha(cif.read_bytes()) == source_hashes.get("cif"),
                "SCIENTIFIC_PROTEIN_H_REGISTERED_SOURCE_CHANGED")
        flags = value.get("state_flags", {})
        require(value.get("status") == "computed_diagnostic_pending_human_review" and
                flags.get("computed_diagnostic") is True and
                flags.get("formal_scientific_approval") is False,
                "SCIENTIFIC_PROTEIN_H_DIAGNOSTIC_INVALID")
        hydrogens = value.get("protein_hydrogens")
        require(isinstance(hydrogens, list), "SCIENTIFIC_PROTEIN_H_ARRAY")
        review = policy.get("decisions", {}).get("protein_hydrogen_review")
        reviewed = bool(isinstance(review, dict)
                        and review.get("evidence_ref") == ref
                        and isinstance(review.get("decision_id"), str)
                        and review.get("decision_id")
                        and isinstance(review.get("actor"), dict)
                        and isinstance(review["actor"].get("id"), str))
        accepted = bool(reviewed and review.get("accepted") is True)
        authenticated_binding = None
        if reviewed:
            authenticated_binding = {
                "platform_authenticated": True,
                "project_id": self.project,
                "job_id": job_id,
                "evidence_id": row["id"],
                "evidence_ref": _copy(ref),
                "decision_id": review["decision_id"],
                "actor_id": review["actor"]["id"],
                "disposition": "accepted" if accepted else "rejected",
            }
        receipt = {
            "requires_review": True,
            "evidence_ref": ref,
            "evidence_id": row["id"],
            "original_state_flags": _copy(flags),
            "effective_human_review_status": (
                "accepted" if accepted else "rejected" if reviewed else "pending"
            ),
            "human_review": _copy(review) if reviewed else None,
            "platform_authenticated_human_review_binding": authenticated_binding,
            "optimization_auto_approved": False,
            "diagnostic_is_human_approval": False,
        }
        return {"id": row["id"], "ref": ref, "hydrogens": _copy(hydrogens),
                "receipt": receipt}

    def import_protein_hydrogens(self, job_id, evidence_json_path, session_token):
        evidence_path = Path(evidence_json_path).expanduser().resolve()
        require(evidence_path.is_file() and not evidence_path.is_symlink(),
                "SCIENTIFIC_PROTEIN_H_EVIDENCE_FILE")
        artifact_root = evidence_path.parent
        require(evidence_path == artifact_root / "protein_hydrogen_evidence.json",
                "SCIENTIFIC_PROTEIN_H_EVIDENCE_FILENAME")
        submitted_raw = evidence_path.read_bytes()
        submitted = parse_json(submitted_raw)
        require(isinstance(submitted, dict), "SCIENTIFIC_PROTEIN_H_EVIDENCE_JSON")
        receipt = submitted.get("receipt")
        require(isinstance(receipt, dict), "SCIENTIFIC_PROTEIN_H_RECEIPT")
        input_path = Path(receipt.get("input_path", "")).expanduser().resolve()
        cif_candidates = [
            path for path in artifact_root.rglob("*")
            if path.is_file() and not path.is_symlink()
            and path.suffix.lower() in {".cif", ".mmcif"}
            and _sha(path.read_bytes()) == submitted.get("sha256", {}).get("source_cif_file")
        ]
        require(len(cif_candidates) == 1, "SCIENTIFIC_PROTEIN_H_SOURCE_CIF")
        validated = ProteinHydrogenEvidence(
            input_path, cif_candidates[0], receipt, artifact_root
        ).validate()
        require(encoded(validated) == encoded(submitted),
                "SCIENTIFIC_PROTEIN_H_EVIDENCE_NOT_REPRODUCIBLE")
        source_pdb, source_cif = self._protein_source_paths()
        require(_sha(source_pdb.read_bytes()) == validated["sha256"]["source_pdb_file"] and
                _sha(source_cif.read_bytes()) == validated["sha256"]["source_cif_file"],
                "SCIENTIFIC_PROTEIN_H_SOURCE_HASH")

        manifest = validated.get("artifact_manifest")
        require(isinstance(manifest, list) and manifest,
                "SCIENTIFIC_PROTEIN_H_ARTIFACT_MANIFEST")
        manifest_paths = {
            item.get("path") for item in manifest if isinstance(item, dict)
        }
        require(len(manifest_paths) == len(manifest) and
                all(isinstance(item, str) and item for item in manifest_paths),
                "SCIENTIFIC_PROTEIN_H_ARTIFACT_MANIFEST")
        actual_paths = {
            path.relative_to(artifact_root).as_posix()
            for path in artifact_root.rglob("*")
            if path.is_file() and not path.is_symlink()
            and path != evidence_path
        }
        require(actual_paths == manifest_paths,
                "SCIENTIFIC_PROTEIN_H_ARTIFACT_TREE_MISMATCH")
        paths = []
        try:
            with self.store.db() as db:
                db.execute("BEGIN IMMEDIATE")
                actor = self._authenticate(session_token, db)
                _, result, binding, _ = self._verified_result(db, job_id, require_current=True)
                submitted_evidence_sha256 = _sha(submitted_raw)
                prior = db.execute("""SELECT id,ref,binding FROM scientific_evidence
                    WHERE project=? AND job_id=? AND kind='protein_h'
                    AND identity_key=? AND superseded=0 ORDER BY rowid DESC LIMIT 1""",
                    (self.project, job_id, binding["parent_id"])).fetchone()
                if prior is not None:
                    prior_binding = json.loads(prior["binding"])
                    if (prior_binding.get("submitted_evidence_sha256") ==
                            submitted_evidence_sha256):
                        return {
                            "id": prior["id"],
                            "kind": "protein_h",
                            "ref": json.loads(prior["ref"]),
                            "binding": prior_binding,
                            "created": False,
                        }
                registered = {}
                for item in manifest:
                    relative = item.get("path") if isinstance(item, dict) else None
                    require(isinstance(relative, str) and relative and not Path(relative).is_absolute(),
                            "SCIENTIFIC_PROTEIN_H_PORTABLE_PATH")
                    path = (artifact_root / relative).resolve(strict=True)
                    require(artifact_root == path.parent or artifact_root in path.parents,
                            "SCIENTIFIC_PROTEIN_H_PATH_ESCAPE")
                    raw = path.read_bytes()
                    require(_sha(raw) == item.get("sha256"),
                            "SCIENTIFIC_PROTEIN_H_ARTIFACT_HASH")
                    suffix = path.suffix.lower()
                    media = ("chemical/x-pdb" if suffix == ".pdb" else
                             "chemical/x-pqr" if suffix == ".pqr" else
                             "chemical/x-mmcif" if suffix in {".cif", ".mmcif"} else
                             "application/json" if suffix == ".json" else
                             "text/plain" if suffix in {".log", ".txt", ".out", ".err"} else
                             "application/octet-stream")
                    provenance = "source" if relative.startswith("source/") else "computed"
                    registered[relative] = self._write(db, paths, raw, media, provenance)
                submitted_ref = self._write(db, paths, submitted_raw, "application/json", "source")
                portable = _copy(validated)
                path_fields = ("input_path", "output_pqr_path", "output_pdb_path",
                               "stdout_path", "stderr_path")
                for field in path_fields:
                    raw_path = portable["receipt"].get(field)
                    if isinstance(raw_path, str):
                        resolved = Path(raw_path).expanduser().resolve()
                        portable["receipt"][field] = resolved.relative_to(artifact_root).as_posix()
                portable["receipt"]["executable_path"] = None

                def portable_path(raw_path):
                    if not isinstance(raw_path, str) or not raw_path:
                        return raw_path
                    candidate = Path(raw_path).expanduser()
                    if not candidate.is_absolute():
                        normalized = candidate.as_posix()
                        if "/" in raw_path or "\\" in raw_path:
                            require(normalized in manifest_paths,
                                    "SCIENTIFIC_PROTEIN_H_COMMAND_PATH")
                        return normalized if normalized in manifest_paths else raw_path
                    resolved = candidate.resolve(strict=True)
                    require(resolved.is_file() and not resolved.is_symlink(),
                            "SCIENTIFIC_PROTEIN_H_COMMAND_PATH")
                    try:
                        relative = resolved.relative_to(artifact_root).as_posix()
                    except ValueError as error:
                        raise ContractError(
                            "SCIENTIFIC_PROTEIN_H_COMMAND_PATH"
                        ) from error
                    require(relative in manifest_paths,
                            "SCIENTIFIC_PROTEIN_H_COMMAND_PATH")
                    return relative

                command = portable["receipt"].get("command")
                if isinstance(command, list) and command:
                    redacted = [Path(command[0]).name]
                    for argument in command[1:]:
                        if isinstance(argument, str) and argument.startswith("--") and "=" in argument:
                            option, raw_value = argument.split("=", 1)
                            redacted.append(option + "=" + str(portable_path(raw_value)))
                        else:
                            redacted.append(portable_path(argument))
                    portable["receipt"]["command"] = redacted

                entry_groups = []
                for key in ("actual_pka_entries", "actual_pKa_entries"):
                    entries = portable.get(key)
                    if isinstance(entries, list):
                        entry_groups.append(entries)
                for key in ("pKa_entries_from_actual_output",
                            "pka_entries_from_actual_output"):
                    entries = portable["receipt"].get(key)
                    if isinstance(entries, list):
                        entry_groups.append(entries)
                for entries in entry_groups:
                    for entry in entries:
                        if isinstance(entry, dict) and isinstance(entry.get("source_path"), str):
                            converted = portable_path(entry["source_path"])
                            entry["source_path"] = (
                                converted["raw_reference"]
                                if isinstance(converted, dict) else converted
                            )
                portable["registered_artifacts"] = registered
                portable["submitted_evidence_ref"] = submitted_ref
                portable["imported_by"] = actor
                evidence_binding = {
                    **binding,
                    "submitted_evidence_sha256": submitted_evidence_sha256,
                    "source_sha256": {
                        "pdb": _sha(source_pdb.read_bytes()),
                        "cif": _sha(source_cif.read_bytes()),
                    },
                    "measurement_binding": {
                        "kind": "validated_fixed_frame_protein_hydrogen_diagnostic",
                        "module_version": validated.get("format"),
                        "module_fingerprint": validated.get("fingerprint"),
                    },
                }
                return self._register_evidence(
                    db, paths, job_id, "protein_h", binding["parent_id"],
                    evidence_binding, portable
                )
        except BaseException:
            for path in paths:
                path.unlink(missing_ok=True)
            raise

    def _snapshot_for_compute(self, job_id, session_token, *,
                              allow_archived_runtime: bool = False):
        with self.store.db() as db:
            actor = self._authenticate(session_token, db)
            _, result, binding, source = self._verified_reinspection_snapshot(
                db, job_id,
                allow_archived_runtime=allow_archived_runtime,
            )
            policy = self._policy(db, job_id, result)
            protein_h = self._protein_h_for_compute(db, job_id, binding, policy)
            provenance = self._reinspection_report_provenance(source, actor)
            return actor, _copy(result), binding, policy, protein_h, provenance

    @staticmethod
    def _mapless_isomeric_smiles(mol) -> str:
        value = Chem.RemoveHs(Chem.Mol(mol))
        for atom in value.GetAtoms():
            atom.SetAtomMapNum(0)
        return Chem.MolToSmiles(value, canonical=True, isomericSmiles=True)

    @staticmethod
    def _pdbqt_poses(raw: bytes) -> list:
        try:
            import re
            from meeko import PDBQTMolecule, RDKitMolCreate

            text = raw.decode("utf-8", errors="strict")
            blocks = []
            current = []
            model_id = None
            for line in text.splitlines(keepends=True):
                record = line.rstrip("\r\n")
                model_match = re.fullmatch(
                    r"[ \t]*MODEL[ \t]+([0-9]+)[ \t]*", record
                )
                end_match = re.fullmatch(r"[ \t]*ENDMDL[ \t]*", record)
                if model_id is None:
                    if not record.strip():
                        continue
                    if model_match is None:
                        raise ValueError("PDBQT_MODEL_RECORD")
                    model_id = int(model_match.group(1))
                    if model_id != len(blocks) + 1:
                        raise ValueError("PDBQT_MODEL_SEQUENCE")
                    current = [line]
                    continue

                if model_match is not None or record.lstrip().startswith("MODEL"):
                    raise ValueError("PDBQT_NESTED_MODEL")
                current.append(line)
                if end_match is not None:
                    block = "".join(current)
                    if not block.endswith(("\n", "\r")):
                        block += "\n"
                    blocks.append(block)
                    current = []
                    model_id = None
                elif record.lstrip().startswith("ENDMDL"):
                    raise ValueError("PDBQT_ENDMDL_RECORD")

            if model_id is not None or current or not blocks:
                raise ValueError("PDBQT_MODEL_PAIRING")

            restored = []
            for block in blocks:
                converted = RDKitMolCreate.from_pdbqt_mol(
                    PDBQTMolecule(block, skip_typing=True)
                )
                if (not isinstance(converted, list) or len(converted) != 1
                        or converted[0] is None
                        or converted[0].GetNumConformers() != 1):
                    raise ValueError("PDBQT_MODEL_RESTORATION")
                restored.append(converted[0])
            return restored
        except Exception as error:
            raise ContractError("SCIENTIFIC_INTERACTION_PDBQT_RESTORATION") from error

    def _repair_actual_pose_maps(self, raw_pose, restored_pose, expected):
        require(raw_pose is not None and raw_pose.GetNumConformers() == 1
                and restored_pose is not None
                and restored_pose.GetNumConformers() == 1,
                "SCIENTIFIC_INTERACTION_POSE_FRAME")
        raw = Chem.RemoveHs(Chem.Mol(raw_pose))
        restored = Chem.RemoveHs(Chem.Mol(restored_pose))
        expected = Chem.RemoveHs(Chem.Mol(expected))
        require(raw.GetNumAtoms() == restored.GetNumAtoms() == expected.GetNumAtoms(),
                "SCIENTIFIC_INTERACTION_POSE_GRAPH")

        expected_maps = [atom.GetAtomMapNum() for atom in expected.GetAtoms()]
        restored_maps = [atom.GetAtomMapNum() for atom in restored.GetAtoms()]
        require(all(type(value) is int and value > 0 for value in expected_maps)
                and len(expected_maps) == len(set(expected_maps))
                and all(type(value) is int and value > 0 for value in restored_maps)
                and len(restored_maps) == len(set(restored_maps))
                and set(restored_maps) == set(expected_maps),
                "SCIENTIFIC_INTERACTION_PDBQT_MAP_FRAME")
        require(Chem.MolToSmiles(restored, canonical=True, isomericSmiles=True) ==
                Chem.MolToSmiles(expected, canonical=True, isomericSmiles=True),
                "SCIENTIFIC_INTERACTION_PDBQT_GRAPH")
        require(self._mapless_isomeric_smiles(raw) ==
                self._mapless_isomeric_smiles(restored),
                "SCIENTIFIC_INTERACTION_POSE_GRAPH")

        raw_query = Chem.Mol(raw)
        restored_target = Chem.Mol(restored)
        for molecule in (raw_query, restored_target):
            for atom in molecule.GetAtoms():
                atom.SetAtomMapNum(0)
        matches = restored_target.GetSubstructMatches(
            raw_query, uniquify=False, useChirality=True
        )
        raw_conf = raw.GetConformer()
        restored_conf = restored.GetConformer()
        accepted = []
        tolerance = 1.0e-4
        for match in matches:
            if len(match) != raw.GetNumAtoms():
                continue
            mapping = []
            valid = True
            for raw_index, restored_index in enumerate(match):
                raw_position = raw_conf.GetAtomPosition(raw_index)
                restored_position = restored_conf.GetAtomPosition(restored_index)
                if raw_position.Distance(restored_position) > tolerance:
                    valid = False
                    break
                raw_map = raw.GetAtomWithIdx(raw_index).GetAtomMapNum()
                full_map = restored.GetAtomWithIdx(restored_index).GetAtomMapNum()
                truncated = int(str(full_map)[:3]) if full_map >= 1000 else None
                if raw_map not in (0, full_map) and raw_map != truncated:
                    valid = False
                    break
                mapping.append((raw_index, full_map, raw_map))
            if valid:
                accepted.append(mapping)
        require(len(accepted) == 1,
                "SCIENTIFIC_INTERACTION_MAP_RESTORATION_AMBIGUOUS")

        normalized = Chem.Mol(raw_pose)
        heavy_indices = [atom.GetIdx() for atom in normalized.GetAtoms()
                         if atom.GetAtomicNum() > 1]
        require(len(heavy_indices) == len(accepted[0]),
                "SCIENTIFIC_INTERACTION_POSE_GRAPH")
        repaired = []
        for (heavy_index, full_map, raw_map), atom_index in zip(
                accepted[0], heavy_indices):
            require(heavy_index == len(repaired),
                    "SCIENTIFIC_INTERACTION_MAP_RESTORATION_AMBIGUOUS")
            normalized.GetAtomWithIdx(atom_index).SetAtomMapNum(full_map)
            repaired.append({
                "heavy_atom_index_zero_based": heavy_index,
                "raw_sdf_atom_map": raw_map,
                "restored_atom_map": full_map,
                "restoration": (
                    "missing" if raw_map == 0 else
                    "v2000_first_three_digits" if raw_map != full_map else
                    "unchanged"
                ),
            })
        require(Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(normalized)),
                                 canonical=True, isomericSmiles=True) ==
                Chem.MolToSmiles(expected, canonical=True,
                                 isomericSmiles=True),
                "SCIENTIFIC_INTERACTION_POSE_GRAPH")
        return normalized, {
            "status": "explicit_atom_map_restoration",
            "scientific_approval": False,
            "coordinate_preservation": {
                "alignment_applied": False,
                "optimization_applied": False,
                "all_heavy_atom_coordinates_compared": True,
                "maximum_allowed_distance_angstrom": tolerance,
            },
            "full_atom_map_mapping": repaired,
        }

    def _actual_passing_pose(self, diagnostic_source: Any, poses: list,
                             analog: dict, pdbqt_loader=None):
        require(isinstance(diagnostic_source, dict),
                "SCIENTIFIC_INTERACTION_DIAGNOSTIC_SOURCE")
        require(diagnostic_source.get("real_docking") is True,
                "SCIENTIFIC_INTERACTION_NOT_ACTUAL_DOCKING")
        results = diagnostic_source.get("results")
        require(isinstance(results, dict),
                "SCIENTIFIC_INTERACTION_DIAGNOSTIC_RESULTS")
        preservation = results.get("pose_preservation")
        require(isinstance(preservation, dict)
                and preservation.get("status") == "pass"
                and preservation.get("docking_pose_preserved") is True,
                "SCIENTIFIC_INTERACTION_PRESERVATION_STATUS")
        diagnostics = preservation.get("all_poses_diagnostics")
        require(isinstance(diagnostics, list),
                "SCIENTIFIC_INTERACTION_DIAGNOSTICS")
        require(type(results.get("pose_count")) is int
                and results["pose_count"] == len(poses),
                "SCIENTIFIC_INTERACTION_POSE_COUNT_MISMATCH")

        declared_poses = results.get("poses")
        require(isinstance(declared_poses, list)
                and len(declared_poses) == len(poses),
                "SCIENTIFIC_INTERACTION_POSE_MAPPING")

        chosen = None
        for row in diagnostics:
            if not isinstance(row, dict) or row.get("status") != "pass":
                continue
            require(row.get("docking_pose_preserved") is True,
                    "SCIENTIFIC_INTERACTION_PASS_MAPPING")
            require(row.get("alignment_applied") is False,
                    "SCIENTIFIC_INTERACTION_ALIGNED_POSE_REJECTED")
            require(row.get("pose_origin") ==
                    "caller supplied and not inferred or verified",
                    "SCIENTIFIC_INTERACTION_POSE_CALLER")
            index = row.get("pose_index_zero_based")
            require(type(index) is int and index >= 0,
                    "SCIENTIFIC_INTERACTION_POSE_INDEX")
            require(index < len(poses) and poses[index] is not None,
                    "SCIENTIFIC_DECLARED_POSE_MISSING")
            chosen = (index, row)
            break
        require(chosen is not None, "SCIENTIFIC_INTERACTION_NO_PASSING_POSE")
        index, row = chosen

        # Docking mode is separate producer metadata. It is never used to select
        # an SDF record; this consistency check only audits the explicit adapter
        # between one-based modes and immutable zero-based SDF record indices.
        declared = declared_poses[index]
        require(isinstance(declared, dict)
                and type(declared.get("mode")) is int
                and declared["mode"] == index + 1
                and declared.get("mapping_ambiguous") is False
                and type(declared.get("mapping_count")) is int
                and declared["mapping_count"] == 1
                and isinstance(declared.get("selected_mapping"), str)
                and bool(declared["selected_mapping"]),
                "SCIENTIFIC_INTERACTION_POSE_MAPPING")

        pose = poses[index]
        require(pose.GetNumConformers() == 1,
                "SCIENTIFIC_INTERACTION_POSE_FRAME")
        expected = Chem.MolFromSmiles(analog.get("mapped_smiles", ""))
        require(expected is not None, "SCIENTIFIC_INTERACTION_ANALOG_GRAPH")
        expected = Chem.RemoveHs(expected)
        actual = Chem.RemoveHs(Chem.Mol(pose))
        expected_maps = [atom.GetAtomMapNum() for atom in expected.GetAtoms()
                         if atom.GetAtomicNum() > 1]
        actual_maps = [atom.GetAtomMapNum() for atom in actual.GetAtoms()
                       if atom.GetAtomicNum() > 1]
        require(all(type(value) is int and value > 0 for value in expected_maps)
                and len(expected_maps) == len(set(expected_maps)),
                "SCIENTIFIC_INTERACTION_ANALOG_GRAPH")
        maps_valid = (
            all(type(value) is int and value > 0 for value in actual_maps)
            and len(actual_maps) == len(set(actual_maps))
            and set(actual_maps) == set(expected_maps)
        )
        map_repair_receipt = None
        if not maps_valid:
            require(callable(pdbqt_loader),
                    "SCIENTIFIC_INTERACTION_PDBQT_SOURCE_REQUIRED")
            pdbqt_raw = pdbqt_loader()
            require(isinstance(pdbqt_raw, bytes) and bool(pdbqt_raw),
                    "SCIENTIFIC_INTERACTION_PDBQT_SOURCE_REQUIRED")
            restored_poses = self._pdbqt_poses(pdbqt_raw)
            require(len(restored_poses) == len(poses)
                    and all(value is not None for value in restored_poses),
                    "SCIENTIFIC_INTERACTION_PDBQT_POSE_COUNT_MISMATCH")
            pose, map_repair_receipt = self._repair_actual_pose_maps(
                pose, restored_poses[index], expected
            )
            actual = Chem.RemoveHs(Chem.Mol(pose))
            actual_maps = [atom.GetAtomMapNum() for atom in actual.GetAtoms()
                           if atom.GetAtomicNum() > 1]
        require(all(type(value) is int and value > 0 for value in actual_maps)
                and len(actual_maps) == len(set(actual_maps))
                and set(actual_maps) == set(expected_maps),
                "SCIENTIFIC_INTERACTION_POSE_MAP_FRAME")
        require(Chem.MolToSmiles(actual, canonical=True, isomericSmiles=True) ==
                Chem.MolToSmiles(expected, canonical=True, isomericSmiles=True),
                "SCIENTIFIC_INTERACTION_POSE_GRAPH")
        common_maps = row.get("common_atom_maps")
        require(isinstance(common_maps, list)
                and all(type(value) is int and value > 0 for value in common_maps)
                and set(common_maps) <= set(actual_maps),
                "SCIENTIFIC_INTERACTION_DIAGNOSTIC_MAPS")
        return index, row, declared, pose, map_repair_receipt

    def _replay_actual_ternary(self, db, paths, job_id: str, result: dict,
                               binding: dict, policy: dict,
                               report_provenance: dict) -> list:
        """Reinterpret registered actual measurements under the current policy.

        This operates only on immutable, platform-registered import records. It
        neither reads producer paths nor reruns or estimates a GPU measurement.
        """
        policy_digest = self._policy_digest(policy)
        geometry = policy.get("decisions", {}).get(
            "ternary_geometry_criterion"
        )
        candidates = self._candidate_map(result)
        rows = db.execute("""SELECT id,identity_key,binding,ref FROM scientific_evidence
            WHERE project=? AND job_id=? AND kind='ternary' ORDER BY rowid""",
            (self.project, job_id)).fetchall()
        authoritative = {}
        from scripts.run_novel_ternary import _candidate as qualified_candidate

        for row in rows:
            stored = json.loads(row["binding"])
            ref = json.loads(row["ref"])
            require(_ref(ref), "SCIENTIFIC_TERNARY_REPLAY_SOURCE_REF")
            value = self.port.json(ref)
            require(_digest(value) == stored.get("evidence_sha256"),
                    "SCIENTIFIC_TERNARY_REPLAY_SOURCE_HASH")
            self._verify_refs(db, value)
            require(stored.get("project", self.project) == self.project and
                    stored.get("job_id", job_id) == job_id and
                    stored.get("parent_id") == binding.get("parent_id") and
                    stored.get("input_sha256") == binding.get("input_sha256") and
                    stored.get("result_sha256") == binding.get("result_sha256"),
                    "SCIENTIFIC_TERNARY_REPLAY_SOURCE_BINDING")
            if value.get("derived_interpretation") is True:
                continue

            measurement = stored.get("measurement_binding")
            protocol = value.get("measurement_protocol")
            producer = value.get("producer_policy_binding")
            source_evaluation = value.get("evaluation_policy")
            require(isinstance(measurement, dict) and
                    isinstance(protocol, dict) and
                    isinstance(producer, dict) and
                    isinstance(source_evaluation, dict) and
                    stored.get("evaluation_policy") == source_evaluation and
                    stored.get("evaluation_policy_digest") ==
                    source_evaluation.get("digest") and
                    measurement.get("producer_policy") == producer and
                    measurement.get("plan_digest") == protocol.get("plan_digest") and
                    measurement.get("candidate_graph_sha256") ==
                    protocol.get("candidate_graph_sha256"),
                    "SCIENTIFIC_TERNARY_REPLAY_METADATA")
            require(type(value.get("actual_computation")) is bool and
                    type(value.get("execution_success")) is bool and
                    isinstance(value.get("inspection"), dict) and
                    isinstance(value.get("inspection_summary"), dict) and
                    isinstance(value.get("registered_artifacts"), dict) and
                    bool(value["registered_artifacts"]) and
                    _ref(value.get("plan_ref")) and
                    _ref(value.get("aggregate_receipt_ref")),
                    "SCIENTIFIC_TERNARY_REPLAY_AUTHORITATIVE")

            candidate_id = value.get("candidate_id")
            e3_type = value.get("e3_type")
            seed = value.get("seed")
            require(isinstance(candidate_id, str) and candidate_id in candidates and
                    isinstance(e3_type, str) and e3_type and
                    type(seed) is int and protocol.get("seed") == seed,
                    "SCIENTIFIC_TERNARY_REPLAY_IDENTITY")
            candidate = candidates[candidate_id]
            require(candidate.get("e3_type") == e3_type,
                    "SCIENTIFIC_TERNARY_REPLAY_E3_BINDING")
            graph = novel_ternary._mapped_graph(
                qualified_candidate(result, candidate_id)
            )
            require(graph.get("graph_sha256") ==
                    protocol.get("candidate_graph_sha256"),
                    "SCIENTIFIC_TERNARY_REPLAY_GRAPH_BINDING")

            plan = self.port.json(value["plan_ref"])
            job_binding = plan.get("bindings", {}).get("job", {})
            plan_graph = plan.get("candidate_graph", {})
            require(plan.get("plan_digest") == protocol.get("plan_digest") and
                    job_binding.get("project") == self.project and
                    job_binding.get("job_id") == job_id and
                    job_binding.get("input_sha256") == binding.get("input_sha256") and
                    job_binding.get("result_sha256") == binding.get("result_sha256") and
                    job_binding.get("candidate_graph_sha256") == graph["graph_sha256"] and
                    plan_graph.get("candidate_id") == candidate_id and
                    plan_graph.get("e3_type") == e3_type and
                    plan_graph.get("graph_sha256") == graph["graph_sha256"],
                    "SCIENTIFIC_TERNARY_REPLAY_PLAN_BINDING")
            identity = f"{candidate_id}|{e3_type}|{seed}"
            require(row["identity_key"] == identity,
                    "SCIENTIFIC_TERNARY_REPLAY_IDENTITY")
            authoritative[identity] = (ref, stored, value)

        created = []
        for identity, (source_ref, source_binding, source) in authoritative.items():
            measured = self._flatten_metrics(
                source["inspection"].get("descriptive_metrics", {})
            )
            metric = geometry.get("metric") if isinstance(geometry, dict) else None
            threshold = geometry.get("threshold") if isinstance(geometry, dict) else None
            comparison = geometry.get("comparison") if isinstance(geometry, dict) else None
            observed = measured.get(metric) if isinstance(metric, str) else None
            valid_observed = (type(observed) in (int, float) and
                              math.isfinite(float(observed)))
            valid_threshold = (type(threshold) in (int, float) and
                               math.isfinite(float(threshold)))
            valid_criterion = (isinstance(metric, str) and bool(metric) and
                               valid_threshold and comparison in {
                                   "lte", "max_lte", "gte", "min_gte"
                               })
            success = (source.get("actual_computation") is True and
                       source.get("execution_success") is True and
                       source.get("failure") is None)
            compatible = None
            reason = None
            if not success:
                reason = "original registered execution did not succeed"
            elif not valid_criterion:
                reason = "current geometry criterion is absent or invalid"
            elif not valid_observed:
                reason = "current geometry metric is absent from measured output"
            elif comparison in {"lte", "max_lte"}:
                compatible = observed <= threshold
            else:
                compatible = observed >= threshold

            derived = _copy(source)
            derived.update({
                "derived_interpretation": True,
                "original_measurement_ref": _copy(source_ref),
                "original_inspection_summary": _copy(source["inspection_summary"]),
                "original_evaluation_policy": _copy(source.get("evaluation_policy")),
                "evaluation_policy": {
                    "revision": policy["revision"], "digest": policy_digest,
                },
                "inspection_summary": {
                    "compatibility": ("compatible" if compatible is True else
                                      "incompatible" if compatible is False else
                                      "unknown"),
                    "metrics": measured,
                    "criterion_metric": metric,
                    "criterion_observed": observed if valid_observed else None,
                    "criterion_threshold": threshold if valid_threshold else None,
                    "criterion_comparison": comparison if valid_criterion else None,
                },
                "policy_incompatibilities": ([] if reason is None else
                    list(dict.fromkeys([
                        *source.get("policy_incompatibilities", []), reason
                    ]))),
                "source_chain_refs": [
                    _copy(source_ref),
                    *[_copy(ref) for ref in source.get("source_chain_refs", [])
                      if _ref(ref)],
                ],
                "scientific_approved": False,
                "efficacy_approved": False,
                "hydrogen_approved": False,
            })
            derived["report_provenance"] = _copy(report_provenance)
            replay_binding = {
                **binding,
                "evaluation_policy_digest": policy_digest,
                "evaluation_policy": {
                    "revision": policy["revision"], "digest": policy_digest,
                },
                "measurement_binding": _copy(source_binding["measurement_binding"]),
                "interpretation_binding": {
                    "kind": "stored_actual_ternary_policy_re_evaluation",
                    "original_evidence_sha256": source_binding["evidence_sha256"],
                    "original_measurement_ref": _copy(source_ref),
                },
                "reinspection": _copy(report_provenance),
            }
            created.append(self._register_evidence(
                db, paths, job_id, "ternary", identity, replay_binding, derived
            ))
        return created

    def _absent_protein_hydrogen_receipt(self, report: dict) -> dict:
        pending = {
            "requires_review": True,
            "effective_human_review_status": "pending",
            "reason": "no_current_bound_validated_protein_hydrogen_evidence",
            "diagnostic_is_human_approval": False,
        }
        requirements = report.get("requirements")
        if not isinstance(requirements, list) or not requirements:
            return pending
        if any(not isinstance(row, dict) for row in requirements):
            return pending
        required = [row for row in requirements if row.get("required") is True]
        directional = [row for row in required
                       if row.get("kind") == "directional_hbond"]
        if not required or not directional:
            return pending

        def finite_measurement(value: Any) -> bool:
            return (type(value) in (int, float)
                    and math.isfinite(float(value)))

        def ligand_donor(match: Any, requirement: dict) -> bool:
            ligand_maps = requirement.get("ligand_maps")
            protein_ids = requirement.get("protein_atom_ids")
            if (not isinstance(ligand_maps, list) or not ligand_maps
                    or any(type(value) is not int or value <= 0
                           for value in ligand_maps)
                    or not isinstance(protein_ids, list) or not protein_ids
                    or any(not isinstance(value, str) or not value
                           for value in protein_ids)):
                return False
            participants = match.get("participants") if isinstance(match, dict) else None
            return (
                isinstance(match, dict)
                and match.get("kind") == "directional_hbond"
                and match.get("requires_review") is False
                and finite_measurement(match.get("distance_HA_A"))
                and finite_measurement(match.get("angle_DHA_deg"))
                and isinstance(participants, list)
                and len(participants) == 2
                and participants[0] in {f"L:{value}" for value in ligand_maps}
                and participants[1] in set(protein_ids)
            )

        resolved = all(
            row.get("status") == "preserved"
            and row.get("computed_pass") is True
            and row.get("pending_missing_protein_hydrogen") is not True
            and row.get("needs_expert") is not True
            for row in required
        )
        resolved = resolved and all(
            isinstance(row.get("before_matches"), list)
            and isinstance(row.get("after_matches"), list)
            and any(ligand_donor(match, row)
                    for match in row["before_matches"])
            and any(ligand_donor(match, row)
                    for match in row["after_matches"])
            for row in directional
        )
        resolved = resolved and not report.get("failures") and (
            report.get("protein_hydrogen_status") ==
            "not_required_or_not_implicated")
        if not resolved:
            return pending
        return {
            "requires_review": False,
            "effective_human_review_status": "not_required",
            "reason": "direction_specific_ligand_donor_to_protein_acceptor_checks_computed_pass_without_required_protein_donor_hydrogen",
            "diagnostic_is_human_approval": False,
        }

    def compute(self, job_id, session_token, candidate_ids=None, analog_ids=None,
                *, allow_archived_runtime: bool = False):
        request = {"candidate_ids": candidate_ids, "analog_ids": analog_ids}
        self._validate(request, "compute")
        actor, result, binding, policy, protein_h, report_provenance = (
            self._snapshot_for_compute(
                job_id, session_token,
                allow_archived_runtime=allow_archived_runtime,
            )
        )
        policy_digest = self._policy_digest(policy)
        candidates = self._candidate_map(result)
        chosen = policy.get("decisions", {}).get("synthesis_criterion", {}).get(
            "chosen_candidate_ids")
        expected_candidates = sorted(chosen if isinstance(chosen, list) and chosen else candidates)
        if candidate_ids is not None:
            require(set(candidate_ids) == set(expected_candidates),
                    "SCIENTIFIC_SYNTHESIS_COVERAGE_REQUIRED")
        expected_analogs = sorted(
            row["id"] for row in result.get("analogs", [])
            if isinstance(row, dict) and row.get("selected") is True
            and row.get("pipeline_status") == "qualified"
        )
        if analog_ids is not None:
            require(set(analog_ids) == set(expected_analogs),
                    "SCIENTIFIC_INTERACTION_COVERAGE_REQUIRED")

        computed = []
        resolved = policy.get("_resolved_synthesis", {})
        for candidate_id in expected_candidates:
            candidate = candidates[candidate_id]
            try:
                value = synthesis_review.assess_synthesis(
                    candidate,
                    exact_route_records=resolved.get("exact_route_records", {}).get(candidate_id, []),
                    constraints=resolved.get("constraints", {}).get(candidate_id),
                )
            except Exception as error:
                value = {
                    "format": synthesis_review.VERSION,
                    "candidate_id": candidate_id,
                    "canonical_smiles": candidate.get("canonical_smiles"),
                    "overall_status": "blocked",
                    "documented_route_for_graph": False,
                    "failure_type": type(error).__name__,
                    "failure": str(error),
                    "failure_preserved": True,
                    "scientific_approved": False,
                }
            value["report_provenance"] = _copy(report_provenance)
            computed.append(("synthesis", candidate_id, value))

        requirements = {
            row["analog_id"]: row["requirements"]
            for row in policy["decisions"].get("interaction_requirements", [])
        }
        choices = policy["decisions"].get("microstate_expert_decision", {})
        ph = choices.get("pH_conditions", {}) if isinstance(choices, dict) else {}
        ph_value = ph.get("pH", 7.4) if isinstance(ph, dict) else 7.4
        parent = self._parent_mol(binding["parent_id"])
        protein = self._protein_atoms()
        analog_map = self._analog_map(result)
        for analog_id in expected_analogs:
            analog = analog_map[analog_id]
            diagnostics = []
            try:
                docking = analog.get("docking", {})
                pose_ref = docking.get("files", {}).get("poses_sdf")
                pdbqt_ref = docking.get("files", {}).get("poses_pdbqt")
                source_ref = docking.get("source_ref")
                require(_ref(pose_ref) and _ref(source_ref),
                        "SCIENTIFIC_INTERACTION_SOURCE_REFS")
                diagnostic_source = self.port.json(source_ref)
                diagnostics = diagnostic_source.get("results", {}).get(
                    "pose_preservation", {}).get("all_poses_diagnostics", [])
                poses = list(Chem.ForwardSDMolSupplier(
                    io.BytesIO(self.port.read(pose_ref)), removeHs=False))
                def load_pdbqt():
                    require(_ref(pdbqt_ref),
                            "SCIENTIFIC_INTERACTION_PDBQT_SOURCE_REQUIRED")
                    return self.port.read(pdbqt_ref)

                index, chosen_pose, declared_pose, normalized_pose, map_receipt = (
                    self._actual_passing_pose(
                        diagnostic_source, poses, analog,
                        pdbqt_loader=load_pdbqt,
                    )
                )
                if map_receipt is not None:
                    map_receipt.update({
                        "actual_raw_sdf_source_ref": _copy(pose_ref),
                        "actual_raw_pdbqt_source_ref": _copy(pdbqt_ref),
                        "actual_raw_sdf_sha256": pose_ref["sha256"],
                        "actual_raw_pdbqt_sha256": pdbqt_ref["sha256"],
                    })
                actual_hydrogens = protein_h["hydrogens"] if protein_h is not None else None
                value = interaction_review.assess_interactions(
                    parent, normalized_pose, protein,
                    protein_hydrogens=actual_hydrogens,
                    requirements=requirements.get(analog_id, []), pH=ph_value)
                value.update({
                    "report_id": "interaction-" + analog_id,
                    "analog_id": analog_id,
                    "selected_pose_index_zero_based": index,
                    "selected_pose_mode": declared_pose["mode"],
                    "selected_pose_policy": (
                        "first diagnostic pass using immutable zero-based SDF record index; "
                        "one-based mode audited separately and never used for selection"
                    ),
                    "selected_pose_diagnostic": _copy(chosen_pose),
                    "all_pose_diagnostics": diagnostics,
                    "pose_source_ref": pose_ref,
                    "diagnostic_source_ref": source_ref,
                    "map_repair_receipt": map_receipt,
                    "pH_conditions": ph,
                    "protein_hydrogen_evidence_ref": (
                        protein_h["ref"] if protein_h is not None else None
                    ),
                    "protein_hydrogen_receipt": (
                        protein_h["receipt"] if protein_h is not None
                        else self._absent_protein_hydrogen_receipt(value)
                    ),
                    "scientific_approved": False,
                    "scientific_review_status": "pending_authenticated_expert_review",
                })
            except Exception as error:
                value = {
                    "format": interaction_review.VERSION,
                    "report_id": "interaction-" + analog_id,
                    "analog_id": analog_id,
                    "requirements": [],
                    "state_alternatives": [],
                    "all_pose_diagnostics": diagnostics,
                    "failures": [{
                        "status": "failed",
                        "reason": str(error),
                        "failure_type": type(error).__name__,
                        "preserved_in_failure_ledger": True,
                    }],
                    "policy_status": "registered_computation_failure",
                    "scientific_approved": False,
                }
            value["report_provenance"] = _copy(report_provenance)
            computed.append(("interaction", analog_id, value))

        paths = []
        try:
            with self.store.db() as db:
                db.execute("BEGIN IMMEDIATE")
                current_actor = self._authenticate(session_token, db)
                _, current_result, current_binding, current_source = (
                    self._verified_reinspection_snapshot(
                        db, job_id,
                        allow_archived_runtime=allow_archived_runtime,
                    )
                )
                current_provenance = self._reinspection_report_provenance(
                    current_source, current_actor
                )
                current_policy = self._policy(db, job_id, current_result)
                current_protein_h = self._protein_h_for_compute(
                    db, job_id, current_binding, current_policy
                )
                require(current_actor["id"] == actor["id"] and
                        current_binding == binding and
                        current_provenance == report_provenance and
                        self._policy_digest(current_policy) == policy_digest and
                        ((protein_h is None and current_protein_h is None) or
                         (protein_h is not None and current_protein_h is not None and
                          protein_h["id"] == current_protein_h["id"] and
                          protein_h["receipt"] == current_protein_h["receipt"])),
                        "SCIENTIFIC_COMPUTE_INPUT_CHANGED")
                evidence_binding = {
                    **binding,
                    "evaluation_policy_digest": policy_digest,
                    "evaluation_policy": {"revision": policy["revision"],
                                          "digest": policy_digest},
                    "measurement_binding": {
                        "engine_fingerprint": engine_fingerprint(),
                        "kind": "cpu_scientific_review",
                    },
                    "reinspection": _copy(report_provenance),
                }
                created = [self._register_evidence(
                    db, paths, job_id, kind, identity, evidence_binding, value)
                    for kind, identity, value in computed]
                created.extend(self._replay_actual_ternary(
                    db, paths, job_id, current_result, current_binding,
                    current_policy, report_provenance
                ))
                return {"job_id": job_id, "policy_revision": policy["revision"],
                        "policy_digest": policy_digest, "evidence": created}
        except BaseException:
            for path in paths:
                path.unlink(missing_ok=True)
            raise

    def _safe_path(self, root: Path, value: Any, *, directory: bool = False) -> Path:
        require(isinstance(value, str) and value, "SCIENTIFIC_TERNARY_PATH")
        root = root.resolve(strict=True)
        raw = Path(value)
        candidates = [raw] if raw.is_absolute() else [Path.cwd() / raw, root / raw]
        path = next((candidate for candidate in candidates if candidate.exists()), candidates[-1])
        unresolved = path.absolute()
        for ancestor in (unresolved, *unresolved.parents):
            if ancestor.exists():
                require(not ancestor.is_symlink(), "SCIENTIFIC_TERNARY_SYMLINK")
            if ancestor == root:
                break
        path = path.resolve(strict=True)
        require(path == root or root in path.parents, "SCIENTIFIC_TERNARY_PATH_ESCAPE")
        require(path.is_dir() if directory else path.is_file(),
                "SCIENTIFIC_TERNARY_DIRECTORY" if directory else "SCIENTIFIC_TERNARY_FILE")
        require(not path.is_symlink(), "SCIENTIFIC_TERNARY_SYMLINK")
        return path

    def _flatten_metrics(self, value: Any, prefix="") -> dict:
        result = {}
        if isinstance(value, dict):
            for key, item in value.items():
                name = f"{prefix}.{key}" if prefix else str(key)
                result.update(self._flatten_metrics(item, name))
        elif type(value) in (int, float):
            result[prefix] = float(value)
        return result

    def import_calibration_distribution(self, job_id, additional_receipt_paths,
                                        session_token, *,
                                        allow_archived_runtime: bool = False):
        with self.store.db() as db:
            actor = self._authenticate(session_token, db)
            _, result, binding, source = self._verified_reinspection_snapshot(
                db, job_id,
                allow_archived_runtime=allow_archived_runtime,
            )
            policy = self._policy(db, job_id, result)
            report_provenance = self._reinspection_report_provenance(
                source, actor
            )
        policy_digest = self._policy_digest(policy)

        verified = calibration_import.verify_distribution(
            list(additional_receipt_paths))
        paths = []
        try:
            with self.store.db() as db:
                db.execute("BEGIN IMMEDIATE")
                current_actor = self._authenticate(session_token, db)
                _, current_result, current_binding, current_source = (
                    self._verified_reinspection_snapshot(
                        db, job_id,
                        allow_archived_runtime=allow_archived_runtime,
                    )
                )
                current_policy = self._policy(db, job_id, current_result)
                current_provenance = self._reinspection_report_provenance(
                    current_source, current_actor
                )
                require(current_actor["id"] == actor["id"] and
                        current_binding == binding and
                        current_provenance == report_provenance and
                        self._policy_digest(current_policy) == policy_digest,
                        "SCIENTIFIC_CALIBRATION_INPUT_CHANGED")

                artifact_refs = {}
                for name, source in verified["artifacts"]:
                    source = Path(source)
                    require(source.is_file() and not source.is_symlink(),
                            "SCIENTIFIC_CALIBRATION_ARTIFACT_FILE")
                    suffix = source.suffix.lower()
                    media = ("application/json" if suffix == ".json" else
                             "text/yaml" if suffix in {".yaml", ".yml"} else
                             "text/plain" if suffix in {".log", ".txt"} else
                             "chemical/x-mmcif" if suffix in {".cif", ".mmcif"} else
                             "application/octet-stream")
                    artifact_refs[name] = self._write(
                        db, paths, source.read_bytes(), media, "computed")

                value = {
                    "format": "tpd-source-bound-calibration-report/1",
                    "id": "CRBN_6BOY",
                    "kind": "source_independent_benchmark",
                    "candidate_id": None,
                    "source": {
                        "case": "6BOY",
                        "e3_type": "CRBN",
                        "base_manifest_sha256": verified["base_manifest_sha256"],
                    },
                    "reinspection": _copy(report_provenance),
                    "distribution": verified["distribution"],
                    "registered_artifacts": artifact_refs,
                    "scientific_approval": False,
                    "expert_threshold_present": False,
                }
                value["report_provenance"] = _copy(report_provenance)
                evidence_binding = {
                    **binding,
                    "evaluation_policy_digest": policy_digest,
                    "evaluation_policy": {
                        "revision": policy["revision"],
                        "digest": policy_digest,
                    },
                    "measurement_binding": {
                        "source_id": "CRBN_6BOY",
                        "module_version": verified["measurement_version"],
                        "module_fingerprint": verified["measurement_fingerprint"],
                        "base_manifest_sha256": verified["base_manifest_sha256"],
                    },
                }
                evidence = self._register_evidence(
                    db, paths, job_id, "calibration_distribution",
                    "CRBN_6BOY", evidence_binding, value)
                return {
                    "job_id": job_id,
                    "source_id": "CRBN_6BOY",
                    "technical_success_count":
                        verified["distribution"]["technical_success_count"],
                    "scientific_model_quality_success_count": None,
                    "evidence": evidence,
                }
        except BaseException:
            for path in paths:
                path.unlink(missing_ok=True)
            raise

    def _matching_ternary_inspection(self, run_root: Path, rerun: Any,
                                     recorded: Any) -> bool:
        if not isinstance(rerun, dict) or not isinstance(recorded, dict):
            return False

        normalized_rerun = dict(rerun)
        normalized_recorded = dict(recorded)
        allowed_suffixes = {
            "prediction": {".cif", ".mmcif"},
            "confidence_file": {".json"},
        }
        try:
            for key, suffixes in allowed_suffixes.items():
                if key not in rerun or key not in recorded:
                    return False
                rerun_path = self._safe_path(run_root, rerun[key])
                recorded_path = self._safe_path(run_root, recorded[key])
                if rerun_path.suffix.lower() not in suffixes:
                    return False
                if recorded_path != rerun_path:
                    return False
                normalized_rerun[key] = rerun_path
                normalized_recorded[key] = recorded_path
        except Exception:
            return False

        return normalized_rerun == normalized_recorded

    def import_ternary(self, job_id, plan_path, run_directory, session_token,
                       *, allow_archived_runtime: bool = False):
        plan_path = Path(plan_path).resolve()
        run_root = Path(run_directory).resolve()
        require(plan_path.is_file() and not plan_path.is_symlink(),
                "SCIENTIFIC_TERNARY_PLAN_FILE")
        require(run_root.is_dir() and not run_root.is_symlink(),
                "SCIENTIFIC_TERNARY_RUN_DIRECTORY")
        plan_raw = plan_path.read_bytes()
        plan = parse_json(plan_raw)
        novel_ternary.verify_plan(plan)

        with self.store.db() as db:
            actor = self._authenticate(session_token, db)
            _, result, binding, source = self._verified_reinspection_snapshot(
                db, job_id,
                allow_archived_runtime=allow_archived_runtime,
            )
            policy = self._policy(db, job_id, result)
            report_provenance = self._reinspection_report_provenance(
                source, actor
            )
        policy_digest = self._policy_digest(policy)
        job_binding = plan["bindings"]["job"]
        require(job_binding["project"] == self.project and job_binding["job_id"] == job_id,
                "SCIENTIFIC_TERNARY_JOB_BINDING")
        require(job_binding["input_sha256"] == binding["input_sha256"]
                and job_binding["result_sha256"] == binding["result_sha256"],
                "SCIENTIFIC_TERNARY_SOURCE_BINDING")
        from scripts.run_novel_ternary import _candidate as qualified_candidate
        candidates = self._candidate_map(result)
        candidate_id = plan["candidate_graph"]["candidate_id"]
        require(candidate_id in candidates, "SCIENTIFIC_TERNARY_CANDIDATE")
        actual_graph = novel_ternary._mapped_graph(qualified_candidate(result, candidate_id))
        require(actual_graph["graph_sha256"] == job_binding["candidate_graph_sha256"]
                == plan["candidate_graph"]["graph_sha256"],
                "SCIENTIFIC_TERNARY_GRAPH_BINDING")

        aggregate_path = self._safe_path(run_root, "receipt.json")
        aggregate = parse_json(aggregate_path.read_bytes())
        require(aggregate.get("plan_digest") == plan["plan_digest"],
                "SCIENTIFIC_TERNARY_RECEIPT_PLAN")
        require(aggregate.get("candidate_id") == candidate_id,
                "SCIENTIFIC_TERNARY_RECEIPT_CANDIDATE")
        require(aggregate.get("bindings") == {
                    "project": self.project,
                    "job_id": job_id,
                    "policy_digest": plan["bindings"]["policy"]["digest"],
                }, "SCIENTIFIC_TERNARY_AGGREGATE_BINDING")
        receipt_paths = aggregate.get("all_seed_receipts")
        require(isinstance(receipt_paths, list) and len(receipt_paths) == len(plan["seeds"]),
                "SCIENTIFIC_TERNARY_RECEIPT_COVERAGE")

        receipts = []
        seen_seeds = set()
        for receipt_name in receipt_paths:
            receipt_path = self._safe_path(run_root, receipt_name)
            receipt = parse_json(receipt_path.read_bytes())
            seed = receipt.get("seed")
            require(seed in plan["seeds"] and seed not in seen_seeds,
                    "SCIENTIFIC_TERNARY_RECEIPT_SEED")
            seen_seeds.add(seed)
            require(receipt.get("plan_digest") == plan["plan_digest"]
                    and receipt.get("candidate_id") == candidate_id
                    and receipt.get("e3_type") == plan["candidate_graph"]["e3_type"]
                    and receipt.get("bindings") == aggregate.get("bindings"),
                    "SCIENTIFIC_TERNARY_RECEIPT_BINDING")
            artifacts = receipt.get("artifacts", {})
            copied = {"receipt": receipt_path}
            for key in ("input_yaml", "log"):
                if artifacts.get(key):
                    copied[key] = self._safe_path(run_root, artifacts[key])
            hashes = receipt.get("hashes", {})
            require(copied.get("input_yaml") is not None and
                    isinstance(hashes.get("input_yaml_sha256"), str) and
                    _sha(copied["input_yaml"].read_bytes()) == hashes["input_yaml_sha256"] ==
                    plan["boltz_input"]["sha256"], "SCIENTIFIC_TERNARY_YAML_HASH")
            require(copied.get("log") is not None and
                    isinstance(hashes.get("inference_log_sha256"), str) and
                    _sha(copied["log"].read_bytes()) == hashes["inference_log_sha256"],
                    "SCIENTIFIC_TERNARY_LOG_HASH")

            success = (receipt.get("status") == "completed"
                       and receipt.get("process_state") == "completed"
                       and receipt.get("exit_code") == 0
                       and receipt.get("actual_computation") is True)
            inspection = receipt.get("inspection", {})
            if success:
                output_value = artifacts.get("output_directory")
                require(isinstance(output_value, str), "SCIENTIFIC_TERNARY_OUTPUT_DIRECTORY")
                output_dir = self._safe_path(run_root, output_value, directory=True)
                recorded_tree = hashes.get("output_files_sha256")
                require(isinstance(recorded_tree, dict) and recorded_tree,
                        "SCIENTIFIC_TERNARY_OUTPUT_HASH_TREE")
                require(novel_ternary.worker._hash_tree(output_dir) == recorded_tree,
                        "SCIENTIFIC_TERNARY_OUTPUT_HASH_TREE")
                rerun = novel_ternary.assess_novel_outputs(output_dir, plan)
                require(rerun.get("status") == "computed_hypothesis",
                        "SCIENTIFIC_TERNARY_REASSESSMENT_FAILED")
                require(inspection.get("status") == "computed_hypothesis" and
                        inspection.get("actual_computation") is True,
                        "SCIENTIFIC_TERNARY_RECORDED_INSPECTION")
                require(self._matching_ternary_inspection(run_root, rerun, inspection),
                        "SCIENTIFIC_TERNARY_INSPECTION_MISMATCH")
                inspection = rerun
                copied["prediction"] = self._safe_path(run_root, rerun["prediction"])
                copied["confidence"] = self._safe_path(run_root, rerun["confidence_file"])
            metrics = self._flatten_metrics(inspection.get("descriptive_metrics", {}))
            geometry = policy.get("decisions", {}).get("ternary_geometry_criterion", {})
            metric = geometry.get("metric")
            threshold = geometry.get("threshold")
            mode = geometry.get("comparison")
            observed = metrics.get(metric) if isinstance(metric, str) else None
            compatible = None
            if type(observed) in (int, float) and type(threshold) in (int, float):
                compatible = observed <= threshold if mode in {"lte", "max_lte"} else observed >= threshold
            value = {
                "format": novel_ternary.RECEIPT_FORMAT,
                "candidate_id": candidate_id,
                "e3_type": plan["candidate_graph"]["e3_type"],
                "seed": seed,
                "reference_free": True,
                "actual_computation": success,
                "execution_success": success,
                "failure": None if success else {
                    "status": receipt.get("status"),
                    "process_state": receipt.get("process_state"),
                    "exit_code": receipt.get("exit_code"),
                    "reason": receipt.get("failure_reason"),
                },
                "inspection": inspection,
                "inspection_summary": {
                    "compatibility": ("compatible" if compatible is True else
                                      "incompatible" if compatible is False else "unknown"),
                    "metrics": metrics,
                    "criterion_metric": metric,
                    "criterion_observed": observed,
                },
                "producer_policy_binding": plan["bindings"]["policy"],
                "measurement_protocol": {
                    "plan_digest": plan["plan_digest"],
                    "candidate_graph_sha256": plan["candidate_graph"]["graph_sha256"],
                    "msa_mode": plan["msa_mode"],
                    "boltz_input_sha256": plan["boltz_input"]["sha256"],
                    "seed": seed,
                },
                "evaluation_policy": {"revision": policy["revision"],
                                      "digest": policy_digest},
                "policy_incompatibilities": [] if compatible is not None else [
                    "current geometry metric is absent from measured output"
                ],
                "scientific_status": "computed_hypothesis" if success else "not_computed",
            }
            value["report_provenance"] = _copy(report_provenance)
            receipts.append((seed, value, copied))
        require(seen_seeds == set(plan["seeds"]), "SCIENTIFIC_TERNARY_SEED_COVERAGE")

        paths = []
        try:
            with self.store.db() as db:
                db.execute("BEGIN IMMEDIATE")
                current_actor = self._authenticate(session_token, db)
                _, current_result, current_binding, current_source = (
                    self._verified_reinspection_snapshot(
                        db, job_id,
                        allow_archived_runtime=allow_archived_runtime,
                    )
                )
                current_policy = self._policy(db, job_id, current_result)
                current_provenance = self._reinspection_report_provenance(
                    current_source, current_actor
                )
                require(current_actor["id"] == actor["id"] and
                        current_binding == binding and
                        current_provenance == report_provenance and
                        self._policy_digest(current_policy) == policy_digest,
                        "SCIENTIFIC_TERNARY_INPUT_CHANGED")
                plan_ref = self._write(db, paths, plan_raw, provenance="computed")
                aggregate_ref = self._write(db, paths, aggregate_path.read_bytes(), provenance="computed")
                created = []
                for seed, value, copied in receipts:
                    artifact_refs = {}
                    for key, path in copied.items():
                        media = ("application/json" if path.suffix == ".json" else
                                 "text/yaml" if path.suffix in {".yaml", ".yml"} else
                                 "text/plain" if path.suffix in {".log", ".txt"} else
                                 "chemical/x-mmcif" if path.suffix in {".cif", ".mmcif"} else
                                 "application/octet-stream")
                        artifact_refs[key] = self._write(
                            db, paths, path.read_bytes(), media, "computed")
                    value["registered_artifacts"] = artifact_refs
                    value["plan_ref"] = plan_ref
                    value["aggregate_receipt_ref"] = aggregate_ref
                    evidence_binding = {
                        **binding,
                        "evaluation_policy_digest": policy_digest,
                        "evaluation_policy": {"revision": policy["revision"],
                                              "digest": policy_digest},
                        "measurement_binding": {
                            "producer_policy": plan["bindings"]["policy"],
                            "plan_digest": plan["plan_digest"],
                            "candidate_graph_sha256": plan["candidate_graph"]["graph_sha256"],
                            "msa_mode": plan["msa_mode"],
                        },
                        "reinspection": _copy(report_provenance),
                    }
                    identity = f"{candidate_id}|{plan['candidate_graph']['e3_type']}|{seed}"
                    created.append(self._register_evidence(
                        db, paths, job_id, "ternary", identity, evidence_binding, value))
                return {"job_id": job_id, "candidate_id": candidate_id,
                        "plan_digest": plan["plan_digest"], "evidence": created}
        except BaseException:
            for path in paths:
                path.unlink(missing_ok=True)
            raise

    def export(self, assessment_id, session_token):
        paths = []
        try:
            with self.store.db() as db:
                db.execute("BEGIN IMMEDIATE")
                actor = self._authenticate(session_token, db)
                assessment = self._effective(db, self._assessment(db, assessment_id))
                job, result, _, source = self._verified_result(
                    db, assessment["job_id"], require_current=False)
                policy = self._policy(db, assessment["job_id"], result)
                packet = {
                    "format": "scientific-acceptance-export/2.0",
                    "created_at": now(),
                    "created_by": actor,
                    "assessment": assessment,
                    "current_implementation": assessment["current_implementation"],
                    "current_job_runtime": assessment["current_job_runtime"],
                    "current_source_inputs": source["source_inputs_current"],
                    "current_policy": self._public_policy(policy),
                    "pending_scientific_questions": assessment.get("expert_questions", []),
                    "pending_actions": assessment.get("pending_actions", []),
                    "scientific_gates": assessment.get("criteria", []),
                    "evidence_links": assessment.get("evidence_refs", []),
                    "job": {
                        "id": job["id"], "state": job["state"],
                        "input_digest": job["input_digest"],
                        "runtime": job.get("runtime"), "outputs": job.get("outputs", []),
                    },
                    "development_status": {
                        "design_status": result.get("status"),
                        "summary": result.get("summary", {}),
                        "candidate_progress": [{
                            "candidate_id": row.get("candidate_id"),
                            "e3_type": row.get("e3_type"),
                            "status": row.get("status", row.get("assembly_status")),
                            "files": row.get("files", {}),
                        } for row in result.get("protac_candidates", [])
                           if isinstance(row, dict)],
                    },
                }
                json_ref = self._write(db, paths, packet)
                lines = [scientific_assessment.report_markdown(assessment),
                         "\n## Current implementation and source state\n",
                         f"- Current implementation: {assessment['current_implementation']}",
                         f"- Current job runtime: {assessment['current_job_runtime']}",
                         f"- Current source inputs: {assessment['source_inputs_current']}",
                         f"- Current policy: {assessment['current_policy']}",
                         "\n## Evidence links\n",
                         json.dumps(packet["evidence_links"], ensure_ascii=False, indent=2),
                         "\n## Scoped policy choices\n",
                         json.dumps(packet["current_policy"], ensure_ascii=False, indent=2)]
                markdown_ref = self._write(
                    db, paths, "\n".join(lines).encode("utf-8"), "text/markdown")
                return {"assessment_id": assessment_id, "json_ref": json_ref,
                        "report_ref": markdown_ref}
        except BaseException:
            for path in paths:
                path.unlink(missing_ok=True)
            raise


__all__ = [
    "ScientificAcceptanceService", "VERSION", "engine_fingerprint",
    "MAX_RESULT_BYTES", "MAX_MICROSTATE_SOURCE_BYTES", "POLICY_KINDS",
]
