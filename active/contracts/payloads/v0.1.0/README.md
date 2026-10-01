> 경로 정리: 2026-09-22. 아래는 기존 계약/시제품의 작성·검사 이력이다. 목표 제품과 최신 분업은 [INTERFACES](<../../../../docs/INTERFACES.md>), 실제 구현 범위는 [STATUS](<../../../../docs/STATUS.md>)를 따른다. 과거 A/B 배분·public_explore·다음 작업은 현재 배정/공개 범위/개발 재개를 뜻하지 않는다.

# I1: 표적·문헌 payload 계약 v0.1.0

상태: **I1 JSON Schema·오프라인 검증 완료. 후속 I2에서 두 작업의 로컬 실제 연결도 검증했다.** [13 계약 설계](<../../../../../Bio_JumpAI_보관자료/본선_설계변경이력/2026-09-11_12_Codex_원본/20260911/13_I1_표적문헌_데이터계약.md>), [현재 I2 인터페이스](<../../../docs/i2_interfaces.md>)를 함께 읽는다. 공개 서비스·생물학적 검증은 미완료다. D001·D006·D007에 따라 새로 작성했으며 Study/JH 자산을 이식하지 않았다.

## 계약 파일과 버전

- [domain.schema.json](<domain.schema.json>): JSON Schema 2020-12. 38개 정의, 7개 최상위 payload. 모든 객체는 미등록 필드를 거부한다.
- [build_schemas.py](<build_schemas.py>): 위 JSON의 원본 정의. 변경 후 생성하며, 검사기가 생성물과의 일치를 확인한다.
- [validation.py](<validation.py>): schema + ID/참조/coverage 정합성 + 저장 바이트 hash + JobSpec/결과 연결 검사. 네트워크 조회 없이 두 schema를 메모리 registry에 등록한다.
- [examples/artifact-index.json](<examples/artifact-index.json>): 독립적으로 만든 가상 파일 16개의 ID·버전·실제 SHA-256과 상대 경로.
- [examples/target.job.json](<examples/target.job.json>), [examples/target.result.json](<examples/target.result.json>): 표적 식별 요청/결과 봉투.
- [examples/literature.job.json](<examples/literature.job.json>), [examples/literature.result.json](<examples/literature.result.json>): 검색 원문 1건과 추출 실패 PDF 1건을 포함한 부분 결과.
- [build_fixtures.py](<build_fixtures.py>): 예시 재생성. 실제 표적·논문·조회·추출을 사용하지 않는다.
- [validate_examples.py](<validate_examples.py>), [verification.json](<verification.json>): 실행 가능한 계약 검사와 결과 기록.

경계 봉투 버전은 기존 `contracts/v0.1.0/`의 `contract_version=0.1.0`을 유지한다. 이 폴더는 별도로 `payload_version=0.1.0`을 사용한다. 두 버전이 같은 숫자라도 같은 계약은 아니다. 기존 I0 예시의 임시 hash를 검증 가능한 실제 데이터로 승격한 것이 아니라, 새로운 I1 예시를 작성했다.

Schema 식별자:

```text
urn:tpd-navigator:boundary:0.1.0
urn:tpd-navigator:payloads:0.1.0
urn:tpd-navigator:payloads:0.1.0#/$defs/TargetQuery
urn:tpd-navigator:payloads:0.1.0#/$defs/TargetResolution
urn:tpd-navigator:payloads:0.1.0#/$defs/LiteratureCollectionRequest
urn:tpd-navigator:payloads:0.1.0#/$defs/EvidenceBundle
```

URN은 문서 식별자이며 접속할 URL이 아니다. 다른 언어의 소비자도 두 schema를 위 URN으로 등록한다. Python 검사기는 이를 처리한다. [작업 목록](<../../v0.1.0/operations.json>)에서 I2의 2개 operation만 local-single-user-i2 실행을 허용하며 공개 dispatch·나머지 12개는 비활성이다.

## 개발자가 사용할 진입점

| 데이터 | 생산자 → 소비자 | 주요 의미 |
|---|---|---|
| TargetQuery | M2 → M4 | query/kind/taxon, isoform·mutation·construct, 선택 연구 맥락 |
| TargetResolution | M4 → M2/M3 | 후보들, resolved/ambiguous/not_found, 선택 근거, 제품 범위 상태 |
| LiteratureCollectionRequest | M2 → M3 | 식별된 표적 참조, 목적별 질문/질의, 선택 PDF 참조, 수집 정책 참조 |
| EvidenceBundle | M3 → M2/M6/M4 | 원문·추출·읽기 상태, locator, 원문 관측, 화합물 서술, coverage와 부분 실패 |
| InputManifest | M2 → 각 모듈 | 작업 입력·정책·인자·모드·revision의 불변 묶음. 원본 바이트 hash 사용 |
| CollectionPolicy | M2 → M3 | Europe PMC, 레코드/문서/구간/원문 바이트 상한. 예시값은 운영 기본값 아님 |
| OperationParameters | M2 → 각 모듈 | 구현 profile_id. I2는 local-cpu-i2-v1을 사용하며 실행 권한 프로필과 별도로 검사 |

원문을 찾았지만 분석 전이면 `EvidenceBundle.evidence_records=[]`를 반환한다. 이미 표준화된 원문 관측을 수집할 수 있는 경로는 EvidenceRecord를 채울 수 있다. LLM 해석·요약·PolicyProposal을 수행하는 `assess_literature`의 요청/결과는 I3의 별도 계약이며 이번에 구현하지 않았다.

## 값과 상태의 규칙

1. 알려진 문자열은 `{state:known, value:..., reason:null}`이다. 미상/미제공/미보고/해당 없음은 value=null과 이유를 함께 보낸다. 빈 문자열·0으로 대체하지 않는다.
2. Measurement는 scalar(값·부등호·단위·원문), interval(경계·포함 여부), qualitative, missing의 네 종류다. `not_detected`는 관측된 정성 결과로 표현하고 `not_measured`와 구분한다. 수치가 있지만 단위가 미상인 경우 unit의 state로 명시하며 자동 정규화하지 않는다.
3. SourceLocator는 원문과 추출 snapshot을 모두 참조한다. XML은 block_id+절대 XPath, PDF는 1부터 시작하는 page+segment_id, 초록은 record/abstract/segment, metadata는 record/field다. PDF bbox는 0~1 정규화 좌표이며 좌상단 원점이다. 문서 화면상의 인쇄 페이지 번호와 파일 페이지 번호를 섞지 않는다.
4. 원문 확보, 유효 전문 접근 수준, 추출 완료, LLM 읽기, 관측 검토는 별개다. 형식/내용을 확인 못한 파일은 access_level=unverified다. 예시의 잘린 PDF는 받은 파일이지만 읽을 수 있는 전문으로 집계하지 않는다.
5. `basis_kind=experiment`와 `producer.kind=llm`은 함께 존재할 수 있다. 원문 실험을 LLM이 추출했다는 뜻이다. 실험의 진위·추출 정확성이 schema 검증으로 입증되는 것은 아니다.
6. Compound의 `provided_unverified`는 구조 파일이 있다는 의미다. 분자 타당성·원자 매핑 완료를 뜻하지 않는다. AttachmentDescription에는 원문 R1 서술만 넣을 수 있고 제품 atom_map_id는 허용하지 않는다. F03의 실제 매핑·비교는 M4 후속 계약이다.
7. 검색에서 나온 다른 표적의 관측도 보존한다. EvidenceRecord의 표적을 현재 연구 표적으로 덮어쓰거나 같다고 강제하지 않는다. 현재 표적에 적용할 수 있는지는 F03 등 후속 판정의 책임이다.
8. 논문별 레코드 수와 독립 실험 수를 분리한다. 원 실험 ID를 확인하지 못한 실험 레코드가 있으면 독립 실험 수는 unknown이다. 서로 다른 출처에 기록된 같은 origin_id는 한 번만 센다.

## 호출과 검사

이 폴더에서:

```text
python validate_examples.py
```

확인한 환경: Python 3.12, jsonschema 4.26.0, referencing 0.37.0. 이번 검사를 위해 패키지를 설치하지 않았다. 제품 런타임 전체의 의존성 버전을 이 환경으로 확정한 것은 아니다.

구현 모듈은 다음 포트를 사용할 수 있다.

```python
from validation import ArtifactReader, validate_exchange

# read_bytes(ref)는 M2 ArtifactStore가 제공할 함수다.
# 프로젝트 접근권한·등록 메타데이터를 검사한 뒤 정확한 bytes를 돌려준다.
reader = ArtifactReader(read_bytes)
check = validate_exchange(job_spec, reader, module_result)
assert check["contract_valid"]
assert check["dispatch_authorized"] is False
```

위 import는 현재 폴더가 Python 경로에 있을 때의 예시다. active 런타임은 `packages.contracts`를 통해 같은 검사기를 사용한다. `validate_payload(kind,value)`는 구조/내부 검사, `validate_exchange`는 저장 파일·요청/결과 대조이며 후자가 호출돼도 작업을 실행하지 않는다.

`InputManifest`에는 approval_ref·attempt·limits·reservation_id를 넣지 않는다. 승인 artifact가 입력 digest를 참조해도 서로의 hash를 재귀적으로 요구하지 않도록 한다. 이 값들은 JobSpec에 남고, M2가 승인 범위·실행 한도·예약·재시도를 검사한다. 원문/정책/파라미터/실행 모드 변경은 manifest와 input_revision을 갱신한다.

## 검사 범위와 남은 작업

검사하는 것: schema, 참조 ID/버전/추출기 일치, locator 형식, 유한 수치·구간 순서, 중복 ID, 표적 선택 일관성·요청 종 일치, 근거/화합물 연결, coverage 계산, LLM 읽기 기록, 원 실험 중복, 파일 hash, 요청·결과 상관 ID, 입력 선언, 일부 수집 한도, fixture의 real 모드 혼입.

별도 구현이 필요한 것: 실제 UniProt 검색과 ID의 생물학적 유효성, HTTP 검색·PDF/XML 파서, XPath/페이지/추출 내용의 원문 대조, 화학 구조·assay 비교, LLM 추출 정확성, 접근권한·공개 권리의 진위, 비용/승인/재시도/타임아웃 통제. `read_bytes`의 접근 통제나 provenance의 진위를 이 검사기가 대신 보장하지 않는다. 모든 원시 파일의 임의 내부 형식까지 재귀적으로 해석하지도 않는다.

특히 source/typed artifact를 읽는 콜백은 신뢰할 수 있는 등록 정보를 사용해야 한다. 호출자가 provenance나 schema_id를 바꿔 전송했다고 저장소의 등록 정보를 바꿔주는 구현은 허용되지 않는다.

예시 PDF는 의도적으로 잘린 **오류 fixture**이고 XML의 표적·수치도 가상이다. 실제 평가·학습·제출 재생 데이터로 사용하지 않는다. 후속 I2는 이 fixture를 실제 데이터로 승격하지 않고 별도로 공개 원문을 수집했다. 과학 검증은 여전히 별도다.
