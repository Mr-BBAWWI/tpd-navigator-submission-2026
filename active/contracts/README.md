> 경로 정리: 2026-09-22. 아래는 기존 계약/시제품의 작성·검사 이력이다. 목표 제품과 최신 분업은 [INTERFACES](<../../docs/INTERFACES.md>), 실제 구현 범위는 [STATUS](<../../docs/STATUS.md>)를 따른다. 과거 A/B 배분·public_explore·다음 작업은 현재 배정/공개 범위/개발 재개를 뜻하지 않는다.

2026-09-24 상태 구분: `drafts/`라는 경로명만으로 실행 여부를 판단하지 않는다. B의 `b_review_material`/`b_candidate_evidence`/`b_literature_exchange`는 실제 B 소비자와 A 인계 어댑터가 사용한다. A의 `candidate_dossier`/`a_review_handoff`는 서비스의 검토 자료, `literature_agent`는 문헌 응답, `agent_review`는 계층형 서비스 응답, `claim_review`는 비교 실험 응답 검사에 사용된다. 초안 규격 사용은 공식 Phase1/GatePacket/ApprovalRecord나 나머지 I0 작업 실행 활성화를 뜻하지 않는다. 아래 과거 “drafts는 실행 schema가 아님” 문구는 당시 시점의 기록이다.

2026-09-23 A 연결: I3 0.1.0의 분석/전달/검토 계약을 내부 `local-literature-review` 서비스에서 사용한다. `drafts/literature_agent.schema.json`은 해당 LLM 응답 규격이며 I1 관측 정의를 참조한다. GatePacket/ApprovalRecord 실서비스와 I0 작업 등록 통합은 미완료다. [연결 명세](../docs/a_literature_service.md)를 따른다.


# 개발자 간 공통 계약

2026-09-23 에이전트 추가: [agent_review schema](drafts/agent_review.schema.json)는 내부 실험의 계획·자료 읽기 요청·해석·비판·보고 순서를 검증하는 실제 로컬 검사 규격이다. `drafts/` 중 이 파일과 B 결과 schema는 각각 독립 시제품이 사용하지만 **M2 승인/작업 실행 계약으로 활성화된 것은 아니다.** [실행 범위](../docs/agent_review_experiment.md)를 확인한다.

2026-09-23 B 착수 추가: [b_case_preparation 초안](drafts/b_case_preparation.schema.json)은 공개 원자료 기반 CPU 결과 파일을 검사하는 B 전용 인계 초안이다. 기존 I0/I1/I3 schema·작업 registry를 바꾸지 않았고 M2 통합 계약으로 채택되지 않았다. [실행/인계 안내](../docs/b_smarca2_start.md)를 참고한다.

작업 기준: [분업·모듈·통합 명세](<../../../Bio_JumpAI_보관자료/본선_설계변경이력/2026-09-11_12_Codex_원본/20260911/11_분업을위한_모듈경계와통합계약.md>), 현재 개발 기준 D007. 새로 작성한 계약이며 Study/JH 코드·데이터를 가져오지 않았다.

`v0.1.0/`은 I0 **경계 봉투**, [payloads/v0.1.0](<payloads/v0.1.0/README.md>)은 I1 **표적·문헌 payload**다. I1은 38 정의·가상 파일 16개·요청/결과 2쌍·정상 40/오류 41 검사를 포함한다. 후속 [I2](<../docs/i2_interfaces.md>)에서 실제 로컬 표적/수집기·웹을 연결했다. 실제 담당 배정·전체 도메인·공개 서버 완성을 뜻하지 않는다.

후속 F03 조사로 [부착점 SAR 도메인 명세](<../../docs/design/notes/attachment_sar.md>)와 [draft 수용 사례](<drafts/README.md>)를 작성했다. `drafts/`는 실행 schema가 아니며 현재 작업 등록·실행 가능 상태를 변경하지 않는다.

후속 **I3-1**: [review/v0.1.0](<review/v0.1.0/README.md>)에 분석·정책·후보 비교·검토·G1의 별도 schema와 오프라인 연결 검사를 추가했다. 25 정의·14 진입점·정상 25·거부 45. I1 schema와 실행 registry는 변경하지 않았으며, 새 계약이 실제 LLM/CPU/승인 서비스에 연결된 것은 아니다. [I3 분업 명세](<../docs/i3_review_interfaces.md>)를 함께 읽는다.

- `module_boundary.schema.json`: JSON Schema 2020-12. JobSpec, ModuleResult, ArtifactRef, Issue의 형식·일부 구조 조건.
- `operations.json`: `resolve_target`, `collect_literature`만 schema_validated + live_execution_ready=true, 허용 프로필은 local-single-user-i2다. 실제 I2 검증 기록을 참조한다. 공개 dispatch와 나머지 12개는 비활성이다. schema 유효성과 실행 허가는 별개다.
- `examples/`: `data_mode=test_fixture`인 문헌 수집 요청·부분 결과 예시. `a` 64자의 hash와 참조 ID는 형식 검사용 가상 값이다. 실제 논문·실측·산출물 파일이나 사용량 측정 증거가 아니다.
- `validate_examples.py`: 기존 Python 환경의 `jsonschema` 패키지로 schema·예시·금지 형태·작업 ID 일치를 검사한다. 패키지 설치나 API 호출은 하지 않는다.

검사 명령:

```text
python validate_examples.py
```

위 명령은 이 파일 아래의 `v0.1.0/` 디렉터리에서 실행한다. Python 환경에 `jsonschema`가 있어야 한다. 앱의 프레임워크·운영 의존성 버전을 이 검사 환경으로 확정하는 것은 아니다.

I0 검증 범위: 봉투의 필드/타입, 허용 작업 ID, 결과 상태, 공개 GPU 요청의 일부 구조 조건, 예시의 요청-결과 상관 ID. I1 검사기는 추가로 두 operation의 payload·참조·coverage와 저장 바이트 hash 및 요청/결과의 일치를 검사한다. 실제 접근 통제·실행 주체 인증·승인 범위·비용 원장·중복 실행 방지·원문 내용과 과학적 타당성은 별도 구현/검증이다. 유효한 JSON은 실행 허가가 아니다.

입력 digest는 Orchestrator가 저장한 **입력 manifest 원본 바이트의 SHA-256**으로 정의한다. 전달된 `input_manifest.sha256`와 같아야 한다. 워커는 manifest 바이트를 읽고 hash 및 내부의 operation/input/policy/parameters 참조 일치를 검사한다. 각 언어가 JSON을 임의로 재직렬화한 결과를 따로 hash하지 않는다. manifest 실제 schema와 검사기는 I1에 있다. approval_ref·attempt·limits·reservation_id는 JobSpec에 두며 manifest와 승인 기록의 hash 순환을 만들지 않는다. 이 정보의 권한·범위·집행은 M2의 별도 책임이다.

기존 `../../../Bio_JumpAI_보관자료/본선_설계변경이력/2026-09-11_12_Codex_원본/20260911/agent_contracts.v0.1.json`은 **LLM 역할·행동 설계 목록**이다. 이 폴더의 모듈 간 작업 계약과 목적이 다르며 파일 하나가 다른 하나를 대체하지 않는다. LLM 행동을 바로 GPU JobSpec으로 전달하지 않는다.

공통 계약 변경은 생산자와 소비자를 함께 검토한다. 필드 이름·타입·필수 여부·단위·enum·null 의미 변경은 호환성 변경이다. 새 버전과 예시·이행 작업을 함께 준비하며, 서로 다른 작업 브랜치가 같은 계약을 각자 수정하지 않도록 담당자를 정한다.
