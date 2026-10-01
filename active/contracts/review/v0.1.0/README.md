> 경로 정리: 2026-09-22. 아래는 기존 계약/시제품의 작성·검사 이력이다. 목표 제품과 최신 분업은 [INTERFACES](<../../../../docs/INTERFACES.md>), 실제 구현 범위는 [STATUS](<../../../../docs/STATUS.md>)를 따른다. 과거 A/B 배분·public_explore·다음 작업은 현재 배정/공개 범위/개발 재개를 뜻하지 않는다.

# I3-1 — 문헌 분석·정책 변경·후보 미리보기·G1 계약

2026-09-11. **25개 schema 정의·14개 진입점과 오프라인 연결 검사를 완료했다.** 실제 LLM, 화학 계산, 승인 DB/API/UI는 구현하지 않았다. I2의 실행 registry·서비스·저장 결과는 변경하지 않았다.

이 폴더의 버전은 `urn:tpd-navigator:review:0.1.0`이다. I1 `payloads:0.1.0`을 덮어쓰지 않고, 기존 EvidenceBundle/ArtifactRef를 참조한다. 같은 버전 숫자가 같은 schema라는 뜻은 아니다. 공개/로컬 실행 허가는 별도다.

## 파일

- [build_schema.py](<build_schema.py>): schema 원본. [workflow.schema.json](<workflow.schema.json>): 생성물.
- [validation.py](<validation.py>): 구조·파일 hash·인용·정책 차이·후보 비교·검토 대상·G1 바인딩 검사. 과학적 정답·인증·비용 집행을 대신하지 않는다.
- [build_fixtures.py](<build_fixtures.py>), [예시 목록](<examples/artifact-index.json>): 새 가상 artifact 24개 + 기존 새 구현 I1 fixture 16개 참조. 실제 논문·분자·모델 출력·승인 기록이 아니다.
- [validate_examples.py](<validate_examples.py>), [검증 기록](<verification.json>): 정상 25개·거부 45개·40파일 hash 확인.
- [분업·연결 명세](<../../../docs/i3_review_interfaces.md>): 생산자/소비자, 상태·권한·승인 저장의 후속 구현 요구.

```powershell
# active/에서 실행
& '.venv/Scripts/python.exe' contracts/review/v0.1.0/build_schema.py
& '.venv/Scripts/python.exe' contracts/review/v0.1.0/validate_examples.py
```

두 번째 명령은 가상 예시를 재생성하고 검사한다. 네트워크·모델·GPU·실험 데이터 없이 실행하며 테스트용 승인 객체는 `test_fixture`, `dispatch_authorized=false`다.

## 진입점과 소유자

| 객체 | 작성 주체 → 소비자 | 핵심 의미 |
|---|---|---|
| PolicySnapshot | M2 → M3/M4/M6 | 표적, 필드별 값·단위/endpoint·변경 허용·범위·필수 여부. 임의 과학 기본값 없음 |
| LiteratureAssessmentRequest | M2 → M3 | I1 원문/표적/연구 맥락, 현재 정책, 허용 변경 필드, 요청 구간, 출력 한도 |
| DeliveryReceipt | M2 LLM 통로 → M3/M6 | 실제 전달 구간, 호출 완료/불완전/실패/거절, 모델·프롬프트 버전·원응답 참조 |
| LiteratureAssessment | M3 → M2/M4/M6 | 원문 관측 bundle의 새 버전, claims·요약·상충·변경 제안·미전달 범위 |
| CandidateUniverse | M4 → M2/M4 | 비교에 공통으로 쓸 후보 ID·분자 구조 버전·생성 맥락 |
| Phase1Result | M4 → M2/M6 | 같은 후보 입력에 정책별 CPU 판정, 매핑/화학 상태·근거 종류·미확인. 완성 PROTAC/효능 결과 아님 |
| PreviewRequest | M2 → CPU 작업 조정 | 선택한 제안 IDs, 변경 전후 정책, 공통 후보/계산 맥락, 사전 Critic 검토 |
| CpuPreview | M2 → M1/M6 | M4 전후 결과 참조, 추가/제외/보류 후보 IDs. 활성 정책 변경 없음 |
| ReviewRequest | M2 → M6 | LiteratureAssessment 또는 CpuPreview의 정확한 버전과 필수 검토 claim IDs |
| ReviewRecord | M6 → M2 | Critic 상태·검토 범위·blocking/warning/info·원응답 참조 |
| CodeCheckReport | M2 검증 조정 → M2/M6 | 코드 검사 완료 여부, 항목별 pass/fail/not_run, 지적 |
| GatePacket | M2 → M1/연구자 | G1의 미리보기·검토·선택 후보·revision·후속 계산/비용 범위를 고정 |
| GateDecisionRequest | 연구자/M1 → M2 | 본 packet ref·revision·승인/반려·주의 확인·사유·중복 방지 키. 승인자 필드 없음 |
| ApprovalRecord | M2 → 저장/후속 실행 | 신뢰된 인간 식별·서버 시각·packet hash·결정. 이 객체만으로 실행 허가 없음 |

`PolicyProposal`은 LiteratureAssessment 안에 포함한다. proposal이 assessment를 역참조하지 않아 hash 순환을 만들지 않는다. 같은 이유로 GatePacket은 아직 생성되지 않은 ApprovalRecord를 참조하지 않는다.

## 검사기가 보장하는 것

- 등록 참조 콜백과 실제 bytes/hash, project/run/data_mode 정합. real 그래프에 fixture 포함 금지.
- 인용 segment가 요청/전달 범위에 있고, `text[start:end] == quote`인지 확인. Python Unicode code point 기준, end 미포함이다. JS에서는 `Array.from(text)` 기준으로 맞춘다.
- 기존 원문·검색·관측을 분석 단계가 덮어쓰지 못함. 추가 읽기 기록/새 LLM 관측은 해당 전달 기록과 연결. M3의 자체 human_reviewed 승격 금지.
- 실험/계산 claim은 같은 basis의 EvidenceRecord를 요구함. IC50/Kd 등 endpoint나 단위를 정책 변경 과정에서 몰래 바꾸지 못함. 인용·자료 존재만으로 과학적 지지가 입증되지는 않음.
- 변경은 허용 필드·기존 값·범위와 일치해야 함. 선택되지 않은 변경·hard/locked 변경·타입/단위/endpoint 변경은 거부.
- 미리보기는 같은 후보 집합·분자 참조·근거·계산 설정을 사용. 정책이 같으면 후보 결과도 같아야 함. 후보를 다시 생성해야 하는 변경을 고정 집합 비교로 위장하지 못함.
- 매핑 미해결/화학 invalid 후보를 eligible로 승격하지 못함. 실제 매핑 계산·분자 유효성은 M4가 별도로 입증해야 함.
- 불완전 Critic, 필수 코드 검사 미통과, 필수 정책 미정, 빈/부적합 선택은 승인 불가. 검토용 표시와 반려는 가능.
- 선택 후보의 제한/가설과 사전 검토 warning도 최종 확인에 포함. 반려는 경고 전체를 수락할 필요 없음.
- 본 packet 또는 입력 revision이 현재와 다르면 오래된 결정 거부. 승인자와 시간은 요청 body가 아닌 신뢰된 M2 인자로 받음.

검사하지 않는 것: 원문 의미·약학적 비교/컷오프 정당성, LLM 실제 API 성공/토큰, 분자/원자 매핑 실측, 비용 견적의 진위·만료·잔액, 인증 세션 진위, 영속 트랜잭션·경합·중복 결정 처리. raw artifact 내부 임의 구조까지 재귀 검증하지 않으므로 신뢰된 저장소/생산자의 역할을 없애지 않는다.

## 사용할 함수

```python
reader = ReviewReader(project_scoped_read_bytes)
reader.verify(packet_ref)                  # graph/contract check, no execution
readiness = reader.readiness(packet_ref)    # review_ready / approval_contract_ready
record_candidate = validate_decision(
    decision_request, reader,
    current_packet_ref=current_packet_ref,
    current_revision=current_revision,
    authenticated_actor=server_actor,
    decided_at=server_timestamp,
)
```

이 예시는 해당 폴더를 import할 수 있게 설정한 오프라인 호출이다. 제품 공통 import/JobSpec dispatch는 아직 연결하지 않았다. `record_candidate`는 저장 후보 객체이며 DB에 승인 기록을 남기거나 정책을 바꾸지 않는다. 현재 상태를 읽고 이 함수를 호출한 뒤 별도 트랜잭션으로 저장하면 경합이 생기므로, 후속 M2는 아래 인터페이스 문서의 원자적 compare-and-swap 규칙을 함께 구현해야 한다.
