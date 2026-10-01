> 경로 정리: 2026-09-22. 아래는 기존 계약/시제품의 작성·검사 이력이다. 목표 제품과 최신 분업은 [INTERFACES](<../../../docs/INTERFACES.md>), 실제 구현 범위는 [STATUS](<../../../docs/STATUS.md>)를 따른다. 과거 A/B 배분·public_explore·다음 작업은 현재 배정/공개 범위/개발 재개를 뜻하지 않는다.

# 검토 중인 도메인 계약 자료

이 폴더는 [공통 봉투 v0.1.0](<../README.md>)에 연결할 도메인 명세의 검토 자료다. 실행할 payload schema로 등록된 폴더가 아니다.

- [F03 부착점 SAR 해결안](<../../../docs/design/notes/attachment_sar.md>): `assess_exit_vectors`의 입력/출력 필드, 문헌→과학→검토 연결, 원자 매핑·근거 범위·승인 규칙.
- [f03_acceptance_cases.json](<f03_acceptance_cases.json>): 위 설계를 검토할 가상 사례 17개. JSON Schema, 실제 입력/출력 payload, 실행 테스트, 실측 데이터가 아니다. `given`은 가정이고 `must`/`must_not`은 향후 구현의 기대 행동이다.

이 자료만으로 `assess_exit_vectors`의 `payload_schema_status=pending`, `live_execution_ready=false`를 바꾸지 않는다. 후속 [I1](<../payloads/v0.1.0/README.md>)에 표적·문헌 관측/원문 위치의 실제 schema를 작성했지만 부착점 평가 전체 schema와 분자 매핑은 남아 있다. I1 참조·관측 필드에 맞춰 ExitVectorRequest/Assessment를 구현하고 원문/구조 검수를 따로 수행해야 한다.

모든 사례는 새로 작성한 설계 fixture다. Study/JH 자산에서 이식하지 않았으며, 실제 평가·학습·제출 재생 묶음에 넣지 않는다. A/B는 권장 담당 트랙이고 사람 배정은 미정이다.
