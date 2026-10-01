# 공개 라이브 실행 데모

이 디렉터리는 동결된 공개 영수증을 근거로 세 가지 제한된 질의를 실제 모델에 실행하는 데모이다. 저장된 정적 재생 결과는 계속 별도로 제공된다.

## 허용 범위

`POST /api/run`은 다음 작업만 허용한다.

- `candidate_compare`
- `failure_analysis`
- `acceptance_plan`

요청 형식:

```json
{
  "operation": "candidate_compare | failure_analysis | acceptance_plan",
  "filters": {
    "parent_id": "선택 사항",
    "e3_type": "CRBN 또는 VHL",
    "candidate_ids": ["기존 후보 ID, 1~3개"]
  }
}
```

임의 프롬프트, 임의 표적, 임의 SMILES, PDF 업로드, 신규 웹 검색, 신규 GPU 작업, 신규 도킹은 허용하지 않는다. 사람 승인, M2 기록, 상태 승격도 수행하지 않는다.

## 로컬 실행

요구 사항은 Fetch API를 포함한 Node.js 22 이상이다. `DACON_API_KEY`와 선택적인 `ALLOWED_ORIGINS`는 사용자의 비밀 관리 수단을 통해 실행 프로세스 환경에 이미 설정되어 있어야 한다. 이 실행기는 `.env` 파일을 읽지 않는다.

공개 Git 저장소 체크아웃 루트에서 다음을 실행한다. `active/submission/live-demo`를 로컬 프로젝트 루트로 직접 체크아웃한 Sites 사용자는 첫 `cd`를 생략할 수 있다.

```bash
cd active/submission/live-demo
npm install
npm test
node build.mjs --self-test
npm run build
npm run dev
```

개발 서버 주소는 `http://127.0.0.1:8826`이다. Windows에서는 심볼릭 링크 관련 자체 검사가 건너뛰어질 수 있다. 빌더는 알려진 비밀·민감 자료의 포함 여부를 점검하지만 완전한 탐지를 보장하지 않는다.

## 실행 가능한 요청 예시

상태 확인:

```bash
curl -sS http://127.0.0.1:8826/api/health
```

공개 카탈로그 확인:

```bash
curl -sS http://127.0.0.1:8826/api/catalog
```

세 허용 작업 실행:

```bash
curl -sS -X POST http://127.0.0.1:8826/api/run -H 'Origin: http://127.0.0.1:8826' -H 'Content-Type: application/json' --data '{"operation":"candidate_compare","filters":{}}'
curl -sS -X POST http://127.0.0.1:8826/api/run -H 'Origin: http://127.0.0.1:8826' -H 'Content-Type: application/json' --data '{"operation":"failure_analysis","filters":{}}'
curl -sS -X POST http://127.0.0.1:8826/api/run -H 'Origin: http://127.0.0.1:8826' -H 'Content-Type: application/json' --data '{"operation":"acceptance_plan","filters":{}}'
```

`Origin` 없는 POST는 거부된다. `ALLOWED_ORIGINS`는 정확한 출처를 쉼표로 구분한 환경값이며 동일 출처는 허용된다. CORS 출처 검사는 인증이 아니다.

## 고정 실행 조건

- 모델: `gpt-5.6-sol`
- 공급자: 고정 DACON Responses 엔드포인트
- `store: false`
- 최대 출력: 1,800
- 컨텍스트: 최대 15,000바이트
- 요청 본문: 최대 4,096바이트
- 런타임 제한: 85초
- 브라우저 제한: 180초
- 자동 재시도 및 답변 캐시 없음

API 키는 서버 환경에만 존재하며 클라이언트에 전달하지 않는다. 브라우저의 수동 취소는 대기 중단일 뿐이며, 공급자 작업이 이미 실행되었거나 과금되었을 수 있어 환불을 보장하지 않는다.

IP당 분당 4회와 동시 실행 2회 제한은 격리 인스턴스별 메모리 기반 최선 노력 제한이다. 전역 하드 쿼터나 강한 남용 방지가 아니다. 대규모 공개 트래픽에는 영속적인 전역 예산·속도 제어가 추가로 필요하다.

## 데이터와 표시 의미

카탈로그는 공개 후보 324개와 부모 9개를 제공한다. `candidate_compare`의 산술과 필터 집계는 결정론적이며, 표시는 최대 3개이다. 선택은 사전식 순서를 사용하고 가능하면 같은 부모의 CRBN·VHL을 함께 보여 주며 효능 순위가 아니다.

`failure_analysis`의 동결 집계는 시도 234, 완료 231, 실패 3, 포즈 제외 91이다. E3 또는 후보 필터는 부모 집합만 선택하며 E3별 실패 원인을 귀속하지 않는다. 누락된 실패 원인은 생성하지 않는다.

`acceptance_plan`은 형식 기준 14개에 대해 통과 3, 실패 4, 보류 7인 전역 스냅샷을 사용한다. 필터는 후보 문맥 수만 바꾸며 기준 상태를 바꾸지 않는다. 실제 수령 자료에서 포함된 동일 6BOY 홀드아웃 3/5는 하드코딩 대체값이 아니며 형식 승격을 일으키지 않는다.

응답에는 정규 도구 출력, 스냅샷 SHA, 도구 추적, 새 UTC 시각, 요청 ID, 요청·반환 모델, 응답 ID, 사용량, 답변 근거가 포함된다. 엄격한 파싱에 실패한 답변은 실패 처리한다. 알려진 실패 사용량은 보존하고, 사용량이 미제공·확인 불가이면 JSON에서 생략한 채 0 또는 추정치로 채우지 않는다.

## 배포 및 공개 자료 상태

Sites 프로젝트 ID는 `appgprj_6abe82b4c5788191b7216f9bf0c1b3ea`, 예상 출처는 `https://tpd-navigator-live.wlgudyun.chatgpt.site`이며, 이 정보만으로 현재 라이브 상태를 주장하지 않는다. 원본 `.openai/hosting.json`은 Sites 정체성 정보를 그대로 보존하며, 빌드 결과는 `dist/server/index.js`, `dist/client/*`, `dist/.openai/hosting.json`이고 배포용 `dist/.openai/hosting.json`에는 지원되는 `d1`/`r2` 키와 각각의 `null` 값만 포함한다. 문서 작성과 구현 검사만으로 공개 배포 성공을 주장하지 않는다. 실제 배포·예제 3개 실행 여부와 현행 파일 해시는 별도의 최신 `live-release-manifest.json` 및 운영 검증 기록을 확인한다.

공개 저장소는 `https://github.com/Mr-BBAWWI/tpd-navigator-submission-2026`이다. 경로는 이 공개 저장소 기준이다. 기존 `source-export-manifest.json`은 원래 V1 동결 페이로드 792개만 다루며 새 라이브 파일을 포함한다고 주장할 수 없다. 새 검토 대상 라이브 페이로드는 루트에서 기계적으로 생성되는 `live-release-manifest.json`으로 검증한다. V1 매니페스트는 수정하지 않는다.

저장 결과 공개 주소로 계획된 위치는 `https://mr-bbawwi.github.io/tpd-navigator-submission-2026/stored-results.html`이다. 이 주소의 현재 공개 여부는 최신 `live-release-manifest.json` 및 운영 검증 기록으로 판단한다. 키 유효성이나 공급자 자원은 종료될 수 있으며, 그 경우 라이브 서비스는 실패할 수 있지만 저장된 재생 결과는 독립적으로 유지된다.
