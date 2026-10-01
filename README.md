# TPD Navigator 공개 라이브 비교실

`index.html`은 이제 저장 결과만 재생하는 정적 화면이 아니라 라이브 프런트엔드입니다. 다음 예시 질의를 선택해 실행할 수 있습니다.

1. CRBN·VHL 후보와 linker를 비교해 주세요.
2. 도킹 실패와 pose 제외 원인을 분석해 주세요.
3. 14개 기준을 확인하고 후속 실행 계획을 제안해 주세요.

예시 카드를 선택하고 Parent·E3 필터를 선택한 뒤 실행 버튼을 누릅니다. 브라우저는 `public-config.json`의 `apiBase`인 `https://tpd-navigator-live.wlgudyun.chatgpt.site`로 POST 요청을 보냅니다. 클라우드 서버는 저장된 근거를 도구로 분석하고 새 LLM 응답을 생성하며, 화면에는 요청·사용량·trace가 표시됩니다.

브라우저에는 API 키가 없습니다. 비공개 provider 자격 증명은 서버에만 있으며 응답이나 다운로드에 포함되지 않습니다. 라이브 실행에는 인터넷 연결과 구성되어 작동 중인 클라우드 서비스가 필요합니다.

권장 접속 주소는 <https://mr-bbawwi.github.io/tpd-navigator-submission-2026/>입니다. `python -m http.server`로 로컬 프런트엔드를 제공할 수는 있지만, 백엔드 origin allowlist가 승인되지 않은 로컬 origin을 거부할 수 있습니다. 임의의 `file://` 실행이나 모든 로컬 환경에서의 동작을 보장하지 않습니다.

## 저장 결과 재생

[`stored-results.html`](stored-results.html)은 이전 정적 공개본의 저장 결과를 재생합니다. 직접 열거나 오프라인으로 볼 수 있으며 새 LLM·GPU·외부 서비스를 호출하지 않습니다.

라이브와 저장 결과 모두 새로운 GPU 구조 예측을 수행하지 않으며 과학적 승인 또는 인간 전문가 승인을 새로 부여하지 않습니다. 보충 holdout은 같은 6BOY 개발 구조의 fresh-seed 검증으로, formal 14개 기준을 변경하거나 새 target·독립 구조·training holdout·통계적 일반화를 입증하지 않습니다.

후보 키와 job/assessment ID는 정제된 저장 receipt 안의 불투명 공개 식별자입니다. 공개 receipt hash는 원본 canonical digest가 아니라 정제된 파일의 hash입니다.

원본 source index SHA256: `86928b0a375896e4c84185b1574148c3ae71872d45c96fcba0e218b687039029`  
V1 공개 정적 `index.html` SHA256 (현재 `stored-results.html`과 동일): `0f0cb90c707675496287a3fb23af8c160e111936e5c46eef567e19e944e15b19`  
원본 science summary SHA256: `fd8bb27fd7170c8a39164a84fb6433140b5e2d370cbca7a1aad920a52d479750`

이 provenance·hash·개인정보 설명은 정제되어 저장된 receipt 범위에만 적용됩니다. 개인 메타데이터, 원시 API 요청·응답, 자격 증명, 내부 경로와 private asset은 저장 공개본에 포함하지 않습니다. 새 라이브 브라우저 응답과 request·usage·trace는 화면에 표시하고 다운로드할 수 있지만 비공개 provider 자격 증명은 포함하지 않습니다.

저장소 루트의 `live-release-manifest.json`은 새 라이브 asset manifest입니다. 저장소 루트의 `source-export-manifest.json`은 V1 source history를 보존하며, `submission/public-demo/manifest.json`은 V1 static history를 보존합니다. V1 static manifest는 새 라이브 `index.html`을 검증하지 않습니다.
