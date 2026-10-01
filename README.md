# TPD Navigator — public submission

## 실행형 공개 데모

호스팅된 루트 <https://mr-bbawwi.github.io/tpd-navigator-submission-2026/>는 현재 실행형 UI입니다. 다음 예시 질의를 선택할 수 있습니다.

1. **CRBN·VHL 후보와 linker를 비교해 주세요.**
2. **도킹 실패와 pose 제외 원인을 분석해 주세요.**
3. **14개 기준을 확인하고 후속 실행 계획을 제안해 주세요.**

사용자는 예시와 선택적 필터를 고른 뒤 **실행**을 누릅니다. 서버 도구가 저장된 근거를 먼저 분석하고, 모델은 그 결과를 바탕으로 새로운 LLM 답변을 생성합니다. 실제 사용 토큰과 trace가 화면에 표시됩니다. 브라우저 API 키는 필요 없지만 인터넷 연결과 정상 작동하는 클라우드 API가 필요합니다. 새 GPU 예측이나 과학적·인적 승인은 수행하지 않습니다.

[실행형 데모의 로컬 Node.js 실행 안내](active/submission/live-demo/README.md) (Node.js 22와 서버 프로세스 환경 변수 키가 필요하며 자세한 내용은 링크된 README를 참조하십시오.)

## 이전 정적 결과 재생

호스팅된 <https://mr-bbawwi.github.io/tpd-navigator-submission-2026/stored-results.html>은 이전 정적 결과를 재생하는 화면입니다. Python 3.12 로컬 HTTP 서버와 launcher도 이 정적 재생만 제공합니다.

```bash
python -m http.server 8000 --bind 127.0.0.1 --directory active/submission/public-demo
```

그다음 `http://127.0.0.1:8000/stored-results.html`을 엽니다. 저장소 루트에서는 다음 launcher도 사용할 수 있습니다.

```bash
python launch_public_demo.py
python launch_public_demo.py --port 8080
```

Windows에서는 `launch_public_demo.cmd`를 사용할 수 있습니다. 이 정적 모드의 서버는 `127.0.0.1`에만 바인딩되며 API나 GPU를 사용하지 않습니다. 데모 HTML은 실행 가능한 외부 스크립트를 요구하지 않는 저장 결과물입니다.

## 소스 환경

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# POSIX: source .venv/bin/activate
python -m pip install -r active/requirements.txt
```

선택적 과학 계산 의존성이 필요한 경우:

```bash
python -m pip install -r active/requirements-science.txt
```

DACON 연동이 필요한 수동 실행에서만 환경 변수로 키를 설정합니다. 키를 파일이나 명령행 옵션에 넣지 마십시오.

```powershell
$env:DACON_API_KEY = '<your-key>'
```

```bash
export DACON_API_KEY='<your-key>'
```

전체 내부 재현 실행에는 공개 스냅샷에 포함되지 않은 동결된 비공개 근거 자료와 소유자 접근 권한이 필요합니다.

## 결과 해석 범위

- 14개 검토 기준의 현재 저장 상태는 **통과 3, 실패 4, 보류 7**이며 전체 과학 검토 완료를 뜻하지 않습니다.
- 324개 조립 결과는 탐색 집합입니다. 엄격 필터의 3개 후보와 native-context 비교 집합 15개를 같은 의미로 해석하면 안 됩니다.
- pKa 3개 값은 분리된 모델 예측이며 결합 상태 측정값이나 실험값이 아닙니다.
- holdout은 동일 BRD4–CRBN 6BOY 구조에서 새 seed 5개·25개 raw 모델을 계산하고, 새 seed 실행 전에 고정한 complex_iplddt 최고 선택에서 3/5가 기준을 충족한 점검입니다. 규칙은 앞선 개발 seed 결과에서 도출했으며 독립 표적 검증이나 공식 14개 기준 통과가 아닙니다.
- 후속 기술 검토에는 scope-review, native context 및 M2 관련 차이가 남아 있습니다.
- 내부 기준선에서 1202개 검사가 사용되었지만 비공개 동결 데이터가 필요하므로, 공개 clone에서 동일 결과가 보장된다고 주장하지 않습니다. 공개 export의 컴파일·해시·비밀정보 폐쇄 검사는 별도 범위입니다.

## 저장소와 무결성

공개 저장소: <https://github.com/Mr-BBAWWI/tpd-navigator-submission-2026>

`source-export-manifest.json`과 `active/submission/public-demo/manifest.json`은 현재 실행형 파일 목록이 아니라 원본 V1 소스 export와 정적 공개 화면의 이력을 기록합니다. 원본 V1의 `source-export-manifest.json`은 각 payload의 원본/공개 SHA-256, 허용된 표시명 치환 및 개인정보 메타데이터 제거 사유를 기록합니다. `source_git_commit`은 당시 작업 트리의 기준 revision이며, 실제로 export된 작업 파일은 파일별 `source_sha256`으로 식별합니다. manifest는 자기 자신의 해시 폐쇄에서 명시적으로 제외됩니다. V1의 `index.html` 해시는 현재 `stored-results.html`로 보존된 바이트와 같습니다.

루트 `live-release-manifest.json`은 새 실행형 데모와 이 README 및 launcher를 포함한 live 변경을 추적합니다. 따라서 원본 V1의 해시 폐쇄 검사가 수정된 현재 파일 전체를 검증한다고 주장하지 않습니다.

원본 V1 exporter의 보안 제한도 그대로 적용됩니다. 선택적 PDF/PPTX는 저장소 루트에서 최종 검토한 출력물만 승인된 파일명으로 exporter에 전달해야 합니다. PPTX는 압축 member를 메모리에서 제한적으로 펼쳐 비밀정보를 검사하지만, PDF는 raw byte 검사와 알려진 생성 출처에 대한 수동 검토에 의존하며 압축된 PDF 텍스트를 포괄적으로 검사한다고 주장하지 않습니다. V1에서 export된 README와 PDF도 공개 전 저장소 루트에서 최종 검토해야 합니다. 생성된 V1 export 디렉터리에 파일을 임의로 복사하면 해당 스냅샷의 독립 검증이 실패합니다.

## 제출 문서

- [한국어 본선 발표 PDF](submission-documents/presentation.pdf)
- [상세 기술서 PDF](submission-documents/technical-report.pdf)
- [편집용 PPTX](submission-documents/presentation.pptx)
