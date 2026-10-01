"""Disposable PDF text worker. No OCR, figure interpretation, LLM, or network."""
import io
import json
import sys

from pypdf import PdfReader


def main():
    raw = sys.stdin.buffer.read(6_000_001)
    if len(raw) > 6_000_000:
        raise ValueError("PDF_BYTE_LIMIT")
    reader = PdfReader(io.BytesIO(raw), strict=False)
    if reader.is_encrypted:
        raise ValueError("PDF_ENCRYPTED")
    limit = min(int(sys.argv[1]), 60)
    chunks = []
    truncated = len(reader.pages) > limit
    for index, page in enumerate(reader.pages[:limit]):
        text = (page.extract_text() or "").strip()
        if len(text) > 12000:
            text = text[:12000]
            truncated = True
        if text:
            chunks.append({"page": index + 1, "text": text, "content_kind": "text"})
    sys.stdout.buffer.write(json.dumps({"chunks": chunks, "truncated": truncated, "page_count": len(reader.pages)}, ensure_ascii=False).encode("utf-8"))


if __name__ == "__main__":
    main()
