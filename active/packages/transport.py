"""Bounded reads from the two explicitly configured public data providers."""
import time
from dataclasses import dataclass
from urllib.parse import urlsplit

import httpx


class FetchError(Exception):
    def __init__(self, code, message):
        self.code = code
        super().__init__(message)


@dataclass
class Fetch:
    body: bytes
    url: str
    status: int
    headers: dict


class PublicTransport:
    hosts = {"rest.uniprot.org", "www.ebi.ac.uk"}

    def __init__(self, seconds=150):
        self.deadline = time.monotonic() + seconds
        self.calls = 0
        self.bytes_received = 0

    def get(self, url, params=None, max_bytes=2_000_000):
        parsed = urlsplit(url)
        if parsed.scheme != "https" or parsed.hostname not in self.hosts or parsed.port not in {None, 443}:
            raise FetchError("SOURCE_NOT_ALLOWED", "등록되지 않은 자료 출처입니다.")
        left = self.deadline - time.monotonic()
        if left <= 0:
            raise FetchError("TIME_LIMIT", "자료 수집 시간 한도에 도달했습니다.")
        self.calls += 1
        try:
            with httpx.Client(timeout=min(18, left), follow_redirects=False) as client:
                with client.stream("GET", url, params=params, headers={"User-Agent": "TPD-Navigator-I2-local/0.1"}) as response:
                    if response.status_code != 200:
                        raise FetchError("HTTP_" + str(response.status_code), f"자료 제공 서버 응답: HTTP {response.status_code}")
                    chunks, size = [], 0
                    for chunk in response.iter_bytes():
                        size += len(chunk)
                        self.bytes_received += len(chunk)
                        if size > max_bytes or self.bytes_received > 24_000_000:
                            raise FetchError("BYTE_LIMIT", "자료 크기가 이번 실행의 한도를 초과했습니다.")
                        if time.monotonic() > self.deadline:
                            raise FetchError("TIME_LIMIT", "자료 수집 시간 한도에 도달했습니다.")
                        chunks.append(chunk)
                    return Fetch(b"".join(chunks), str(response.url), response.status_code, dict(response.headers))
        except httpx.HTTPError as exc:
            raise FetchError("NETWORK_ERROR", "공개 자료 서버에 연결하지 못했습니다.") from exc
