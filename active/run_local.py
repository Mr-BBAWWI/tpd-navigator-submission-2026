"""Start the I2 local app; intentionally binds only the loopback interface."""
import argparse

import uvicorn

from apps.api.main import create_app


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8767)
    parser.add_argument("--enable-llm", action="store_true", help="Enable explicit local literature analysis using the competition API")
    args = parser.parse_args()
    factory = None
    if args.enable_llm:
        from packages.agents.provider import DaconProvider, load_key
        load_key()
        factory = lambda: DaconProvider(load_key())
    uvicorn.run(create_app(provider_factory=factory), host="127.0.0.1", port=args.port, log_level="info")
