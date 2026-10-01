"""Offline baseline by default; source-backed science requires --suite science."""
import argparse
import hashlib
import importlib.util
import json
import platform
import shutil
import subprocess
import sys
import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=["base", "science"], default="base")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or ROOT / ".localdata/verification" / (args.suite + "-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + ".json")
    if output.exists():
        parser.error("Output exists; historical verification is never overwritten.")
    if args.suite == "science":
        missing = [name for name in ("rdkit", "gemmi", "numpy", "yaml") if importlib.util.find_spec(name) is None]
        missing += [name for name in ("6HAZ.cif", "6HAY.cif", "6HAX.cif", "FX5.cif", "FX8.cif", "FWZ.cif", "article.xml")
                    if not (ROOT / ".localdata/smarca2/raw" / name).is_file()]
        if missing:
            parser.error("Science prerequisites missing: " + ", ".join(missing) + "; see docs/b_smarca2_start.md")
    folder = "tests" if args.suite == "base" else "tests_science"
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.discover(str(ROOT / folder)))
    checks = []
    if args.suite == "base":
        with tempfile.TemporaryDirectory(prefix="tpd-contracts-") as temp:
            scratch = Path(temp)
            for name in ("contracts", "packages", "verification"):
                shutil.copytree(ROOT / name, scratch / name, ignore=shutil.ignore_patterns("__pycache__"))
            for script in ("contracts/v0.1.0/validate_examples.py", "contracts/payloads/v0.1.0/validate_examples.py", "contracts/review/v0.1.0/validate_examples.py"):
                run = subprocess.run([sys.executable, str(scratch / script)], capture_output=True, text=True, cwd=scratch)
                checks.append({"script": script, "exit_code": run.returncode, "output": run.stdout.strip(), "stderr": run.stderr.strip()})
                print(run.stdout.strip() or run.stderr.strip())
    paths = [p for name in ("packages", "apps", "tests", "tests_science", "scripts") for p in (ROOT / name).rglob("*")
             if p.is_file() and p.suffix in {".py", ".js", ".css", ".html"}]
    passed = result.wasSuccessful() and not result.skipped and all(c["exit_code"] == 0 for c in checks)
    report = {"verified_at": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(), "suite": args.suite,
              "passed": passed, "unittest": {"run": result.testsRun, "failures": len(result.failures), "errors": len(result.errors), "skipped": len(result.skipped)},
              "contracts": checks, "live_llm_calls": 0,
              "code_sha256": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Suite {args.suite}: {result.testsRun} tests; passed={passed}; {output}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
