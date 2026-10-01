#!/usr/bin/env python3
"""Render a self-contained, offline research-campaign dashboard."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_json(path: Path) -> tuple[bytes, object]:
    raw = path.read_bytes()
    return raw, json.loads(raw.decode("utf-8"))


def json_script(value: object) -> str:
    text = json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return text.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")


def require_evidence_root(value: object) -> dict:
    if not isinstance(value, dict):
        raise ValueError("evidence.json root must be an object")
    required = {
        "candidates": list,
        "analogs": list,
        "parents": list,
        "campaign_counts_recomputed_from_records": dict,
        "assessment_snapshot": dict,
    }
    for key, expected in required.items():
        if key not in value or not isinstance(value[key], expected):
            raise ValueError(f"evidence.{key} must be a {expected.__name__}")
    unexpected = {"candidate_examples", "analog_examples"} & set(value)
    if unexpected:
        raise ValueError(f"unexpected prompt-illustration evidence keys: {sorted(unexpected)}")
    if any(not isinstance(row, dict) for row in value["candidates"]):
        raise ValueError("evidence.candidates entries must be objects")
    if any(not isinstance(row, dict) for row in value["analogs"]):
        raise ValueError("evidence.analogs entries must be objects")
    return value


def render_graphs(evidence: dict, image_dir: Path) -> tuple[dict, list]:
    graph_map, files, smiles_files = {}, [], {}
    candidates = evidence["candidates"]
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
    except Exception:
        for row in candidates:
            if isinstance(row, dict) and isinstance(row.get("candidate_key"), str):
                graph_map[row["candidate_key"]] = {"file": None, "reason": "RDKit 사용 불가"}
        return graph_map, files

    image_dir.mkdir()
    unique_count = 0
    for row in candidates:
        if not isinstance(row, dict) or not isinstance(row.get("candidate_key"), str):
            continue
        key, smiles = row["candidate_key"], row.get("canonical_smiles")
        if not isinstance(smiles, str) or not smiles:
            graph_map[key] = {"file": None, "reason": "canonical SMILES 없음"}
            continue
        if smiles in smiles_files:
            graph_map[key] = {"file": smiles_files[smiles], "reason": None}
            continue
        if unique_count >= 400:
            graph_map[key] = {"file": None, "reason": "400개 고유 그래프 렌더 한도 초과"}
            continue
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("invalid smiles")
            rdDepictor.Compute2DCoords(mol)
            drawer = rdMolDraw2D.MolDraw2DSVG(560, 260)
            opts = drawer.drawOptions()
            opts.clearBackground = False
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText().encode("utf-8")
            name = sha256((key + "\0" + smiles).encode("utf-8"))[:24] + ".svg"
            (image_dir / name).write_bytes(svg)
            rel = "images/" + name
            smiles_files[smiles] = rel
            graph_map[key] = {"file": rel, "reason": None}
            files.append({"path": rel, "sha256": sha256(svg), "bytes": len(svg),
                          "kind": "RDKit 2D identity graph", "canonical_smiles": smiles})
            unique_count += 1
        except Exception:
            graph_map[key] = {"file": None, "reason": "canonical SMILES 2D 렌더 실패"}
    return graph_map, files


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    campaign = args.campaign_dir.resolve()
    output = args.output_dir.resolve()
    required = {
        "evidence": campaign / "evidence.json",
        "rounds": campaign / "rounds" / "result.json",
        "tokens": campaign / "token-report.json",
    }
    for path in required.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    if output.exists():
        raise FileExistsError(str(output))
    output.mkdir(parents=True)
    raw_dir = output / "raw"
    raw_dir.mkdir()

    parsed, inputs = {}, []
    raw_names = {"evidence": "evidence.json", "rounds": "rounds-result.json", "tokens": "token-report.json"}
    for name, source in required.items():
        raw, value = load_json(source)
        parsed[name] = value
        target = raw_dir / raw_names[name]
        target.write_bytes(raw)
        inputs.append({"name": name, "source_leaf": source.name, "output": "raw/" + target.name,
                       "sha256": sha256(raw), "bytes": len(raw)})

    parsed["evidence"] = require_evidence_root(parsed["evidence"])
    if not isinstance(parsed["rounds"], dict):
        raise ValueError("rounds result root must be an object")
    if not isinstance(parsed["tokens"], dict):
        raise ValueError("token report root must be an object")

    repo = Path(__file__).resolve().parents[1]
    asset_dir = repo / "apps" / "web"
    sources = {}
    source_paths = [Path(__file__).resolve()] + [
        asset_dir / name for name in
        ("research-campaign.template.html", "research-campaign.js", "research-campaign.css")
    ]
    for path in source_paths:
        data = path.read_bytes()
        name = "scripts/render_research_campaign.py" if path == Path(__file__).resolve() else path.name
        sources[name] = {"sha256": sha256(data), "bytes": len(data)}
    css = (asset_dir / "research-campaign.css").read_text(encoding="utf-8")
    js = (asset_dir / "research-campaign.js").read_text(encoding="utf-8")
    template = (asset_dir / "research-campaign.template.html").read_text(encoding="utf-8")

    graphs, image_files = render_graphs(parsed["evidence"], output / "images")
    payload = {**parsed, "graphs": graphs, "rendered_at": datetime.now(timezone.utc).isoformat()}
    html = template.replace("{{CSS}}", css).replace("{{DATA}}", json_script(payload)).replace("{{JS}}", js)
    index_bytes = html.encode("utf-8")
    (output / "index.html").write_bytes(index_bytes)
    manifest = {
        "format": "research-dashboard-manifest/1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "inputs": inputs,
        "renderer_sources": sources,
        "outputs": ([{"path": "index.html", "sha256": sha256(index_bytes), "bytes": len(index_bytes)}]
                    + [{"path": item["output"], "sha256": item["sha256"], "bytes": item["bytes"],
                        "kind": "raw input copy"} for item in inputs]
                    + image_files),
        "notes": ["Images are deterministic RDKit 2D identity graphs, not experimental or steric poses.",
                  "Only fixed campaign input files were read."],
    }
    manifest_bytes = (json.dumps(manifest, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    (output / "manifest.json").write_bytes(manifest_bytes)
    print(json.dumps({"ok": True, "output": str(output), "inputs": len(inputs),
                      "graphs": len(image_files), "manifest_sha256": sha256(manifest_bytes)}, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({"ok": False, "error": type(exc).__name__, "message": str(exc)}, ensure_ascii=False), file=sys.stderr)
        raise SystemExit(1)
