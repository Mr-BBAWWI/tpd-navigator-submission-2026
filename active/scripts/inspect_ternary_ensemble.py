"""Reference-free diagnostics for complete six-candidate ternary ensembles."""
from __future__ import annotations

import argparse
import heapq
import html
import json
import os
import sys
from collections import Counter
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
from rdkit import Chem

from packages.science import boltz_worker as worker
from packages.science import novel_ternary as novel

SEEDS = [23, 41, 61, 79, 97]
PRIORITY = ["D-99b12e64986a", "D-b39273b7a53b"]
THRESHOLDS = [0.3, 0.5, 0.7]
OUTPUT_FILES = {"summary.json", "report.html", "output-manifest.json"}


class InspectionError(ValueError):
    """An ensemble input, path, hash, or structural invariant failed."""


def _check(condition: bool, code: str) -> None:
    if not condition:
        raise InspectionError(code)


def _sha(path: Path) -> str:
    return worker.sha256_file(path)


def _json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)


def _regular(path: Path, code: str) -> Path:
    _check(path.is_file() and not path.is_symlink(), code)
    return path


def _within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except (OSError, ValueError):
        return False


def _reject_symlink_chain(path: Path, root: Path) -> None:
    _check(_within(path, root), "PATH_TRAVERSAL_REJECTED")
    current = path
    while True:
        _check(not current.is_symlink(), "SYMLINK_REJECTED")
        if current.resolve() == root.resolve():
            break
        current = current.parent


def safe_output_key(value: Any) -> str:
    """Validate a receipt output-tree key and return normalized POSIX text."""
    _check(isinstance(value, str) and value != "", "OUTPUT_KEY_INVALID")
    _check("\\" not in value and "\x00" not in value, "OUTPUT_KEY_INVALID")
    posix, windows = PurePosixPath(value), PureWindowsPath(value)
    _check(not posix.is_absolute() and not windows.is_absolute(), "OUTPUT_KEY_ABSOLUTE")
    _check(
        not windows.drive and all(part not in {"", ".", ".."} for part in posix.parts),
        "OUTPUT_KEY_TRAVERSAL",
    )
    return posix.as_posix()


def _load_json(path: Path) -> dict:
    _regular(path, "JSON_REGULAR_FILE_REQUIRED")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise InspectionError("JSON_INVALID") from error
    _check(isinstance(value, dict), "JSON_OBJECT_REQUIRED")
    return value


def _stored_source_path(value: Any) -> Path:
    """Resolve a stored source path from the repository root, not the plan directory."""
    _check(isinstance(value, str) and value != "", "SOURCE_PATH_INVALID")
    supplied = Path(value)
    resolved = supplied if supplied.is_absolute() else ROOT / supplied
    return _regular(resolved, "SOURCE_FILE_MISSING")


def _xyz(row: dict) -> np.ndarray:
    return np.asarray(row["xyz"], dtype=float)


def _heavy(row: dict) -> bool:
    symbol = str(row.get("element", row.get("type_symbol", ""))).strip().upper()
    return symbol not in {"H", "D"}


def protein_contact_signature(proteins: dict, cutoff: float = 5.0) -> dict:
    """Return target/E3 heavy-atom residue contacts, including the cutoff boundary."""
    target = [(key, row) for key, row in proteins.items() if key[0] == "target" and _heavy(row)]
    e3 = [(key, row) for key, row in proteins.items() if key[0] == "e3" and _heavy(row)]
    minima: dict[tuple[Any, Any], float] = {}
    if e3:
        e3_xyz = np.asarray([_xyz(row) for _, row in e3])
        for tkey, trow in target:
            distances = np.linalg.norm(e3_xyz - _xyz(trow), axis=1)
            for index in np.flatnonzero(distances <= cutoff):
                pair = (tkey[1], e3[int(index)][0][1])
                minima[pair] = min(minima.get(pair, float("inf")), float(distances[index]))
    ordered = sorted(minima, key=lambda pair: (str(pair[0]), str(pair[1])))
    contacts = [
        {
            "label_seq_target": pair[0],
            "label_seq_e3": pair[1],
            "minimum_distance_A": minima[pair],
        }
        for pair in ordered
    ]
    keys = [f"{item['label_seq_target']}|{item['label_seq_e3']}" for item in contacts]
    return {
        "definition": f"protein-protein heavy-atom residue pair distance <= {cutoff} A",
        "count": len(keys),
        "keys": keys,
        "full_contact_array": contacts,
    }


def site_signature(proteins: dict, ligand: dict, graph: dict, role: str) -> dict:
    """Return the role's own protein-residue/ligand-map contacts."""
    group = "warhead_maps" if role == "target" else "recruiter_maps"
    map_to_name = graph["atom_map_to_predicted_atom_name"]
    mapped = [
        (int(atom_map), map_to_name[str(atom_map)])
        for atom_map in graph["atom_roles"][group]
    ]
    minima: dict[tuple[Any, int, str], float] = {}
    atoms = [(key, row) for key, row in proteins.items() if key[0] == role and _heavy(row)]
    for atom_map, name in mapped:
        ligand_point = _xyz(ligand[name])
        for key, row in atoms:
            distance = float(np.linalg.norm(_xyz(row) - ligand_point))
            if distance <= 4.0:
                item = (key[1], atom_map, name)
                minima[item] = min(minima.get(item, float("inf")), distance)
    ordered = sorted(minima, key=lambda item: (str(item[0]), item[1], item[2]))
    contacts = [
        {
            "protein_role": role,
            "label_seq_id": item[0],
            "ligand_atom_map": item[1],
            "ligand_atom_name": item[2],
            "minimum_distance_A": minima[item],
        }
        for item in ordered
    ]
    keys = [f"{role}|{item[0]}|{item[1]}|{item[2]}" for item in ordered]
    return {
        "definition": (
            f"{role} protein residue x {group.removesuffix('_maps')} ligand atom "
            "distance <= 4.0 A"
        ),
        "count": len(keys),
        "keys": keys,
        "full_contact_array": contacts,
    }


def interface_clashes(proteins: dict) -> dict:
    """Return descriptive target/E3 clash counts and an explicitly global denominator."""
    target = [(key, row) for key, row in proteins.items() if key[0] == "target" and _heavy(row)]
    e3 = [(key, row) for key, row in proteins.items() if key[0] == "e3" and _heavy(row)]
    table = Chem.GetPeriodicTable()

    def radius(row: dict) -> float:
        symbol = str(row.get("element", row.get("type_symbol", "C"))).strip().title()
        try:
            number = table.GetAtomicNumber(symbol)
            return float(table.GetRvdw(number)) if number else 1.7
        except Exception:
            return 1.7

    e3_xyz = np.asarray([_xyz(row) for _, row in e3]) if e3 else np.empty((0, 3))
    e3_radii = np.asarray([radius(row) for _, row in e3])
    vdw_count = 0
    raw_count = 0
    closest: list[tuple[float, int, dict]] = []
    serial = 0
    for tkey, trow in target:
        distances = np.linalg.norm(e3_xyz - _xyz(trow), axis=1)
        thresholds = 0.75 * (radius(trow) + e3_radii)
        vdw_count += int(np.sum(distances < thresholds))
        raw_count += int(np.sum(distances < 2.0))
        nearest_count = min(20, len(distances))
        nearest = (
            np.argpartition(distances, nearest_count - 1)[:nearest_count]
            if nearest_count
            else []
        )
        for index in nearest:
            ekey = e3[int(index)][0]
            item = {
                "distance_A": float(distances[index]),
                "target": {"label_seq_id": tkey[1], "atom_name": tkey[2]},
                "e3": {"label_seq_id": ekey[1], "atom_name": ekey[2]},
                "vdw_0.75_sum_threshold_A": float(thresholds[index]),
            }
            entry = (-float(distances[index]), serial, item)
            serial += 1
            if len(closest) < 20:
                heapq.heappush(closest, entry)
            elif entry > closest[0]:
                heapq.heapreplace(closest, entry)
    top = [entry[2] for entry in sorted(closest, key=lambda item: (-item[0], item[1]))]
    total_heavy_atoms = len(target) + len(e3)
    return {
        "vdw_definition": "distance < 0.75 times van der Waals radii sum; descriptive",
        "raw_definition": "distance < 2.0 A; descriptive",
        "vdw_pair_count": vdw_count,
        "raw_under_2A_pair_count": raw_count,
        "total_protein_heavy_atom_count": total_heavy_atoms,
        "vdw_pairs_per_total_heavy_atom": (
            vdw_count / total_heavy_atoms if total_heavy_atoms else None
        ),
        "closest_20_atom_pairs": top,
    }


def jaccard(first: Iterable[str], second: Iterable[str]) -> float | None:
    a, b = set(first), set(second)
    union = a | b
    return len(a & b) / len(union) if union else None


def complete_link_clusters(items: dict[Any, Iterable[str]], threshold: float) -> list[list[Any]]:
    """Deterministic grouping in which every pair in each cluster must pass."""
    keys = sorted(items, key=lambda value: str(value))
    sets = {key: set(items[key]) for key in keys}
    clusters: list[list[Any]] = []
    for key in keys:
        for cluster in clusters:
            scores = [jaccard(sets[key], sets[member]) for member in cluster]
            if all(score is not None and score >= threshold for score in scores):
                cluster.append(key)
                break
        else:
            clusters.append([key])
    return clusters


def _verified_tree(output_dir: Path, recorded: dict, root: Path) -> tuple[dict, Path]:
    _check(isinstance(recorded, dict), "OUTPUT_HASH_TABLE_INVALID")
    _check(output_dir.is_dir() and not output_dir.is_symlink(), "OUTPUT_DIRECTORY_INVALID")
    _reject_symlink_chain(output_dir, root)
    expected = {safe_output_key(key): value for key, value in recorded.items()}
    _check(len(expected) == len(recorded), "OUTPUT_KEY_DUPLICATE")
    actual: dict[str, Path] = {}
    for current, directories, files in os.walk(output_dir, followlinks=False):
        base = Path(current)
        for name in directories:
            _check(not (base / name).is_symlink(), "OUTPUT_SYMLINK_REJECTED")
        for name in files:
            path = base / name
            _regular(path, "OUTPUT_REGULAR_FILE_REQUIRED")
            _reject_symlink_chain(path, root)
            actual[path.relative_to(output_dir).as_posix()] = path
    _check(set(actual) == set(expected), "OUTPUT_FILE_SET_MISMATCH")
    hashes = {}
    for key in sorted(expected):
        _check(
            isinstance(expected[key], str) and _sha(actual[key]) == expected[key],
            "OUTPUT_HASH_MISMATCH",
        )
        hashes[key] = expected[key]
    predictions = [actual[key] for key in sorted(actual) if key.lower().endswith(".cif")]
    _check(len(predictions) == 1, "PREDICTION_CIF_COUNT_INVALID")
    return hashes, predictions[0]


def _preflight(batch_root: Path) -> list[dict]:
    _check(batch_root.is_dir() and not batch_root.is_symlink(), "BATCH_ROOT_INVALID")
    plans = sorted(batch_root.rglob("*.plan.json"))
    _check(len(plans) == 6, "PLAN_COUNT_NOT_SIX")
    prepared = []
    identities = set()
    for plan_path in plans:
        _reject_symlink_chain(plan_path, batch_root)
        plan = _load_json(plan_path)
        novel.verify_plan(plan)
        _check(plan.get("seeds") == SEEDS, "PLAN_SEEDS_NOT_EXACT")
        identity = (
            plan["candidate_graph"]["candidate_id"],
            plan["candidate_graph"]["e3_type"],
        )
        _check(identity not in identities, "CANDIDATE_DUPLICATE")
        identities.add(identity)
        run_dir = plan_path.with_name(plan_path.name[: -len(".plan.json")])
        _check(run_dir.is_dir() and not run_dir.is_symlink(), "RUN_DIRECTORY_INVALID")
        seed_data = []
        for seed in SEEDS:
            seed_dir = run_dir / f"seed-{seed}"
            receipt_path = seed_dir / "receipt.json"
            input_path = seed_dir / "novel_ternary.yaml"
            _reject_symlink_chain(receipt_path, batch_root)
            receipt = _load_json(receipt_path)
            _check(receipt.get("plan_digest") == plan["plan_digest"], "RECEIPT_PLAN_MISMATCH")
            _check(receipt.get("candidate_id") == identity[0], "RECEIPT_CANDIDATE_MISMATCH")
            _check(receipt.get("e3_type") == identity[1], "RECEIPT_E3_MISMATCH")
            _check(receipt.get("seed") == seed, "RECEIPT_SEED_MISMATCH")
            _check(
                receipt.get("status") == "completed"
                and receipt.get("exit_code") == 0
                and receipt.get("process_state") == "completed"
                and receipt.get("actual_computation") is True,
                "RECEIPT_NOT_COMPLETED",
            )
            hashes = receipt.get("hashes", {})
            _regular(input_path, "INPUT_YAML_MISSING")
            _check(
                _sha(input_path)
                == plan["boltz_input"]["sha256"]
                == hashes.get("input_yaml_sha256"),
                "INPUT_YAML_HASH_MISMATCH",
            )
            output_dir = seed_dir / "boltz_output"
            output_hashes, prediction = _verified_tree(
                output_dir,
                hashes.get("output_files_sha256", {}),
                batch_root,
            )
            seed_data.append(
                {
                    "seed": seed,
                    "receipt": receipt,
                    "receipt_path": receipt_path,
                    "input_path": input_path,
                    "prediction": prediction,
                    "output_hashes": output_hashes,
                }
            )
        prepared.append(
            {
                "plan": plan,
                "plan_path": plan_path,
                "run_dir": run_dir,
                "identity": identity,
                "seeds": seed_data,
            }
        )
    return prepared


def _frequency(seeds: list[dict], signature: str) -> dict:
    counts = Counter(key for seed in seeds for key in seed[signature]["keys"])
    denominator = len(SEEDS)
    return {
        "all_contacts": [
            {
                "key": key,
                "seed_count": counts[key],
                "frequency": counts[key] / denominator,
            }
            for key in sorted(counts)
        ],
        "shared_at_least_4_of_5": [
            {"key": key, "seed_count": counts[key]}
            for key in sorted(counts)
            if counts[key] >= 4
        ],
    }


def _analyze_candidate(item: dict) -> dict:
    plan = item["plan"]
    graph = plan["candidate_graph"]
    analyzed = []
    consistency_receipts = []
    for seed_item in item["seeds"]:
        block = novel._source_block(seed_item["prediction"])
        chains = novel._prediction_chain_mapping(block, plan)
        _, _, ligand = novel._ligand_atoms(block, graph)
        proteins = {}
        proteins.update(worker._protein_atoms(block, chains["target"], "target"))
        proteins.update(worker._protein_atoms(block, chains["e3"], "e3"))
        analyzed.append(
            {
                "seed": seed_item["seed"],
                "protein_protein": protein_contact_signature(proteins),
                "target_warhead_site": site_signature(proteins, ligand, graph, "target"),
                "e3_recruiter_site": site_signature(proteins, ligand, graph, "e3"),
                "interface_clashes": interface_clashes(proteins),
                "hashes": {
                    "receipt_sha256": _sha(seed_item["receipt_path"]),
                    "input_yaml_sha256": _sha(seed_item["input_path"]),
                    "prediction_cif_sha256": _sha(seed_item["prediction"]),
                    "confirmed_output_files_sha256": seed_item["output_hashes"],
                    "recorded_runtime_hashes": seed_item["receipt"].get("hashes", {}),
                },
            }
        )
        consistency_receipts.append(
            {
                "seed": seed_item["seed"],
                "inspection": {"prediction": str(seed_item["prediction"])},
            }
        )

    pairs = []
    for index, first in enumerate(analyzed):
        for second in analyzed[index + 1 :]:
            pairs.append(
                {
                    "seed_a": first["seed"],
                    "seed_b": second["seed"],
                    "protein_protein_contact_jaccard": jaccard(
                        first["protein_protein"]["keys"],
                        second["protein_protein"]["keys"],
                    ),
                    "target_warhead_site_jaccard": jaccard(
                        first["target_warhead_site"]["keys"],
                        second["target_warhead_site"]["keys"],
                    ),
                    "e3_recruiter_site_jaccard": jaccard(
                        first["e3_recruiter_site"]["keys"],
                        second["e3_recruiter_site"]["keys"],
                    ),
                }
            )

    consistency = []
    consistency_error = None
    try:
        consistency = novel._pairwise_consistency(consistency_receipts, plan)
    except Exception as error:
        consistency_error = {
            "type": type(error).__name__,
            "reason": str(error) or type(error).__name__,
        }

    sets = {seed["seed"]: seed["protein_protein"]["keys"] for seed in analyzed}
    clusters = [
        {
            "threshold": threshold,
            "memberships": complete_link_clusters(sets, threshold),
        }
        for threshold in THRESHOLDS
    ]
    sizes = [seed["protein_protein"]["count"] for seed in analyzed]
    ambiguity = (
        "Protein-protein contact signature is empty or small in at least one seed; "
        "contact-based grouping is ambiguous."
        if min(sizes) < 5
        else None
    )
    sources = {}
    for role in ("target", "e3"):
        stored_path = plan["sources"][role]["path"]
        sources[role] = {
            "path": stored_path,
            "sha256": _sha(_stored_source_path(stored_path)),
        }
    return {
        "candidate_id": graph["candidate_id"],
        "e3_type": graph["e3_type"],
        "priority_diagnostic": graph["candidate_id"] in PRIORITY,
        "plan": {
            "path": str(item["plan_path"]),
            "sha256": _sha(item["plan_path"]),
            "plan_digest": plan["plan_digest"],
        },
        "source_hashes": sources,
        "seeds": analyzed,
        "pairwise_contact_jaccard": pairs,
        "reference_free_pairwise_coordinate_consistency": consistency,
        "coordinate_consistency_error": consistency_error,
        "contact_frequencies": {
            name: _frequency(analyzed, name)
            for name in (
                "protein_protein",
                "target_warhead_site",
                "e3_recruiter_site",
            )
        },
        "complete_link_contact_jaccard_sensitivity": clusters,
        "contact_signature_ambiguity": ambiguity,
    }


def _escape(value: Any) -> str:
    return html.escape(str(value), quote=True)


def _render_html(summary: dict) -> str:
    overview_rows = []
    pair_sections = []
    clash_sections = []
    clustering_sections = []
    for candidate in summary["candidates"]:
        candidate_id = _escape(candidate["candidate_id"])
        seed_counts = ", ".join(
            f"{seed['seed']}: {seed['protein_protein']['count']}"
            for seed in candidate["seeds"]
        )
        shared = len(
            candidate["contact_frequencies"]["protein_protein"]["shared_at_least_4_of_5"]
        )
        overview_rows.append(
            "<tr>"
            f"<td>{candidate_id}</td>"
            f"<td>{_escape(candidate['e3_type'])}</td>"
            f"<td>{_escape(seed_counts)}</td>"
            f"<td>{_escape(shared)}</td>"
            "</tr>"
        )

        pair_rows = []
        for pair in candidate["pairwise_contact_jaccard"]:
            pair_rows.append(
                "<tr>"
                f"<td>{_escape(pair['seed_a'])}</td>"
                f"<td>{_escape(pair['seed_b'])}</td>"
                f"<td>{_escape(pair['protein_protein_contact_jaccard'])}</td>"
                f"<td>{_escape(pair['target_warhead_site_jaccard'])}</td>"
                f"<td>{_escape(pair['e3_recruiter_site_jaccard'])}</td>"
                "</tr>"
            )
        pair_sections.append(
            f"<h3>{candidate_id}</h3>"
            "<table><thead><tr><th>Seed A</th><th>Seed B</th>"
            "<th>Protein-protein</th><th>Target-warhead</th><th>E3-recruiter</th>"
            "</tr></thead><tbody>"
            + "".join(pair_rows)
            + "</tbody></table>"
        )

        clash_rows = []
        for seed in candidate["seeds"]:
            clash = seed["interface_clashes"]
            clash_rows.append(
                "<tr>"
                f"<td>{_escape(seed['seed'])}</td>"
                f"<td>{_escape(clash['vdw_pair_count'])}</td>"
                f"<td>{_escape(clash['raw_under_2A_pair_count'])}</td>"
                f"<td>{_escape(clash['total_protein_heavy_atom_count'])}</td>"
                f"<td>{_escape(clash['vdw_pairs_per_total_heavy_atom'])}</td>"
                "</tr>"
            )
        clash_sections.append(
            f"<h3>{candidate_id}</h3>"
            "<table><thead><tr><th>Seed</th><th>VDW clash pairs</th>"
            "<th>Pairs under 2 Å</th><th>Total protein heavy atoms</th>"
            "<th>VDW pairs per total heavy atom</th></tr></thead><tbody>"
            + "".join(clash_rows)
            + "</tbody></table>"
        )

        cluster_rows = []
        for sensitivity in candidate["complete_link_contact_jaccard_sensitivity"]:
            memberships = json.dumps(sensitivity["memberships"], ensure_ascii=False)
            cluster_rows.append(
                "<tr>"
                f"<td>{_escape(sensitivity['threshold'])}</td>"
                f"<td>{_escape(memberships)}</td>"
                "</tr>"
            )
        clustering_sections.append(
            f"<h3>{candidate_id}</h3>"
            "<table><thead><tr><th>Jaccard threshold</th>"
            "<th>Complete-link memberships</th></tr></thead><tbody>"
            + "".join(cluster_rows)
            + "</tbody></table>"
        )

    limitations = "".join(
        f"<li>{_escape(value)}</li>" for value in summary["limitations"]
    )
    style = (
        "body{font-family:sans-serif;margin:2rem}"
        "table{border-collapse:collapse;margin-bottom:1rem}"
        "th,td{border:1px solid #888;padding:.35rem;text-align:left;vertical-align:top}"
        "th{background:#eee}code{white-space:pre-wrap}"
    )
    return (
        "<!doctype html><html><head><meta charset=\"utf-8\">"
        "<title>Ternary ensemble diagnostic</title>"
        f"<style>{style}</style></head><body>"
        "<h1>Reference-free ternary ensemble diagnostic</h1>"
        "<p>Technical descriptive results only; no topology is selected.</p>"
        "<table><thead><tr><th>Scientific approved</th>"
        "<th>Human review performed</th><th>Raw seed results preserved</th>"
        "</tr></thead><tbody><tr>"
        f"<td>{_escape(summary['scientific_approved'])}</td>"
        f"<td>{_escape(summary['human_review_performed'])}</td>"
        f"<td>{_escape(summary['raw_seed_result_count'])}</td>"
        "</tr></tbody></table>"
        "<h2>Candidate overview</h2>"
        "<table><thead><tr><th>Candidate</th><th>E3</th>"
        "<th>Protein-protein contacts by seed</th><th>Shared ≥4/5</th>"
        "</tr></thead><tbody>"
        + "".join(overview_rows)
        + "</tbody></table>"
        "<h2>Pairwise contact Jaccard metrics</h2>"
        + "".join(pair_sections)
        + "<h2>Clash sensitivity</h2>"
        + "".join(clash_sections)
        + "<h2>Complete-link clustering sensitivity</h2>"
        + "".join(clustering_sections)
        + "<h2>Limitations</h2><ul>"
        + limitations
        + "</ul></body></html>"
    )


def inspect(batch_root: Path, output: Path) -> dict:
    batch_root, output = Path(batch_root), Path(output)
    _check(not output.exists() and not output.is_symlink(), "OUTPUT_ALREADY_EXISTS")
    prepared = _preflight(batch_root)
    candidates = [_analyze_candidate(item) for item in prepared]
    candidates.sort(
        key=lambda item: (
            PRIORITY.index(item["candidate_id"])
            if item["candidate_id"] in PRIORITY
            else len(PRIORITY),
            item["candidate_id"],
            item["e3_type"],
        )
    )
    raw_seed_result_count = sum(len(candidate["seeds"]) for candidate in candidates)
    _check(
        len(candidates) == 6
        and raw_seed_result_count == 30
        and all([seed["seed"] for seed in candidate["seeds"]] == SEEDS for candidate in candidates),
        "RAW_ENSEMBLE_NOT_COMPLETE",
    )
    summary = {
        "format": "tpd-ternary-ensemble-diagnostic/1.0",
        "reference_free": True,
        "actual_computation": True,
        "scientific_approved": False,
        "human_review_performed": False,
        "candidate_count": 6,
        "raw_seed_result_count": raw_seed_result_count,
        "seeds": SEEDS,
        "prioritized_diagnostics": [
            {
                "candidate_id": candidate,
                "present": any(item["candidate_id"] == candidate for item in candidates),
            }
            for candidate in PRIORITY
        ],
        "candidates": candidates,
        "clustering": {
            "method": "deterministic complete-link; every pair must meet each tested threshold",
            "thresholds": THRESHOLDS,
            "selected_optimum": None,
            "interpretation": "Sensitivity diagnostic only; no topology is selected or endorsed.",
        },
        "limitations": [
            "Reference-free prediction consistency is not reference accuracy.",
            "Contact and clash thresholds are descriptive geometric definitions, not scientific acceptance cutoffs.",
            "No binding, degradation, efficacy, ubiquitination, or topology acceptance is inferred.",
            "Empty contact unions produce null Jaccard values.",
        ],
    }
    output.mkdir(parents=True, exist_ok=False)
    summary_path = output / "summary.json"
    html_path = output / "report.html"
    summary_path.write_text(_json(summary) + "\n", encoding="utf-8", newline="\n")
    html_path.write_text(_render_html(summary), encoding="utf-8", newline="\n")
    manifest = {
        "format": "tpd-ternary-ensemble-output-manifest/1.0",
        "exact_output_file_set": sorted(OUTPUT_FILES),
        "hashed_files": {
            "summary.json": _sha(summary_path),
            "report.html": _sha(html_path),
        },
        "self_hash_omitted": True,
    }
    manifest_path = output / "output-manifest.json"
    manifest_path.write_text(_json(manifest) + "\n", encoding="utf-8", newline="\n")
    _check(
        {path.name for path in output.iterdir()} == OUTPUT_FILES,
        "GENERATED_OUTPUT_SET_INVALID",
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    inspect(args.batch_root, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
