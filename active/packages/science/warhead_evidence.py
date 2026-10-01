"""Collect auditable co-crystal ligand evidence from the official RCSB APIs.

The collector deliberately makes only structural observations. Protein/ligand
heavy-atom contacts are not binding measurements, SAR, or docking evidence.
Every downloaded response is retained as a content-addressed immutable source
snapshot and is hash-verified by :func:`load_catalog`.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urlparse

import gemmi
import httpx
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from .molecules import canonical, identity, read_ccd, rows, write_sdf
from .structures import atom_sites, ligand_coordinates

SCHEMA_VERSION = "warhead-evidence-v1"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ALLOWED_HOSTS = {"search.rcsb.org", "files.rcsb.org"}
MAX_RESPONSE_BYTES = 64 * 1024 * 1024

TARGET_REGISTRY = {
    "SMARCA2": {
        "name": "SMARCA2",
        "uniprot_accession": "P51531",
        "search_attribute": (
            "rcsb_polymer_entity_container_identifiers."
            "reference_sequence_identifiers.database_accession"
        ),
    },
    "P51531": {
        "name": "SMARCA2",
        "uniprot_accession": "P51531",
        "search_attribute": (
            "rcsb_polymer_entity_container_identifiers."
            "reference_sequence_identifiers.database_accession"
        ),
    },
}

# Components which can pass the requested size filters but are normally
# crystallization additives rather than useful target-ligand evidence.
_ADDITIVE_CCD = {
    "12P", "15P", "BOG", "BME", "CAP", "CIT", "DMS", "EDO", "EPE",
    "GOL", "HEP", "MES", "MPD", "P6G", "PEG", "PG4", "SO4", "TRS",
}
_ADDITIVE_WORDS = {
    "buffer", "detergent", "glycerol", "polyethylene glycol", "solvent",
    "tris(hydroxymethyl)", "crystallization additive",
}
_PROTAC_WORDS = {
    "protac", "proteolysis targeting chimera", "protein degrader",
    "heterobifunctional degrader", "bifunctional degrader",
}

# Natural metabolites and cofactors are deposited target-associated molecules,
# but are not suitable warhead evidence without separate, curated evidence.
_NATURAL_COFACTOR_CCD = {
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "CTP", "CDP", "CMP",
    "UTP", "UDP", "UMP", "NAD", "NAP", "NAI", "FAD", "FMN", "SAM",
    "SAH", "COA", "ACO", "PLP", "PMP", "TPP", "THM", "BTN", "B12",
    "MTA", "HEM", "HEC", "UQ1",
}

# Skip CCD downloads which predictably do not yield a useful organic ligand.
# Unknown CCD failures remain errors and are not silently converted to exclusions.
_PREFETCH_SKIP_CCD = _ADDITIVE_CCD | {
    "HOH", "WAT", "DOD", "NA", "K", "CL", "BR", "IOD", "F", "LI",
    "MG", "CA", "MN", "FE", "FE2", "FE3", "CO", "NI", "CU", "CU1",
    "CU2", "ZN", "CD", "HG", "SR", "CS", "BA", "AL", "NH4", "NO3",
    "PO4",
}

_BROMODOMAIN_REFERENCE = {
    "label": "bromodomain",
    "basis": "curated 6HAZ construct overlap",
    "reference_pdb": "6HAZ",
    "reference_accession": "P51531-2",
    "uniprot_sequence_begin": 1373,
    "uniprot_sequence_end": 1493,
    "minimum_overlap_residues": 60,
}


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
            ensure_ascii=False) + "\n").encode("utf-8")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_name(url: str) -> str:
    name = Path(urlparse(url).path).name or "response"
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


def _safe_local_path(root: Path, relative_path: str | Path, *, must_exist: bool) -> Path:
    """Resolve a catalog path without permitting absolute paths, traversal, or symlinks."""
    root = root.resolve()
    relative = Path(relative_path)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise ValueError(f"Unsafe relative path: {relative_path}")
    candidate = root / relative
    current = root
    for part in relative.parts:
        current = current / part
        if current.exists() and current.is_symlink():
            raise ValueError(f"Symlinks are not permitted in evidence paths: {current}")
    if must_exist and not candidate.is_file():
        raise ValueError(f"Evidence file is missing: {candidate}")
    resolved = candidate.resolve(strict=must_exist)
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"Evidence path escapes its root: {relative_path}") from exc
    return candidate


class SnapshotStore:
    """Content-addressed response store with a mutable URL lookup index.

    Source bytes are never overwritten. The index only points URLs to immutable
    hashes and is not itself treated as source evidence.
    """

    def __init__(self, root: str | Path, timeout: float = 30.0,
                 max_bytes: int = MAX_RESPONSE_BYTES):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.timeout = float(timeout)
        self.max_bytes = int(max_bytes)
        self.index_path = self.root / "url_index.json"
        if self.index_path.exists():
            self.index = json.loads(self.index_path.read_text(encoding="utf-8"))
        else:
            self.index = {}
        self._client = httpx.Client(timeout=self.timeout, follow_redirects=False)

    def close(self) -> None:
        self._client.close()

    def _save_index(self) -> None:
        temporary = self.index_path.with_suffix(".tmp")
        temporary.write_bytes(_json_bytes(self.index))
        os.replace(temporary, self.index_path)

    def fetch(self, url: str, *, params: dict[str, str] | None = None,
              reuse: bool = True) -> tuple[Path, dict[str, Any]]:
        parsed = urlparse(url)
        if parsed.scheme != "https" or parsed.hostname not in ALLOWED_HOSTS:
            raise ValueError(f"URL is outside the RCSB HTTPS allowlist: {url}")
        key = url
        if params:
            key += "?" + "&".join(f"{k}={params[k]}" for k in sorted(params))
        cached = self.index.get(key)
        if reuse and cached:
            if not isinstance(cached, dict) or not isinstance(cached.get("sha256"), str):
                raise ValueError("Cached source metadata is invalid")
            path = _safe_local_path(self.root, cached.get("relative_path", ""), must_exist=True)
            if _file_sha256(path) != cached["sha256"]:
                raise ValueError(f"Cached source hash verification failed: {path}")
            return path, {**cached, "retrieval": "immutable_cache_reuse", "url": url}

        with self._client.stream("GET", url, params=params) as response:
            response.raise_for_status()
            length = response.headers.get("content-length")
            if length and int(length) > self.max_bytes:
                raise ValueError(f"RCSB response exceeds size bound: {url}")
            chunks: list[bytes] = []
            total = 0
            for chunk in response.iter_bytes():
                total += len(chunk)
                if total > self.max_bytes:
                    raise ValueError(f"RCSB response exceeds size bound: {url}")
                chunks.append(chunk)
        data = b"".join(chunks)
        digest = _sha256(data)
        relative = Path(digest) / _safe_name(url)
        destination = _safe_local_path(self.root, relative, must_exist=False)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination = _safe_local_path(self.root, relative, must_exist=False)
        if destination.exists():
            if _file_sha256(destination) != digest:
                raise ValueError(f"Immutable snapshot collision: {destination}")
        else:
            with destination.open("xb") as stream:
                stream.write(data)
        metadata = {
            "sha256": digest,
            "bytes": len(data),
            "relative_path": relative.as_posix(),
            "url": url,
        }
        self.index[key] = metadata
        self._save_index()
        return destination, {**metadata, "retrieval": "fresh_fetch"}


def _search_query(accession: str, max_entries: int) -> dict[str, Any]:
    return {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": TARGET_REGISTRY[accession]["search_attribute"],
                "operator": "exact_match",
                "value": accession,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": max_entries},
            "results_content_type": ["experimental"],
        },
    }


def _clean(value: Any) -> str:
    if value in (None, False, "?", "."):
        return ""
    return str(value).strip()


def _references(block: gemmi.cif.Block) -> list[dict[str, str]]:
    result = []
    for row in rows(block, "_struct_ref."):
        result.append({key: _clean(value) for key, value in sorted(row.items())})
    return sorted(result, key=lambda item: json.dumps(item, sort_keys=True))


def _accession_matches(value: Any, accession: str) -> bool:
    cleaned = _clean(value).upper()
    accession = accession.upper()
    return cleaned == accession or cleaned.startswith(accession + "-")


def _target_chains(block: gemmi.cif.Block, accession: str) -> tuple[list[str], list[dict[str, str]]]:
    references = _references(block)
    entity_ids = {
        row.get("entity_id", "")
        for row in references
        if _accession_matches(row.get("pdbx_db_accession", ""), accession)
        or _accession_matches(row.get("db_code", ""), accession)
    }
    polymer_entities = {
        _clean(row.get("entity_id"))
        for row in rows(block, "_entity_poly.")
        if _clean(row.get("entity_id")) in entity_ids
        and _clean(row.get("type")) == "polypeptide(L)"
    }
    chains = sorted(
        _clean(row.get("id"))
        for row in rows(block, "_struct_asym.")
        if _clean(row.get("entity_id")) in polymer_entities
    )
    return [chain for chain in chains if chain], references


def _citation_dois(block: gemmi.cif.Block) -> list[str]:
    values = {
        _clean(row.get("pdbx_database_id_DOI"))
        for row in rows(block, "_citation.")
    }
    return sorted(value for value in values if value)


def _target_construct_mapping(block: gemmi.cif.Block, accession: str) -> list[dict[str, Any]]:
    """Return deposited PDB-to-UniProt construct mappings without guessing a domain."""
    reference_ids = {
        _clean(row.get("id"))
        for row in rows(block, "_struct_ref.")
        if _accession_matches(row.get("pdbx_db_accession"), accession)
        or _accession_matches(row.get("db_code"), accession)
    }
    result = []
    for row in rows(block, "_struct_ref_seq."):
        if _clean(row.get("ref_id")) not in reference_ids:
            continue
        mapping: dict[str, Any] = {
            "alignment_id": _clean(row.get("align_id")),
            "pdb_chain_ids": sorted(filter(None, (
                value.strip() for value in _clean(row.get("pdbx_strand_id")).split(",")
            ))),
            "domain_annotation": "unidentified",
            "domain_annotation_requires_review": True,
        }
        for source_key, destination_key in (
            ("seq_align_beg", "pdb_sequence_begin"),
            ("seq_align_end", "pdb_sequence_end"),
            ("db_align_beg", "uniprot_sequence_begin"),
            ("db_align_end", "uniprot_sequence_end"),
        ):
            value = _clean(row.get(source_key))
            try:
                mapping[destination_key] = int(value)
            except ValueError:
                mapping[destination_key] = None
        result.append(mapping)
    return sorted(result, key=lambda item: json.dumps(item, sort_keys=True))


def _label_to_auth_chains(block: gemmi.cif.Block) -> dict[str, list[str]]:
    """Map atom-site label chain IDs to deposited author chain IDs."""
    mapping: dict[str, set[str]] = {}
    for row in rows(block, "_atom_site."):
        label = _clean(row.get("label_asym_id"))
        auth = _clean(row.get("auth_asym_id"))
        if label and auth:
            mapping.setdefault(label, set()).add(auth)
    return {label: sorted(values) for label, values in sorted(mapping.items())}


def classify_construct_domain(
    construct_mapping: list[dict[str, Any]],
    contacting_label_chains: list[str],
    label_to_auth_chains: dict[str, list[str]],
    domain: str = "bromodomain",
) -> dict[str, Any]:
    """Qualify contacted chains using deposited mappings, not inferred domains."""
    if domain not in {"bromodomain", "all_domains"}:
        raise ValueError("domain must be 'bromodomain' or 'all_domains'")
    contacted_auth = {
        auth
        for label in contacting_label_chains
        for auth in label_to_auth_chains.get(label, [])
    }
    mapped = []
    for item in construct_mapping:
        if not contacted_auth.intersection(item.get("pdb_chain_ids", [])):
            continue
        begin = item.get("uniprot_sequence_begin")
        end = item.get("uniprot_sequence_end")
        if not isinstance(begin, int) or not isinstance(end, int) or begin > end:
            continue
        mapped.append(item)
    if not mapped:
        return {
            "qualified_for_cards": False,
            "domain_annotation": "unidentified",
            "domain_partition": "unidentified",
            "domain_annotation_requires_review": True,
            "reason": "contacted_chain_has_no_unambiguous_struct_ref_seq_range",
            "reference_basis": _BROMODOMAIN_REFERENCE,
            "contacting_auth_asym_ids": sorted(contacted_auth),
            "matched_construct_mappings": [],
        }

    if domain == "all_domains":
        ranges = sorted({
            (item["uniprot_sequence_begin"], item["uniprot_sequence_end"])
            for item in mapped
        })
        partition = "+".join(f"P51531:{begin}-{end}" for begin, end in ranges)
        return {
            "qualified_for_cards": True,
            "domain_annotation": "construct_range_partition",
            "domain_partition": partition,
            "domain_annotation_requires_review": False,
            "reason": None,
            "reference_basis": "deposited _struct_ref_seq range; no cross-domain comparison",
            "contacting_auth_asym_ids": sorted(contacted_auth),
            "matched_construct_mappings": mapped,
        }

    reference_begin = _BROMODOMAIN_REFERENCE["uniprot_sequence_begin"]
    reference_end = _BROMODOMAIN_REFERENCE["uniprot_sequence_end"]
    overlaps = [
        max(0, min(item["uniprot_sequence_end"], reference_end)
            - max(item["uniprot_sequence_begin"], reference_begin) + 1)
        for item in mapped
    ]
    maximum_overlap = max(overlaps)
    qualified = maximum_overlap >= _BROMODOMAIN_REFERENCE["minimum_overlap_residues"]
    return {
        "qualified_for_cards": qualified,
        "domain_annotation": "bromodomain" if qualified else "different_construct_range",
        "domain_partition": "bromodomain" if qualified else "outside_bromodomain_reference",
        "domain_annotation_requires_review": not qualified,
        "reason": None if qualified else "contacted_construct_does_not_overlap_bromodomain_reference_by_60_residues",
        "reference_basis": _BROMODOMAIN_REFERENCE,
        "maximum_reference_overlap_residues": maximum_overlap,
        "contacting_auth_asym_ids": sorted(contacted_auth),
        "matched_construct_mappings": mapped,
    }


def _component_names(block: gemmi.cif.Block) -> dict[str, str]:
    return {
        _clean(row.get("id")): _clean(row.get("name"))
        for row in rows(block, "_chem_comp.")
    }


def _nonpolymer_instances(block: gemmi.cif.Block) -> list[dict[str, str]]:
    entity_to_comp = {
        _clean(row.get("entity_id")): _clean(row.get("comp_id"))
        for row in rows(block, "_pdbx_entity_nonpoly.")
    }
    names = _component_names(block)
    result = []
    for row in rows(block, "_struct_asym."):
        entity = _clean(row.get("entity_id"))
        comp = entity_to_comp.get(entity)
        if comp:
            result.append({
                "label_asym_id": _clean(row.get("id")),
                "entity_id": entity,
                "ccd": comp,
                "name": names.get(comp, ""),
            })
    return sorted(result, key=lambda item: (item["ccd"], item["label_asym_id"]))


def _crosslinked_nonpolymers(block: gemmi.cif.Block,
                             nonpolymer_chains: set[str]) -> set[str]:
    excluded: set[str] = set()
    for row in rows(block, "_struct_conn."):
        first = _clean(row.get("ptnr1_label_asym_id"))
        second = _clean(row.get("ptnr2_label_asym_id"))
        kind = _clean(row.get("conn_type_id")).lower()
        if first in nonpolymer_chains and second in nonpolymer_chains and first != second:
            if kind not in {"hydrog", "metalc"}:
                excluded.update((first, second))
    return excluded


def evaluate_candidate(candidate: dict[str, Any], accession: str = "P51531") -> str | None:
    """Return a deterministic exclusion reason, or ``None`` if eligible."""
    if accession not in set(candidate.get("target_accessions", [])):
        return "non_target_chain"
    heavy = int(candidate.get("heavy_atom_count", 0))
    if heavy < 10 or heavy > 50:
        return "heavy_atom_count_outside_10_50"
    mw = float(candidate.get("molecular_weight_Da", math.inf))
    if not math.isfinite(mw) or mw > 600:
        return "molecular_weight_over_600_Da"
    if not candidate.get("contains_carbon", False):
        return "not_an_organic_ligand"
    ccd = str(candidate.get("ccd", "")).upper()
    name = str(candidate.get("name", "")).lower()
    if ccd in _NATURAL_COFACTOR_CCD:
        return "natural_metabolite_or_cofactor"
    if ccd in _ADDITIVE_CCD or any(word in name for word in _ADDITIVE_WORDS):
        return "solvent_buffer_or_additive"
    if candidate.get("multiligand_crosslinked", False):
        return "multiligand_crosslinked"
    if candidate.get("complete_protac", False) or any(word in name for word in _PROTAC_WORDS):
        return "complete_PROTAC_not_a_warhead"
    if int(candidate.get("contacting_ligand_heavy_atoms", 0)) < 3:
        return "fewer_than_3_contacting_ligand_heavy_atoms"
    if int(candidate.get("contacting_protein_heavy_atoms", 0)) < 3:
        return "fewer_than_3_contacting_protein_heavy_atoms"
    return None


def _contact_summary(cif_path: Path, ligand: Chem.Mol, target_chains: list[str],
                     cutoff: float = 4.5) -> dict[str, Any]:
    ligand_xyz = np.asarray(ligand.GetConformer().GetPositions(), dtype=float)
    contacting_ligand: set[int] = set()
    contacting_protein: set[tuple[str, str]] = set()
    by_chain = []
    for chain in target_chains:
        protein = atom_sites(cif_path, [chain])
        protein_xyz = np.asarray([row["xyz"] for row in protein], dtype=float)
        distances = np.linalg.norm(ligand_xyz[:, None, :] - protein_xyz[None, :, :], axis=2)
        pairs = np.argwhere(distances <= cutoff)
        chain_ligand = {int(pair[0]) for pair in pairs}
        chain_protein = {int(pair[1]) for pair in pairs}
        contacting_ligand.update(chain_ligand)
        contacting_protein.update((chain, protein[index]["id"]) for index in chain_protein)
        if len(pairs) > 0:
            by_chain.append({
                "target_label_asym_id": chain,
                "contact_pair_count": int(len(pairs)),
                "contacting_ligand_heavy_atom_count": len(chain_ligand),
                "contacting_protein_heavy_atom_count": len(chain_protein),
            })
    return {
        "cutoff_A": cutoff,
        "contacting_ligand_heavy_atoms": len(contacting_ligand),
        "contacting_protein_heavy_atoms": len(contacting_protein),
        "contacting_target_label_asym_ids": [item["target_label_asym_id"] for item in by_chain],
        "by_exact_target_chain": by_chain,
        "interpretation": (
            "Heavy-atom proximity only; it is not a binding measurement, "
            "hydrogen-bond assignment, or SAR observation."
        ),
    }


def _connectivity_key(mol: Chem.Mol) -> str:
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    # Charge and stereochemistry are intentionally preserved.
    return Chem.MolToSmiles(copy, isomericSmiles=True, canonical=True)


def deduplicate_records(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """Group identical charged stereochemical parents while retaining all structures."""
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for record in records:
        key = record["parent_connectivity_key"]
        partition = record.get("domain_partition", "")
        groups.setdefault((key, partition), []).append(record)
    result = []
    for key, partition in sorted(groups):
        versions = sorted(groups[(key, partition)], key=lambda r: (r["pdb"], r["ligand_label_asym_id"], r["ccd"]))
        representative = min(
            versions,
            key=lambda r: (r["ccd"] != "FX5", r["pdb"] != "6HAZ",
                           r["pdb"], r["ligand_label_asym_id"]),
        )
        result.append({
            "group_id": "ligand-" + hashlib.sha256((key + "\0" + partition).encode()).hexdigest()[:16],
            "parent_connectivity_key": key,
            "domain_partition": partition or "unspecified",
            "canonical_isomeric_smiles": representative["identity"]["canonical_isomeric_smiles"],
            "representative_record_id": representative["record_id"],
            "representative_mol": representative.get("_mol"),
            "experimental_versions": [
                {k: v for k, v in version.items() if k != "_mol"}
                for version in versions
            ],
        })
    return result


def _fingerprint(group: dict[str, Any]):
    mol = group.get("representative_mol")
    if mol is None:
        mol = Chem.MolFromSmiles(group["canonical_isomeric_smiles"])
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def select_diverse(groups: list[dict[str, Any]], limit: int = 10) -> list[dict[str, Any]]:
    """Deterministic max-min Morgan selection, with the curated FX5 source first."""
    if not groups or limit <= 0:
        return []
    ordered = sorted(groups, key=lambda group: group["group_id"])
    fx5 = [group for group in ordered if any(
        version["ccd"] == "FX5" and version["pdb"] == "6HAZ"
        for version in group["experimental_versions"]
    )]
    first = fx5[0] if fx5 else ordered[0]
    selected = [first]
    remaining = [group for group in ordered if group is not first]
    fingerprints = {group["group_id"]: _fingerprint(group) for group in ordered}
    while remaining and len(selected) < limit:
        def diversity(group: dict[str, Any]) -> tuple[float, str]:
            similarities = [DataStructs.TanimotoSimilarity(
                fingerprints[group["group_id"]], fingerprints[item["group_id"]]
            ) for item in selected]
            return min(1.0 - value for value in similarities), group["group_id"]
        best = max(remaining, key=diversity)
        selected.append(best)
        remaining.remove(best)
    return selected


def _load_curated(path: str | Path | None, accession: str) -> list[dict[str, Any]]:
    if path is None:
        return []
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    entries = data.get("records", data) if isinstance(data, dict) else data
    if not isinstance(entries, list):
        raise ValueError("Curated evidence must be a list or contain a records list")
    result = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(f"Curated record {index} must be an object")
        missing = [field for field in ("source", "locator", "target_condition") if not entry.get(field)]
        if missing:
            raise ValueError(f"Curated record {index} lacks provenance: {', '.join(missing)}")
        if entry.get("target_accession") != accession:
            continue
        result.append(dict(entry))
    return result


def _canonical_curated_structure(entry: dict[str, Any]) -> str | None:
    value = entry.get("canonical_isomeric_smiles")
    if value is None and isinstance(entry.get("identity"), dict):
        value = entry["identity"].get("canonical_isomeric_smiles")
    if not value:
        return None
    molecule = Chem.MolFromSmiles(str(value))
    if molecule is None:
        raise ValueError("Curated canonical_isomeric_smiles is invalid")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def _matching_curated(entries: list[dict[str, Any]], pdb: str, ccd: str,
                      canonical_isomeric_smiles: str) -> list[dict[str, Any]]:
    """Match only an exact PDB+CCD pair or an exact charged stereochemical structure."""
    molecule = Chem.MolFromSmiles(canonical_isomeric_smiles)
    if molecule is None:
        raise ValueError("Collected ligand identity is not a valid SMILES")
    structure_key = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    result = []
    for entry in entries:
        pair_match = bool(entry.get("pdb") and entry.get("ccd")) and (
            str(entry["pdb"]).upper() == pdb.upper()
            and str(entry["ccd"]).upper() == ccd.upper()
        )
        curated_key = _canonical_curated_structure(entry)
        if pair_match or (curated_key is not None and curated_key == structure_key):
            result.append(entry)
    return result


def collect_evidence(target: str, output_dir: str | Path, *, max_entries: int = 32,
                     max_final: int = 10, timeout: float = 30.0,
                     curated_path: str | Path | None = None,
                     reuse_cache: bool = True, domain: str = "bromodomain",
                     cache_root: str | Path | None = None) -> dict[str, Any]:
    """Run the target-input evidence funnel and write ``catalog.json`` plus SDFs."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    if (output / "catalog.json").exists():
        raise FileExistsError(f"Existing evidence catalog will not be overwritten: {output / 'catalog.json'}")
    registry = TARGET_REGISTRY.get(target.upper())
    if registry is None:
        catalog = {
            "schema_version": SCHEMA_VERSION,
            "status": "unsupported",
            "requested_target": target,
            "supported_targets": ["SMARCA2", "P51531"],
            "cards": [], "all_records": [], "exclusions": [], "errors": [],
            "sources": [],
        }
        (output / "catalog.json").write_bytes(_json_bytes(catalog))
        return catalog
    if not 1 <= max_entries <= 32:
        raise ValueError("max_entries must be between 1 and 32")
    if not 1 <= max_final <= 10:
        raise ValueError("max_final must be between 1 and 10")
    if domain not in {"bromodomain", "all_domains"}:
        raise ValueError("domain must be 'bromodomain' or 'all_domains'")

    accession = registry["uniprot_accession"]
    curated = _load_curated(curated_path, accession)
    source_root = output / "source_snapshots"
    store_root = Path(cache_root) if cache_root is not None else source_root
    store = SnapshotStore(store_root, timeout=timeout)
    sources: dict[str, dict[str, Any]] = {}

    def retain_source(path: Path, metadata: dict[str, Any]) -> tuple[Path, dict[str, Any]]:
        """Copy externally cached immutable bytes into the portable output catalog."""
        digest = metadata["sha256"]
        if _file_sha256(path) != digest:
            raise ValueError(f"Source hash verification failed before retention: {path}")
        relative = Path(digest) / path.name
        destination = _safe_local_path(source_root, relative, must_exist=False)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if _file_sha256(destination) != digest:
                raise ValueError(f"Immutable snapshot collision: {destination}")
        else:
            with destination.open("xb") as stream:
                stream.write(path.read_bytes())
        retained = {**metadata, "relative_path": relative.as_posix()}
        sources[digest] = retained
        return destination, retained
    artifacts: list[dict[str, Any]] = []
    if curated_path is not None:
        curated_bytes = Path(curated_path).read_bytes()
        curated_hash = _sha256(curated_bytes)
        curated_relative = Path("curated_snapshots") / curated_hash / Path(curated_path).name
        curated_copy = _safe_local_path(output, curated_relative, must_exist=False)
        curated_copy.parent.mkdir(parents=True, exist_ok=True)
        curated_copy = _safe_local_path(output, curated_relative, must_exist=False)
        with curated_copy.open("xb") as stream:
            stream.write(curated_bytes)
        artifacts.append({
            "kind": "curated_evidence_json",
            "relative_path": curated_relative.as_posix(),
            "sha256": curated_hash,
        })
    records: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    reference_identifiers: dict[str, list[dict[str, str]]] = {}
    try:
        query = _search_query(accession, max_entries)
        query_bytes = _json_bytes(query)
        query_path = output / "search_query.json"
        query_path.write_bytes(query_bytes)
        artifacts.append({
            "kind": "RCSB_search_query",
            "relative_path": "search_query.json",
            "sha256": _sha256(query_bytes),
        })
        pdb_ids: list[str] = []
        try:
            search_path, search_meta = store.fetch(
                SEARCH_URL,
                params={"json": json.dumps(query, separators=(",", ":"))},
                reuse=reuse_cache,
            )
            search_path, search_meta = retain_source(search_path, search_meta)
            search = json.loads(search_path.read_text(encoding="utf-8"))
            pdb_ids = sorted({str(item["identifier"]).upper()
                              for item in search.get("result_set", [])})[:max_entries]
        except Exception as exc:
            errors.append({"stage": "search", "error": f"{type(exc).__name__}: {exc}"})

        for pdb_id in pdb_ids:
            try:
                cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
                cif_path, cif_meta = store.fetch(cif_url, reuse=reuse_cache)
                cif_path, cif_meta = retain_source(cif_path, cif_meta)
                block = gemmi.cif.read_file(str(cif_path)).sole_block()
                target_chains, references = _target_chains(block, accession)
                construct_mapping = _target_construct_mapping(block, accession)
                label_to_auth = _label_to_auth_chains(block)
                reference_identifiers[pdb_id] = references
                if not target_chains:
                    exclusions.append({"pdb": pdb_id, "reason": "no_exact_target_reference_chain"})
                    continue
                instances = _nonpolymer_instances(block)
                nonpolymer_chains = {item["label_asym_id"] for item in instances}
                crosslinked = _crosslinked_nonpolymers(block, nonpolymer_chains)
                dois = _citation_dois(block)
                for instance in instances:
                    ccd = instance["ccd"]
                    asym = instance["label_asym_id"]
                    base = {"pdb": pdb_id, "ccd": ccd, "ligand_label_asym_id": asym}
                    if ccd in _NATURAL_COFACTOR_CCD:
                        exclusions.append({**base, "reason": "natural_metabolite_or_cofactor"})
                        continue
                    if ccd in _PREFETCH_SKIP_CCD:
                        exclusions.append({**base, "reason": "known_ion_water_or_additive"})
                        continue
                    try:
                        ccd_url = f"https://files.rcsb.org/ligands/download/{ccd}.cif"
                        ccd_path, ccd_meta = store.fetch(ccd_url, reuse=reuse_cache)
                        ccd_path, ccd_meta = retain_source(ccd_path, ccd_meta)
                        ccd_mol = read_ccd(ccd_path)
                        ligand, mapping = ligand_coordinates(cif_path, ccd_mol, asym, ccd)
                        chemical_identity = identity(ligand)
                        contacts = _contact_summary(cif_path, ligand, target_chains)
                        candidate = {
                            **base,
                            "name": instance["name"],
                            "target_accessions": [accession],
                            "heavy_atom_count": chemical_identity["heavy_atom_count"],
                            "molecular_weight_Da": chemical_identity["properties"]["molecular_weight_Da"],
                            "contains_carbon": any(atom.GetAtomicNum() == 6 for atom in ligand.GetAtoms()),
                            "multiligand_crosslinked": asym in crosslinked,
                            "contacting_ligand_heavy_atoms": contacts["contacting_ligand_heavy_atoms"],
                            "contacting_protein_heavy_atoms": contacts["contacting_protein_heavy_atoms"],
                        }
                        reason = evaluate_candidate(candidate, accession)
                        if reason:
                            exclusions.append({**base, "reason": reason})
                            continue
                        record_id = f"{pdb_id}:{asym}:{ccd}"
                        domain_qualification = classify_construct_domain(
                            construct_mapping,
                            contacts["contacting_target_label_asym_ids"],
                            label_to_auth,
                            domain,
                        )
                        if not domain_qualification["qualified_for_cards"]:
                            exclusions.append({
                                **base,
                                "reason": domain_qualification["reason"],
                                "scope": "qualified_cards_only; retained_in_all_records",
                            })
                        curated_matches = _matching_curated(
                            curated, pdb_id, ccd,
                            chemical_identity["canonical_isomeric_smiles"],
                        )
                        records.append({
                            "record_id": record_id,
                            **base,
                            "name": instance["name"],
                            "target": registry["name"],
                            "target_accession": accession,
                            "exact_target_label_asym_ids": target_chains,
                            "target_construct_mapping": construct_mapping,
                            "label_to_auth_chain_mapping": label_to_auth,
                            "domain_annotation": domain_qualification["domain_annotation"],
                            "domain_partition": domain_qualification["domain_partition"],
                            "domain_annotation_requires_review": domain_qualification["domain_annotation_requires_review"],
                            "domain_qualification": domain_qualification,
                            "qualified_for_cards": domain_qualification["qualified_for_cards"],
                            "identity": chemical_identity,
                            "parent_connectivity_key": _connectivity_key(ligand),
                            "atom_mapping": mapping,
                            "contacts": contacts,
                            "source_pdb_citation_dois": dois,
                            "source_hashes": {
                                "pdb_mmcif_sha256": cif_meta["sha256"],
                                "ccd_mmcif_sha256": ccd_meta["sha256"],
                            },
                            "evidence_tiers": {
                                "co_crystal": "co_crystal_bound_ligand",
                                "binding": "curated_source_only" if curated_matches else "unknown",
                                "PROTAC": "unknown",
                                "database": "RCSB_PDB",
                            },
                            "curated_evidence": curated_matches,
                            "ligand_relevance_review": {
                                "required": True,
                                "reason": "Proximity establishes a bound deposited ligand, not target-specific affinity or warhead suitability.",
                            },
                            "measured_Kd": None,
                            "measured_Ki": None,
                            "measured_IC50": None,
                            "SAR": "unknown",
                            "exit_vector": "unknown",
                            "docking": "unknown",
                            "claim": "co_crystal_bound_ligand",
                            "claim_limit": "No binding measurement is inferred from contacts.",
                            "_mol": ligand,
                        })
                    except Exception as exc:
                        errors.append({**base, "stage": "ligand", "error": f"{type(exc).__name__}: {exc}"})
            except Exception as exc:
                errors.append({"pdb": pdb_id, "stage": "pdb", "error": f"{type(exc).__name__}: {exc}"})
    finally:
        store.close()

    records.sort(key=lambda record: (record["pdb"], record["ligand_label_asym_id"], record["ccd"]))
    qualified_records = [record for record in records if record.get("qualified_for_cards")]
    groups = deduplicate_records(qualified_records)
    if domain == "all_domains":
        # Select within deposited construct partitions so fingerprints from broad
        # ATPase-region constructs are never compared with bromodomain ligands.
        selected = []
        partitions = sorted({group["domain_partition"] for group in groups})
        partition_groups = {
            partition: [group for group in groups if group["domain_partition"] == partition]
            for partition in partitions
        }
        while len(selected) < max_final:
            changed = False
            for partition in partitions:
                already = sum(group["domain_partition"] == partition for group in selected)
                ranked = select_diverse(partition_groups[partition], already + 1)
                if len(ranked) > already:
                    selected.append(ranked[-1])
                    changed = True
                    if len(selected) == max_final:
                        break
            if not changed:
                break
    else:
        selected = select_diverse(groups, max_final)
    ligand_dir = output / "bound_ligands"
    ligand_dir.mkdir(exist_ok=True)
    cards = []
    for rank, group in enumerate(selected, start=1):
        mol = group.pop("representative_mol", None)
        sdf_name = f"{rank:02d}_{group['group_id']}.sdf"
        sdf_path = ligand_dir / sdf_name
        write_sdf(mol, sdf_path)
        sdf_hash = _file_sha256(sdf_path)
        artifacts.append({
            "kind": "bound_ligand_sdf",
            "relative_path": sdf_path.relative_to(output).as_posix(),
            "sha256": sdf_hash,
        })
        cards.append({
            **group,
            "selection_rank": rank,
            "selection_method": "deterministic max-min Morgan radius-2 diversity; 6HAZ/FX5 first when eligible",
            "bound_sdf": {"relative_path": sdf_path.relative_to(output).as_posix(), "sha256": sdf_hash},
            "evidence_tiers": {
                "co_crystal": "co_crystal_bound_ligand",
                "binding": "see experimental_versions; otherwise unknown",
                "PROTAC": "unknown",
                "database": "RCSB_PDB",
            },
            "missing_evidence": {"SAR": "unknown", "exit_vector": "unknown", "docking": "unknown"},
        })

    clean_records = [{k: v for k, v in record.items() if k != "_mol"} for record in records]
    minimum_cards_met = len(cards) >= 5
    catalog = {
        "schema_version": SCHEMA_VERSION,
        "status": ("complete" if not errors else "partial") if minimum_cards_met else "insufficient_evidence",
        "requested_target": target,
        "target": registry,
        "collection_policy": {
            "official_RCSB_HTTPS_only": True,
            "timeout_seconds": timeout,
            "retries": 0,
            "maximum_response_bytes": MAX_RESPONSE_BYTES,
            "maximum_entries": max_entries,
            "maximum_final_cards": max_final,
            "selected_domain": domain,
            "domain_basis": _BROMODOMAIN_REFERENCE if domain == "bromodomain" else "deposited construct-range partitions",
            "minimum_cards_for_acceptance": 5,
            "contact_cutoff_A": 4.5,
            "minimum_contacting_ligand_heavy_atoms": 3,
            "minimum_contacting_protein_heavy_atoms": 3,
            "charge_policy": "formal charge preserved during identity and deduplication",
            "measurement_policy": "Kd, Ki, and IC50 are never inferred, normalized, or compared",
        },
        "acceptance_summary": {
            "minimum_required_cards": 5,
            "qualified_card_count": len(cards),
            "minimum_5_cards_met": minimum_cards_met,
            "natural_cofactors_are_not_substitutes": True,
        },
        "cards": cards,
        "all_records": clean_records,
        "exclusions": sorted(exclusions, key=lambda item: json.dumps(item, sort_keys=True)),
        "errors": sorted(errors, key=lambda item: json.dumps(item, sort_keys=True)),
        "retrieved_reference_identifiers": reference_identifiers,
        "sources_root": "source_snapshots",
        "sources": sorted(sources.values(), key=lambda item: (item["sha256"], item["url"])),
        "artifacts": artifacts,
    }
    (output / "catalog.json").write_bytes(_json_bytes(catalog))
    return catalog


def load_catalog(path: str | Path, verify_hashes: bool = True) -> dict[str, Any]:
    """Load a catalog and verify all immutable source and SDF hashes."""
    catalog_path = Path(path)
    if catalog_path.is_dir():
        catalog_path = catalog_path / "catalog.json"
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    if catalog.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported warhead evidence catalog schema")
    if not verify_hashes:
        return catalog
    root = catalog_path.parent
    source_root_relative = Path(catalog.get("sources_root", "source_snapshots"))
    source_root = _safe_local_path(root, source_root_relative, must_exist=False)
    for source in catalog.get("sources", []):
        snapshot = _safe_local_path(source_root, source.get("relative_path", ""), must_exist=True)
        if _file_sha256(snapshot) != source.get("sha256"):
            raise ValueError(f"Source hash verification failed: {snapshot}")
    for artifact in catalog.get("artifacts", []):
        artifact_path = _safe_local_path(root, artifact.get("relative_path", ""), must_exist=True)
        if _file_sha256(artifact_path) != artifact.get("sha256"):
            raise ValueError(f"Artifact hash verification failed: {artifact_path}")
    return catalog
