"""AutoDock Vina 1.2.7 adapter for the fixed SMARCA2 docking case."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate

_METHOD = "Vina1.2.7"
_LIMITATIONS = [
    "none_of_this_establishes_efficacy",
    "no_pH_selection",
    "no_explicit_water",
]
_VINA_RESULT = re.compile(
    r"REMARK\s+VINA\s+RESULT:\s*"
    r"([-+]?\d+(?:\.\d+)?)"
    r"(?:\s+([-+]?\d+(?:\.\d+)?))?"
    r"(?:\s+([-+]?\d+(?:\.\d+)?))?"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _text_hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _short_windows_path(path: Path) -> str:
    value = str(path.resolve())
    if os.name != "nt":
        return value
    import ctypes

    size = ctypes.windll.kernel32.GetShortPathNameW(value, None, 0)
    if not size:
        return value
    buffer = ctypes.create_unicode_buffer(size)
    if not ctypes.windll.kernel32.GetShortPathNameW(value, buffer, size):
        return value
    return buffer.value


def _json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    if hasattr(value, "item"):
        return _json_safe(value.item())
    if hasattr(value, "tolist"):
        return _json_safe(value.tolist())
    return str(value)


def _base_result(seed: int, hashes: dict[str, str]) -> dict[str, Any]:
    return {
        "status": "failed",
        "method": _METHOD,
        "real_docking": False,
        "parent_redocking": False,
        "seed": seed,
        "results": {},
        "sourcehashes": hashes,
        "files": {
            "ligand_pdbqt": "ligand.pdbqt",
            "receptor_pdbqt": "receptor.pdbqt",
            "poses_pdbqt": "poses.pdbqt",
            "poses_sdf": "poses.sdf",
            "stdout": "vina.stdout.txt",
            "stderr": "vina.stderr.txt",
        },
        "limitations": list(_LIMITATIONS),
    }


def _validate_inputs(mol: Chem.Mol, parent: Chem.Mol, protected_maps: Any) -> set[int]:
    if not isinstance(mol, Chem.Mol) or not isinstance(parent, Chem.Mol):
        raise TypeError("mol and parent must be RDKit molecules")
    if mol.GetNumConformers():
        raise ValueError("mol must not contain a 3D conformer")
    if parent.GetNumConformers() < 1 or not parent.GetConformer().Is3D():
        raise ValueError("parent must contain an experimental 3D conformer")
    Chem.SanitizeMol(mol)
    Chem.SanitizeMol(parent)
    if Chem.GetFormalCharge(mol) != 0:
        raise ValueError("mol must be neutral")

    clean = Chem.Mol(mol)
    for atom in clean.GetAtoms(): atom.SetAtomMapNum(0)
    if any(str(info.specified) == 'Unspecified' for info in Chem.FindPotentialStereo(clean)):
        raise ValueError('DOCKING_UNSPECIFIED_STEREOCHEMISTRY')

    maps = [atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not maps or any(number <= 0 for number in maps) or len(maps) != len(set(maps)):
        raise ValueError("all ligand heavy atoms must have unique positive atom maps")

    protected = {int(number) for number in protected_maps}
    mol_maps = set(maps)
    parent_maps = {
        atom.GetAtomMapNum()
        for atom in parent.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom.GetAtomMapNum() > 0
    }
    if not protected <= mol_maps or not protected <= parent_maps:
        raise ValueError("protected_maps must be present in both mol and parent")
    return protected


def _box(parent: Chem.Mol) -> tuple[list[float], list[float]]:
    conf = parent.GetConformer()
    points = [conf.GetAtomPosition(i) for i in range(parent.GetNumAtoms())]
    if not points:
        raise ValueError("parent has no atoms")
    center = [
        sum(getattr(point, axis) for point in points) / len(points)
        for axis in ("x", "y", "z")
    ]
    size = []
    for axis in ("x", "y", "z"):
        values = [getattr(point, axis) for point in points]
        size.append(max(20.0, min(32.0, max(values) - min(values) + 12.0)))
    return center, size


def _canonical_heavy(molecule: Chem.Mol) -> str:
    heavy = Chem.RemoveHs(Chem.Mol(molecule), sanitize=True)
    for atom in heavy.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(heavy, canonical=True, isomericSmiles=True)


def _canonical_mapped_isomeric(molecule: Chem.Mol) -> str:
    copy = Chem.Mol(molecule)
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)


def _atom_identity(molecule: Chem.Mol) -> list[tuple[Any, ...]]:
    return [
        (
            atom.GetAtomicNum(),
            atom.GetIsotope(),
            atom.GetFormalCharge(),
            atom.GetIsAromatic(),
            int(atom.GetChiralTag()),
            atom.GetAtomMapNum(),
        )
        for atom in molecule.GetAtoms()
    ]


def _coordinates_match(actual: Chem.Mol, expected: Chem.Mol, tolerance: float = 0.002) -> bool:
    import numpy as np

    if actual.GetNumConformers() != 1 or expected.GetNumConformers() != 1:
        return False
    actual_xyz = np.asarray(actual.GetConformer().GetPositions(), dtype=float)
    expected_xyz = np.asarray(expected.GetConformer().GetPositions(), dtype=float)
    if actual_xyz.shape != expected_xyz.shape:
        return False
    if actual_xyz.shape != (actual.GetNumAtoms(), 3):
        return False
    if not np.isfinite(actual_xyz).all() or not np.isfinite(expected_xyz).all():
        return False
    return bool(np.allclose(actual_xyz, expected_xyz, atol=tolerance, rtol=0.0))


def _validated_parent_context(
    parent: Chem.Mol,
    protein: Any,
    selected_parent_id: str,
    module_root: Path,
) -> tuple[Chem.Mol, list[dict[str, Any]], dict[str, Any], dict[str, str]]:
    import numpy as np

    from packages.science.structures import atom_sites

    source = module_root / "cases" / "design_sources"
    fixed_cif = source / "6HAZ.cif"
    fixed_protein = atom_sites(fixed_cif, ["A"])
    protein_xyz = np.asarray([row["xyz"] for row in protein], dtype=float)
    fixed_xyz = np.asarray([row["xyz"] for row in fixed_protein], dtype=float)
    fields = (
        "label_asym_id",
        "auth_seq_id",
        "label_comp_id",
        "label_atom_id",
        "type_symbol",
    )
    if (
        len(protein) != len(fixed_protein)
        or protein_xyz.shape != fixed_xyz.shape
        or protein_xyz.shape != (len(fixed_protein), 3)
        or not np.isfinite(protein_xyz).all()
        or not np.isfinite(fixed_xyz).all()
        or not np.allclose(protein_xyz, fixed_xyz, atol=0.002, rtol=0.0)
    ):
        raise ValueError("Fixed SMARCA2 receptor/parent coordinate context required")
    if any(
        tuple(actual[key] for key in fields) != tuple(expected[key] for key in fields)
        for actual, expected in zip(protein, fixed_protein)
    ):
        raise ValueError("DOCKING_PROTEIN_IDENTITY_MISMATCH")

    context_hashes = {
        "6HAZ.cif": _sha256(fixed_cif),
    }
    source_file_hashes: dict[str, str] = {
        "fixed_6HAZ_mmcif_sha256": context_hashes["6HAZ.cif"],
    }
    if selected_parent_id == "SMARCA2-FX5":
        parent_source = source / "SMARCA2-neutral-design.sdf"
        fixed_parent = Chem.MolFromMolBlock(parent_source.read_text())
        if fixed_parent is None:
            raise ValueError("Fixed SMARCA2 parent is unreadable")
        context_hashes["SMARCA2-neutral-design.sdf"] = _sha256(parent_source)
        source_file_hashes["parent_sdf_sha256"] = context_hashes[
            "SMARCA2-neutral-design.sdf"
        ]
        matches_parent = (
            _canonical_heavy(parent) == _canonical_heavy(fixed_parent)
            and [atom.GetAtomMapNum() for atom in parent.GetAtoms()]
            == [atom.GetAtomMapNum() for atom in fixed_parent.GetAtoms()]
            and _coordinates_match(parent, fixed_parent)
        )
        source_kind = "original_neutral_default"
    else:
        from packages.science.dual_e3 import REFERENCE_SOURCE
        from packages.science.reference_parents import (
            load_reference_parent,
            reference_source_binding,
        )

        fixed_parent, metadata = load_reference_parent(
            selected_parent_id, REFERENCE_SOURCE
        )
        binding = reference_source_binding(REFERENCE_SOURCE)
        if binding.get("status") != "configured_hash_verified":
            raise ValueError("Reference-parent registry is not hash verified")
        required_binding_names = (
            "reference_parents/manifest.json",
            "reference_parents/index.json",
            f"reference_parents/parents/{selected_parent_id}.sdf",
            f"reference_parents/parents/{selected_parent_id}.metadata.json",
        )
        for name in required_binding_names:
            record = binding.get("files", {}).get(name)
            if not isinstance(record, dict) or record.get("status") != "configured_hash_verified":
                raise ValueError("Reference-parent source files are not fully registered")
            source_file_hashes[name] = record["sha256"]
            context_hashes[name] = record["sha256"]
        for name, value in metadata.get("source_input_hashes", {}).items():
            if isinstance(value, str):
                source_file_hashes[str(name)] = value
        matches_parent = (
            _canonical_mapped_isomeric(parent)
            == _canonical_mapped_isomeric(fixed_parent)
            and _atom_identity(parent) == _atom_identity(fixed_parent)
            and _coordinates_match(parent, fixed_parent)
        )
        source_kind = "hash_verified_aligned_literature_parent"

    if not matches_parent:
        raise ValueError("Fixed SMARCA2 receptor/parent coordinate context required")

    receipt = {
        "original_source_parent_id": selected_parent_id,
        "receptor_frame": "6HAZ chain A",
        "source_kind": source_kind,
        "source_file_hashes": source_file_hashes,
        "technical_role": "exploratory_pose_comparison_not_scientific_approval",
    }
    return fixed_parent, fixed_protein, receipt, context_hashes


def _restore_maps(
    pose: Chem.Mol,
    source: Chem.Mol,
    parent: Chem.Mol,
    protected: set[int],
) -> tuple[Chem.Mol, int, float | None]:
    target = Chem.RemoveHs(Chem.Mol(pose), sanitize=True)
    query = Chem.RemoveHs(Chem.Mol(source), sanitize=True)
    if _canonical_heavy(target) != _canonical_heavy(query):
        raise ValueError("restored pose heavy-atom graph differs from input ligand")

    if Chem.MolToSmiles(target, isomericSmiles=True) != Chem.MolToSmiles(query, isomericSmiles=True):
        raise ValueError('DOCKING_ATOM_PROVENANCE_MISMATCH')
    if not protected:
        return target, 1, None
    pidx={a.GetAtomMapNum():a.GetIdx() for a in parent.GetAtoms()}
    tidx={a.GetAtomMapNum():a.GetIdx() for a in target.GetAtoms()}
    squares=[target.GetConformer().GetAtomPosition(tidx[m]).Distance(parent.GetConformer().GetAtomPosition(pidx[m]))**2 for m in protected]
    return target, 1, math.sqrt(sum(squares)/len(squares))


def _kill(process: subprocess.Popen[Any]) -> None:
    if process.poll() is None:
        process.kill()
    process.wait()


def _failed(
    result: dict[str, Any], status: str, error: Exception | str
) -> dict[str, Any]:
    result["status"] = status
    result["results"] = {"error": str(error)}
    return _json_safe(result)


def dock(
    mol,
    parent,
    protein,
    protected_maps,
    workdir,
    check_active=lambda: None,
    exhaustiveness=8,
    seed=23,
    selected_parent_id="SMARCA2-FX5",
):
    module_root = Path(__file__).resolve().parents[2]
    vina = module_root / "vendor" / "vina" / "vina_1.2.7_win.exe"
    receptor_source = (
        module_root
        / "cases"
        / "design_sources"
        / "SMARCA2-receptor.pdbqt"
    )
    hashes: dict[str, str] = {}
    result = _base_result(int(seed), hashes)

    try:
        if not vina.is_file() or not receptor_source.is_file():
            raise FileNotFoundError("pinned Vina binary or fixed SMARCA2 receptor is missing")
        hashes["vina_1.2.7_win.exe"] = _sha256(vina)
        hashes["SMARCA2-receptor.pdbqt"] = _sha256(receptor_source)
        protected = _validate_inputs(mol, parent, protected_maps)
        fixed_parent, fixed_protein, parent_receipt, context_hashes = (
            _validated_parent_context(
                parent,
                protein,
                selected_parent_id,
                module_root,
            )
        )
        hashes.update(context_hashes)
        parent_receipt["source_file_hashes"]["fixed_receptor_pdbqt_sha256"] = hashes[
            "SMARCA2-receptor.pdbqt"
        ]
        result["source_parent_receipt"] = parent_receipt
        parent = fixed_parent
        protein = fixed_protein
        hashes["ligand_smiles"] = _text_hash(
            Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        )
        hashes["parent_smiles"] = _text_hash(
            Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        )

        run_dir = Path(workdir).resolve()
        run_dir.mkdir(parents=True, exist_ok=True)
        receptor = run_dir / "receptor.pdbqt"
        ligand_path = run_dir / "ligand.pdbqt"
        poses_path = run_dir / "poses.pdbqt"
        sdf_path = run_dir / "poses.sdf"
        shutil.copyfile(receptor_source, receptor)

        ligand = Chem.AddHs(Chem.Mol(mol))
        params = AllChem.ETKDGv3()
        params.randomSeed = int(seed)
        params.clearConfs = True
        if AllChem.EmbedMolecule(ligand, params) != 0:
            raise RuntimeError("ETKDGv3 ligand embedding failed")
        if not AllChem.MMFFHasAllMoleculeParams(ligand):
            raise RuntimeError("MMFF parameters are unavailable for the ligand")
        AllChem.MMFFOptimizeMolecule(ligand, maxIters=100)

        setups = MoleculePreparation().prepare(ligand)
        if not setups:
            raise RuntimeError("Meeko produced no ligand setup")
        written = PDBQTWriterLegacy.write_string(setups[0])
        pdbqt, success, message = written
        if not success:
            raise RuntimeError(f"Meeko PDBQT writing failed: {message}")
        ligand_path.write_text(pdbqt, encoding="ascii")

        center, size = _box(parent)
        command = [
            _short_windows_path(vina),
            "--receptor", "receptor.pdbqt",
            "--ligand", "ligand.pdbqt",
            "--center_x", f"{center[0]:.6f}",
            "--center_y", f"{center[1]:.6f}",
            "--center_z", f"{center[2]:.6f}",
            "--size_x", f"{size[0]:.6f}",
            "--size_y", f"{size[1]:.6f}",
            "--size_z", f"{size[2]:.6f}",
            "--exhaustiveness", str(int(exhaustiveness)),
            "--seed", str(int(seed)),
            "--cpu", "2",
            "--num_modes", "5",
            "--out", "poses.pdbqt",
        ]
        creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        check_active()
        with (run_dir / "vina.stdout.txt").open("wb") as stdout, (
            run_dir / "vina.stderr.txt"
        ).open("wb") as stderr:
            process = subprocess.Popen(
                command,
                cwd=str(run_dir),
                stdout=stdout,
                stderr=stderr,
                shell=False,
                creationflags=creationflags,
            )
            result['real_docking'] = True
            started = time.monotonic()
            try:
                while process.poll() is None:
                    check_active()
                    if time.monotonic() - started > 180.0:
                        _kill(process)
                        raise TimeoutError("Vina exceeded the 180 second timeout")
                    time.sleep(0.2)
            except BaseException:
                _kill(process)
                raise
            if process.returncode != 0:
                raise RuntimeError(f"Vina failed with exit code {process.returncode}")

        if not poses_path.is_file() or poses_path.stat().st_size == 0:
            raise RuntimeError("Vina produced no pose output")
        output_text = poses_path.read_text(encoding="utf-8", errors="replace")
        scores = [
            {
                "affinity": float(match.group(1)),
                "rmsd_lb": float(match.group(2)) if match.group(2) else None,
                "rmsd_ub": float(match.group(3)) if match.group(3) else None,
            }
            for match in _VINA_RESULT.finditer(output_text)
        ]
        restored = []
        mode_blocks = re.findall(r'MODEL\s+\d+.*?ENDMDL', output_text, flags=re.S)
        if not mode_blocks or len(mode_blocks) != len(scores):
            raise ValueError('DOCKING_MODEL_SCORE_PAIRING')
        for block in mode_blocks:
            converted = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(block+'\n', skip_typing=True))
            if len(converted) != 1 or converted[0] is None or converted[0].GetNumConformers() != 1:
                raise ValueError('DOCKING_MODEL_RESTORATION')
            restored.append(converted[0])
        if not restored:
            raise RuntimeError("Meeko could not restore any valid docked pose")

        pose_records, summaries = [], []
        sdf_stream = sdf_path.open('w', encoding='utf-8')
        writer = Chem.SDWriter(sdf_stream)
        try:
            for index, raw_pose in enumerate(restored):
                pose, mapping_count, core_rmsd = _restore_maps(
                    raw_pose, mol, parent, protected
                )
                score = scores[index] if len(scores) == len(restored) else None
                pose.SetIntProp("vina_mode", index + 1)
                pose.SetIntProp("atom_mapping_count", mapping_count)
                pose.SetProp("atom_mapping_selected", "retained_Meeko_SMILES_atom_map_provenance")
                if core_rmsd is not None:
                    pose.SetDoubleProp("parent_protected_core_RMSD", core_rmsd)
                if score is not None:
                    pose.SetDoubleProp("vina_affinity", score["affinity"])
                writer.write(pose)
                pose_records.append({
                    "mol": pose,
                    "score": None if score is None else score["affinity"],
                })
                summaries.append({
                    "mode": index + 1,
                    "score": score,
                    "mapping_count": mapping_count,
                    "mapping_ambiguous": mapping_count > 1,
                    "selected_mapping": "retained_Meeko_SMILES_atom_map_provenance",
                    "parent_protected_core_rmsd": core_rmsd,
                })
        finally:
            writer.close()
            sdf_stream.close()

        if protected:
            from packages.science.warhead_sites import pose_preservation

            preservation = pose_preservation(
                parent, pose_records, protein, protected_maps
            )
        else:
            preservation = {
                "status": "review_protected_core_not_defined",
                "docking_pose_preserved": "review",
                "protected_maps": [],
                "all_poses_diagnostics": [
                    {
                        "mode": index + 1,
                        "status": "review",
                        "reason": "PROTECTED_CORE_NOT_DEFINED",
                        "core_RMSD_A_in_receptor_frame": None,
                        "raw_pose_available_for_comparison": True,
                    }
                    for index in range(len(pose_records))
                ],
                "interpretation": (
                    "Vina poses were produced for exploratory comparison, but no "
                    "scientific protected map was defined and no preservation pass "
                    "is assigned."
                ),
            }
        result["status"] = "completed_with_limits"
        result["results"] = {
            "poses": summaries,
            "pose_count": len(summaries),
            "pose_preservation": _json_safe(preservation),
            "box_center": center,
            "box_size": size,
        }
        return _json_safe(result)
    except TimeoutError as exc:
        return _failed(result, "failed_timeout", exc)
    except Exception as exc:
        from packages.agents.provider import AgentError
        if isinstance(exc, AgentError): raise
        return _failed(result, "failed_docking", exc)
