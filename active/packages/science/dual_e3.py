from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem
from rdkit.Chem import Descriptors

from . import linker_design
from .acceptance_evidence import load_acceptance_evidence, source_binding
from .analog_generation import RULE_CATALOG
from .linker_assessment import expanded_library
from .molecules import canonical, identity, join_fragments

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "cases" / "design_sources"
REFERENCE_SOURCE = ROOT / "cases" / "reference_parents"
DEFAULT_PARENT_ID = "SMARCA2-FX5"

_DEFAULT_LINKER_IDS = (
    "alkyl_c6",
    "peg3",
    "peg_alkyl",
    "piperazine",
    "triazole",
    "aryl_para",
)
_LINKER_IDS = _DEFAULT_LINKER_IDS + (
    "peg_short", "peg_long", "alkyl_short", "alkyl_long",
    "peg_alkyl_short", "peg_alkyl_long", "piperazine_short", "piperazine_long",
    "triazole_short", "triazole_long", "aryl_short", "aryl_long",
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_sources():
    manifest_path = SOURCE / 'manifest.json'
    if manifest_path.is_symlink():
        raise ValueError('DESIGN_SOURCE_PATH')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    checked = {}
    for name, expected in manifest['files'].items():
        relative = Path(name)
        path = SOURCE / relative
        if relative.is_absolute() or '..' in relative.parts or path.is_symlink():
            raise ValueError('DESIGN_SOURCE_PATH')
        resolved = path.resolve()
        if not resolved.is_relative_to(SOURCE.resolve()):
            raise ValueError('DESIGN_SOURCE_PATH')
        if _sha256_file(resolved) != expected:
            raise ValueError('DESIGN_SOURCE_HASH_MISMATCH')
        checked['design_sources/' + name] = {'status': 'configured_hash_verified', 'sha256': expected}
    checked['acceptance_sources'] = source_binding()
    reference_module = ROOT / "packages" / "science" / "reference_parents.py"
    checked['code/reference_parents.py'] = {
        'status': 'configured_hash_verified',
        'sha256': _sha256_file(reference_module),
    }
    from .reference_parents import reference_source_binding
    checked['reference_parents'] = reference_source_binding(REFERENCE_SOURCE)
    return checked


def catalog():
    verify_sources()
    from .reference_parents import reference_parent_catalog
    return {'recruiters': json.loads((SOURCE/'recruiters.json').read_text(encoding='utf-8')),
            'linkers': expanded_library(),
            'CRBNbenchmark': json.loads((SOURCE/'crbn_benchmark.json').read_text(encoding='utf-8')),
            'warhead': json.loads((SOURCE/'warhead.json').read_text(encoding='utf-8')),
            'reference_parents': reference_parent_catalog(REFERENCE_SOURCE),
            'acceptance_evidence': load_acceptance_evidence(),
            'rule_catalog': RULE_CATALOG}


def validate(parameters: dict[str, Any]) -> dict[str, Any]:
    if type(parameters) is not dict:
        raise TypeError("parameters must be a dict")

    allowed = {
        "exploratory",
        "dock",
        "use_api",
        "panel_size",
        "linker_ids",
        "target",
        "parent_id",
    }
    unknown = sorted(set(parameters) - allowed)
    if unknown:
        raise ValueError(f"Unknown parameters: {unknown}")

    result = {
        "exploratory": False,
        "dock": True,
        "use_api": False,
        "panel_size": 12,
        "linker_ids": list(_DEFAULT_LINKER_IDS),
        "target": "SMARCA2",
    }

    for key in ("exploratory", "dock", "use_api"):
        if key in parameters:
            if type(parameters[key]) is not bool:
                raise TypeError(f"{key} must be a bool")
            result[key] = parameters[key]

    if "panel_size" in parameters:
        value = parameters["panel_size"]
        if type(value) is not int:
            raise TypeError("panel_size must be an int, not bool or another numeric type")
        if not 10 <= value <= 20:
            raise ValueError("panel_size must be between 10 and 20")
        result["panel_size"] = value

    if "linker_ids" in parameters:
        value = parameters["linker_ids"]
        if type(value) is not list:
            raise TypeError("linker_ids must be a list")
        if not value:
            raise ValueError("linker_ids must not be empty")
        if len(value) > 18:
            raise ValueError("linker_ids may contain at most eighteen IDs")
        if any(type(item) is not str for item in value):
            raise TypeError("Every linker ID must be a string")
        if len(set(value)) != len(value):
            raise ValueError("linker_ids must be unique")
        unknown_ids = sorted(set(value) - set(_LINKER_IDS))
        if unknown_ids:
            raise ValueError(f"Unknown linker IDs: {unknown_ids}")
        result["linker_ids"] = list(value)

    if "target" in parameters:
        value = parameters["target"]
        if type(value) is not str:
            raise TypeError("target must be a string")
        if value != "SMARCA2":
            raise ValueError("Only target 'SMARCA2' is supported")
        result["target"] = value

    if "parent_id" in parameters:
        value = parameters["parent_id"]
        if type(value) is not str:
            raise TypeError("parent_id must be a string")
        from .reference_parents import reference_parent_catalog
        available = {
            row["id"] for row in reference_parent_catalog(REFERENCE_SOURCE)["available"]
            if isinstance(row, dict) and row.get("ready") is True
        }
        if value not in available:
            raise ValueError("Unknown or unavailable parent_id")
        result["parent_id"] = value

    return result


validate_parameters = validate


def _validated_artifact_ref(value):
    if type(value) is str:
        try:value=json.loads(value)
        except (TypeError,ValueError):raise ValueError('SCIENTIFIC_POLICY_SOURCE_REF')
    required={'artifact_id','version','sha256','media_type','schema_id','provenance'}
    if type(value) is not dict or set(value)!=required:
        raise ValueError('SCIENTIFIC_POLICY_SOURCE_REF')
    if type(value['artifact_id']) is not str or not value['artifact_id']:
        raise ValueError('SCIENTIFIC_POLICY_SOURCE_REF')
    if type(value['version']) is not int or type(value['version']) is bool or value['version']<1:
        raise ValueError('SCIENTIFIC_POLICY_SOURCE_REF')
    sha=value['sha256']
    if type(sha) is not str or len(sha)!=64 or any(character not in '0123456789abcdefABCDEF' for character in sha):
        raise ValueError('SCIENTIFIC_POLICY_SOURCE_REF')
    if type(value['media_type']) is not str or not value['media_type']:
        raise ValueError('SCIENTIFIC_POLICY_SOURCE_REF')
    if type(value['schema_id']) is not str or not value['schema_id']:
        raise ValueError('SCIENTIFIC_POLICY_SOURCE_REF')
    if value['provenance'] not in {'computed','source'}:
        raise ValueError('SCIENTIFIC_POLICY_SOURCE_REF')
    return json.loads(json.dumps(value))


def validate_scientific_policy(policy,project_id=None,policy_id=None,parent_id=None):
    if type(policy) is not dict or set(policy)!={'scope','site_policy','parent_funnel'}:
        raise ValueError('SCIENTIFIC_POLICY')
    policy=json.loads(json.dumps(policy))
    scope=policy.get('scope')
    if type(scope) is not dict:
        raise ValueError('SCIENTIFIC_POLICY_SCOPE')
    required={'project_id','job_id','parent_id','assessment_id','policy_revision','policy_digest'}
    if set(scope)!=required:
        raise ValueError('SCIENTIFIC_POLICY_SCOPE')
    if project_id is not None and scope['project_id']!=project_id:
        raise ValueError('SCIENTIFIC_POLICY_PROJECT_SCOPE_MISMATCH')
    if policy_id is not None and scope['assessment_id']!=policy_id:
        raise ValueError('SCIENTIFIC_POLICY_ASSESSMENT_SCOPE_MISMATCH')
    if parent_id is not None and scope['parent_id']!=parent_id:
        raise ValueError('SCIENTIFIC_POLICY_PARENT_SCOPE_MISMATCH')
    if type(scope['policy_revision']) is not int or type(scope['policy_revision']) is bool or scope['policy_revision']<1:
        raise ValueError('SCIENTIFIC_POLICY_REVISION')
    if type(scope['policy_digest']) is not str or not scope['policy_digest']:
        raise ValueError('SCIENTIFIC_POLICY_DIGEST')
    for key in ('site_policy','parent_funnel'):
        if type(policy.get(key)) is not list:
            raise ValueError('SCIENTIFIC_POLICY_'+key.upper())
    from packages.science.analog_generation import _RULE_BY_ID
    allowed_scopes=set(_RULE_BY_ID)
    allowed_scopes.update(rule['transformation_class'] for rule in _RULE_BY_ID.values())
    for declaration in policy['site_policy']:
        if type(declaration) is not dict:
            raise ValueError('SCIENTIFIC_POLICY_SITE_DECLARATION')
        required={'atom_map','state','allowed_transforms','rationale','source_ref','allow_release_protected'}
        if set(declaration)!=required:
            raise ValueError('SCIENTIFIC_POLICY_SITE_DECLARATION')
        transforms=declaration['allowed_transforms']
        if type(transforms) is not list or any(type(value) is not str or value not in allowed_scopes for value in transforms):
            raise ValueError('SCIENTIFIC_POLICY_ALLOWED_TRANSFORMS')
        declaration['source_ref']=_validated_artifact_ref(declaration['source_ref'])
    return policy


def apply_scientific_policy(sites,policy,source_protected_maps=()):
    policy=validate_scientific_policy(policy)
    if type(sites) is not list:
        raise TypeError('sites must be a list')
    output=json.loads(json.dumps(sites))
    by_map={row.get('atom_map'):row for row in output}
    if len(by_map)!=len(output) or any(type(value) is not int or type(value) is bool for value in by_map):
        raise ValueError('SCIENTIFIC_POLICY_SITE_MAPS')
    declarations={}
    for declaration in policy['site_policy']:
        atom_map=declaration['atom_map']
        if type(atom_map) is not int or type(atom_map) is bool or atom_map not in by_map or atom_map in declarations:
            raise ValueError('SCIENTIFIC_POLICY_SITE_MAP_SCOPE')
        if declaration['state'] not in {'PROTECTED','MODIFIABLE','UNKNOWN'}:
            raise ValueError('SCIENTIFIC_POLICY_SITE_STATE')
        if type(declaration['rationale']) is not str:
            raise ValueError('SCIENTIFIC_POLICY_PROVENANCE')
        if type(declaration['allow_release_protected']) is not bool:
            raise ValueError('SCIENTIFIC_POLICY_RELEASE_FLAG')
        declarations[atom_map]=declaration
    source_protected=set(source_protected_maps)
    if any(type(atom_map) is not int or type(atom_map) is bool or atom_map not in by_map for atom_map in source_protected):
        raise ValueError('SCIENTIFIC_POLICY_PROTECTED_MAP_SCOPE')
    for atom_map,row in by_map.items():
        original=row.get('state','UNKNOWN')
        row['original_state']=original
        declaration=declarations.get(atom_map)
        if declaration is None:
            if original=='UNKNOWN':row['scientific_pending']=True
            continue
        provenance={'assessment_id':policy['scope']['assessment_id'],
            'policy_revision':policy['scope']['policy_revision'],'policy_digest':policy['scope']['policy_digest'],
            'rationale':declaration['rationale'],'source_ref':declaration['source_ref']}
        row.setdefault('evidence',{})['scientific_policy']=provenance
        requested=declaration['state']
        release=(atom_map in source_protected or original=='PROTECTED') and requested!='PROTECTED'
        if release and not (declaration['allow_release_protected'] and declaration['rationale']):
            requested='PROTECTED'
        if original=='UNKNOWN' and requested=='MODIFIABLE' and not declaration['rationale']:
            requested='UNKNOWN'
        row['state']=requested
        row['protected']=requested=='PROTECTED'
        row['scientific_pending']=requested=='UNKNOWN'
        sar=row.setdefault('evidence',{}).setdefault('SAR',{})
        sar['allowed_transformations']=list(declaration['allowed_transforms']) if requested=='MODIFIABLE' else []
        sar['scientific_policy_source_ref']=declaration['source_ref']
    return output


def _mol_from_smiles(smiles: str, label: str) -> Chem.Mol:
    if type(smiles) is not str or not smiles:
        raise ValueError(f"{label} must provide a non-empty SMILES string")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid {label} SMILES")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def _mapped_atom(mol: Chem.Mol, atom_map: int) -> Chem.Atom:
    matches = [
        atom
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() > 0 and atom.GetAtomMapNum() == atom_map
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one mapped heavy atom for map {atom_map}, found {len(matches)}"
        )
    return matches[0]


def _mapped_heavy_maps(mol: Chem.Mol) -> list[int]:
    maps = [
        atom.GetAtomMapNum()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() > 0
    ]
    if any(atom_map <= 0 for atom_map in maps):
        raise ValueError("All heavy atoms must have positive atom maps")
    if len(maps) != len(set(maps)):
        raise ValueError("Heavy-atom maps must be unique")
    return sorted(maps)


def _protected_maps(analog: dict[str, Any]) -> set[int]:
    raw = analog.get("protected_atom_maps")
    if raw is None and isinstance(analog.get("metadata"), dict):
        raw = analog["metadata"].get("protected_atom_maps")
    if raw is None:
        raw = []
    if type(raw) is not list or any(type(value) is not int for value in raw):
        raise TypeError("protected_atom_maps metadata must be a list of integers")
    return set(raw)


def _is_amide_nitrogen(atom: Chem.Atom) -> bool:
    if atom.GetAtomicNum() != 7:
        return False
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() != 6:
            continue
        for bond in neighbor.GetBonds():
            other = bond.GetOtherAtom(neighbor)
            if (
                other.GetIdx() != atom.GetIdx()
                and other.GetAtomicNum() in (8, 16)
                and bond.GetBondType() == Chem.BondType.DOUBLE
            ):
                return True
    return False


def _handle_record(atom: Chem.Atom, source: str) -> dict[str, Any]:
    return {
        "atom_map": atom.GetAtomMapNum(),
        "element": atom.GetSymbol(),
        "source": source,
        "formal_charge": atom.GetFormalCharge(),
        "hydrogens_before_attachment": int(atom.GetTotalNumHs()),
        "heavy_atom_degree": sum(
            1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 0
        ),
        "selection_basis": (
            "new terminal neutral N-H/O-H handle"
            if source == "added_atom_maps"
            else "configured original map-19 neutral N-H fallback"
        ),
    }


def attachment_options(analog: dict[str, Any]) -> list[dict[str, Any]]:
    if type(analog) is not dict:
        raise TypeError("analog must be a dict")
    mol = _mol_from_smiles(analog.get("mapped_smiles"), "analog mapped")
    protected = _protected_maps(analog)

    added = analog.get("added_atom_maps", [])
    if type(added) is not list or any(type(value) is not int for value in added):
        raise TypeError("added_atom_maps must be a list of integers")

    options: list[dict[str, Any]] = []
    for atom_map in sorted(set(added)):
        if atom_map in protected:
            continue
        atom = _mapped_atom(mol, atom_map)
        if atom.GetAtomicNum() not in (7, 8):
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
            continue
        if atom.GetTotalNumHs() < 1:
            continue
        heavy_degree = sum(
            1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 0
        )
        if heavy_degree != 1:
            continue
        if atom.GetAtomicNum() == 7 and _is_amide_nitrogen(atom):
            continue
        options.append(_handle_record(atom, "added_atom_maps"))

    if not options and 19 not in protected:
        try:
            atom = _mapped_atom(mol, 19)
        except ValueError:
            atom = None
        if (
            atom is not None
            and atom.GetAtomicNum() == 7
            and atom.GetFormalCharge() == 0
            and not atom.GetIsAromatic()
            and atom.GetTotalNumHs() >= 1
            and not _is_amide_nitrogen(atom)
        ):
            options.append(_handle_record(atom, "map19_fallback"))

    return options[:2]


def _attach_dummy(
    mol: Chem.Mol,
    attachment_map: int,
    dummy_map: int,
) -> tuple[Chem.Mol, int, int]:
    editable = Chem.RWMol(mol)
    atom = _mapped_atom(editable, attachment_map)
    if atom.GetAtomicNum() not in (7, 8):
        raise ValueError("Attachment atom must be nitrogen or oxygen")
    if atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
        raise ValueError("Attachment atom must be neutral and non-aromatic")
    if atom.GetAtomicNum() == 7 and _is_amide_nitrogen(atom):
        raise ValueError("Amide nitrogen attachment is prohibited")

    hydrogens_before = int(atom.GetTotalNumHs())
    if hydrogens_before < 1:
        raise ValueError("Attachment requires an N-H or O-H hydrogen to consume")

    atom.SetNumExplicitHs(hydrogens_before - 1)
    atom.SetNoImplicit(True)

    dummy = Chem.Atom(0)
    dummy.SetAtomMapNum(dummy_map)
    dummy_idx = editable.AddAtom(dummy)
    editable.AddBond(atom.GetIdx(), dummy_idx, Chem.BondType.SINGLE)

    result = editable.GetMol()
    Chem.SanitizeMol(result)
    attached = _mapped_atom(result, attachment_map)
    hydrogens_after = int(attached.GetTotalNumHs())
    if hydrogens_after != hydrogens_before - 1:
        raise ValueError("Attachment did not consume exactly one hydrogen")
    return result, hydrogens_before, hydrogens_after


def _prepare_linker(template: dict[str, Any], orientation: str) -> tuple[Chem.Mol, list[int]]:
    if type(template) is not dict:
        raise TypeError("linker template must be a dict")
    linker_id = template.get("id")
    if linker_id not in _LINKER_IDS:
        raise ValueError(f"Unknown linker template ID: {linker_id!r}")

    smiles = template.get("smiles", template.get("mapped_smiles"))
    mol = _mol_from_smiles(smiles, "linker")
    dummy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    dummy_maps = sorted(atom.GetAtomMapNum() for atom in dummy_atoms)
    if dummy_maps != [1001, 1002]:
        raise ValueError("Linker must contain exactly terminal dummies 1001 and 1002")

    editable = Chem.RWMol(mol)
    if orientation == "reverse":
        for atom in editable.GetAtoms():
            if atom.GetAtomicNum() == 0:
                if atom.GetAtomMapNum() == 1001:
                    atom.SetAtomMapNum(1002)
                elif atom.GetAtomMapNum() == 1002:
                    atom.SetAtomMapNum(1001)

    next_map = 2001
    linker_maps: list[int] = []
    for atom in editable.GetAtoms():
        if atom.GetAtomicNum() == 0:
            continue
        if atom.GetAtomMapNum() != 0:
            raise ValueError("Linker heavy atoms must be unmapped in the source template")
        atom.SetAtomMapNum(next_map)
        linker_maps.append(next_map)
        next_map += 1

    result = editable.GetMol()
    Chem.SanitizeMol(result)
    return result, linker_maps


def _atom_signature(atom: Chem.Atom, include_hydrogens: bool = True) -> tuple[Any, ...]:
    signature: tuple[Any, ...] = (
        atom.GetAtomicNum(),
        atom.GetIsotope(),
        atom.GetFormalCharge(),
        bool(atom.GetIsAromatic()),
        int(atom.GetChiralTag()),
    )
    if include_hydrogens:
        signature += (int(atom.GetTotalNumHs()),)
    return signature


def _bond_signature(bond: Chem.Bond) -> tuple[Any, ...]:
    return (
        str(bond.GetBondType()),
        bool(bond.GetIsAromatic()),
        int(bond.GetStereo()),
    )


def _cip(atom: Chem.Atom) -> str | None:
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None


def _verify_induced_graph(
    source: Chem.Mol,
    product: Chem.Mol,
    maps: set[int],
    *,
    ignore_hydrogen_maps: set[int] | None = None,
) -> None:
    ignore_hydrogen_maps = ignore_hydrogen_maps or set()
    source_maps = set(_mapped_heavy_maps(source))
    product_maps = set(_mapped_heavy_maps(product))
    if not maps <= source_maps or not maps <= product_maps:
        missing = sorted(maps - source_maps) + sorted(maps - product_maps)
        raise ValueError(f"Protected maps are missing from graph comparison: {missing}")

    Chem.AssignStereochemistry(source, cleanIt=True, force=True)
    Chem.AssignStereochemistry(product, cleanIt=True, force=True)

    for atom_map in sorted(maps):
        old = _mapped_atom(source, atom_map)
        new = _mapped_atom(product, atom_map)
        include_h = atom_map not in ignore_hydrogen_maps
        if _atom_signature(old, include_h) != _atom_signature(new, include_h):
            raise ValueError(f"Protected atom properties changed at map {atom_map}")
        if old.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            if _cip(old) != _cip(new):
                raise ValueError(
                    f"Protected stereocenter CIP assignment changed at map {atom_map}; "
                    "rejected conservatively because substituent-priority effects are unresolved"
                )

    for left in sorted(maps):
        old_left = _mapped_atom(source, left)
        new_left = _mapped_atom(product, left)
        for right in sorted(value for value in maps if value > left):
            old_right = _mapped_atom(source, right)
            new_right = _mapped_atom(product, right)
            old_bond = source.GetBondBetweenAtoms(old_left.GetIdx(), old_right.GetIdx())
            new_bond = product.GetBondBetweenAtoms(new_left.GetIdx(), new_right.GetIdx())
            if (old_bond is None) != (new_bond is None):
                raise ValueError(f"Protected bond topology changed between maps {left} and {right}")
            if old_bond is not None and _bond_signature(old_bond) != _bond_signature(new_bond):
                raise ValueError(f"Protected bond properties changed between maps {left} and {right}")


def _clear_maps(mol: Chem.Mol) -> Chem.Mol:
    result = Chem.Mol(mol)
    for atom in result.GetAtoms():
        atom.SetAtomMapNum(0)
    return result


def _as_mol(value: Any) -> Chem.Mol:
    if isinstance(value, Chem.Mol):
        return Chem.Mol(value)
    if isinstance(value, str):
        return _mol_from_smiles(value, "joined candidate")
    raise TypeError("join_fragments returned neither an RDKit molecule nor a SMILES string")


def _json_safe(value: Any) -> Any:
    if value is None or type(value) in (str, int, float, bool):
        return value
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    return str(value)


def _identity_value(mol: Chem.Mol, canonical_smiles: str) -> Any:
    try:
        value = identity(mol)
    except (TypeError, AttributeError):
        value = identity(canonical_smiles)
    return _json_safe(value)


def _actual_attachment(
    product: Chem.Mol,
    attachment_map: int,
    partner_maps: set[int],
) -> dict[str, Any]:
    atom = _mapped_atom(product, attachment_map)
    bonds = []
    for bond in atom.GetBonds():
        neighbor = bond.GetOtherAtom(atom)
        if neighbor.GetAtomMapNum() in partner_maps:
            bonds.append(
                {
                    "attachment_atom_map": attachment_map,
                    "partner_atom_map": neighbor.GetAtomMapNum(),
                    "bond_type": str(bond.GetBondType()),
                }
            )
    if len(bonds) != 1:
        raise ValueError(
            f"Expected one actual linker bond at map {attachment_map}, found {len(bonds)}"
        )
    return bonds[0]


def _risk_flags(mol: Chem.Mol) -> list[str]:
    flags: list[str] = ["linker_geometry_3d_unknown"]
    azide = Chem.MolFromSmarts("[N-]=[N+]=N")
    if azide is not None and mol.HasSubstructMatch(azide):
        flags.append("reactive_azide")
    molecular_weight = float(Descriptors.MolWt(mol))
    if molecular_weight >= 900:
        flags.append("large_molecular_weight")
    if Descriptors.NumRotatableBonds(mol) >= 15:
        flags.append("high_rotatable_bond_count")
    return flags


def assemble(
    analog: dict[str, Any],
    recruiter: dict[str, Any],
    linker_template: dict[str, Any],
    orientation: str = "forward",
    attachment_map: int | None = None,
) -> dict[str, Any]:
    if type(analog) is not dict:
        raise TypeError("analog must be a dict")
    if type(recruiter) is not dict:
        raise TypeError("recruiter must be a dict")
    if type(linker_template) is not dict:
        raise TypeError("linker template must be a dict")
    if orientation not in {"forward", "reverse"}:
        raise ValueError("orientation must be 'forward' or 'reverse'")
    if attachment_map is not None and type(attachment_map) is not int:
        raise TypeError("attachment_map must be an int or None")

    options = attachment_options(analog)
    if not options:
        raise ValueError("No permitted N-H/O-H attachment handles were found")
    option_by_map = {option["atom_map"]: option for option in options}
    if attachment_map is None:
        chosen = options[0]
        attachment_map = chosen["atom_map"]
    elif attachment_map not in option_by_map:
        raise ValueError(
            f"attachment_map {attachment_map} is not a permitted option; "
            f"choose one of {sorted(option_by_map)}"
        )
    else:
        chosen = option_by_map[attachment_map]

    analog_mol = _mol_from_smiles(analog.get("mapped_smiles"), "analog mapped")
    warhead_maps = _mapped_heavy_maps(analog_mol)
    protected_warhead = _protected_maps(analog)
    if attachment_map in protected_warhead:
        raise ValueError("A protected atom cannot be selected as an attachment handle")

    analog_fragment, h_before, h_after = _attach_dummy(
        analog_mol, attachment_map, 1001
    )
    linker_fragment, linker_maps = _prepare_linker(linker_template, orientation)

    recruiter_smiles = recruiter.get("mapped_smiles")
    recruiter_mol = _mol_from_smiles(recruiter_smiles, "recruiter mapped")
    recruiter_maps = _mapped_heavy_maps(recruiter_mol)
    recruiter_dummy = [
        atom for atom in recruiter_mol.GetAtoms()
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == 1002
    ]
    if len(recruiter_dummy) != 1:
        raise ValueError("Recruiter must contain exactly one dummy with map 1002")
    if any(
        atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() != 1002
        for atom in recruiter_mol.GetAtoms()
    ):
        raise ValueError("Recruiter contains an unexpected dummy map")

    recruiter_attachment_map = recruiter.get("attachment_atom_map")
    if type(recruiter_attachment_map) is not int:
        raise TypeError("recruiter attachment_atom_map must be an int")
    recruiter_attachment = _mapped_atom(recruiter_mol, recruiter_attachment_map)
    if recruiter_attachment.GetAtomicNum() not in (7, 8):
        raise ValueError("Recruiter attachment atom must be nitrogen or oxygen")

    fragment_smiles = [
        Chem.MolToSmiles(analog_fragment, isomericSmiles=True),
        Chem.MolToSmiles(linker_fragment, isomericSmiles=True),
        Chem.MolToSmiles(recruiter_mol, isomericSmiles=True),
    ]
    product = _as_mol(join_fragments(fragment_smiles))
    Chem.SanitizeMol(product)
    Chem.AssignStereochemistry(product, cleanIt=True, force=True)

    fragments = Chem.GetMolFrags(product)
    if len(fragments) != 1:
        raise ValueError("Joined candidate is not one connected molecule")
    if any(atom.GetAtomicNum() == 0 for atom in product.GetAtoms()):
        raise ValueError("Joined candidate retains dummy atoms")

    all_role_maps = set(warhead_maps) | set(linker_maps) | set(recruiter_maps)
    product_maps = set(_mapped_heavy_maps(product))
    if product_maps != all_role_maps:
        raise ValueError(
            f"Final atom maps differ from expected roles: "
            f"missing={sorted(all_role_maps - product_maps)}, "
            f"unexpected={sorted(product_maps - all_role_maps)}"
        )

    _verify_induced_graph(
        analog_mol,
        product,
        protected_warhead,
        ignore_hydrogen_maps={attachment_map},
    )
    recruiter_protected_raw = recruiter.get("protected_core_maps")
    if (
        type(recruiter_protected_raw) is not list
        or any(type(value) is not int for value in recruiter_protected_raw)
    ):
        raise TypeError("recruiter protected_core_maps must be a list of integers")
    recruiter_protected = set(recruiter_protected_raw)
    _verify_induced_graph(recruiter_mol, product, recruiter_protected)

    warhead_bond = _actual_attachment(product, attachment_map, set(linker_maps))
    recruiter_bond = _actual_attachment(
        product, recruiter_attachment_map, set(linker_maps)
    )

    mapped_smiles = Chem.MolToSmiles(product, isomericSmiles=True, canonical=True)
    unmapped = _clear_maps(product)
    canonical_smiles = canonical(unmapped)
    if not isinstance(canonical_smiles, str):
        canonical_smiles = str(canonical_smiles)
    canonical_hash = hashlib.sha256(canonical_smiles.encode("utf-8")).hexdigest()

    e3_type = recruiter.get("e3_type")
    if e3_type == "VHL":
        calibration_status = (
            "VHL_known_saved_benchmark_not_transferable_to_new_combinations"
        )
    elif e3_type == "CRBN":
        calibration_status = "CRBN_reference_pending"
    else:
        calibration_status = "e3_specific_reference_unavailable"

    analog_id = analog.get("id", analog.get("analog_id"))
    if not isinstance(analog_id, str) or not analog_id:
        analog_id = "sha256:" + hashlib.sha256(
            analog["mapped_smiles"].encode("utf-8")
        ).hexdigest()

    result = {
        "candidate_id": f"D-{canonical_hash[:12]}",
        "e3_type": e3_type,
        "recruiter_id": recruiter.get("id"),
        "parent_warhead": analog.get("parent_warhead", analog.get("parent_id")),
        "warhead_analog_id": analog_id,
        "canonical_smiles": canonical_smiles,
        "mapped_smiles": mapped_smiles,
        "identity": _identity_value(unmapped, canonical_smiles),
        "linker_id": linker_template.get("id"),
        "orientation": orientation,
        "attachment_metadata": {
            "selected_attachment_map": attachment_map,
            "selection": chosen,
            "available_options": options,
            "hydrogens_before_attachment": h_before,
            "hydrogens_after_attachment": h_after,
            "warhead_linker_bond": warhead_bond,
            "recruiter_linker_bond": recruiter_bond,
            "graph_scope": "actual mapped product bonds; not reaction or route proof",
        },
        "protected_graph_preserved": True,
        "status": "hypothesis_pending_review",
        "final_candidate": False,
        "synthetic_route_status": "not_assessed",
        "ternary_structure_status": "not_run",
        "calibration_status": calibration_status,
        "atom_roles": {
            "warhead_maps": warhead_maps,
            "linker_maps": linker_maps,
            "recruiter_maps": recruiter_maps,
        },
        "risk_flags": _risk_flags(product),
        "linker_geometry_status": "3D_unknown",
    }
    return _json_safe(result)