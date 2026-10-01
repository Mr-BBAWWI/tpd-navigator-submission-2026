"""Bounded R5 linker hypotheses. Graph validation does not validate binding or synthesis."""
import copy
import json
from pathlib import Path

from packages.contracts import ContractError

LIBRARY_PATH = Path(__file__).resolve().parents[2] / 'cases/linker_library_20260930.json'


def require(value, message):
    if not value:
        raise ContractError(message)


def library():
    return json.loads(LIBRARY_PATH.read_text(encoding='utf-8'))


def validate_parameters(parameters):
    require(type(parameters) is dict, 'CANDIDATE_PARAMETERS')
    require(type(parameters.get('parent', 'C01')) is str and parameters.get('parent', 'C01') in {'C01', 'C02'}, 'CANDIDATE_PARAMETERS')
    if 'linker_id' in parameters:
        require(set(parameters) <= {'parent', 'linker_id', 'orientation'}, 'CANDIDATE_MIXED_PARAMETERS')
        require(type(parameters['linker_id']) is str and parameters['linker_id'] in {x['id'] for x in library()['templates']}, 'CANDIDATE_LINKER_UNKNOWN')
        require(type(parameters.get('orientation', 'forward')) is str and parameters.get('orientation', 'forward') in {'forward', 'reverse'}, 'CANDIDATE_ORIENTATION')
    else:
        require(set(parameters) <= {'parent', 'linker_extension'}, 'CANDIDATE_PARAMETERS')
        extension = parameters.get('linker_extension', 1)
        require(type(extension) is int and extension in {1, 2}, 'CANDIDATE_PARAMETERS')


def design_class(result):
    # Presentation annotation for immutable pre-R5 results, not a rewritten opinion/decision.
    if result.get('transformation', {}).get('kind') == 'linker_boundary_methylene_extension':
        return 'linker_perturbation_prototype'
    return result.get('design_class', 'unclassified_design_hypothesis')


def classification_label(result):
    return ('초기 linker perturbation 기능 검증용 prototype'
            if design_class(result) == 'linker_perturbation_prototype'
            else 'Linker 설계 가설 · 과학 검토 대기')


def _mapped_atoms(mol):
    pairs = [(a.GetAtomMapNum(), a) for a in mol.GetAtoms() if a.GetAtomicNum()]
    require(all(k > 0 for k, _ in pairs) and len({k for k, _ in pairs}) == len(pairs), 'CANDIDATE_ATOM_MAPS')
    return dict(pairs)


def _bonds(mol, allowed):
    return {frozenset((b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum())):
            (b.GetBondTypeAsDouble(), str(b.GetStereo())) for b in mol.GetBonds()
            if b.GetBeginAtom().GetAtomMapNum() in allowed and b.GetEndAtom().GetAtomMapNum() in allowed}


def _stereo_by_map(mol):
    from rdkit import Chem
    copy_mol = Chem.Mol(mol)
    Chem.AssignStereochemistry(copy_mol, cleanIt=True, force=True)
    return {copy_mol.GetAtomWithIdx(i).GetAtomMapNum(): label
            for i, label in Chem.FindMolChiralCenters(copy_mol, includeUnassigned=True)}


def linker_metrics(mol):
    from rdkit import Chem
    ends = sorted((a.GetAtomMapNum(), a.GetNeighbors()[0].GetIdx()) for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
    require(len(ends) == 2, 'CANDIDATE_LINKER_BOUNDARIES')
    path = Chem.GetShortestPath(mol, ends[0][1], ends[1][1])
    return {'heavy_atom_count': mol.GetNumHeavyAtoms(),
            'heteroatom_count': sum(a.GetAtomicNum() not in (0, 1, 6) for a in mol.GetAtoms()),
            'ring_count': mol.GetRingInfo().NumRings(),
            'terminal_to_terminal_shortest_path_bonds': len(path) - 1,
            'length_definition': 'Graph path between terminal linker heavy atoms; excludes ligand attachments, not a 3D distance.'}


def replace_linker(parts, parameters):
    from rdkit import Chem
    from packages.science.molecules import join_fragments, canonical, identity
    validate_parameters(parameters)
    source = library()
    template = next(x for x in source['templates'] if x['id'] == parameters['linker_id'])
    require(sorted(p['role'] for p in parts) == ['linker', 'recruiter', 'warhead'], 'CANDIDATE_FRAGMENT_ROLES')
    original = [p['mapped_smiles'] for p in parts]
    parent = join_fragments(original)
    index = next(i for i, p in enumerate(parts) if p['role'] == 'linker')
    old_linker = Chem.MolFromSmiles(original[index])
    protected = set()
    endpoints = []
    for p in parts:
        if p['role'] == 'linker':
            continue
        fragment = Chem.MolFromSmiles(p['mapped_smiles'])
        protected.update(_mapped_atoms(fragment))
        dummies = [a for a in fragment.GetAtoms() if not a.GetAtomicNum()]
        require(len(dummies) == 1 and dummies[0].GetDegree() == 1, 'CANDIDATE_ATTACHMENT')
        dummy = dummies[0]
        neighbor = dummy.GetNeighbors()[0]
        expected = (1001, 'N') if p['role'] == 'warhead' else (1002, 'O')
        require((dummy.GetAtomMapNum(), neighbor.GetSymbol()) == expected, 'CANDIDATE_ATTACHMENT_CONTEXT')
        endpoints.append({'role': p['role'], 'attachment_dummy_map': expected[0],
                          'retained_atom_map': neighbor.GetAtomMapNum(), 'element': neighbor.GetSymbol(),
                          'selection': 'inherited_reference_attachment_not_new_exit_vector_prediction',
                          'bond_kind': 'N-C' if p['role'] == 'warhead' else 'O-C'})
    changed = Chem.MolFromSmiles(template['smiles'])
    require(changed is not None and len(Chem.GetMolFrags(changed)) == 1, 'CANDIDATE_TEMPLATE_INVALID')
    new_maps = []
    for atom in changed.GetAtoms():
        if atom.GetAtomicNum() == 0:
            require(atom.GetAtomMapNum() in (1001, 1002) and atom.GetDegree() == 1, 'CANDIDATE_TEMPLATE_ATTACHMENT')
            require(atom.GetNeighbors()[0].GetAtomicNum() == 6 and not atom.GetNeighbors()[0].GetIsAromatic(), 'CANDIDATE_TERMINAL_CARBON_REQUIRED')
            if parameters.get('orientation', 'forward') == 'reverse':
                atom.SetAtomMapNum(2003 - atom.GetAtomMapNum())
        else:
            mapping = 2001 + len(new_maps)
            require(mapping not in _mapped_atoms(parent), 'CANDIDATE_MAP_COLLISION')
            atom.SetAtomMapNum(mapping)
            new_maps.append(mapping)
    fragments = original.copy()
    fragments[index] = Chem.MolToSmiles(changed, isomericSmiles=True)
    molecule = join_fragments(fragments)
    require(canonical(molecule) != canonical(parent), 'CANDIDATE_IDENTICAL_TO_PARENT')
    old_atoms, new_atoms = _mapped_atoms(parent), _mapped_atoms(molecule)
    removed_maps = sorted(set(old_atoms) - protected)
    require(set(new_atoms) == protected | set(new_maps), 'CANDIDATE_REPLACEMENT_MAPS')
    for mapping in protected:
        a, b = old_atoms[mapping], new_atoms[mapping]
        features = lambda x: (x.GetAtomicNum(), x.GetIsotope(), x.GetFormalCharge(), x.GetIsAromatic(), x.GetTotalNumHs())
        require(features(a) == features(b), 'CANDIDATE_PROTECTED_ATOM_CHANGED')
    require(_bonds(parent, protected) == _bonds(molecule, protected), 'CANDIDATE_PROTECTED_BOND_CHANGED')
    old_stereo, new_stereo = _stereo_by_map(parent), _stereo_by_map(molecule)
    require({m:s for m,s in old_stereo.items() if m in protected} ==
            {m:s for m,s in new_stereo.items() if m in protected}, 'CANDIDATE_PROTECTED_STEREO_CHANGED')
    for endpoint in endpoints:
        atom = new_atoms[endpoint['retained_atom_map']]
        added_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomMapNum() in new_maps]
        require(len(added_neighbors) == 1, 'CANDIDATE_ATTACHMENT_CHANGED')
        neighbor = added_neighbors[0]
        require(neighbor.GetAtomicNum() == 6 and molecule.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE,
                'CANDIDATE_ATTACHMENT_BOND')
        endpoint['new_linker_atom_map'] = neighbor.GetAtomMapNum()
    parent_info, info = identity(parent), identity(molecule)
    comparison = {'parent': linker_metrics(old_linker), 'hypothesis': linker_metrics(changed),
                  'whole_molecule_property_delta': {k: info['properties'][k] - v for k, v in parent_info['properties'].items()},
                  'meaning': 'Calculated graph descriptors; not measured permeability, linker strain, binding or efficacy.'}
    return {'molecule': molecule, 'parent': parent, 'fragments': fragments, 'identity': info,
            'new_maps': new_maps, 'removed_maps': removed_maps, 'protected_maps': sorted(protected),
            'template': copy.deepcopy(template), 'library_version': source['version'],
            'sources': [s for s in source['sources'] if s['id'] in template['source_ids']],
            'attachments': endpoints, 'comparison': comparison}
