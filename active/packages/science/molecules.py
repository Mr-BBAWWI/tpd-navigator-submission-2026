"""Fresh CCD chemistry and explicit reconstruction of reference compounds.

No efficacy scoring, implicit neutralization, or automatic attachment selection.
CCD ideal coordinates establish stereo only; they are never prediction results.
"""
from pathlib import Path

import gemmi
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Draw, rdDepictor


def rows(block, category):
    values = block.get_mmcif_category(category)
    return [dict(zip(values, row)) for row in zip(*values.values())] if values else []


def canonical(mol):
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    # Atom maps distinguish otherwise equivalent substituents during stereo assignment.
    # Remove that bookkeeping distinction before comparing chemical identity.
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    return Chem.MolToSmiles(copy, isomericSmiles=True)


def atom_index(mol, name):
    matches = [a.GetIdx() for a in mol.GetAtoms()
               if a.HasProp('ccd_atom_id') and a.GetProp('ccd_atom_id') == name]
    if len(matches) != 1:
        raise ValueError(f'Atom name must identify exactly one atom: {name}')
    return matches[0]


def read_ccd(path):
    block = gemmi.cif.read_file(str(path)).sole_block()
    atoms = rows(block, '_chem_comp_atom.')
    bonds = rows(block, '_chem_comp_bond.')
    if not atoms or not bonds:
        raise ValueError('CCD atom/bond tables are required')
    rw = Chem.RWMol()
    names = {}
    conf = Chem.Conformer(len(atoms))
    conf.Set3D(True)
    for row in atoms:
        name = row['atom_id']
        if name in names:
            raise ValueError(f'Duplicate CCD atom: {name}')
        atom = Chem.Atom(row['type_symbol'].title())
        atom.SetFormalCharge(int(row['charge']))
        atom.SetProp('ccd_atom_id', name)
        idx = rw.AddAtom(atom)
        names[name] = idx
        conf.SetAtomPosition(idx, tuple(float(row['pdbx_model_Cartn_'+axis+'_ideal'])
                                       for axis in 'xyz'))
    orders = {'SING': Chem.BondType.SINGLE, 'DOUB': Chem.BondType.DOUBLE,
              'TRIP': Chem.BondType.TRIPLE, 'AROM': Chem.BondType.AROMATIC}
    for row in bonds:
        rw.AddBond(names[row['atom_id_1']], names[row['atom_id_2']], orders[row['value_order']])
    mol = rw.GetMol()
    mol.AddConformer(conf)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    mol = Chem.RemoveHs(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for row in atoms:
        expected = row['pdbx_stereo_config']
        if expected in ('R', 'S'):
            atom = mol.GetAtomWithIdx(atom_index(mol, row['atom_id']))
            if not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode') != expected:
                raise ValueError(f'CCD stereochemistry disagreement: {row["atom_id"]}')
    descriptors = rows(block, '_pdbx_chem_comp_descriptor.')
    stereo_smiles = [r['descriptor'] for r in descriptors if r['type'] == 'SMILES_CANONICAL']
    if not stereo_smiles or any(canonical(Chem.MolFromSmiles(s)) != canonical(mol)
                               for s in stereo_smiles):
        raise ValueError('CCD graph disagrees with independent canonical SMILES descriptors')
    expected_keys = {r['descriptor'] for r in descriptors if r['type'] == 'InChIKey'}
    if expected_keys != {Chem.MolToInchiKey(mol)}:
        raise ValueError('CCD InChIKey mismatch')
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError('Disconnected CCD compound')
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    mol.SetProp('_Name', block.name)
    return mol


def identity(mol):
    clean = Chem.MolFromSmiles(canonical(mol))
    return {'canonical_isomeric_smiles': canonical(mol),
            'inchikey': Chem.MolToInchiKey(clean),
            'formula': rdMolDescriptors.CalcMolFormula(clean),
            'formal_charge': Chem.GetFormalCharge(clean), 'heavy_atom_count': clean.GetNumHeavyAtoms(),
            'stereocenters': Chem.FindMolChiralCenters(clean, includeUnassigned=True),
            'properties': {'molecular_weight_Da': Descriptors.MolWt(clean),
                           'TPSA_A2': rdMolDescriptors.CalcTPSA(clean),
                           'cLogP_RDKit': Crippen.MolLogP(clean),
                           'HBD': rdMolDescriptors.CalcNumHBD(clean),
                           'HBA': rdMolDescriptors.CalcNumHBA(clean),
                           'rotatable_bonds_RDKit_strict': rdMolDescriptors.CalcNumRotatableBonds(clean)}}


def join_fragments(smiles):
    """Join pairs of mapped dummy atoms; persistent heavy-atom maps must be unique."""
    combined = None
    dummies, heavy_maps = {}, set()
    for value in smiles:
        mol = Chem.MolFromSmiles(value)
        if mol is None or len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError('Each component must be one valid connected molecule')
        for atom in mol.GetAtoms():
            key = atom.GetAtomMapNum()
            if key <= 0:
                raise ValueError('Every fragment atom needs an explicit positive map')
            if atom.GetAtomicNum() == 0:
                if atom.GetDegree() != 1:
                    raise ValueError('Attachment dummy must have one neighbor')
                dummies[key] = dummies.get(key, 0) + 1
            else:
                if key in heavy_maps:
                    raise ValueError('Duplicate heavy-atom map')
                heavy_maps.add(key)
        combined = mol if combined is None else Chem.CombineMols(combined, mol)
    if not dummies or any(n != 2 for n in dummies.values()) or set(dummies) & heavy_maps:
        raise ValueError('Each attachment map must occur on exactly two dummies only')
    result = Chem.molzip(combined)
    Chem.SanitizeMol(result)
    if len(Chem.GetMolFrags(result)) != 1 or any(a.GetAtomicNum() == 0 for a in result.GetAtoms()):
        raise ValueError('Assembly left disconnected or unpaired atoms')
    if {a.GetAtomMapNum() for a in result.GetAtoms()} != heavy_maps:
        raise ValueError('Assembly changed atom identity')
    return result


def reconstruct(mol, cuts):
    bond_ids = []
    for pair in cuts:
        bond = mol.GetBondBetweenAtoms(*(atom_index(mol, n) for n in pair))
        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
            raise ValueError(f'Cut must be a non-ring single bond: {pair}')
        bond_ids.append(bond.GetIdx())
    if len(set(bond_ids)) != len(bond_ids):
        raise ValueError('Duplicate cut')
    fragment = Chem.FragmentOnBonds(mol, bond_ids, addDummies=True,
                                   dummyLabels=[(1001+i, 1001+i) for i in range(len(cuts))])
    for atom in fragment.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomMapNum(atom.GetIsotope())
            atom.SetIsotope(0)
    parts = Chem.GetMolFrags(fragment, asMols=True)
    if len(parts) != len(cuts) + 1:
        raise ValueError('Cuts did not produce the expected number of parts')
    # Serialize and independently reparse: reconstruction cannot retain reference coordinates.
    part_smiles = [Chem.MolToSmiles(p, isomericSmiles=True) for p in parts]
    joined = join_fragments(part_smiles)
    if canonical(joined) != canonical(mol):
        raise ValueError('Reconstructed candidate differs from reference chemistry/stereo')
    labels = {a.GetAtomMapNum(): a.GetProp('ccd_atom_id') for a in mol.GetAtoms()}
    for a in joined.GetAtoms():
        a.SetProp('ccd_atom_id', labels[a.GetAtomMapNum()])
    return joined, part_smiles


def map_warhead(start, candidate, start_attachment, candidate_attachment):
    query = Chem.Mol(start)
    atom = query.GetAtomWithIdx(atom_index(query, start_attachment))
    before = {'formal_charge': atom.GetFormalCharge(), 'total_H': atom.GetTotalNumHs()}
    # Case-specific, recorded +1 piperazinium -> neutral substituted amine comparison.
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1 or atom.GetDegree() != 2:
        raise ValueError('Expected a secondary piperazinium attachment atom')
    atom.SetFormalCharge(0)
    atom.SetNumExplicitHs(0)
    atom.SetNoImplicit(False)
    Chem.SanitizeMol(query)
    qi = atom_index(query, start_attachment)
    ci = atom_index(candidate, candidate_attachment)
    matches = [m for m in candidate.GetSubstructMatches(query, uniquify=False, useChirality=True)
               if m[qi] == ci]
    if not matches:
        raise ValueError('Starting warhead is absent at the configured attachment atom')
    chosen = min(matches)
    external = [(i, n.GetIdx()) for i in chosen for n in candidate.GetAtomWithIdx(i).GetNeighbors()
                if n.GetIdx() not in chosen]
    if len(external) != 1 or external[0][0] != ci:
        raise ValueError('Warhead has unexpected external substitutions')
    return {'status': 'computed_pending_human_review', 'equivalent_mapping_count': len(matches),
            'mapping_choice': 'lexicographically first; symmetry-equivalent ring atom maps are not unique',
            'normalization_for_matching_only': {'atom': start_attachment, 'before': before,
                                               'after_formal_charge': 0},
            'atoms': [{'source_ccd_atom': a.GetProp('ccd_atom_id'),
                       'candidate_ccd_atom': candidate.GetAtomWithIdx(chosen[a.GetIdx()]).GetProp('ccd_atom_id'),
                       'candidate_map': candidate.GetAtomWithIdx(chosen[a.GetIdx()]).GetAtomMapNum()}
                      for a in query.GetAtoms()]}


def write_sdf(mol, path):
    # Python opens Unicode Windows paths; RDKit's narrow filename API may fail.
    with Path(path).open('w', encoding='utf-8', newline='\n') as stream:
        with Chem.SDWriter(stream) as writer:
            writer.write(mol)


def write_molecule(mol, stem):
    stem = Path(stem)
    copy = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(copy)
    copy.SetProp('coordinate_kind', '2D depiction, not an experimental or predicted conformation')
    copy.SetProp('ccd_atom_names_in_sdf_order', ','.join(a.GetProp('ccd_atom_id') for a in copy.GetAtoms()))
    write_sdf(copy, stem.with_suffix('.sdf'))
    stem.with_suffix('.smi').write_text(canonical(mol)+'\n', encoding='utf-8')
    for a in copy.GetAtoms():
        a.SetAtomMapNum(0)
    Draw.MolToFile(copy, str(stem.with_suffix('.svg')), size=(1100, 550))
