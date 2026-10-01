"""Explicit mmCIF chain/atom selection and descriptive CPU geometry."""
import math

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFreeSASA

from .molecules import atom_index, rows, canonical


def atom_sites(path, asym_ids, altloc='A'):
    block = gemmi.cif.read_file(str(path)).sole_block()
    selected, seen = [], set()
    for r in rows(block, '_atom_site.'):
        if r['label_asym_id'] not in asym_ids or r['pdbx_PDB_model_num'] != '1':
            continue
        if r['type_symbol'] in ('H', 'D') or r['label_alt_id'] not in (False, None, '', altloc):
            continue
        if float(r['occupancy']) <= 0:
            continue
        key = (r['label_asym_id'], r['auth_seq_id'], r['pdbx_PDB_ins_code'], r['label_atom_id'])
        if key in seen:
            raise ValueError(f'Ambiguous selected atom: {key}')
        seen.add(key)
        r['xyz'] = [float(r['Cartn_'+a]) for a in 'xyz']
        selected.append(r)
    if {r['label_asym_id'] for r in selected} != set(asym_ids):
        raise ValueError('A selected chain has no eligible heavy atoms')
    return selected


def ligand_coordinates(path, mol, asym_id, comp_id):
    sites = atom_sites(path, [asym_id])
    if {r['label_comp_id'] for r in sites} != {comp_id}:
        raise ValueError('Selected ligand chain has unexpected chemical components')
    by_name = {r['label_atom_id']: r for r in sites}
    if len(by_name) != len(sites) or set(by_name) != {a.GetProp('ccd_atom_id') for a in mol.GetAtoms()}:
        raise ValueError('Missing, extra or duplicate observed ligand heavy atoms')
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    mapping = []
    for atom in mol.GetAtoms():
        r = by_name[atom.GetProp('ccd_atom_id')]
        if r['type_symbol'].title() != atom.GetSymbol():
            raise ValueError('CCD/coordinate element mismatch')
        conf.SetAtomPosition(atom.GetIdx(), r['xyz'])
        mapping.append({'atom_map': atom.GetAtomMapNum(), 'rdkit_index_zero_based': atom.GetIdx(),
                        'ccd_atom_id': r['label_atom_id'], 'atom_site_id': r['id'],
                        'label_asym_id': asym_id, 'auth_asym_id': r['auth_asym_id'],
                        'auth_seq_id': r['auth_seq_id'], 'occupancy': float(r['occupancy']),
                        'altloc': r['label_alt_id'] or None})
    positioned = Chem.Mol(mol)
    positioned.RemoveAllConformers()
    positioned.AddConformer(conf)
    observed_stereo = Chem.Mol(positioned)
    Chem.AssignStereochemistryFrom3D(observed_stereo, replaceExistingTags=True)
    if canonical(observed_stereo) != canonical(mol):
        raise ValueError('Observed ligand coordinates disagree with CCD stereochemistry')
    # Check observed bond lengths without assigning an efficacy/pose score.
    lengths = [conf.GetAtomPosition(b.GetBeginAtomIdx()).Distance(conf.GetAtomPosition(b.GetEndAtomIdx()))
               for b in mol.GetBonds()]
    if any(not math.isfinite(x) or x < 0.5 or x > 2.5 for x in lengths):
        raise ValueError('Implausible observed covalent bond distance; inspect structure')
    return positioned, mapping


def assembly_proteins(path, assembly_id, chain_ids):
    block = gemmi.cif.read_file(str(path)).sole_block()
    assemblies = [r for r in rows(block, '_pdbx_struct_assembly_gen.') if r['assembly_id'] == assembly_id]
    if len(assemblies) != 1 or assemblies[0]['oper_expression'] != '1':
        raise ValueError('This adapter only supports an explicitly selected identity-operation assembly')
    assembly_chains = assemblies[0]['asym_id_list'].split(',')
    if not set(chain_ids) <= set(assembly_chains):
        raise ValueError('Selected protein chains are outside the chosen assembly')
    entities = {r['entity_id']: r for r in rows(block, '_entity_poly.')}
    asym = {r['id']: r['entity_id'] for r in rows(block, '_struct_asym.')}
    references = rows(block, '_struct_ref.')
    alignments = rows(block, '_struct_ref_seq.')
    result = []
    for chain in chain_ids:
        row = entities[asym[chain]]
        if row['type'] != 'polypeptide(L)':
            raise ValueError('Only L-polypeptide inputs are supported')
        sequence = ''.join(row['pdbx_seq_one_letter_code_can'].split())
        if not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('Unsupported or ambiguous protein sequence')
        sites = atom_sites(path, [chain])
        auth_chains = sorted({r['auth_asym_id'] for r in sites})
        result.append({'label_asym_id': chain, 'entity_id': asym[chain],
                       'auth_asym_ids': auth_chains,
                       'database_references': [r for r in references if r['entity_id'] == asym[chain]],
                       'database_alignments': [r for r in alignments if r['pdbx_strand_id'] in auth_chains],
                       'sequence': sequence, 'length': len(sequence),
                       'modeled_residue_count': len({r['label_seq_id'] for r in sites}),
                       'sequence_source': '_entity_poly.pdbx_seq_one_letter_code_can; includes deposited construct tags'})
    return result, assembly_chains


def attachment_geometry(path, ligand, ligand_asym, protein_asym, attachment, cutoff=4.0):
    if cutoff <= 0 or not math.isfinite(cutoff):
        raise ValueError('Contact reporting radius must be finite and positive')
    protein = atom_sites(path, [protein_asym])
    coord = np.array(ligand.GetConformer().GetPositions())
    idx = atom_index(ligand, attachment)
    distances = np.linalg.norm(np.array([r['xyz'] for r in protein]) - coord[idx], axis=1)
    contacts = [{'label_asym_id': r['label_asym_id'], 'auth_asym_id': r['auth_asym_id'],
                 'residue': r['label_comp_id'], 'auth_seq_id': r['auth_seq_id'],
                 'atom': r['label_atom_id'], 'distance_A': round(float(d), 4)}
                for r, d in zip(protein, distances) if d <= cutoff]
    rw = Chem.RWMol()
    xyz = list(coord) + [r['xyz'] for r in protein]
    elements = [a.GetSymbol() for a in ligand.GetAtoms()] + [r['type_symbol'] for r in protein]
    conf = Chem.Conformer(len(elements))
    radii = []
    for i, (element, point) in enumerate(zip(elements, xyz)):
        rw.AddAtom(Chem.Atom(element))
        conf.SetAtomPosition(i, tuple(float(x) for x in point))
        radii.append(Chem.GetPeriodicTable().GetRvdw(element))
    combined = rw.GetMol()
    combined.AddConformer(conf)
    options = rdFreeSASA.SASAOpts()
    options.probeRadius = 1.4
    rdFreeSASA.CalcSASA(combined, radii, opts=options)
    isolated = Chem.Mol(ligand)
    rdFreeSASA.CalcSASA(isolated, radii[:len(coord)], opts=options)
    return {'kind': 'computed_descriptive_geometry', 'attachment_ccd_atom': attachment,
            'ligand_label_asym_id': ligand_asym, 'protein_label_asym_id': protein_asym,
            'SASA_A2': combined.GetAtomWithIdx(idx).GetDoubleProp('SASA'),
            'isolated_ligand_same_pose_SASA_A2': isolated.GetAtomWithIdx(idx).GetDoubleProp('SASA'),
            'nearest_protein_heavy_atom_A': float(distances.min()), 'contacts': contacts,
            'method': {'tool': 'RDKit rdFreeSASA', 'algorithm': str(options.algorithm),
                       'radii': 'RDKit periodic-table van der Waals radii, explicit for every atom',
                       'probe_radius_A': 1.4, 'contact_reporting_radius_A': cutoff,
                       'hydrogens': 'excluded', 'waters_ions_other_chains_crystal_mates': 'excluded',
                       'model': 1, 'altloc': 'blank or A'},
            'interpretation': 'Descriptive only; no exposure cutoff, hydrogen-bond claim, efficacy or attachment approval'}


def aligned_rmsd(reference_xyz, moving_xyz):
    """Kabsch helper for future prediction comparisons; caller must supply validated atom pairs."""
    ref, mov = np.asarray(reference_xyz, dtype=float), np.asarray(moving_xyz, dtype=float)
    if ref.shape != mov.shape or ref.ndim != 2 or ref.shape[1] != 3 or len(ref) < 3:
        raise ValueError('At least three paired 3D coordinates required')
    if not np.isfinite(ref).all() or not np.isfinite(mov).all():
        raise ValueError('Nonfinite coordinates')
    x, y = mov - mov.mean(0), ref - ref.mean(0)
    if np.linalg.matrix_rank(x) < 2 or np.linalg.matrix_rank(y) < 2:
        raise ValueError('Collinear points cannot establish an unambiguous alignment')
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.eye(3)
    correction[2, 2] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    translation = ref.mean(0) - mov.mean(0) @ rotation
    fitted = mov @ rotation + translation
    return {'rmsd_A': float(np.sqrt(np.mean(np.sum((fitted-ref)**2, axis=1)))),
            'rotation_row_vectors': rotation.tolist(), 'translation_A': translation.tolist()}
