"""Read-only overlays using the exact experimental bytes and stored target transform."""
import json
import math
from pathlib import Path

from packages.platform.dossiers import require,digest
from packages.platform.saved_results import BASE


def ensure_table(store):
    with store.db() as db:
        db.execute('''CREATE TABLE IF NOT EXISTS structure_sources(
            project TEXT NOT NULL,result_id TEXT NOT NULL,pdb TEXT NOT NULL,ref TEXT NOT NULL,
            PRIMARY KEY(project,result_id,pdb))''')


def register_references(saved,result_id,directory):
    ensure_table(saved.store);data=saved.view(result_id);refs={}
    for c in data['candidates']:
        reports=[saved.port.json(s['reports']['comparison'])['comparison'] for s in c['samples']]
        names={r['reference_pdb'] for r in reports};hashes={r['reference_sha256'] for r in reports}
        require(len(names)==len(hashes)==1,'OVERLAY_REFERENCE_VERSIONS')
        name=next(iter(names));raw=(Path(directory)/(name+'.cif')).read_bytes()
        require(digest(raw) in hashes,'OVERLAY_REFERENCE_HASH')
        with saved.store.db() as db:
            previous=db.execute('SELECT ref FROM structure_sources WHERE project=? AND result_id=? AND pdb=?',(saved.project,result_id,name)).fetchone()
        if previous:
            ref=json.loads(previous[0]);require(saved.port.read(ref)==raw,'OVERLAY_REFERENCE_CONFLICT')
        else:
            ref=saved.port.put_raw(raw,'chemical/x-mmcif','source')
            with saved.store.db() as db:db.execute('INSERT INTO structure_sources VALUES(?,?,?,?)',(saved.project,result_id,name,json.dumps(ref)))
        refs[name]=ref
    return refs


def structure_view(saved,result_id,candidate_id,sample_index):
    import gemmi
    import numpy as np
    from rdkit import Chem
    from packages.science.molecules import rows
    require(candidate_id in {'C01','C02'} and type(sample_index) is int and sample_index>=0,'OVERLAY_SELECTION')
    data=saved.view(result_id);c=next((c for c in data['candidates'] if c['compound_id']==candidate_id),None)
    require(c is not None and sample_index<len(c['samples']),'OVERLAY_SAMPLE')
    sample=c['samples'][sample_index];report=saved.port.json(sample['reports']['comparison'])['comparison']
    ensure_table(saved.store)
    with saved.store.db() as db:r=db.execute('SELECT ref FROM structure_sources WHERE project=? AND result_id=? AND pdb=?',(saved.project,result_id,report['reference_pdb'])).fetchone()
    require(r is not None,'OVERLAY_REFERENCE_NOT_REGISTERED')
    reference_ref=json.loads(r[0]);reference=saved.port.read(reference_ref)
    require(reference_ref['sha256']==report['reference_sha256'],'OVERLAY_REFERENCE_HASH')
    files=saved.files(result_id);matches=[ref for name,ref in files.items() if name.endswith('.cif') and ref['sha256']==sample['binding']['structure_sha256']]
    require(bool(matches),'OVERLAY_PREDICTION_MISSING');prediction=saved.port.read(matches[0])
    config=saved.port.json(files[BASE+'cpu/case-input.json'])
    candidate=next(c for c in config['candidates'] if c['id']==candidate_id)
    parts=saved.port.json(files[BASE+'cpu/'+candidate_id+'-parts.json'])['parts']
    part_maps={a.GetAtomMapNum():p['role'] for p in parts for a in Chem.MolFromSmiles(p['mapped_smiles']).GetAtoms() if a.GetAtomicNum()>0}
    # Mapping zero is an explicit equivalent label correspondence, never the best-scoring sample.
    mapping=report['ligand']['all_mappings'][0]['pairs']
    predicted_parts={a['boltz_atom_name']:part_maps[a['atom_map']] for a in mapping}
    reference_parts={a['ccd_atom_id']:part_maps[a['atom_map']] for a in mapping}
    transform=report['target_alignment'];rotation=np.array(transform['rotation_row_vectors']);translation=np.array(transform['translation_A'])
    require(rotation.shape==(3,3) and translation.shape==(3,) and np.isfinite(rotation).all() and np.isfinite(translation).all(),'OVERLAY_TRANSFORM')
    def points(raw,predicted):
        block=gemmi.cif.read_string(raw.decode()).sole_block();points=[]
        for row in rows(block,'_atom_site.'):
            chain=row['label_asym_id'];name=row['label_atom_id'];is_protein=chain in report['proteins']
            is_ligand=(name in predicted_parts and not is_protein) if predicted else chain==candidate['ligand_label_asym_id']
            if not is_protein and not is_ligand:continue
            if row['pdbx_PDB_model_num']!='1' or row['label_alt_id'] not in (False,None,'','A') or float(row['occupancy'])<=0 or row['type_symbol'].upper() in ('H','D'):continue
            # Protein C-alpha traces keep a lightweight, auditable overview; ligand uses all heavy atoms.
            if is_protein and name!='CA':continue
            xyz=np.array([float(row['Cartn_'+axis]) for axis in 'xyz'])
            if predicted:xyz=xyz@rotation+translation
            require(np.isfinite(xyz).all(),'OVERLAY_NONFINITE')
            role=('target' if chain==report['target_chain'] else 'vhl' if chain==report['vhl_chain'] else 'complex_partner') if is_protein else (predicted_parts if predicted else reference_parts)[name]
            points.append({'serial':len(points),'elem':row['type_symbol'],'atom':name,'resn':row['label_comp_id'],
                'resi':int(row['label_seq_id']) if is_protein else 1,'chain':chain,'x':float(xyz[0]),'y':float(xyz[1]),'z':float(xyz[2]),
                'part':role,'hetflag':not is_protein,'bonds':[],'bondOrder':[]})
        # Display bonds: exact source molecular graph for ligand; adjacent observed C-alpha trace segments for proteins.
        ligand_lookup={p['atom']:i for i,p in enumerate(points) if p['hetflag']}
        source_pairs={a['atom_map']:a['boltz_atom_name' if predicted else 'ccd_atom_id'] for a in mapping}
        import io
        source=next(Chem.ForwardSDMolSupplier(io.BytesIO(saved.port.read(files[BASE+'cpu/'+candidate_id+'.sdf']))))
        for bond in source.GetBonds():
            names=[source_pairs[a.GetAtomMapNum()] for a in (bond.GetBeginAtom(),bond.GetEndAtom())]
            require(all(n in ligand_lookup for n in names),'OVERLAY_LIGAND_COVERAGE')
            i,j=[ligand_lookup[n] for n in names]
            for a,b in ((i,j),(j,i)):points[a]['bonds'].append(b);points[a]['bondOrder'].append(bond.GetBondTypeAsDouble())
        previous={}
        for i,p in enumerate(points):
            if p['hetflag']:continue
            old=previous.get(p['chain']);previous[p['chain']]=i
            if old is not None and p['resi']==points[old]['resi']+1:
                distance=math.dist([p[k] for k in 'xyz'],[points[old][k] for k in 'xyz'])
                if distance<5:
                    for a,b in ((i,old),(old,i)):points[a]['bonds'].append(b);points[a]['bondOrder'].append(1)
        return points
    return {'candidate_id':candidate_id,'sample_index':sample_index,'run_id':sample['run_id'],
        'reference_pdb':report['reference_pdb'],'reference_ref':reference_ref,'prediction_ref':matches[0],
        'comparison_ref':sample['reports']['comparison'],'transform':transform,'method':report['method'],
        'mapping_index':0,'mapping_count':report['ligand']['mapping_count'],
        'reference_atoms':points(reference,False),'predicted_atoms':points(prediction,True),
        'summary':sample['summary'],'parts':sample['parts_rmsd_ranges_A'],
        'notice':'Protein C-alpha traces and ligand heavy atoms. One target transform for all parts; no VHL/ligand refit. Equivalent atom mapping 0; not a unique chemical assignment or efficacy prediction.'}
