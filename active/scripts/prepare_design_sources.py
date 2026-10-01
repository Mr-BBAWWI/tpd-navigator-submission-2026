"""Rebuild reviewed design metadata from explicit, newly retrieved PDB/CCD inputs."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rdkit import Chem
from packages.science.molecules import read_ccd, identity, rows
from packages.science.structures import ligand_coordinates, atom_sites
import gemmi

def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')

def main():
    source = ROOT/'cases/design_sources'
    parent = read_ccd(source/'FX5.cif')
    parent, mapping = ligand_coordinates(source/'6HAZ.cif', parent, 'C', 'FX5')
    # Preserve the deposited +1 species, and declare the separate neutral design state.
    neutral = Chem.Mol(parent)
    nitrogen = next(a for a in neutral.GetAtoms() if a.GetAtomMapNum() == 19)
    assert nitrogen.GetFormalCharge() == 1 and nitrogen.GetTotalNumHs() == 2
    nitrogen.SetFormalCharge(0); nitrogen.SetNumExplicitHs(1)
    Chem.SanitizeMol(neutral)
    for name, mol in [('SMARCA2-parent', parent), ('SMARCA2-neutral-design', neutral)]:
        (source/(name+'.sdf')).write_text(Chem.MolToMolBlock(mol)+'\n$$$$\n', encoding='utf-8')
    protected = [a.GetAtomMapNum() for a in neutral.GetAtoms() if a.GetAtomMapNum() not in {8,9,11,12,19}]
    meta = {'id':'SMARCA2-FX5','target':'SMARCA2','pdb':'6HAZ','ccd':'FX5',
        'target_chain':'A','ligand_chain':'C','evidence_tier':'co_crystal',
        'parent_identity':identity(parent),'design_identity':identity(neutral),'atom_mapping':mapping,
        'protected_maps':protected,'known_attachment_map':19,
        'protected_basis':{'source':'https://doi.org/10.1038/s41589-019-0294-6',
            'locator':'Results: Identification of a partial SMARCA2/4 degrader; Figure 1a',
            'key_interactions':['phenol / Y1421','aminopyridazine / N1464'],
            'scope':'The whole aminopyridazine/phenol binding scaffold and its connecting N are conservatively protected; not every atom is claimed essential.'},
        'chemical_state':{'source_charge':1,'design_charge':0,'modified_map':19,
            'operation':'explicit deprotonation of terminal piperazinium to secondary amine',
            'status':'declared_design_microstate_not_pH_prediction','source_unchanged':True,
            'rationale':'Graph derivatization requires an explicit free-base model. Binding and population of this state remain unvalidated.'},
        'sar':{'19':{'status':'known_attachment_precedent',
            'source':'https://doi.org/10.1038/s41589-019-0294-6',
            'locator':'Results / Identification of a partial SMARCA2/4 degrader; Fig. 1a; compounds 1 and 2',
            'allowed_transformations':['linker_handle_introduction'],
            'limitation':'Known N-substitution precedent is not proof that all substitutions retain activity.'}},
        'unknown_sites':[8,9,11,12], 'new_activity_values':None,
        'limitations':['Only one atom-specific attachment precedent is curated. Two SAR-qualified sites have not been established.',
            'No target-wide 5–10 experimentally supported parent ligands are claimed from this single registered case.']}
    dump(source/'warhead.json',meta)
    parts=json.loads((ROOT/'outputs/smarca2_20260923/C01-parts.json').read_text())['parts']
    vhl=next(p for p in parts if p['role']=='recruiter')['mapped_smiles']
    crbn=read_ccd(source/'RN6.cif')
    # dBET6 C15-O bond defines the experimentally used phenoxy recruiter vector.
    bond=crbn.GetBondBetweenAtoms(2,9)
    assert {crbn.GetAtomWithIdx(i).GetProp('ccd_atom_id') for i in (2,9)}=={'O','C15'}
    fragmented=Chem.FragmentOnBonds(crbn,[bond.GetIdx()],addDummies=True,dummyLabels=[(0,0)])
    fragments=Chem.GetMolFrags(fragmented,asMols=True,sanitizeFrags=True)
    core=next(m for m in fragments if any(a.GetAtomMapNum()==3 for a in m.GetAtoms()))
    dummies=[a for a in core.GetAtoms() if a.GetAtomicNum()==0];assert len(dummies)==1
    dummies[0].SetAtomMapNum(1002);dummies[0].SetIsotope(0)
    recruiters=[]
    for name,e3,mol,attach,source_id,description in [
        ('VHL-FX8-O37','VHL',Chem.MolFromSmiles(vhl),37,'6HAY/FX8','VHL ligand fragment from PROTAC 1, reference phenoxy exit vector'),
        ('CRBN-RN6-O','CRBN',core,3,'6BOY/RN6','4-oxy-thalidomide-derived recruiter from dBET6; deposited S glutarimide stereochemistry')]:
        # Shared namespace: warhead <1000; linker 2000+; recruiter 4000+.
        mapping={}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum():
                old=atom.GetAtomMapNum();mapping[old]=4000+old;atom.SetAtomMapNum(4000+old)
        Chem.SanitizeMol(mol)
        recruiters.append({'id':name,'e3_type':e3,'mapped_smiles':Chem.MolToSmiles(mol),
            'attachment_atom_map':mapping[attach],'dummy_map':1002,'attachment_chemistry':'O-C single bond',
            'source_id':source_id,'source_url':'https://www.rcsb.org/structure/'+source_id.split('/')[0],
            'description':description,'protected_core_maps':sorted(mapping.values()),
            'attachment_precedent':'observed covalent graph in the deposited PROTAC; not a proven route for new combinations',
            'stereochemistry':'retained from deposited CCD; no epimerization/population prediction',
            'e3_binding_new_combinations':'not_tested'})
    dump(source/'recruiters.json',{'version':'dual-e3/20260930.1','recruiters':recruiters})
    benchmark={'id':'CRBN-dBET6-6BOY','target':'BRD4_BD1','e3_type':'CRBN','recruiter_id':'CRBN-RN6-O',
        'pdb':'6BOY','ccd':'RN6','reference_target_chain':'C','reference_e3_chain':'B','reference_ligand_chain':'E',
        'source_url':'https://www.rcsb.org/structure/6BOY','doi':'10.1038/s41589-018-0055-y',
        'role':'independently retrieved CRBN calibration reference; not legacy BRD4 data or model reuse',
        'prediction_status':'not_run','calibration_status':'reference_ready_prediction_pending',
        'comparison_metrics':['target_CA_RMSD_A','e3_CA_RMSD_after_target_alignment_A','ligand_RMSD_after_target_alignment_A','seed_consistency'],
        'prohibited_inference':'CRBN versus VHL efficacy or raw-confidence ranking'}
    dump(source/'crbn_benchmark.json',benchmark)
    source_urls={n:'https://files.rcsb.org/'+('download/' if n.startswith(('6BOY','6HAZ')) else 'ligands/download/')+n for n in ['6BOY.cif','RN6.cif','6HAZ.cif','FX5.cif']}
    dump(source/'manifest.json',{'version':'design-sources/20260930.1',
        'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.iterdir()) if p.is_file() and p.name!='manifest.json'},
        'official_sources':source_urls,'expert_documents':json.loads((ROOT/'.localdata/expert-design-20260930/sources.json').read_text(encoding='utf-8')),
        'receptor_preparation':{'tool':'Meeko 0.8.0','input':'6HAZ chain A heavy atoms','output':'SMARCA2-receptor.pdbqt',
            'charges':'Gasteiger template','residues_deleted':False,'pH_or_tautomer_population':'not_predicted',
            'waters_ions_crystal_neighbors':'excluded','status':'fixed receptor model for preliminary docking'}})
    print('Prepared source-bound warhead, 2 recruiter branches, CRBN reference and hash manifest.')

if __name__=='__main__':main()
