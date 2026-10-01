"""M2 adapters invoking science functions on registered immutable molecule sources."""
import io
import json

from packages.platform.dossiers import require,digest
from packages.contracts import encoded
from packages.platform.saved_results import AUTHORITY,BASE


def run_tool(operation,saved,result_id,parameters,put,check_active):
    if operation=='candidate':return design_candidate(saved,result_id,parameters,put,check_active)
    if operation=='chemistry':return hydrogen_trial(saved,result_id,parameters,put,check_active)
    require(False,'SCIENCE_TOOL_NOT_ALLOWED')


def design_candidate(saved,result_id,parameters,put,check_active):
    from packages.science.linker_design import validate_parameters
    validate_parameters(parameters)
    if 'linker_id' in parameters:
        return design_linker_candidate(saved,result_id,parameters,put,check_active)
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import rdDepictor,Draw
    from packages.science.molecules import join_fragments,identity,canonical
    parent=parameters.get('parent','C01');extension=parameters.get('linker_extension',1)
    require(parent in {'C01','C02'} and type(extension) is int and extension in {1,2},'CANDIDATE_PARAMETERS')
    files=saved.files(result_id);part_ref=files[BASE+'cpu/'+parent+'-parts.json']
    parts=saved.port.json(part_ref)['parts'];original=[p['mapped_smiles'] for p in parts]
    parent_mol=join_fragments(original)
    index=next(i for i,p in enumerate(parts) if p['role']=='linker')
    linker=Chem.MolFromSmiles(original[index]);rw=Chem.RWMol(linker)
    dummies=sorted((a.GetAtomMapNum(),a.GetIdx()) for a in linker.GetAtoms() if a.GetAtomicNum()==0)
    require(len(dummies)==2,'CANDIDATE_LINKER_BOUNDARIES')
    dummy_map,dummy=dummies[0];neighbor=linker.GetAtomWithIdx(dummy).GetNeighbors()[0].GetIdx()
    bond=linker.GetBondBetweenAtoms(dummy,neighbor)
    require(bond.GetBondType()==Chem.BondType.SINGLE,'CANDIDATE_LINKER_BOND')
    rw.RemoveBond(dummy,neighbor);previous=neighbor;new_maps=[]
    used={a.GetAtomMapNum() for a in parent_mol.GetAtoms()}
    for offset in range(extension):
        atom=Chem.Atom(6);mapping=2001+offset;require(mapping not in used,'CANDIDATE_MAP_COLLISION')
        atom.SetAtomMapNum(mapping);new_maps.append(mapping);i=rw.AddAtom(atom)
        rw.AddBond(previous,i,Chem.BondType.SINGLE);previous=i
    rw.AddBond(previous,dummy,Chem.BondType.SINGLE)
    changed=rw.GetMol();Chem.SanitizeMol(changed)
    fragments=original.copy();fragments[index]=Chem.MolToSmiles(changed,isomericSmiles=True)
    molecule=join_fragments(fragments);require(canonical(molecule)!=canonical(parent_mol),'CANDIDATE_UNCHANGED')
    require(molecule.GetNumHeavyAtoms()==parent_mol.GetNumHeavyAtoms()+extension,'CANDIDATE_ATOM_COUNT')
    info=identity(molecule);prototype_id='C03-'+digest(info['canonical_isomeric_smiles'].encode())[:8]
    # Protect all non-linker atoms and all original bonds, except the chosen boundary bond.
    old_atoms={a.GetAtomMapNum():a for a in parent_mol.GetAtoms()};new_atoms={a.GetAtomMapNum():a for a in molecule.GetAtoms()}
    for mapping,atom in old_atoms.items():
        other=new_atoms[mapping]
        require((atom.GetAtomicNum(),atom.GetFormalCharge(),atom.GetChiralTag())==
                (other.GetAtomicNum(),other.GetFormalCharge(),other.GetChiralTag()),'CANDIDATE_PROTECTED_ATOM_CHANGED')
    boundary_map=linker.GetAtomWithIdx(neighbor).GetAtomMapNum()
    other_boundary=[]
    for i,smiles in enumerate(original):
        if i==index:continue
        for atom in Chem.MolFromSmiles(smiles).GetAtoms():
            if atom.GetAtomicNum()==0 and atom.GetAtomMapNum()==dummy_map:
                other_boundary.append(atom.GetNeighbors()[0].GetAtomMapNum())
    require(len(other_boundary)==1,'CANDIDATE_BOUNDARY_PARTNER')
    removed=frozenset([boundary_map,other_boundary[0]])
    def bonds(mol):
        return {frozenset([b.GetBeginAtom().GetAtomMapNum(),b.GetEndAtom().GetAtomMapNum()]):
                (b.GetBondTypeAsDouble(),str(b.GetStereo())) for b in mol.GetBonds()}
    old_bonds=bonds(parent_mol);new_bonds=bonds(molecule)
    require(removed in old_bonds and removed not in new_bonds,'CANDIDATE_BOUNDARY_REPLACEMENT')
    require(all(new_bonds.get(key)==value for key,value in old_bonds.items() if key!=removed),'CANDIDATE_PROTECTED_BOND_CHANGED')
    require(len(new_bonds)==len(old_bonds)+extension,'CANDIDATE_BOND_COUNT')
    check_active();molecule.SetProp('_Name',prototype_id+' / design hypothesis / not synthesized')
    rdDepictor.Compute2DCoords(molecule)
    sdf=put(prototype_id+'.sdf',(Chem.MolToMolBlock(molecule)+'\n$$$$\n').encode(),'chemical/x-mdl-sdfile')
    smi=put(prototype_id+'.smi',(info['canonical_isomeric_smiles']+'\n').encode(),'text/plain')
    picture=Chem.Mol(molecule)
    for atom in picture.GetAtoms():atom.SetAtomMapNum(0)
    image=Draw.MolToImage(picture,size=(1000,520),highlightAtoms=[a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomMapNum() in new_maps])
    stream=io.BytesIO();image.save(stream,format='PNG');png=put(prototype_id+'.png',stream.getvalue(),'image/png')
    return {'format':'tpd-design-hypothesis/0.1.0','candidate_id':prototype_id,'status':'hypothesis_pending_review',
        'design_class':'linker_perturbation_prototype','final_candidate':False,
        'parent_candidate':parent,'source_parts_ref':part_ref,'parent_identity':identity(parent_mol),'identity':info,
        'transformation':{'kind':'linker_boundary_methylene_extension','attachment_dummy_map':dummy_map,'added_methylene_count':extension,
            'new_atom_maps':new_maps,'protected':'Original mapped atoms, formal charges and stereochemical tags retained; warhead and recruiter fragments unchanged.'},
        'rationale':'Explore a small linker-length perturbation while preserving the supplied binding fragments. This is a geometric design hypothesis, not evidence of improved efficacy.',
        'fragments':[{'role':p['role'],'mapped_smiles':fragments[i]} for i,p in enumerate(parts)],
        'files':{'sdf':sdf,'smiles':smi,'image':png},'coordinates_kind':'2d_depiction_not_bound_pose',
        'chemical_validity':{'rdkit_sanitized':True,'one_connected_molecule':True,'rdkit_version':rdkit.__version__},
        'synthesis':{'status':'requires_route_review','verified_route':None,'note':'Graph cuts/joining are not synthetic reagents or a reaction feasibility proof.'},
        'remaining':['변경 linker 도입 반응과 경로 근거 검토','수소/화학 상태와 결합 부위 보존 검토','실제 구조 예측·실험 검증'],
        'authority':AUTHORITY.copy()}


def design_linker_candidate(saved,result_id,parameters,put,check_active):
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import rdDepictor,Draw
    from packages.science.linker_design import replace_linker,library
    from packages.science.molecules import identity
    check_active()
    parent=parameters.get('parent','C01')
    part_ref=saved.files(result_id)[BASE+'cpu/'+parent+'-parts.json']
    parts=saved.port.json(part_ref)['parts']
    design=replace_linker(parts,parameters)
    molecule=design['molecule'];info=design['identity']
    identifier='C03-'+digest(info['canonical_isomeric_smiles'].encode())[:8]
    molecule.SetProp('_Name',identifier+' / linker design hypothesis / not a validated candidate')
    rdDepictor.Compute2DCoords(molecule)
    picture=Chem.Mol(molecule)
    for atom in picture.GetAtoms():atom.SetAtomMapNum(0)
    img=Draw.MolToImage(picture,size=(1000,520),highlightAtoms=[a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomMapNum() in design['new_maps']])
    stream=io.BytesIO();img.save(stream,format='PNG')
    check_active()
    library_ref=put('linker-library.json',encoded(library()),'application/json')
    sdf=put(identifier+'.sdf',(Chem.MolToMolBlock(molecule)+'\n$$$$\n').encode(),'chemical/x-mdl-sdfile')
    smi=put(identifier+'.smi',(info['canonical_isomeric_smiles']+'\n').encode(),'text/plain')
    png=put(identifier+'.png',stream.getvalue(),'image/png')
    return {'format':'tpd-design-hypothesis/0.2.0','candidate_id':identifier,'status':'hypothesis_pending_review',
        'design_class':'medchem_linker_hypothesis','final_candidate':False,
        'parent_candidate':parent,'source_parts_ref':part_ref,'parent_identity':identity(design['parent']),'identity':info,
        'transformation':{'kind':'whole_linker_replacement','template_id':parameters['linker_id'],
            'label':design['template']['label'],'family':design['template']['family'],
            'orientation':parameters.get('orientation','forward'),'new_atom_maps':design['new_maps'],
            'removed_linker_atom_maps':design['removed_maps'],'protected_atom_maps':design['protected_maps'],
            'protected':'Warhead/recruiter atoms, internal bonds, charge, hydrogens and mapped stereocenters checked; reference attachment sites retained.'},
        'rationale':design['template']['rationale'],'linker_comparison':design['comparison'],
        'attachments':design['attachments'],'chemical_state_review':design['template']['state_review'],
        'design_evidence':{'library_ref':library_ref,'library_version':design['library_version'],
            'template_origin':'developer_defined_motif_not_exact_paper_compound','sources':design['sources'],
            'expert_reply_sha256':library()['review_source_sha256'],'novelty':'not_assessed'},
        'fragments':[{'role':p['role'],'mapped_smiles':design['fragments'][i]} for i,p in enumerate(parts)],
        'files':{'sdf':sdf,'smiles':smi,'image':png},'coordinates_kind':'2d_depiction_not_bound_pose',
        'chemical_validity':{'rdkit_sanitized':True,'one_connected_molecule':True,'rdkit_version':rdkit.__version__,
            'protected_fragments_verified':True,'reference_attachment_sites_verified':True},
        'synthesis':{'status':'not_assessed','route_search_performed':False,'verified_route':None,
            'attachment_chemistry':design['attachments'],
            'parent_si_status':'route_not_curated',
            'review_checklist':['starting material','coupling and linker introduction','reaction compatibility and protecting groups',
                                'yield','purification','structure confirmation'],
            'note':'Graph attachment is not a reaction template match. No route search was run; this is not a finding of synthetic infeasibility.'},
        'promotion_checks':{'chemical_graph':'passed','protected_fragments':'passed',
            'attachment_geometry':'not_run','bound_structure_prediction':'not_run','full_complex_contacts':'not_run',
            'reaction_or_retrosynthesis_review':'not_assessed','expert_candidate_decision':'pending','efficacy':'not_established'},
        'remaining':['C01/C02 SI의 실제 합성 경로와 변경 linker 반응 적용성 확인',
                     '화학 상태와 exit-vector/부착 기하 검토','실제 GPU 결합 구조·복합체 접촉 평가',
                     '부모 대비 결합 부위 보존과 전문가 검토; 최종 후보 승격 전 별도 확인'],
        'authority':AUTHORITY.copy()}


def hydrogen_trial(saved,result_id,parameters,put,check_active):
    import rdkit
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem,Lipinski
    cid=parameters.get('candidate_id','C01');index=parameters.get('sample_index',0)
    data=saved.view(result_id);candidate=next((c for c in data['candidates'] if c['compound_id']==cid),None)
    require(candidate is not None and type(index) is int and 0<=index<len(candidate['samples']),'CHEMISTRY_SAMPLE')
    sample=candidate['samples'][index];preparation=saved.port.json(sample['reports']['preparation'])
    expected=preparation['prepared_structure']['sha256']
    refs=[ref for name,ref in sorted(saved.files(result_id).items()) if name.endswith('.sdf') and ref['sha256']==expected]
    require(bool(refs),'CHEMISTRY_PREPARED_SDF_MISSING');raw=saved.port.read(refs[0])
    molecules=list(Chem.ForwardSDMolSupplier(io.BytesIO(raw),removeHs=False));require(len(molecules)==1 and molecules[0] is not None,'CHEMISTRY_SDF')
    mol=molecules[0];require(mol.GetNumConformers()==1,'CHEMISTRY_CONFORMER')
    heavy=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()>1];before=np.array(mol.GetConformer().GetPositions())
    require(AllChem.MMFFHasAllMoleculeParams(mol),'CHEMISTRY_MMFF_UNSUPPORTED')
    properties=AllChem.MMFFGetMoleculeProperties(mol,mmffVariant='MMFF94s')
    ff=AllChem.MMFFGetMoleculeForceField(mol,properties)
    for i in heavy:ff.AddFixedPoint(i)
    ff.Initialize();initial=ff.CalcEnergy();converged=False
    for iteration in range(10):
        check_active()
        if ff.Minimize(maxIts=100)==0:converged=True;break
    final=ff.CalcEnergy();after=np.array(mol.GetConformer().GetPositions())
    displacement=float(np.linalg.norm(after[heavy]-before[heavy],axis=1).max())
    require(displacement<1e-7,'CHEMISTRY_HEAVY_ATOMS_MOVED')
    source_ref=put('ligand-before.sdf',raw,'chemical/x-mdl-sdfile')
    prepared_ref=put('ligand-hydrogen-trial.sdf',(Chem.MolToMolBlock(mol)+'\n$$$$\n').encode(),'chemical/x-mdl-sdfile')
    acceptors={i for match in Lipinski._HAcceptors(mol) for i in match};donors={i for match in Lipinski._HDonors(mol) for i in match}
    atoms=[{'atom_map':a.GetAtomMapNum(),'index':a.GetIdx(),'element':a.GetSymbol(),'formal_charge':a.GetFormalCharge(),
        'attached_hydrogens':sum(n.GetAtomicNum()==1 for n in a.GetNeighbors()),'rdkit_donor':a.GetIdx() in donors,
        'rdkit_acceptor':a.GetIdx() in acceptors,
        'bonds':[{'neighbor_map':b.GetOtherAtom(a).GetAtomMapNum(),'order':b.GetBondTypeAsDouble()} for b in a.GetBonds()]}
        for a in mol.GetAtoms() if a.GetAtomicNum()>1 and a.GetSymbol()=='N']
    return {'format':'tpd-hydrogen-trial/0.1.0','status':'converged_with_limits' if converged else 'not_converged',
        'candidate_id':cid,'sample_run_id':sample['run_id'],'binding':sample['binding'],'source_ref':refs[0],
        'before_ref':source_ref,'prepared_ref':prepared_ref,'original_preparation_ref':sample['reports']['preparation'],
        'method':'RDKit MMFF94s, ligand-only, every heavy atom fixed; at most 1000 minimization iterations',
        'runtime':{'rdkit':rdkit.__version__},'heavy_atom_max_displacement_A':displacement,
        'ligand_forcefield_energy_kcal_mol':{'before':initial,'after':final,'meaning':'Local optimizer diagnostic, not binding energy, linker strain or efficacy.'},
        'nitrogen_inventory':atoms,'contact_recalculation':'not_run',
        'limitations':['입력 protonation/tautomer를 유지한 국소 시험이며 실험 pH 상태를 예측하지 않음',
                      '단백질·용매·결합 환경을 힘장에 포함하지 않음','전체 복합체 수소 방향/Probe 접촉 재계산과 별도',
                      '효능/합성/최종 화학 상태 판정에 사용하지 않음'], 'authority':AUTHORITY.copy()}
