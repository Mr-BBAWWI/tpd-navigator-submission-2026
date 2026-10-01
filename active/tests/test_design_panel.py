"""Scientific boundary and execution regression tests; fixtures are not SAR evidence."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from rdkit import Chem
from fastapi.testclient import TestClient
from apps.api.main import create_app,PROJECT
from packages.platform.store import Store
from packages.platform.workbench import WorkbenchService
from packages.platform.design_panel import _stereoisomers
from packages.science.analog_generation import generate,select_diverse,RULE_CATALOG
from packages.science.dual_e3 import SOURCE,catalog,assemble,validate,verify_sources
from packages.science.structures import atom_sites
from packages.science.warhead_sites import analyze_sites,pose_preservation,funnel
from packages.science.e3_benchmark import compare_paired,reference_packet
from packages.science.molecules import read_ccd,canonical
from packages.contracts import ContractError

class ScienceDesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog=catalog();cls.meta=cls.catalog['warhead']
        cls.parent=Chem.MolFromMolBlock((SOURCE/'SMARCA2-neutral-design.sdf').read_text())
        cls.protein=atom_sites(SOURCE/'6HAZ.cif',['A'])
        cls.sites=analyze_sites(cls.parent,cls.protein,cls.meta['sar'],cls.meta['protected_maps'])

    def test_source_hashes_and_crbn_halogen_stereo(self):
        self.assertGreater(len(verify_sources()),8)
        molecule=read_ccd(SOURCE/'RN6.cif')
        self.assertTrue(any(a.GetSymbol()=='Cl' for a in molecule.GetAtoms()))
        self.assertEqual(reference_packet()['prediction_status'],'not_run')
        self.assertFalse(reference_packet().get('calibration_complete',False))

    def test_known_exit_vector_recovered_missing_sar_does_not_pass(self):
        atoms={a['atom_map']:a for a in self.sites['atoms']}
        self.assertEqual(atoms[19]['state'],'MODIFIABLE')
        self.assertEqual(atoms[20]['state'],'PROTECTED')
        absent=analyze_sites(self.parent,self.protein,{},self.meta['protected_maps'])
        self.assertEqual(next(a['state'] for a in absent['atoms'] if a['atom_map']==19),'UNKNOWN')
        self.assertEqual(sum(a['state']=='MODIFIABLE' for a in self.sites['atoms']),1)

    def test_missing_scope_and_wrong_maps_fail_closed(self):
        sites=copy.deepcopy(self.sites['atoms'])
        next(a for a in sites if a['atom_map']==19)['evidence']['SAR']['allowed_transformations']=[]
        self.assertEqual(generate(self.parent,sites,False)['analogs'],[])
        sites[0]['atom_map']=999
        with self.assertRaisesRegex(ValueError,'COVERAGE'):generate(self.parent,sites,True)
        with self.assertRaisesRegex(ValueError,'HEAVY'):generate(Chem.AddHs(self.parent),self.sites['atoms'],True)

    def test_real_graph_classes_and_parent_preservation(self):
        before=Chem.MolToMolBlock(self.parent)
        result=generate(self.parent,self.sites['atoms'],True)
        self.assertGreaterEqual(len(result['analogs']),20)
        self.assertGreaterEqual(len({a['transformation_class'] for a in result['analogs']}),5)
        protected={a['atom_map'] for a in self.sites['atoms'] if a['state']=='PROTECTED'}
        for a in result['analogs']:
            self.assertTrue(a['protected_graph_preserved'])
            self.assertFalse(protected & set(a['modified_atom_maps']+a['removed_atom_maps']))
            m=Chem.MolFromSmiles(a['mapped_smiles'])
            self.assertEqual(len(Chem.GetMolFrags(m)),1)
            self.assertEqual(a['pharmacophore_preserved'],'review')
        self.assertTrue(any(a['reason']=='PROTECTED_TOUCHED' for a in result['rejections']))
        self.assertEqual(Chem.MolToMolBlock(self.parent),before)

    def test_strict_generation_never_changes_unknown(self):
        result=generate(self.parent,self.sites['atoms'],False)
        self.assertGreater(len(result['analogs']),0)
        for a in result['analogs']:self.assertTrue(all(s['state']=='MODIFIABLE' for s in a['site_status']))

    def test_stereo_enumeration_and_cluster_caps(self):
        raw=generate(self.parent,self.sites['atoms'],True)['analogs'];records=_stereoisomers(raw)
        self.assertGreater(len(records),len(raw))
        for r in records:
            clean=Chem.MolFromSmiles(r['canonical_smiles'])
            unresolved=any(str(s.specified)=='Unspecified' for s in Chem.FindPotentialStereo(clean))
            self.assertEqual(unresolved,r['stereochemistry']['status']=='requires_review')
        selected=select_diverse(records,16)
        from collections import Counter
        self.assertTrue(all(n<=3 for n in Counter(r['cluster_id'] for r in selected).values()))
        self.assertGreater(len({r['transformation_class'] for r in selected}),1)

    def test_bioisostere_and_exit_relocation_are_real_edits(self):
        for smi,required in [('[CH3:1][C:2](=[O:3])[CH2:4][c:5]1[cH:6][c:7]([OH:11])[cH:8][cH:9][cH:10]1',
                             {'bioisosteric_replacement','exit_vector_relocation','hbond_rewiring','heteroatom_swap'})]:
            mol=Chem.MolFromSmiles(smi)
            sites=[{'atom_map':a.GetAtomMapNum(),'state':'UNKNOWN','evidence':{'kind':'synthetic_test_fixture'}} for a in mol.GetAtoms()]
            result=generate(mol,sites,True)
            self.assertTrue(required<={a['transformation_class'] for a in result['analogs']})

    def test_dual_e3_assembly_is_different_chemistry_with_preserved_core(self):
        generated=generate(self.parent,self.sites['atoms'],False)['analogs']

        def attachment_options(record):
            mol=Chem.MolFromSmiles(record['mapped_smiles'])
            added=set(record['added_atom_maps'])
            return [atom.GetAtomMapNum() for atom in mol.GetAtoms()
                    if atom.GetAtomMapNum() in added and atom.GetAtomicNum() in {7,8}
                    and atom.GetTotalNumHs(includeNeighbors=True)>0]

        eligible=[record for record in generated if attachment_options(record)]
        no_handle=[record for record in generated if not attachment_options(record)]
        self.assertTrue(eligible,'generation must provide an explicit added NH/OH attachment option')
        self.assertTrue(no_handle,'fixture must exercise rejection when no added NH/OH handle exists')

        a=_stereoisomers([eligible[0]])[0]
        a.update(parent_warhead=self.meta['id'],protected_atom_maps=self.meta['protected_maps'])
        linker=self.catalog['linkers']['templates'][0]
        recruiters=self.catalog['recruiters']['recruiters']
        results=[assemble(a,rec,linker) for rec in recruiters]
        self.assertEqual({r['e3_type'] for r in results},{'VHL','CRBN'})
        self.assertNotEqual(results[0]['canonical_smiles'],results[1]['canonical_smiles'])
        for r in results:
            self.assertFalse(r['final_candidate']);self.assertEqual(r['ternary_structure_status'],'not_run')
            mol=Chem.MolFromSmiles(r['mapped_smiles'])
            self.assertFalse(any(a.GetAtomicNum()==0 for a in mol.GetAtoms()))
            block=Chem.MolToMolBlock(mol,forceV3000=True)
            restored=Chem.MolFromMolBlock(block)
            self.assertEqual(canonical(restored),canonical(mol))
            self.assertEqual({a.GetAtomMapNum() for a in restored.GetAtoms()},{a.GetAtomMapNum() for a in mol.GetAtoms()})

        invalid=_stereoisomers([no_handle[0]])[0]
        invalid.update(parent_warhead=self.meta['id'],protected_atom_maps=self.meta['protected_maps'])
        with self.assertRaises(ValueError):
            assemble(invalid,recruiters[0],linker)

    def test_forbidden_attachment_and_bad_params(self):
        a={'mapped_smiles':Chem.MolToSmiles(self.parent),'added_atom_maps':[],'protected_atom_maps':[19]}
        with self.assertRaises(ValueError):assemble(a,self.catalog['recruiters']['recruiters'][0],self.catalog['linkers']['templates'][0])
        for value in [[],{'panel_size':True},{'exploratory':1},{'target':'BRD4'},{'linker_ids':['bogus']},{'execute':'shell'}]:
            with self.assertRaises((ValueError,TypeError)):validate(value)

    def test_pose_translation_does_not_get_ligand_aligned_away(self):
        moved=Chem.Mol(self.parent);conf=moved.GetConformer()
        for i,p in enumerate(conf.GetPositions()):conf.SetAtomPosition(i,(p+np.array([20,0,0])).tolist())
        result=pose_preservation(self.parent,[{'mol':moved,'score':-1000}],self.protein,self.meta['protected_maps'])
        self.assertIs(result['docking_pose_preserved'],False)
        self.assertFalse(result['all_poses_diagnostics'][0]['alignment_applied'])

    def test_protected_element_mismatch_never_passes(self):
        changed=Chem.RWMol(self.parent);changed.GetAtomWithIdx(0).SetAtomicNum(7)
        out=pose_preservation(self.parent,[changed.GetMol()],self.protein,self.meta['protected_maps'])
        self.assertIs(out['docking_pose_preserved'],False)

    def test_funnel_does_not_compare_kd_ki_or_pad_missing_candidates(self):
        r=funnel([{'evidence_tier':'known_protac','state':'known','Ki':0.1},
                  {'evidence_tier':'co_crystal','state':'known','Kd':9999}],limit=10)
        self.assertEqual(len(r),2);self.assertEqual(r[0]['Kd'],9999)

    def test_benchmark_shared_alignment_not_self_calibration(self):
        xyz=[[0,0,0],[1,0,0],[0,1,0]]
        ref={k:{'ids':['1','2','3'],'xyz':xyz} for k in ['target','e3','ligand']}
        provenance={'model':'synthetic-test','seed':23,'input_sha256':'a'*64,'prediction_sha256':'b'*64,'e3_type':'CRBN'}
        self.assertEqual(compare_paired(ref,ref,provenance)['evaluation_kind'],'self_comparison_not_calibration')
        pred=copy.deepcopy(ref);pred['e3']['xyz']=(np.array(xyz)+[10,0,0]).tolist()
        r=compare_paired(ref,pred,provenance)
        self.assertAlmostEqual(r['metrics']['e3_RMSD_after_target_alignment_A'],10)
        self.assertFalse(r['cross_e3_ranking_allowed'])
        pred['e3']['ids'][0]='other'
        with self.assertRaises(ValueError):compare_paired(ref,pred,provenance)

class DesignJobTests(unittest.TestCase):
    def test_http_design_operation_can_be_submitted(self):
        with tempfile.TemporaryDirectory() as folder:
            app=create_app(data_root=folder,enable_worker=False)
            with TestClient(app) as client:
                response=client.post('/api/lab/jobs',headers={'X-TPD-Local':'1'},json={
                    'result_id':'design:SMARCA2','operation':'design_panel',
                    'request_key':'design-http-regression','parameters':{'dock':False}})
                self.assertEqual(response.status_code,202,response.text)
                job=client.get('/api/lab/jobs/'+response.json()['job_id']).json()
                self.assertEqual(job['operation'],'design_panel')
                self.assertEqual(job['state'],'queued')

    def test_local_job_admission_cancel_api_required_and_read_only(self):
        with tempfile.TemporaryDirectory() as folder:
            svc=WorkbenchService(Store(Path(folder)),PROJECT)
            row,_=svc.create('design:SMARCA2','design_panel','test',{'dock':False})
            self.assertEqual(svc.view(row['id'])['state'],'queued')
            again,fresh=svc.create('design:SMARCA2','design_panel','test',{'dock':False})
            self.assertFalse(fresh);self.assertEqual(again['id'],row['id'])
            svc.cancel(row['id']);svc.execute(row['id']);self.assertEqual(svc.get(row['id'])['state'],'cancelled')
            with self.assertRaises(ContractError):svc.create('design:SMARCA2','design_panel','api',{'use_api':True})
            app=create_app(data_root=folder,enable_worker=False,read_only=True)
            with TestClient(app) as client:
                self.assertEqual(client.get('/api/design/catalog').status_code,200)
                self.assertEqual(client.get('/design').status_code,200)
                self.assertEqual(client.post('/api/lab/jobs',headers={'X-TPD-Local':'1'},json={'result_id':'design:SMARCA2','operation':'design_panel','request_key':'readonly','parameters':{}}).status_code,403)

if __name__=='__main__':unittest.main()
