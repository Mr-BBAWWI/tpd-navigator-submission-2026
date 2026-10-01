"""R5 graph-level regressions on the registered two first-case inputs, not efficacy tests."""
import io
import json
import unittest
from pathlib import Path
from rdkit import Chem
from packages.contracts import ContractError
from packages.science.linker_design import library,replace_linker,validate_parameters,design_class
from packages.science.molecules import canonical,join_fragments
from packages.platform.cpu_tools import design_candidate
from packages.platform.expert_responses import extract_answers
import test_cpu_tools as cpu_fixture
import test_workbench as workbench_fixture

ROOT=Path(__file__).resolve().parents[1]


class LinkerDesignTests(unittest.TestCase):
    def setUp(self):
        self.parts={cid:json.loads((ROOT/'outputs/smarca2_20260923'/f'{cid}-parts.json').read_text())['parts'] for cid in ('C01','C02')}

    def test_all_templates_both_parents_and_orientations_preserve_binding_fragments(self):
        for cid,parts in self.parts.items():
            before=json.dumps(parts,sort_keys=True)
            for template in library()['templates']:
                for orientation in ('forward','reverse'):
                    with self.subTest(parent=cid,template=template['id'],orientation=orientation):
                        result=replace_linker(parts,{'parent':cid,'linker_id':template['id'],'orientation':orientation})
                        mol=result['molecule']
                        roundtrip=Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
                        self.assertIsNotNone(roundtrip)
                        self.assertEqual(canonical(roundtrip),canonical(mol))
                        self.assertEqual(len(Chem.GetMolFrags(mol)),1)
                        self.assertFalse(any(a.GetAtomicNum()==0 for a in mol.GetAtoms()))
                        self.assertEqual(mol.GetNumHeavyAtoms(),len(result['protected_maps'])+len(result['new_maps']))
                        self.assertFalse(set(result['removed_maps']) & {a.GetAtomMapNum() for a in mol.GetAtoms()})
                        for i,p in enumerate(parts):
                            if p['role']!='linker':self.assertEqual(p['mapped_smiles'],result['fragments'][i])
                        parent_stereo={result['parent'].GetAtomWithIdx(i).GetAtomMapNum():s for i,s in Chem.FindMolChiralCenters(result['parent'])}
                        stereo={mol.GetAtomWithIdx(i).GetAtomMapNum():s for i,s in Chem.FindMolChiralCenters(mol)}
                        self.assertEqual(parent_stereo,stereo)
            self.assertEqual(before,json.dumps(parts,sort_keys=True))

    def test_asymmetric_orientation_changes_structure_and_symmetric_deduplicates(self):
        def design(t,o):return canonical(replace_linker(self.parts['C01'],{'linker_id':t,'orientation':o})['molecule'])
        self.assertNotEqual(design('peg_alkyl','forward'),design('peg_alkyl','reverse'))
        self.assertNotEqual(design('triazole','forward'),design('triazole','reverse'))
        self.assertEqual(design('alkyl_c6','forward'),design('alkyl_c6','reverse'))

    def test_unsupported_or_mixed_inputs_rejected(self):
        for p in ({'linker_id':'x'},{'linker_id':[]},{'parent':[]},{'linker_id':'peg3','orientation':True},
                  {'linker_id':'peg3','linker_extension':1},{'orientation':'reverse'},{'linker_extension':True},
                  {'linker_id':'peg3','smiles':'arbitrary'}):
            with self.subTest(p=p),self.assertRaises(ContractError):validate_parameters(p)

    def test_artifacts_have_source_state_and_promotion_gaps(self):
        fixture=cpu_fixture.CandidateGraphTests();fixture.setUp();self.addCleanup(fixture.doCleanups)
        result=design_candidate(fixture.saved,'synthetic',{'parent':'C02','linker_id':'piperazine'},
                                lambda n,b,m:fixture.port.put_raw(b,m),lambda:None)
        self.assertEqual(design_class(result),'medchem_linker_hypothesis')
        self.assertFalse(result['final_candidate'])
        self.assertFalse(result['authority']['human_approved'])
        self.assertFalse(result['synthesis']['route_search_performed'])
        self.assertEqual(result['promotion_checks']['bound_structure_prediction'],'not_run')
        self.assertEqual(fixture.port.json(result['design_evidence']['library_ref'])['version'],library()['version'])
        mol=next(Chem.ForwardSDMolSupplier(io.BytesIO(fixture.port.read(result['files']['sdf']))))
        self.assertEqual(canonical(mol),result['identity']['canonical_isomeric_smiles'])

    def test_legacy_identity_retained_and_classified(self):
        fixture=cpu_fixture.CandidateGraphTests();fixture.setUp();self.addCleanup(fixture.doCleanups)
        for cid,n,expected in [('C01',1,'C03-a6951643'),('C02',2,'C03-9a437858')]:
            r=design_candidate(fixture.saved,'synthetic',{'parent':cid,'linker_extension':n},lambda n,b,m:fixture.port.put_raw(b,m),lambda:None)
            self.assertEqual(r['candidate_id'],expected)
            self.assertEqual(design_class(r),'linker_perturbation_prototype')
            self.assertFalse(r['final_candidate'])

    def test_r1_r5_doc_stops_before_future_feature_requests(self):
        headings=['2. R1 — 원자형·수소결합 역할 차이','3. R2 — Protonation·tautomer·수소 방향',
                  '4. R3 — RMSD·Boltz confidence·ternary 구조 비교','5. R4 — 문헌 관측값·조건·결측 표현',
                  '6. R5 — 합성 근거와 신규 C03']
        raw=workbench_fixture.DocumentTests().doc([v for i,h in enumerate(headings) for v in (h,'synthetic '+str(i))]+
            ['7. 다음 버전에 추가·수정할 핵심 기능','not part of the R5 opinion'])
        answers=extract_answers(raw)
        self.assertEqual(len(answers),5)
        self.assertEqual(answers['Q05'],'synthetic 4')
