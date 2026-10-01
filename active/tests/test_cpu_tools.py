"""Graph preservation tests on the committed first-case molecular inputs."""
import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from rdkit import Chem
from packages.platform.cpu_tools import design_candidate
from packages.platform.saved_results import BASE,AUTHORITY
from packages.platform.store import Store
from packages.science.molecules import join_fragments,canonical

ROOT=Path(__file__).resolve().parents[1]

class CandidateGraphTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.port=Store(self.temp.name).scope('synthetic-tests')
        self.parts={cid:json.loads((ROOT/'outputs/smarca2_20260923'/f'{cid}-parts.json').read_text()) for cid in ['C01','C02']}
        self.files={BASE+'cpu/'+cid+'-parts.json':self.port.put_json(p) for cid,p in self.parts.items()}
        self.saved=SimpleNamespace(port=self.port,files=lambda _:self.files)
    def test_all_four_bounded_variants_preserve_parent_graph_and_stereo(self):
        for cid in self.parts:
            original=join_fragments([p['mapped_smiles'] for p in self.parts[cid]['parts']]);before=Chem.MolToSmiles(original)
            for extension in (1,2):
                with self.subTest(parent=cid,extension=extension):
                    result=design_candidate(self.saved,'synthetic',{'parent':cid,'linker_extension':extension},
                        lambda n,b,m:self.port.put_raw(b,m),lambda:None)
                    mol=next(Chem.ForwardSDMolSupplier(io.BytesIO(self.port.read(result['files']['sdf']))))
                    self.assertEqual(mol.GetNumHeavyAtoms(),original.GetNumHeavyAtoms()+extension)
                    self.assertEqual(len(Chem.GetMolFrags(mol)),1)
                    self.assertEqual(Chem.FindMolChiralCenters(mol,includeUnassigned=True),Chem.FindMolChiralCenters(original,includeUnassigned=True))
                    self.assertNotEqual(canonical(mol),canonical(original))
                    self.assertEqual(result['authority'],AUTHORITY)
                    self.assertEqual(result['synthesis']['status'],'requires_route_review')
                    self.assertEqual(Chem.MolToSmiles(original),before)
                    self.assertTrue(self.port.read(result['files']['image']).startswith(b'\x89PNG'))
    def test_cancel_check_prevents_any_output(self):
        def cancel():raise RuntimeError('synthetic cancellation')
        with self.assertRaisesRegex(RuntimeError,'cancellation'):
            design_candidate(self.saved,'synthetic',{},lambda *a:self.fail('No output allowed'),cancel)
