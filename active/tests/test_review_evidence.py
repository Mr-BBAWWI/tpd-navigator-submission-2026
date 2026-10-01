"""Binding, missing/failed stages and archive integrity; no scientific execution."""
import copy
from pathlib import Path
import tempfile
import unittest

from packages.science.evidence_common import file_index, seal
from packages.science.handoff import ACTIVE, encoded, sha
from packages.science.ligand_preparation import _report
from packages.science.ligand_contacts import project_contacts
from packages.science.quality_evidence import export_quality_evidence, project
from packages.science.review_evidence import export_review_evidence, markdown, project_review, read_review_evidence
from packages.science.structure_quality import project_request, strict_json
from test_ligand_contacts import fixture as contact_fixture
from test_structure_quality import fixture as quality_fixture, synthetic_evidence


def fixture():
    request, raw = quality_fixture(missing=True)
    quality_report = project_request(request, raw)
    evidence, files = synthetic_evidence(request)
    for candidate in evidence['candidates']:
        candidate.update(literature_observations=[{'kind': 'missing', 'missing_reason': 'synthetic unknown'}],
                         reference_reconstruction={}, computed_properties={'kind': 'computed', 'values': {'MW': 31.06}},
                         review_items=['synthetic pending review'])
    evidence.update(starting_ligand={}, shared_hypothesis={},
                    literature_reconciliation={'status': 'not_provided', 'producer': None, 'records': [], 'meaning': 'synthetic'},
                    authority={'human_review': 'pending', 'approval_record_created': False, 'dispatch_authorized': False,
                               'public_release_ready': False, 'efficacy_claim': 'not_established'})
    quality = project(evidence, files, [('quality/q0000', quality_report)])
    base = {'evidence/' + p: d for p, d in files.items()}
    base.update({'manifest.json': b'synthetic base manifest', 'quality/q0000/report.json': encoded(quality_report)})
    binding = request['binding']
    # Explicitly synthetic inventory to exercise the report contract, not chemical correctness.
    atoms = [{'sdf_atom_index': i, 'element': e, 'atom_map': m, 'source_or_generated_name': n,
              'formal_charge': 0, 'parent_heavy_atom_map': parent, 'xyz_A': [float(i), 0., 0.]}
             for i, (e, m, n, parent) in enumerate([('N', 1, 'N1', None), ('C', 2, 'C2', None), ('H', 0, 'H1', 1)])]
    chemistry = {'canonical_isomeric_smiles': 'CN', 'formula': 'CH5N', 'formal_charge': 0, 'heavy_atom_count': 2,
                 'expected_hydrogen_count': 1, 'generated_hydrogen_count': 1, 'hydrogens_by_heavy_atom_map': {'1': 1, '2': 0},
                 'max_serialized_heavy_displacement_A': 0., 'hydrogen_bond_length_range_A': [1., 1.], 'atom_inventory': atoms}
    prep_request = {'data_mode': 'synthetic_test', 'binding': binding,
                    'atom_mapping': [{'atom_map': 1, 'source_atom_name': 'N1'}, {'atom_map': 2, 'source_atom_name': 'C2'}]}
    prep = _report(prep_request, {'status': 'prepared_for_review', 'error': None, 'chemistry': chemistry,
                                 'runtime': {'rdkit': '2026.03.6', 'gemmi': '0.7.5'}}, {'prepared/ligand.sdf': b'synthetic'})
    context, contact_files, _ = contact_fixture()
    context.update(binding=copy.deepcopy(binding), ligand_preparation_digest=prep['digest'], prior_quality_digest=quality_report['digest'])
    contact = project_contacts(context, contact_files)
    return evidence, quality, base, [('preparation/0000', prep, None)], [('contacts/0000', contact, context)], contact_files


class ReviewProjectionTests(unittest.TestCase):
    def setUp(self):
        self.evidence, self.quality, self.files, self.preps, self.contacts, self.contact_files = fixture()

    def result(self):
        return project_review(self.evidence, self.quality, self.files, self.preps, self.contacts)

    def sample(self):
        return self.result()['candidates'][0]['samples'][0]

    def test_combines_all_stages_without_replacing_original_facts(self):
        result = self.result(); row = result['candidates'][0]['samples'][0]
        self.assertEqual(result['candidates'][0]['evidence'], self.evidence['candidates'][0])
        self.assertEqual(row['prior_quality'], self.quality['candidates'][0]['samples'][0])
        self.assertEqual(row['ligand_preparation']['report']['contact_quality_status'], 'not_assessed')
        self.assertEqual(row['ligand_contacts']['status'], 'completed_with_limits')
        self.assertEqual(row['review_status'], 'review_required')
        self.assertFalse(result['authority']['approval_record_created'])
        self.assertEqual(result['summary']['contacts']['completed_with_limits'], 1)

    def test_missing_stages_are_not_failed_or_zero_clashes(self):
        self.preps = []; self.contacts = []
        row = self.sample()
        self.assertEqual(row['ligand_contacts'], {'status': 'not_provided', 'report_ref': None, 'report': None})
        self.assertIn('LIGAND_PREPARATION_NOT_PROVIDED', [i['code'] for i in row['review_issues']])

    def test_rejected_preparation_keeps_reason(self):
        prep = self.preps[0][1]
        prep.update(status='not_prepared', error='LIGAND_UNSPECIFIED_STEREOCHEMISTRY', chemistry=None, prepared_structure=None, runtime={})
        self.contacts = []
        row = self.sample()
        self.assertEqual(row['ligand_preparation']['status'], 'not_prepared')
        self.assertIn(prep['error'], [i['detail'] for i in row['review_issues']])

    def failed_directions(self, paths):
        for path in paths:
            record = strict_json(self.contact_files[path]); record.update(exit_code=124, error='CONTACT_TOOL_TIMEOUT')
            self.contact_files[path] = encoded(record)
        path, _, context = self.contacts[0]
        self.contacts = [(path, project_contacts(context, self.contact_files), context)]

    def test_all_failed_contacts_keep_null_count(self):
        self.failed_directions([n for n in self.contact_files if n.endswith('execution.json')])
        row = self.sample()
        self.assertEqual(row['ligand_contacts']['status'], 'failed_or_partial')
        self.assertTrue(all(m['bad_overlap_pair_count'] is None for m in row['ligand_contacts']['report']['models']))
        self.assertEqual(sum(i['code'] == 'CONTACT_DIRECTION_FAILED' for i in row['review_issues']), 4)
        self.assertIn('판정 불가', markdown(self.result()).decode())

    def test_partial_contact_keeps_successful_findings_and_failed_direction(self):
        self.failed_directions(['runs/xray/reverse/execution.json'])
        row = self.sample(); model = row['ligand_contacts']['report']['models'][0]
        self.assertEqual(model['bad_overlap_pair_count'], 1)
        self.assertFalse(model['both_directions_completed'])
        issue = next(i for i in row['review_issues'] if i['code'] == 'CONTACT_DIRECTION_FAILED')
        self.assertEqual((issue['mode'], issue['direction']), ('xray', 'reverse'))

    def test_role_findings_deduplicate_directions_but_keep_modes(self):
        for model in self.contacts[0][1]['models']:
            for d in model['directions'].values():
                d['typing']['role_disagreements'] = [{'name': 'N1', 'atom_map': 1, 'rdkit_acceptor': False,
                    'probe_acceptor': True, 'rdkit_donor': True, 'probe_attached_h_donor': True, 'restraint_h_bond_type': 'D'}]
        issues = [i for i in self.sample()['review_issues'] if i['code'] == 'CONTACT_ROLE_DISAGREEMENT']
        self.assertEqual(len(issues), 2)
        self.assertTrue(all(i['detail'].endswith(': 1') for i in issues))

    def test_all_binding_fields_checked(self):
        original = copy.deepcopy(self.preps[0][1]['binding'])
        for field in original:
            with self.subTest(field=field):
                self.preps[0][1]['binding'] = {**original, field: 1 if field == 'model_rank' else 'wrong'}
                with self.assertRaisesRegex(ValueError, 'REVIEW_SAMPLE_BINDING|REVIEW_ORPHAN|REQUIRES_PREPARATION'):
                    self.result()
        self.preps[0][1]['binding'] = original

    def test_real_synthetic_mixing_rejected(self):
        self.preps[0][1]['data_mode'] = 'real'
        with self.assertRaisesRegex(ValueError, 'REVIEW_SAMPLE_ORIGIN'):
            self.result()

    def test_stale_preparation_digest_rejected(self):
        self.contacts[0][2]['ligand_preparation_digest'] = 'f' * 64
        with self.assertRaisesRegex(ValueError, 'REVIEW_STALE_PREPARATION'):
            self.result()

    def test_stale_quality_digest_rejected(self):
        self.contacts[0][1]['prior_quality_digest'] = 'f' * 64
        with self.assertRaisesRegex(ValueError, 'REVIEW_STALE_QUALITY'):
            self.result()

    def test_duplicate_reports_rejected(self):
        for target in (self.preps, self.contacts):
            target.append(target[0])
            with self.assertRaisesRegex(ValueError, 'REVIEW_DUPLICATE'):
                self.result()
            target.pop()

    def test_orphan_report_rejected(self):
        self.contacts[0][1]['binding'] = {**self.contacts[0][1]['binding'], 'run_id': 'unrelated'}
        with self.assertRaisesRegex(ValueError, 'REVIEW_ORPHAN_REPORT'):
            self.result()

    def test_contacts_need_prepared_ligand(self):
        self.preps = []
        with self.assertRaisesRegex(ValueError, 'REVIEW_CONTACT_REQUIRES_PREPARATION'):
            self.result()

    def test_contacts_need_linked_prior_quality(self):
        self.quality['candidates'][0]['samples'][0]['quality_report'] = None
        with self.assertRaisesRegex(ValueError, 'REVIEW_CONTACT_REQUIRES_PRIOR_QUALITY'):
            self.result()

    def test_unparsed_structure_does_not_accept_preparation(self):
        collection = strict_json(self.files['evidence/collection.json'])
        collection['runs'][0]['models'][0]['status'] = 'missing'
        self.files['evidence/collection.json'] = encoded(collection)
        with self.assertRaisesRegex(ValueError, 'REVIEW_UNPARSED_SAMPLE'):
            self.result()

    def test_missing_and_invalid_predictions_preserved_without_stage_reports(self):
        self.preps = []; self.contacts = []
        for state in ('missing', 'invalid'):
            collection = strict_json(self.files['evidence/collection.json'])
            collection['runs'][0]['models'][0]['status'] = state
            self.files['evidence/collection.json'] = encoded(collection)
            quality_row = self.quality['candidates'][0]['samples'][0]
            quality_row.update(structure_status=state, quality_report=None, quality_summary=None, analyses=[], quality_status='not_assessed')
            result = self.result(); row = result['candidates'][0]['samples'][0]
            self.assertIsNone(row['binding'])
            self.assertEqual(result['summary']['structures'][state], 1)
            self.assertIn('STRUCTURE_' + state.upper(), [i['code'] for i in row['review_issues']])

    def test_no_predictions_and_synthesis_remain_unassessed(self):
        candidate = self.result()['candidates'][1]
        self.assertEqual(candidate['samples'], [])
        self.assertEqual(candidate['synthesis_assessment']['status'], 'not_assessed')
        self.assertEqual(candidate['evidence']['literature_observations'][0]['kind'], 'missing')
        self.assertIn('제공된 예측 없음', markdown(self.result()).decode())

    def test_markdown_escapes_untrusted_review_text(self):
        self.evidence['candidates'][0]['review_items'] = ['<script>bad</script>|new\nrow']
        text = markdown(self.result()).decode()
        self.assertNotIn('<script>', text)
        self.assertIn('&lt;script&gt;bad&lt;/script&gt;\\|new row', text)


class ReviewArchiveTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name) / 'quality'; self.out = Path(self.tmp.name) / 'review'
        export_quality_evidence(ACTIVE / 'outputs/b_evidence_20260923', [], self.base)
        self.result = export_review_evidence(self.base, [], [], self.out)

    def reseal(self):
        manifest = strict_json((self.out / 'manifest.json').read_bytes())
        manifest['files'] = file_index({p.relative_to(self.out).as_posix(): p.read_bytes() for p in self.out.rglob('*') if p.is_file() and p != self.out / 'manifest.json'})
        (self.out / 'manifest.json').write_bytes(encoded(seal({k: v for k, v in manifest.items() if k != 'digest'})))

    def test_complete_reader_without_science_imports(self):
        import sys
        before = set(sys.modules)
        self.assertEqual(read_review_evidence(self.out), self.result)
        self.assertFalse({'rdkit', 'gemmi', 'torch'} & (set(sys.modules) - before))
        self.assertEqual(self.result['summary']['sample_count'], 0)

    def test_output_exists_does_not_overwrite(self):
        with self.assertRaisesRegex(ValueError, 'OUTPUT_EXISTS'):
            export_review_evidence(self.base, [], [], self.out)

    def test_changed_source_bytes_rejected(self):
        (self.out / 'source/quality/evidence/README.md').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'REVIEW_FILE_HASH'):
            read_review_evidence(self.out)

    def test_resealed_markdown_cannot_claim_approval(self):
        (self.out / 'README.md').write_bytes(b'approved'); self.reseal()
        with self.assertRaisesRegex(ValueError, 'REVIEW_MARKDOWN_PROJECTION'):
            read_review_evidence(self.out)

    def test_resealed_report_cannot_change_result(self):
        data = strict_json((self.out / 'review-evidence.json').read_bytes())
        data['authority']['approval_record_created'] = True
        (self.out / 'review-evidence.json').write_bytes(encoded(data)); self.reseal()
        with self.assertRaisesRegex(ValueError, 'REVIEW_PROJECTION'):
            read_review_evidence(self.out)

    def test_indexed_unexpected_file_rejected(self):
        (self.out / 'unexpected.json').write_bytes(b'{}'); self.reseal()
        with self.assertRaisesRegex(ValueError, 'REVIEW_UNEXPECTED_FILES'):
            read_review_evidence(self.out)

    def test_unsafe_source_directory_rejected(self):
        manifest = strict_json((self.out / 'manifest.json').read_bytes())
        manifest['preparation_directories'] = ['../outside']
        (self.out / 'manifest.json').write_bytes(encoded(seal({k: v for k, v in manifest.items() if k != 'digest'})))
        with self.assertRaisesRegex(ValueError, 'REVIEW_SOURCE_DIRECTORIES'):
            read_review_evidence(self.out)
