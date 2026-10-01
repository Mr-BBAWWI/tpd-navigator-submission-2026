"""Export saved source excerpts for follow-up, without chemical reinterpretation."""
import argparse
import json
import shlex
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from packages.contracts import encoded
from packages.platform.store import Store
from packages.platform.saved_results import SavedResultsService
from packages.platform.dossiers import digest, require


def atom_rows(text, names):
    # Read the simple atom loop in these source restraint dictionaries only.
    # Semicolon/multiline CIF fields outside this loop are not interpreted.
    headers=[];rows=[];inside=False
    for line in text.splitlines():
        if line.startswith('_chem_comp_atom.'):
            inside=True;headers.append(line.strip());continue
        if inside and line.strip() and not line.startswith('#'):
            if line.startswith(('loop_','data_','_',';')):break
            values=shlex.split(line)
            require(len(values)==len(headers),'FOLLOWUP_DICTIONARY_LAYOUT')
            row=dict(zip(headers,values))
            if row['_chem_comp_atom.atom_id'] in names:rows.append(row)
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True);p.add_argument('--result-id',required=True)
    p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
    require(not a.output_dir.exists(),'FOLLOWUP_OUTPUT_EXISTS');a.output_dir.mkdir(parents=True)
    svc=SavedResultsService(Store(a.data_dir),'local-research');v=svc.view(a.result_id);files=svc.files(a.result_id)
    source_bundle=svc.port.json(files['review-evidence.json']);rows=[];sources={};atom_findings=[]
    for candidate in v['candidates']:
        for sample in candidate['samples']:
            rows.append({k:sample[k] for k in ('run_id','rank','binding','condition_digest','reported_execution','summary','parts_rmsd_ranges_A','mapping_count','alignment_method','bad_overlap_pairs','reports')})
        sample=candidate['samples'][0]
        raw_candidate=next(c for c in source_bundle['candidates'] if c['compound_id']==candidate['compound_id'])
        raw_sample=next(s for s in raw_candidate['samples'] if (s['run_id'],s['rank'])==(sample['run_id'],sample['rank']))
        stage=raw_sample['ligand_preparation'];prep=stage['report'];maps={r['atom_map'] for r in sample['role_disagreements']}
        names={r['name'] for r in sample['role_disagreements']}
        prep_dir=str(Path(stage['report_ref']['path']).parent).replace('\\','/')
        contact_dir=str(Path(raw_sample['ligand_contacts']['report_ref']['path']).parent).replace('\\','/')
        dictionary=contact_dir+'/input/raw/restraints.cif'
        selected={label:path for label,path in [('preparation.json',stage['report_ref']['path']),('contacts.json',raw_sample['ligand_contacts']['report_ref']['path']),('comparison.json',next(name for name,ref in files.items() if ref==sample['reports']['comparison'])),('input-ligand.sdf',prep_dir+'/raw/ligand.sdf'),('prepared-ligand.sdf',prep_dir+'/prepared/ligand.sdf'),('restraints.cif',dictionary)]}
        for label,path in selected.items():
            ref=files[path];raw=svc.port.read(ref);name=candidate['compound_id']+'-'+label
            (a.output_dir/name).write_bytes(raw);sources[name]={'archive_path':path,'ref':ref,'sha256':digest(raw)}
        atom_findings.append({'compound_id':candidate['compound_id'],'run_id':sample['run_id'],
            'reported_atoms':[atom for atom in prep['chemistry']['atom_inventory'] if atom['atom_map'] in maps],
            'reported_hydrogen_counts':{str(m):prep['chemistry']['hydrogens_by_heavy_atom_map'].get(str(m)) for m in maps},
            'dictionary_atom_rows':atom_rows(svc.port.read(files[dictionary]).decode('utf-8'),names),
            'role_disagreements':sample['role_disagreements'],
            'interpretation':'Literal saved source values only; no pKa prediction, chemical-state correction, or atom-role adjudication.'})
    result={'bundle_digest':v['bundle_digest'],'samples':rows,'atom_findings':atom_findings,'source_files':sources,
            'expert_response':v['expert_reviews'],'authority':v['authority']}
    (a.output_dir/'후속검토_근거.json').write_bytes(encoded(result))
    print(json.dumps({'samples':len(rows),'files':len(sources),'atom_candidates':len(atom_findings),'output':str(a.output_dir)},ensure_ascii=False))


if __name__=='__main__':main()
