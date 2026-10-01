"""E3-specific reference intake and explicit paired-coordinate comparison.

This adapter never labels a self-comparison as model calibration.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from .structures import aligned_rmsd, atom_sites

SOURCE = Path(__file__).resolve().parents[2]/'cases/design_sources'

def reference_packet():
    meta=json.loads((SOURCE/'crbn_benchmark.json').read_text())
    groups={role:atom_sites(SOURCE/'6BOY.cif',[meta['reference_'+role+'_chain']])
            for role in ('target','e3','ligand')}
    meta['reference_summary']={role:{'heavy_atom_count':len(atoms),
        'residues':len({(a['label_comp_id'],a['auth_seq_id']) for a in atoms})} for role,atoms in groups.items()}
    meta['reference_sha256']=hashlib.sha256((SOURCE/'6BOY.cif').read_bytes()).hexdigest()
    meta['missing_for_calibration']=['independent model prediction','model/input/MSA/seed provenance','multiple-seed comparison']
    return meta

def compare_paired(reference, prediction, provenance):
    """Each role has coordinates, in the same exact mapped order on both sides.

    Validated sequence/atom pairing is required from the producer; mapping IDs
    are checked here. No sequence-offset guessing or ligand-only fitting.
    """
    required_strings=('model','input_sha256','prediction_sha256','e3_type')
    if not isinstance(provenance,dict) or any(
        not isinstance(provenance.get(k),str) or not provenance[k].strip()
        for k in required_strings
    ):
        raise ValueError('BENCHMARK_PROVENANCE_REQUIRED')
    seed=provenance.get('seed')
    if isinstance(seed,bool) or not isinstance(seed,int) or seed<0:
        raise ValueError('BENCHMARK_PROVENANCE_REQUIRED')
    if provenance['e3_type'] not in {'CRBN','VHL'}:raise ValueError('BENCHMARK_E3')
    arrays={}
    for role in ('target','e3','ligand'):
        a,b=reference[role],prediction[role]
        if a['ids']!=b['ids'] or len(set(a['ids']))!=len(a['ids']):raise ValueError('BENCHMARK_ATOM_PAIRING')
        x,y=np.asarray(a['xyz'],float),np.asarray(b['xyz'],float)
        if x.shape!=y.shape or x.shape!=(len(a['ids']),3) or not np.isfinite(x).all() or not np.isfinite(y).all() or len(x)<3:
            raise ValueError('BENCHMARK_COORDINATES')
        arrays[role]=(x,y)
    fit=aligned_rmsd(*arrays['target']);rotation=np.array(fit['rotation_row_vectors']);shift=np.array(fit['translation_A'])
    metrics={role+'_RMSD_after_target_alignment_A':float(np.sqrt(np.mean(np.sum((mov@rotation+shift-ref)**2,axis=1))))
             for role,(ref,mov) in arrays.items()}
    identical=all(np.allclose(a,b,atol=1e-9) for a,b in arrays.values())
    return {'e3_type':provenance['e3_type'],'provenance':provenance,'metrics':metrics,
        'evaluation_kind':'self_comparison_not_calibration' if identical or provenance.get('origin')=='self_test' else 'prediction_reference_comparison',
        'calibration_complete':False,'cross_e3_ranking_allowed':False,'seed_distribution':'requires_additional_independent_runs',
        'method':'One target Kabsch alignment applied to every part; exact producer-validated mapping',
        'limitations':['Paired atom coverage and model provenance require producer review.','Single case/seed does not establish generalization or E3 superiority.']}

def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--paired',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    result=reference_packet() if a.paired is None else compare_paired(**json.loads(a.paired.read_text(encoding='utf-8')))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)

if __name__=='__main__':main()
