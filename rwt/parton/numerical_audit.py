"""Scan real frozen ParT inputs without repairing records or modifying a fit."""
import json,sys,tarfile,time,warnings
from pathlib import Path
from common import worker,sha,arrays
worker()
import numpy as np
import torch
sys.path.insert(0,str(Path('rwt/transformer/vendor').resolve()))
from ParticleTransformer import pairwise_lv_fts_pp

def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,float) and not np.isfinite(value):return str(value)
    return value

def main():
    ch=sys.argv[1];assert ch in ('zll','zvv');start=time.time();torch.set_num_threads(1)
    with tarfile.open('parton-prepared.tgz') as t:
        assert all((m.name=='parton-data' or m.name.startswith('parton-data/')) and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers())
        t.extractall('.')
    import subprocess
    subprocess.run(['sha256sum','-c','--quiet','SHA256SUMS'],cwd='parton-data',check=True)
    selection=Path(f'production/parton-{ch}-1m-v1.json');m=json.loads(Path('parton-data/audit.json').read_text())
    assert m['selection_sha256']==sha(selection)
    out=dict(schema='ggzh-real-input-numerical-audit-v1',channel=ch,selection_sha256=sha(selection),
             package=Path('VERSION').read_text().strip(),prepared_sha256=sha('parton-prepared.tgz'),
             numpy=np.__version__,torch=str(torch.__version__),device='CPU',samples={},
             interpretation='Diagnosis only; no event/model changes. CPU checks do not certify GPU/optimizer stability.')
    examples={}
    for order in ('lo-cms','nlo'):
        d=arrays(order);n=len(d['weight']);assert n==1000000
        counts=dict(rows=n,active_tokens=0,nonfinite_raw=0,nonfinite_rapidity32_tokens=0,
                    nonfinite_rapidity32_events=0,nonfinite_rapidity64_tokens=0,
                    nonpositive_log_argument32_tokens=0,energy_below_abs_pz=0)
        ex=[]
        for begin in range(0,n,4096):
            end=min(begin+4096,n);v=np.array(d['pf_vectors'][begin:end]);mask=np.asarray(d['pf_mask'][begin:end,0],bool)
            counts['active_tokens']+=int(mask.sum())
            for k in ('pf_vectors','pf_features','globals','weight'):
                counts['nonfinite_raw']+=int((~np.isfinite(d[k][begin:end])).sum())
            e=v[:,3];z=v[:,2]
            with np.errstate(divide='ignore',invalid='ignore',over='ignore'):
                arg=np.float32(1)+(np.float32(2)*z)/np.maximum(e-z,np.float32(1e-20));rap=np.float32(.5)*np.log(arg)
                ee=e.astype(np.float64);zz=z.astype(np.float64)
                rap64=.5*np.log(1+2*zz/np.maximum(ee-zz,1e-20))
                stable=.5*(np.log(np.maximum(ee+zz,1e-20))-np.log(np.maximum(ee-zz,1e-20)))
            bad=mask&~np.isfinite(rap)
            counts['nonfinite_rapidity32_tokens']+=int(bad.sum())
            counts['nonfinite_rapidity32_events']+=int(bad.any(axis=1).sum())
            counts['nonfinite_rapidity64_tokens']+=int((mask&~np.isfinite(rap64)).sum())
            counts['nonpositive_log_argument32_tokens']+=int((mask&(arg<=0)).sum())
            counts['energy_below_abs_pz']+=int((mask&(e<np.abs(z))).sum())
            for local in np.flatnonzero(bad.any(axis=1)):
                if len(ex)>=20:break
                index=begin+int(local);vectors=np.array(d['pf_vectors'][index]);mm=np.array(d['pf_mask'][index,0],bool)
                x=torch.from_numpy(vectors[None,:,:]);pair=pairwise_lv_fts_pp(x.unsqueeze(-1),x.unsqueeze(-2))
                active_pairs=torch.from_numpy((mm[:,None]&mm[None,:])[None,None])
                nonfinite=(~torch.isfinite(pair))&active_pairs
                record=dict(index=index,seed=int(d['seed'][index]),lhe_ordinal=int(d['lhe_ordinal'][index]),
                    table_row=int(d['table_row'][index]),token_indices=np.flatnonzero(bad[local]).tolist(),
                    vectors=vectors.tolist(),mask=mm.tolist(),argument32=arg[local].tolist(),
                    rapidity32=rap[local].tolist(),rapidity64=rap64[local].tolist(),
                    stable_expression_only_not_a_fix=stable[local].tolist(),
                    actual_vendored_pair_nonfinite=int(nonfinite.sum()))
                ex.append(record)
            if begin%200704==0:print('NUMERICAL_SCAN',ch,order,begin,n,counts['nonfinite_rapidity32_events'],flush=True)
        out['samples'][order]=counts;examples[order]=ex
    out['elapsed_seconds']=time.time()-start
    out['status']='NONFINITE_DERIVED_INPUTS_FOUND' if any(x['nonfinite_rapidity32_events'] or x['nonfinite_raw'] for x in out['samples'].values()) else 'NO_NONFINITE_FOUND_IN_SCANNED_EXPRESSIONS'
    Path('numerical-audit.json').write_text(json.dumps(clean(out),indent=2,allow_nan=False)+'\n')
    Path('numerical-examples.json').write_text(json.dumps(clean(examples),indent=2,allow_nan=False)+'\n')
    print(json.dumps(clean(out)),flush=True)
if __name__=='__main__':main()
