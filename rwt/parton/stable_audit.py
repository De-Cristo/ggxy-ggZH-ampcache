"""Complete frozen real-input numerical gate; native shard regression on workers."""
import json,math,sys,tarfile,time,subprocess
from pathlib import Path
from common import worker,sha,arrays
worker()
import numpy as np
import torch
sys.path.insert(0,str(Path('rwt/transformer/vendor').resolve()))
from ParticleTransformer import to_ptrapphim
from stable_pair import rapidity,pairwise,CONVENTION

def scan(d,n):
    counts=dict(rows=n,active_tokens=0,legacy_nonfinite_events=0,corrected_nonfinite_pairs=0,max_finite_legacy_rapidity_difference=0.)
    for start in range(0,n,256):
        f=torch.from_numpy(np.array(d['pf_features'][start:start+256]));v=torch.from_numpy(np.array(d['pf_vectors'][start:start+256]));m=torch.from_numpy(np.array(d['pf_mask'][start:start+256])).bool()
        assert torch.isfinite(f).all() and torch.isfinite(v).all()
        y=rapidity(f,v,m);old=to_ptrapphim(v,False,eps=None)[:,1:2]
        bad=(~torch.isfinite(old))&m;counts['legacy_nonfinite_events']+=int(bad.any(dim=2).sum())
        good=m&torch.isfinite(old)
        if good.any():counts['max_finite_legacy_rapidity_difference']=max(counts['max_finite_legacy_rapidity_difference'],float((old-y)[good].abs().max()))
        counts['active_tokens']+=int(m.sum())
        a=torch.cat((v,y),1);p=pairwise(a.unsqueeze(-1),a.unsqueeze(-2));active=m.unsqueeze(-1)&m.unsqueeze(-2)
        counts['corrected_nonfinite_pairs']+=int(((~torch.isfinite(p))&active).sum())
        if start%200192==0:print('STABLE_SCAN',start,n,flush=True)
    assert counts['corrected_nonfinite_pairs']==0
    return counts

def native():
    import pyarrow.parquet as pq
    assert sha('native.parquet')=='432908ea35d450d73fc5d98699af914c5329330edaf421f188215ce6e569eba5'
    z=np.load('native-inputs.npz',allow_pickle=False)
    assert len(z['seed'])==5000 and (z['seed']==311108).all()
    objects=pq.read_table('native-objects.parquet').to_pylist()
    ys=rapidity(torch.from_numpy(z['pf_features']),torch.from_numpy(z['pf_vectors']),torch.from_numpy(z['pf_mask']).bool()).numpy()[:,0]
    count=0;maximum=0.;example=None
    for batch in pq.ParquetFile('native.parquet').iter_batches(batch_size=16):
        for row in batch.to_pylist():
            i=row['table_row'];ev=row['event'];obj=objects[i]
            assert row['seed']==311108 and row['lhe_ordinal']==int(z['lhe_ordinal'][i])==obj['lhe_ordinal']
            kids=[j for j in ev[row['z_index']]['daughters'] if abs(ev[j]['id']) in (11,12,13,14,15,16)]
            kids.sort(key=lambda j:-math.hypot(*ev[j]['p4'][:2]));kids=kids[:2]
            base=[ev[row['z_index']]['p4'],ev[row['h_index']]['p4']]+[ev[j]['p4'] for j in kids]
            base += [None]*(4-len(base))
            radiation=[ev[j]['p4'] for j in obj['parton_indices']];radiation.sort(key=lambda v:-math.hypot(*v[:2]))
            tokens=base+radiation
            assert len(tokens)<=z['pf_vectors'].shape[2]
            for k,v in enumerate(tokens):
                if v is None:continue
                px,py,pz,e=v;pt=math.hypot(px,py);mass=math.sqrt(max(0,e*e-px*px-py*py-pz*pz))
                y=math.asinh(pz/math.hypot(pt,mass));diff=abs(y-float(ys[i,k]));maximum=max(maximum,diff);count+=1
                if i==1354 and k==3:
                    example=dict(seed=311108,lhe_ordinal=row['lhe_ordinal'],table_row=i,token=k,native_p4=v,native_mass=mass,native_rapidity=y,corrected_rapidity=float(ys[i,k]),absolute_difference=diff)
    assert example and count>0
    report=dict(rows=5000,tokens=count,max_absolute_rapidity_difference=maximum,tolerance=2e-5,example=example,pass_=maximum<=2e-5)
    Path('native-rapidity-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    assert report['pass_'],report
    # Link the native counterexample to the frozen prepared event, exactly.
    d=arrays('lo-cms');ix=np.flatnonzero((d['seed']==311108)&(d['lhe_ordinal']==1354));assert len(ix)==1
    length=d['pf_vectors'].shape[2]
    for key in ('pf_features','pf_vectors','pf_mask'):
        assert np.array_equal(d[key][ix[0]],z[key][1354,:,:length])
    report['prepared_identity_exact']=True
    return report

def main():
    channel=sys.argv[1];torch.set_num_threads(1);start=time.time()
    with tarfile.open('parton-prepared.tgz') as t:
        assert all((m.name=='parton-data' or m.name.startswith('parton-data/')) and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers());t.extractall('.')
    subprocess.run(['sha256sum','-c','--quiet','SHA256SUMS'],cwd='parton-data',check=True)
    selection=sha(f'production/parton-{channel}-1m-v1.json');assert json.loads(Path('parton-data/audit.json').read_text())['selection_sha256']==selection
    out=dict(channel=channel,package=Path('VERSION').read_text().strip(),selection_sha256=selection,prepared_sha256=sha('parton-prepared.tgz'),rapidity_convention=CONVENTION,samples={},status='RUNNING')
    try:
        for order in ('lo-cms','nlo'):
            d=arrays(order);assert len(d['weight'])==1000000
            out['samples'][order]=scan(d,1000000)
            for k in ('weight','globals'):assert np.isfinite(d[k]).all()
        if channel=='zvv':out['native']=native()
        out['status']='PASS'
    finally:
        out['elapsed_seconds']=time.time()-start
        Path('stable-audit.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out),flush=True)
if __name__=='__main__':main()
