"""Prepare real-event scale companions on workers, retaining exact PDF provenance.

The Born ratio is valid only for gg Born events with on-shell quark masses.
It does not reweight shower kernels or reconstruct missing NLO scale weights.
"""
import concurrent.futures
import gzip
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import time
import xml.etree.ElementTree as ET

from common import worker, sha, arrays
worker()
import lhapdf
import numpy as np

POINTS = [(1.,1.),(1.,2.),(1.,.5),(2.,1.),(2.,2.),(2.,.5),(.5,1.),(.5,2.),(.5,.5)]
IDS = [str(i) for i in range(1001,1010)]
PDFS = {325300:'NNPDF31_nnlo_as_0118_mc_hessian_pdfas',93300:'PDF4LHC21_40_pdfas'}
REFERENCE_SHA = '0bfccf9dc9e9aab2d9c5f373dc5226f2ef3d54815639ac2317056ca26821f6d6'

def num(s):
    return float(s.replace('d','e').replace('D','E'))

def header(path):
    lines=[]
    with gzip.open(path,'rt') as f:
        for line in f:
            if re.match(r'\s*<event(?:\s|>)',line):break
            lines.append(line)
            assert len(lines)<20000,'Unexpected oversized LHE header'
    text=''.join(lines)
    defs={m[0]:m[1].strip() for m in re.findall(r'<weight\s+id=[\'\"]([^\'\"]+)[\'\"]\s*>(.*?)</weight>',text,re.S)}
    card={}
    for line in text.splitlines():
        m=re.match(r'^\s*([A-Za-z][A-Za-z0-9_]*)\s+([-+0-9.eEdD]+)(?:\s|$)',line)
        if m:
            try: card[m[1].lower()]=num(m[2])
            except ValueError: pass
    return card,defs

def events(path):
    block=None
    with gzip.open(path,'rt') as f:
        for line in f:
            if re.match(r'\s*<event(?:\s|>)',line):
                assert block is None;block=[]
            if block is None:continue
            end=line.find('</event>');block.append(line if end<0 else line[:end+8])
            if end<0:continue
            node=ET.fromstring(''.join(block));block=None
            rows=[s.split() for s in node.text.splitlines() if s.strip() and not s.lstrip().startswith('#')]
            h=rows[0];n=int(h[0]);assert len(rows)>=n+1
            pp=[dict(id=int(r[0]),status=int(r[1]),v=[num(v) for v in r[6:10]]) for r in rows[1:n+1]]
            yield dict(weight=num(h[2]),alpha=num(h[5]),particles=pp,
                       weights={w.attrib['id']:num(w.text) for w in node.findall('./rwgt/wgt')})
    assert block is None,'Incomplete event'

def check_card(card,order,pdfid,fixed,channel):
    required=dict(lhans1=pdfid,lhans2=pdfid,mtscheme=0,alphas_from_pdf=1,
                  ebeam1=6800,ebeam2=6800,renscfact=1,facscfact=1,
                  zmass=91.1876,hmass=125,topmass=172.5,vdecaymode=11 if channel=='zll' else 12)
    if order=='lo-cms':required.update(bornonly=1,loevents=1)
    else:required.update(bornonly=0,loevents=0)
    if fixed:required['fixedscale']=1
    else:assert card.get('fixedscale',0)!=1 and card.get('dynscale',0)!=3,card
    for k,v in required.items():assert k in card and math.isclose(card[k],v,rel_tol=1e-12),('LHE card mismatch',k,v,card)
    return {k:card[k] for k in required}

def born_ratio(event,base_pdf,pdf,points,fixed):
    incoming=[p for p in event['particles'] if p['status']==-1]
    assert len(incoming)==2 and all(p['id']==21 for p in incoming),'Non-gg Born record'
    outgoing=[p for p in event['particles'] if p['status']==1]
    assert len(outgoing)==3 and sorted(abs(p['id']) for p in outgoing) in ([11,11,25],[13,13,25],[15,15,25],[12,12,25],[14,14,25],[16,16,25]),'Not a ZH Born event'
    for p in incoming:
        assert abs(p['v'][0])+abs(p['v'][1])<1e-8 and math.isclose(abs(p['v'][2]),p['v'][3],rel_tol=1e-8)
    assert incoming[0]['v'][2]*incoming[1]['v'][2]<0
    x1,x2=[p['v'][3]/6800. for p in incoming]
    assert 0<x1<=1 and 0<x2<=1
    if fixed:mu=216.1876
    else:
        v=np.sum([p['v'] for p in incoming],axis=0)
        mu=math.sqrt(v[3]**2-np.dot(v[:3],v[:3]))/2
    a0=base_pdf.alphasQ(mu)
    den=base_pdf.xfxQ(21,x1,mu)*base_pdf.xfxQ(21,x2,mu)
    assert den>0 and math.isfinite(den)
    # Powers and flux/phase-space cancel at fixed kinematics. x*f ratios equal f ratios.
    values=[]
    for r,f in points:
        value=(pdf.alphasQ(r*mu)/a0)**2*pdf.xfxQ(21,x1,f*mu)*pdf.xfxQ(21,x2,f*mu)/den
        assert math.isfinite(value) and value>0
        values.append(value)
    return values

def download(product,path,expected):
    path.parent.mkdir(parents=True,exist_ok=True)
    root=os.environ['PAPER_EOS_ROOT'].rstrip('/')
    for suffix in ('.sha256',''):
        subprocess.run(['xrdcp','--retry','3','--cksum','adler32',root+'/'+product+suffix,str(path)+suffix],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    assert Path(str(path)+'.sha256').read_text().split()[0]==expected
    assert sha(path)==expected,('Input SHA256 mismatch',product)
    return path

def control(pdfs):
    path=Path('reference.lhe.gz');assert sha(path)==REFERENCE_SHA
    card,defs=header(path);check_card(card,'lo-cms',93300,False,'zll')
    worst=0.;n=0
    for e in events(path):
        predicted=born_ratio(e,pdfs[93300],pdfs[93300],POINTS,False)
        predicted+=born_ratio(e,pdfs[93300],pdfs[325300],[(1.,1.)],False)
        exact=np.array([e['weights'][i]/e['weight'] for i in IDS+['2000']])
        predicted=np.array(predicted)
        worst=max(worst,float(np.max(np.abs(predicted-exact)/np.abs(exact))))
        np.testing.assert_allclose(predicted,exact,rtol=1e-4,atol=1e-8)
        n+=1
    assert n==5000
    return dict(events=n,ratios=10*n,max_relative_error=worst,rtol=1e-4,atol=1e-8,
                lhe_sha256=REFERENCE_SHA,card=card,weight_definitions=defs,passed=True)

def unpack_prepared():
    with tarfile.open('parton-prepared.tgz') as t:
        assert all((m.name=='parton-data' or m.name.startswith('parton-data/')) and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers())
        t.extractall('.')
    subprocess.run(['sha256sum','-c','--quiet','SHA256SUMS'],cwd='parton-data',check=True)

def main():
    start=time.time();ch=sys.argv[1];assert ch in ('zll','zvv')
    selection=Path(f'production/parton-{ch}-1m-v1.json')
    m=json.loads(selection.read_text());unpack_prepared()
    oldaudit=json.loads(Path('parton-data/audit.json').read_text());assert oldaudit['selection_sha256']==sha(selection)
    expected_bundle={'zll':'cb07b6b232f147af8ba01cd38bdca82169b085dced41872d4e942e3b4c447e0e',
                     'zvv':'6fc277b7ed8d9788870a29398e618f7fc75dbac31ab9d2de521661be85679111'}[ch]
    assert sha('parton-prepared.tgz')==expected_bundle
    pdfs={i:lhapdf.mkPDF(name,0) for i,name in PDFS.items()}
    out=Path('scale-inputs');assert not out.exists();out.mkdir()
    audit=dict(schema='ggzh-scale-companion-v1',channel=ch,package=Path('VERSION').read_text().strip(),
               selection_sha256=sha(selection),prepared_sha256=expected_bundle,lhapdf=lhapdf.version(),
               scale_ids=IDS,scale_points=POINTS,pdfs=PDFS,samples={},control=control(pdfs),
               limitations=['Fixed shower; no ISR/FSR/tune variations',
                 'Exact NNPDF NLO scale grid absent; no factorized replacement',
                 'Scale-transfer models and real NN resampling ensemble still required'])
    print('REAL_BORN_CONTROL_PASS',json.dumps(audit['control'])[:500],flush=True)
    for order in ('lo-cms','nlo'):
        original=arrays(order);n=len(original['weight']);assert n==1000000
        source=order=='lo-cms';pdfid=325300 if source else 93300
        chunks=sorted(m['samples'][order],key=lambda e:e['seed'])
        data={k:np.asarray(original[k]).copy() for k in ('seed','lhe_ordinal','table_row','source_fold','target_half','in_window')}
        data['nominal_lhe']=np.empty(n,np.float64)
        if source:
            data['scale_nnpdf']=np.empty((n,9),np.float64)
            data['scale_pdf4lhc']=np.empty((n,9),np.float64)
        else:
            data['scale_pdf4lhc']=np.empty((n,9),np.float64)
            data['nominal_nnpdf']=np.empty(n,np.float64)
        suffix='rwt10m-13600' if ch=='zll' else 'vv-rwt10m-13600'
        # Product spelling follows production/paper.sh; never inferred from event arrays.
        sample=f'{order}-{suffix}'
        offset=0;records=[]
        def fetch(entry):
            product=f'lhe/{sample}/events-{entry["seed"]}.lhe.gz'
            return entry,download(product,Path('lhe')/str(entry['seed'])/'events.lhe.gz',entry['lhe_sha256'])
        # At most four active copies; event calculations are single-threaded LHAPDF.
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            futures={i:pool.submit(fetch,e) for i,e in enumerate(chunks[:4])}
            for j in range(len(chunks)):
                entry,path=futures.pop(j).result()
                if j+4<len(chunks):futures[j+4]=pool.submit(fetch,chunks[j+4])
                rows=entry['events'];sl=slice(offset,offset+rows);seed=entry['seed']
                assert np.array_equal(original['seed'][sl],np.full(rows,seed))
                assert np.array_equal(original['table_row'][sl],np.arange(rows))
                ordinal=np.asarray(original['lhe_ordinal'][sl]);assert np.array_equal(np.sort(ordinal),np.arange(rows))
                card,defs=header(path);checked=check_card(card,order,pdfid,source,ch)
                if not source:
                    for k,(r,f) in zip(IDS,POINTS):
                        d=defs[k];assert 'lhapdf=93300' in d
                        assert num(re.search(r'renscfact=(\S+)',d)[1])==r and num(re.search(r'facscfact=(\S+)',d)[1])==f
                    assert 'lhapdf=325300' in defs['2000']
                nominal=[];grid=[];other=[]
                for e in events(path):
                    nominal.append(e['weight'])
                    if source:
                        assert not e['weights'],'Unexpected LO weight block; review schema'
                        grid.append(np.array(born_ratio(e,pdfs[325300],pdfs[325300],POINTS,True))*e['weight'])
                        other.append(np.array(born_ratio(e,pdfs[325300],pdfs[93300],POINTS,True))*e['weight'])
                    else:
                        grid.append([e['weights'][k] for k in IDS]);other.append(e['weights']['2000'])
                assert len(nominal)==rows
                nominal=np.array(nominal)[ordinal];grid=np.array(grid)[ordinal];other=np.array(other)[ordinal]
                np.testing.assert_allclose(nominal if source else other,original['weight'][sl],rtol=1e-6,atol=1e-10)
                if source:np.testing.assert_allclose(grid[:,0],nominal,rtol=1e-12,atol=1e-14)
                for v in (nominal,grid,other):assert np.isfinite(v).all()
                data['nominal_lhe'][sl]=nominal
                if source:
                    data['scale_nnpdf'][sl]=grid;data['scale_pdf4lhc'][sl]=other
                else:
                    data['scale_pdf4lhc'][sl]=grid;data['nominal_nnpdf'][sl]=other
                records.append(dict(seed=seed,rows=rows,start=offset,lhe_sha256=entry['lhe_sha256'],card=checked,
                                    definitions=defs,nominal_sumw=float(nominal.sum())))
                offset+=rows;path.unlink();Path(str(path)+'.sha256').unlink()
                if (j+1)%20==0:print('SCALE_SHARDS',ch,order,j+1,len(chunks),flush=True)
        assert offset==n
        summaries={k:dict(sumw=np.sum(v,axis=0).tolist(),sumw2=np.sum(v*v,axis=0).tolist())
                   for k,v in data.items() if k.startswith(('nominal','scale'))}
        np.savez_compressed(out/(order+'.npz'),**data,scale_points=np.array(POINTS),scale_ids=np.array(IDS))
        audit['samples'][order]=dict(rows=n,joined=True,in_window=int(data['in_window'].sum()),shards=records,
                                     sums=summaries,negative_nominal=int((data['nominal_lhe']<0).sum()))
    audit['elapsed_seconds']=time.time()-start
    audit['status']='INPUT_COMPANIONS_COMPLETE_NOT_UNCERTAINTY_BANDS'
    Path('scale-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    (out/'README.md').write_text('''# Scale companion inputs\n\nJoin to the unchanged prepared bundle by (seed,lhe_ordinal); table_row is retained.\nAll weight arrays are absolute LHE weights in pb. Cross sections use sum(weights)/N_generated;\nno per-variation unit normalization has been applied. Signed target weights are retained.\nThe source scale_nnpdf grid is the genuine Born hard-scale variation at fixed shower.\nThe target scale_pdf4lhc grid is the exact stored PDF4LHC21 grid; nominal_nnpdf is separate.\nDo not label that grid as exact NNPDF3.1 variations. scale_points gives (muR/mu0,muF/mu0).\nThis bundle contains no scale-transfer model or NN uncertainty ensemble. Existing folds\nare cross-validation fits, not resampling replicas. See audit.json for gates/provenance.\n''')
    paths=sorted(out.iterdir());(out/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in paths))
    with tarfile.open('scale-inputs.tgz','w:gz',compresslevel=1) as t:t.add(out,arcname='scale-inputs')
    print(audit['status'],ch,flush=True)

if __name__=='__main__':main()
