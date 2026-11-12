"""Scalar-only bundle of a frozen scalar-table selection (ATLAS-style inputs), CPU worker.

Downloads every selected per-shard event table with its checksum, merges them with ggzh_table into the two scalar
samples that bdt_member.py and ensemble_assess.py read (parton-data/lo-cms.npz = LO+PS source, nlo-cp5.npz = NLO+PS
target; for ATLAS the target nominal is the stored PDF4LHC21 central weight), and writes selection.json, audit.json
and SHA256SUMS.  No particle tensors.
"""
import concurrent.futures,json,os,subprocess,sys,tarfile
from pathlib import Path
from common import worker,configuration,sha
worker()
import numpy as np
sys.path.insert(0,str(Path('rwt').resolve()))
import ggzh_table as Table
from ggzh_reweight_common import Sample

c=configuration();m=json.loads(Path(c['selection_path']).read_text());assert m['schema']=='ggzh-scalar-ensemble-selection-v1'
root=os.environ['PAPER_EOS_ROOT'].rstrip('/');out=Path('parton-data');assert not out.exists();out.mkdir()

def download(item):
    key,seed,f=item;local=Path('staged')/key/f'table-{seed}.csv.gz';local.parent.mkdir(parents=True,exist_ok=True)
    subprocess.run(['xrdcp','--retry','3','--cksum','adler32',root+'/'+f['product'],str(local)],check=True,stdout=subprocess.DEVNULL)
    assert local.stat().st_size==f['bytes'] and sha(local)==f['sha256'],(key,seed)
    return local
items=[(key,x['seed'],next(f for f in x['files'] if f['role']=='legacy_table')) for key,xs in m['samples'].items() for x in xs]
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:list(pool.map(download,items))
audit=dict(study=c['study'],selection_sha256=c['selection_sha256'],tables=len(items),samples={},package=Path('VERSION').read_text().strip(),
           sample_orders=m['sample_orders'],target_nominal='w:nominal')
for key,label in (('lo-cms','lo-cms'),('nlo','nlo-cp5')):
    xs=m['samples'][key];n=sum(x['events'] for x in xs);assert n==1000000
    Table.main(['--label',label,'--tables',f'staged/{key}/table-*.csv.gz','--nominal','nominal','--out',str(out/f'{label}.npz')])
    s=Sample(label,out/f'{label}.npz');assert s.n==n
    seeds=sorted(x['seed'] for x in xs);assert np.array_equal(np.unique(s['run']),seeds)
    for x in xs:assert np.isclose(s.w0[s['run']==x['seed']].sum(),x['sumw'],rtol=1e-6) or s.unweighted
    audit['samples'][key]=dict(rows=n,in_window=int(s.in_window.sum()),sumw=float(s.w0.sum()),shards=len(xs))
    print(key,audit['samples'][key],flush=True)
(out/'selection.json').write_text(json.dumps(m,indent=2)+'\n');(out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
paths=sorted(p for p in out.rglob('*') if p.is_file())
(out/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(out)}\n' for p in paths))
with tarfile.open('parton-prepared.tgz','w:gz',compresslevel=1) as t:t.add(out,arcname='parton-data')
Path('parton-prepare.json').write_text(json.dumps(audit,indent=2)+'\n')
print('PREPARED_SCALAR_1M_EACH',flush=True)
