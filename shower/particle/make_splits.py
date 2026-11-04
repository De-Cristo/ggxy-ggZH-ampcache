#!/usr/bin/env python3
"""Derive exact frozen Zll split prefixes from real cached arrays, on a worker."""
import hashlib,json,os,platform
from pathlib import Path
assert platform.node().split('.')[0]!='cmscon' and os.environ.get('_CONDOR_SCRATCH_DIR')
import numpy as np

def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 return h.hexdigest()

manifest=json.loads(Path('production/rwt10m-available.json').read_text())
out=dict(schema='ggzh-frozen-splits-v1',selection_sha256=manifest['selection_sha256'],
         package=Path('VERSION').read_text().strip(),samples={})
for order,label,expected in [('lo-cms','lo-cms',5000),('nlo','nlo-cp5',2500)]:
 path=Path(f'split-input/{label}.npz')
 with np.load(path,allow_pickle=False) as d:
  run=d['run'];event=d['event'];mass=d['lhe_mll']
 assert len(run)==manifest['expected_events'][label]
 assert np.all(run[1:]>=run[:-1])
 starts=np.r_[0,np.flatnonzero(run[1:]!=run[:-1])+1];ends=np.r_[starts[1:],len(run)]
 sample=dict(rows=len(run),in_window=0,input_sha256=digest(path),seeds={})
 for begin,end in zip(starts,ends):
  seed=int(run[begin]);ids=event[begin:end]
  assert end-begin==expected and np.array_equal(ids,np.arange(expected))
  inside=(mass[begin:end]>15)&(mass[begin:end]<150)
  sample['seeds'][str(seed)]=dict(rows=int(end-begin),row_offset=int(begin),
     window_offset=sample['in_window'],window_rows=int(inside.sum()),
     excluded_rows=ids[~inside].astype(int).tolist())
  sample['in_window']+=int(inside.sum())
 assert sorted(map(int,sample['seeds']))==manifest['included'][label]
 out['samples'][order]=sample
 print(order,'stored',sample['rows'],'window',sample['in_window'],'outside',sample['rows']-sample['in_window'],flush=True)
Path('split-map.json').write_text(json.dumps(out,indent=2)+'\n')
