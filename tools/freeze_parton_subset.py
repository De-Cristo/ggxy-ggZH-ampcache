"""Freeze a real particle-training subset using small manifests/validation only."""
import argparse
import hashlib
import json
from pathlib import Path

def digest(data):
    return hashlib.sha256(data).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    a=p.parse_args()
    assert not a.out.exists()
    raw=a.catalog.read_bytes();catalog=json.loads(raw)
    salt='ggzh-parton-zvv-1m-20261110-v1'
    producer='cff7a446f1044c8d8ec4a9533a35e50aa6dce112'
    selected={};maximum={};rejected=[]
    for order,count,rows in [('lo-cms',200,5000),('nlo',400,2500)]:
        candidates=[(key,x) for key,x in catalog['delivered'].items()
                    if x['channel']=='zvv' and x['order']==order
                    and x['primary_included'] and x['package']==producer]
        candidates.sort(key=lambda item:digest(f'{salt}|{order}|{item[1]["seed"]}'.encode()))
        selected[order]=[];maximum[order]=0
        for key,x in candidates:
            try:
                manifest=Path(x['manifest']);assert manifest.stat().st_size<100000
                b=manifest.read_bytes();assert digest(b)==x['manifest_sha256']
                m=json.loads(b);assert m['primary_included'] and m['events']==rows
                assert m['channel']=='zvv' and m['order']==order and m['seed']==x['seed']
                files=[];validation=None
                for f in m['files']:
                    path=Path(f['path']);assert f['bytes']>0 and path.stat().st_size==f['bytes']
                    product='particle-v1/'+str(path).split('/particle-v1/',1)[1]
                    files.append(dict(role=f['role'],product=product,sha256=f['sha256'],bytes=f['bytes']))
                    if f['role']=='validation':
                        assert path.stat().st_size<100000
                        b=path.read_bytes();assert digest(b)==f['sha256'];validation=json.loads(b)
                assert {f['role'] for f in files}=={'parton_inputs','closure','legacy_table','schema','validation'}
                v=validation;meta=v['metadata']
                assert meta['package'].split()[0]==producer and meta['frozen_manifest'] is None
                assert meta['primary_included'] and meta['channel']=='zvv' and meta['order']==order and meta['seed']==x['seed']
                assert v['events']==v['expected']==rows and not v['missing_ordinals']
                assert v['legacy_comparison']['pass_'] and v['views']['partons-v1']['overflow_events']==0
                assert v['forward_validation']['partons-v1']['finite']
            except (AssertionError,OSError,ValueError,KeyError) as exc:
                rejected.append(dict(key=key,reason=type(exc).__name__));continue
            maximum[order]=max(maximum[order],v['views']['partons-v1']['max_count'])
            selected[order].append(dict(key=key,seed=x['seed'],events=rows,producer=producer,
                primary_included=True,manifest_sha256=x['manifest_sha256'],lhe_sha256=meta['lhe_sha256'],files=files))
            if len(selected[order])==count:break
        assert len(selected[order])==count,(order,len(selected[order]))
    result=dict(schema='ggzh-parton-zvv-1m-selection-v1',study='parton-zvv-1m-v1',channel='zvv',
        selection_rule=f'SHA256 {salt}|order|seed; smallest200 LO/400 NLO eligible primary producerC shards',
        catalog_checked=catalog['checked'],catalog_sha256=digest(raw),parent_selection_sha256=None,
        events_per_sample=1000000,maximum_radiation_tokens=maximum,maximum_tensor_length=4+max(maximum.values()),
        eligibility_rejections=rejected,samples=selected)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(selection_sha256=digest(a.out.read_bytes()),maximum_tensor_length=result['maximum_tensor_length'],
        selected={k:len(v) for k,v in selected.items()},rejected=rejected)))

if __name__=='__main__':main()
