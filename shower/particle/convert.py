#!/usr/bin/env python3
"""Worker-only conversion of native Pythia records to lossless Parquet and ZH tensors."""
import argparse, csv, gzip, hashlib, json, math, os, resource, time
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


def lhe_events(path):
    with gzip.open(path, 'rt') as f:
        for _, node in ET.iterparse(f, events=('end',)):
            if node.tag != 'event':
                continue
            lines=[x.strip() for x in node.text.strip().splitlines() if x.strip() and not x.lstrip().startswith('#')]
            head=lines[0].split(); n=int(head[0]); particles=[]
            for line in lines[1:n+1]:
                x=line.split(); particles.append(dict(id=int(x[0]),status=int(x[1]),mothers=[int(x[2]),int(x[3])],
                    colours=[int(x[4]),int(x[5])],p4=[float(a) for a in x[6:10]],mass=float(x[10]),lifetime=float(x[11]),spin=float(x[12])))
            weights={x.attrib['id']:float(x.text) for x in node.findall('./rwgt/wgt')}
            yield dict(particles=particles,weight=float(head[2]),scalup=float(head[3]),weights=weights)
            node.clear()


def particle(x):
    return dict(id=int(x[0]),status=int(x[1]),p4=x[2:6],mass=x[6],col=int(x[7]),acol=int(x[8]),
                final=bool(x[9]),mothers=x[10],daughters=x[11])


def kin(v):
    x,y,z,e=v;pt=math.hypot(x,y);phi=math.atan2(y,x)
    eta=math.asinh(z/max(pt,1e-30));m=math.sqrt(max(0,e*e-x*x-y*y-z*z))
    return pt,eta,phi,m


def dphi(a,b):return (a-b+math.pi)%(2*math.pi)-math.pi


def encoding(jets):
    jj=jets['jets'];dr=-1.
    if len(jj)>1:
        a,b=kin(jj[0][:4]),kin(jj[1][:4]);dr=math.hypot(a[1]-b[1],dphi(a[2],b[2]))
    d=np.array(jets['dmerge'],dtype=np.float64)
    return [math.log1p(len(jj)),dr,float(dr>=0)]+np.log1p(np.maximum(d,0)).tolist()+(d>=0).astype(float).tolist()


def tensors(rows,view,cap):
    feat=[];vec=[];mask=[];ef=[];ps=[];overflow=[];multiplicity=[]
    for r in rows:
        z=r['boson_tokens'][0][0];h=r['boson_tokens'][1][0];angle=kin(z)[2]
        # Fixed four slots: Z,H,two ancestry-dressed direct Z children; mask absent objects.
        base=r['boson_tokens']
        radiation=[(j[:4],5) for j in r['radiation']['jets']] if view=='jets-v1' else r['truth_tokens']
        radiation.sort(key=lambda q:-kin(q[0])[0]);multiplicity.append(len(radiation));overflow.append(max(0,len(radiation)-cap))
        tokens=base+radiation[:cap];tokens += [(None,0)]*(cap+4-len(tokens))
        ff=np.zeros((7,cap+4),np.float32);vv=np.zeros((4,cap+4),np.float32);mm=np.zeros((1,cap+4),np.float32)
        for k,(v,kind) in enumerate(tokens):
            if v is None:continue
            pt,eta,phi,m=kin(v);ang=dphi(phi,angle);x,y,zz,e=v
            vv[:,k]=[x*math.cos(angle)+y*math.sin(angle),-x*math.sin(angle)+y*math.cos(angle),zz,e]
            ff[:,k]=[math.log1p(pt),eta,math.sin(ang),math.cos(ang),math.log1p(max(e,0)),math.log1p(m),kind];mm[0,k]=1
        feat.append(ff);vec.append(vv);mask.append(mm)
        zz=np.array(z)+np.array(h);k=kin(zz)
        ef.append([math.log1p(k[0]),k[1],math.log1p(k[3]),math.log1p(kin(z)[3]),float(r['channel']=='zvv')])
        ps.append(encoding(r['replacement_kt'] if view=='jets-v1' else r['partons']))
    return dict(pf_features=np.stack(feat),pf_vectors=np.stack(vec),pf_mask=np.stack(mask),
                event_features=np.array(ef,np.float32),ps_encoding_features=np.array(ps,np.float32),
                overflow=np.array(overflow,np.int32),radiation_count=np.array(multiplicity,np.int32))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,required=True);ap.add_argument('--channel',required=True)
    ap.add_argument('--order',required=True);ap.add_argument('--expected',type=int,required=True);a=ap.parse_args()
    start=time.time();out=Path('particle-out');out.mkdir(exist_ok=True)
    metadata=json.loads(Path('particle-context.json').read_text());metadata.update(seed=a.seed,shower_seed=a.seed+100000,channel=a.channel,order=a.order)
    metadata['versions']={'numpy':np.__version__,'pyarrow':pa.__version__}
    metadata['build']=Path('particle-install/BUILD-MANIFEST.txt').read_text()
    metadata['cards']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in Path('shower/cards').glob('*.cmnd')}
    with gzip.open('legacy.csv.gz','rt') as f:reference=list(csv.DictReader(f))
    with gzip.open('particle.table.csv.gz','rt') as f:closure=list(csv.DictReader(f))
    assert len(closure)==a.expected,(len(closure),a.expected)
    assert len(reference)==len(closure),'Legacy row count mismatch'
    numeric=np.array([[float(r[k]) for k in r] for r in closure]);old=np.array([[float(r[k]) for k in r] for r in reference])
    assert list(reference[0])==list(closure[0]),'Legacy schema differs'
    assert np.isfinite(numeric).all(),'Nonfinite closure data'
    assert np.allclose(old,numeric,rtol=2e-5,atol=2e-5),'Replay differs from legacy table'
    lhe=list(lhe_events('input.lhe.gz'));assert len(lhe)==a.expected
    raw_writer=obj_writer=None;rows=[];batch=[];obatch=[];seen=set();stable_check=[]
    with open('particle.particles.jsonl') as f:
        for line in f:
            r=json.loads(line);ordinal=r['lhe_ordinal'];i=r['accepted']
            assert 0<=ordinal<len(lhe) and ordinal not in seen;seen.add(ordinal)
            assert i==len(rows);assert r['z']>0 and r['h']>0
            r['channel']=a.channel
            for key in ('event','process'):
                ev=r[key]
                for p in ev:
                    assert np.isfinite(p[2:7]).all()
                    assert all(0<=x<len(ev) for x in p[10]+p[11]),'Invalid particle ancestry'
            lp=lhe[ordinal];incoming=[p for p in lp['particles'] if p['status']==-1];assert len(incoming)==2
            beam=[float(x) for x in metadata['beams']]
            x=[p['p4'][3]/beam[k] for k,p in enumerate(incoming)]
            w={k:float(v) for k,v in closure[i].items() if k.startswith('w:')}
            lp=dict(lp,weights=[dict(id=k,value=v) for k,v in lp['weights'].items()])
            raw=dict(channel=a.channel,order=a.order,seed=a.seed,lhe_ordinal=ordinal,accepted=i,table_row=i,
                event=[particle(p) for p in r['event']],process=[particle(p) for p in r['process']],
                z_index=r['z'],h_index=r['h'],lhe=lp,lhe_incoming=dict(x=x,ids=[p['id'] for p in incoming]),
                pythia_incoming=r['pythia_incoming'],shower_weights=w)
            batch.append(raw)
            obj=dict(seed=a.seed,lhe_ordinal=ordinal,radiation=r['radiation'],partons=r['partons'],
                     parton_indices=r['parton_indices'],status71_79_indices=r['status71_79_indices'],replacement_kt=r['replacement_kt'])
            obatch.append(obj)
            # Stable analysis-jet definition cross-check with independent native-record clustering.
            js=r['radiation']['jets'];nj=sum(kin(j[:4])[0]>30 for j in js)
            assert nj==int(float(closure[i]['njets30'])),('jets',i,nj,closure[i]['njets30'])
            stable_check.append(nj)
            ds=[j for j in r['event'][r['z']][11] if abs(r['event'][j][0]) in (11,12,13,14,15,16)]
            ds.sort(key=lambda j:-kin(r['event'][j][2:6])[0]);ds=ds[:2]
            bt=[(r['event'][r['z']][2:6],1),(r['event'][r['h']][2:6],2)]+[(r['event'][j][2:6],3 if abs(r['event'][j][0])%2 else 4) for j in ds]
            bt += [(None,0)]*(4-len(bt))
            rows.append(dict(channel=a.channel,lhe_ordinal=ordinal,boson_tokens=bt,
                truth_tokens=[(r['event'][j][2:6],6) for j in r['parton_indices']],
                radiation=r['radiation'],partons=r['partons'],replacement_kt=r['replacement_kt']))
            if len(batch)==16:
                for buf,key in ((batch,'raw'),(obatch,'objects')):
                    t=pa.Table.from_pylist(buf)
                    if key=='raw':
                        if raw_writer is None:raw_writer=pq.ParquetWriter(out/'records.parquet',t.schema,compression='zstd')
                        raw_writer.write_table(t.cast(raw_writer.schema))
                    else:
                        if obj_writer is None:obj_writer=pq.ParquetWriter(out/'objects.parquet',t.schema,compression='zstd')
                        obj_writer.write_table(t.cast(obj_writer.schema))
                    buf.clear()
    for buf,writer in ((batch,raw_writer),(obatch,obj_writer)):
        if buf:writer.write_table(pa.Table.from_pylist(buf,schema=writer.schema))
        writer.close()
    assert len(rows)==len(closure)
    # Explicit mapping remains valid even when Pythia consumes multiple LHE records in one next().
    missing=sorted(set(range(len(lhe)))-seen)
    assert not missing,'Incomplete replay; retained diagnostics, no ready publication'
    sw=np.array([float(r.get('w:2000',r['w:nominal'])) for r in closure])
    cc={k:np.array([float(r[k]) for r in closure]) for k in closure[0]}
    np.savez_compressed(out/'closure.npz',**cc)
    summaries={}
    for view,cap in (('jets-v1',64),('partons-v1',128)):
        t=tensors(rows,view,cap)
        ids=np.array([r['lhe_ordinal'] for r in rows],np.int32)
        # Full shard counts and frozen omitted shards are multiples of 5 and 2.
        # Thus shard-local accepted ordinal reproduces global window folds/target halves.
        assert np.all((cc['mz']>15)&(cc['mz']<150)),'Window-index mapping requires explicit prefix manifest'
        np.savez_compressed(out/f'{view}.npz',**t,lhe_ordinal=ids,table_row=np.arange(len(rows)),seed=np.full(len(rows),a.seed),
            source_fold=np.arange(len(rows))%5,target_half=np.arange(len(rows))%2,weight=sw)
        with np.load(out/f'{view}.npz') as inp:
            for key in ('pf_features','pf_vectors','pf_mask','event_features','ps_encoding_features'):
                assert np.isfinite(inp[key]).all()
            assert np.array_equal(inp['lhe_ordinal'],ids)
        # A real-record forward evaluation through public attention, without fitting or RNG changes to shower.
        # Performed by forward.py after conversion using the public ParT module.
        tail=cc['pth']>450
        summaries[view]=dict(cap=cap,max_count=int(t['radiation_count'].max()),overflow_events=int((t['overflow']>0).sum()),
            overflow_hard_tail=int(((t['overflow']>0)&tail).sum()),hard_tail_events=int(tail.sum()))
    assert pq.read_metadata(out/'records.parquet').num_rows==len(rows)
    assert pq.read_metadata(out/'objects.parquet').num_rows==len(rows)
    report=dict(schema='ggzh-particle-v1',metadata=metadata,events=len(rows),expected=a.expected,missing_ordinals=missing,
        legacy_comparison=dict(pass_=True,max_absolute_difference=float(np.max(abs(old-numeric))),rtol=2e-5,atol=2e-5),
        sumw=float(sw.sum()),sumw2=float(sw@sw),negative_fraction=float(np.mean(sw<0)),
        multiplicity_counts={str(j):stable_check.count(j) for j in sorted(set(stable_check))},views=summaries,
        extraction_seconds=time.time()-start,peak_rss_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        sizes={p.name:p.stat().st_size for p in out.iterdir()})
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('metadata',)},indent=2))

if __name__=='__main__':main()
