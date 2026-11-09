"""Verify the frozen real subset on a CPU worker and create reusable lossless arrays."""
import concurrent.futures,json,os,subprocess,sys,tarfile,time
from pathlib import Path
from common import worker,configuration,sha
worker()
import numpy as np
sys.path.insert(0,str(Path('rwt').resolve()))
import ggzh_table as Table
from ggzh_reweight_common import Sample
c=configuration();m=json.loads(Path('production/parton-zll-1m-v1.json').read_text())
root=os.environ['PAPER_EOS_ROOT'].rstrip('/')
stage=Path('parton-staged');out=Path('parton-data');assert not out.exists();out.mkdir()

def download(item):
    order,seed,entry=item;local=stage/order/str(seed)/entry['role']/Path(entry['product']).name;local.parent.mkdir(parents=True,exist_ok=True)
    remote=root+'/'+entry['product']
    for url,target in ((remote+'.sha256',Path(str(local)+'.sha256')),(remote,local)):
        subprocess.run(['xrdcp','--retry','3','--cksum','adler32',url,str(target)],check=True,stdout=subprocess.DEVNULL)
    assert Path(str(local)+'.sha256').read_text().split()[0]==entry['sha256']
    assert local.stat().st_size==entry['bytes'] and sha(local)==entry['sha256'],str(local)
    return (order,seed,entry['role']),local

items=[(order,x['seed'],f) for order,xs in m['samples'].items() for x in xs for f in x['files']]
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:files=dict(pool.map(download,items))
print('Verified selected artifact files',len(files),flush=True)
audit=dict(study=c['study'],selection_sha256=c['selection_sha256'],artifact_count=len(files),samples={},package=Path('VERSION').read_text().strip())
L=m['maximum_tensor_length'];assert L<=132
for order,selection in m['samples'].items():
    selection=sorted(selection,key=lambda x:x['seed']);n=sum(x['events'] for x in selection);assert n==1000000
    dest=out/order;dest.mkdir();offset=0;maps={};table_dir=Path('tables')/order;table_dir.mkdir(parents=True)
    specs={'pf_features':((7,L),'float32'),'pf_vectors':((4,L),'float32'),'pf_mask':((1,L),'bool'),'globals':((24,),'float32'),
           'weight':((),'float64'),'in_window':((),'bool'),'source_fold':((),'int8'),'target_half':((),'int8'),'seed':((),'uint32'),'lhe_ordinal':((),'int32'),'table_row':((),'int32')}
    for key,(shape,dtype) in specs.items():maps[key]=np.lib.format.open_memmap(dest/(key+'.npy'),mode='w+',dtype=dtype,shape=(n,*shape))
    row_counts=[];sumw=0.;outside=0
    for entry in selection:
        seed=entry['seed'];rows=entry['events'];sl=slice(offset,offset+rows)
        report=json.loads(files[order,seed,'validation'].read_text());meta=report['metadata']
        assert meta['package'].split()[0]==entry['producer'] and meta['lhe_sha256']==entry['lhe_sha256']
        assert meta['channel']=='zll' and meta['order']==order and meta['seed']==seed and meta['primary_included']
        assert meta['frozen_manifest']==m['parent_selection_sha256']
        assert report['events']==report['expected']==rows and not report['missing_ordinals']
        assert report['legacy_comparison']['pass_'] and report['views']['partons-v1']['overflow_events']==0
        assert report['forward_validation']['partons-v1']['finite']
        schema=json.loads(files[order,seed,'schema'].read_text());assert schema['capacities']['partons-v1']==128
        with np.load(files[order,seed,'parton_inputs'],allow_pickle=False) as d,np.load(files[order,seed,'closure'],allow_pickle=False) as cl:
            assert len(d['weight'])==rows and np.array_equal(d['seed'],np.full(rows,seed))
            assert np.array_equal(d['table_row'],np.arange(rows)) and np.array_equal(cl['ievt'],d['table_row'])
            assert np.array_equal(np.sort(d['lhe_ordinal']),np.arange(rows))
            assert not np.any(d['overflow']) and np.array_equal(d['target_half'],np.arange(rows)%2)
            window=(cl['mz']>15)&(cl['mz']<150);assert np.array_equal(window,d['in_window'])
            if order=='lo-cms':assert np.isin(d['source_fold'][window],np.arange(5)).all() and (d['source_fold'][~window]==-1).all()
            else:assert (d['source_fold']==-1).all()
            w=cl['w:2000'] if order=='nlo' else cl['w:nominal'];assert np.array_equal(w,d['weight'])
            assert np.isclose(w.sum(),report['sumw'],rtol=1e-12,atol=1e-12)
            assert np.isclose(np.square(w).sum(),report['sumw2'],rtol=1e-12,atol=1e-12)
            mask=d['pf_mask'];assert np.isin(mask,[0,1]).all() and not mask[:,:,L:].any()
            for key in ('pf_features','pf_vectors'):
                a=d[key];assert np.isfinite(a).all() and not np.any(a[:,:,L:])
                assert not np.any(a*np.logical_not(mask)),'Nonzero padded input'
            for key in ('event_features','ps_encoding_features'):assert np.isfinite(d[key]).all()
            for key in ('pf_features','pf_vectors','pf_mask'):maps[key][sl]=d[key][:,:,:L]
            maps['globals'][sl]=np.concatenate([d['event_features'],d['ps_encoding_features']],axis=1)
            for key in specs.keys()-{'pf_features','pf_vectors','pf_mask','globals'}:maps[key][sl]=d[key]
            sumw+=float(w.sum());outside+=int((~window).sum())
        (table_dir/f'table-{seed}.csv.gz').symlink_to(files[order,seed,'legacy_table'].resolve())
        offset+=rows;row_counts.append(dict(seed=seed,start=offset-rows,rows=rows))
    for a in maps.values():a.flush()
    del maps
    label='lo-cms' if order=='lo-cms' else 'nlo-cp5'
    Table.main(['--label',label,'--tables',str(table_dir/'table-*.csv.gz'),'--nominal','nominal' if order=='lo-cms' else '2000','--out',str(out/(label+'.npz'))])
    with np.load(out/(label+'.npz'),allow_pickle=False) as d:
        cols={k:d[k] for k in d.files}
    seed_array=np.load(dest/'seed.npy',mmap_mode='r');row_array=np.load(dest/'table_row.npy',mmap_mode='r')
    assert np.array_equal(cols['run'],seed_array) and np.array_equal(cols['event'],row_array)
    assert np.array_equal(cols['w0'],np.load(dest/'weight.npy',mmap_mode='r'))
    # Add the original saved splits to the scalar interface without regenerating ranks.
    for key in ('source_fold','target_half'):cols[key]=np.array(np.load(dest/(key+'.npy'),mmap_mode='r'))
    s=Sample(label,cols=cols);assert np.array_equal(s.in_window,np.load(dest/'in_window.npy',mmap_mode='r'))
    np.savez_compressed(out/(label+'.npz'),**cols)
    audit['samples'][order]=dict(rows=n,in_window=n-outside,outside_window=outside,sumw=sumw,shards=row_counts,max_tensor_length=L,zero_overflow=True)
    print(order,audit['samples'][order]['in_window'],'in window',flush=True)
(out/'selection.json').write_text(json.dumps(m,indent=2)+'\n');(out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
paths=sorted(p for p in out.rglob('*') if p.is_file())
(out/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(out)}\n' for p in paths))
with tarfile.open('parton-prepared.tgz','w:gz',compresslevel=1) as t:t.add(out,arcname='parton-data')
Path('parton-prepare.json').write_text(json.dumps(audit,indent=2)+'\n')
print('PREPARED_REAL_1M_EACH',Path('parton-prepared.tgz').stat().st_size,flush=True)
