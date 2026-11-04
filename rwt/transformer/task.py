"""Real full-sample ParT fold jobs and OOF closure merge; worker only."""
import json
import os
from pathlib import Path
import platform
import sys
import tarfile

assert platform.node().split('.')[0]!='cmscon' and os.environ.get('_CONDOR_SCRATCH_DIR')
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ggzh_reweight_bdt as BDT
import ggzh_reweight_common as C
import ggzh_rwt_manifest as M
from ggzh_rwt10m import diagnostics
from train import fit, network, predict_logits, normalized_ratios

task,job=sys.argv[1],int(sys.argv[2])
config=json.loads(Path('rwt/transformer/config.json').read_text())
manifest=M.load(Path(os.environ['RWT_INPUT_MANIFEST']))
M.verify(manifest,Path('rwt-out'))
torch.set_num_threads(4)
assert torch.cuda.is_available(), 'No CPU fallback'
BDT.configure(C.FEATURE_SETS['C'])
src,tgt=C.Sample('source',Path('rwt-out/lo-cms.npz')),C.Sample('target',Path('rwt-out/nlo-cp5.npz'))
xs=BDT.prepare(src).astype('float32')
ix=np.flatnonzero(src.in_window);folds=np.arange(len(ix))%5
assert np.isfinite(xs[ix]).all()
tags=['full','A_to_B','B_to_A']
if task=='fold':
    assert 0<=job<15
    direction,k=divmod(job,5);tag=tags[direction]
    tsel=tgt.in_window.copy()
    if direction:tsel &= np.arange(tgt.n)%2 == (0 if direction==1 else 1)
    xt=BDT.prepare(tgt).astype('float32');assert np.isfinite(xt[tsel]).all()
    selected=ix[folds==k]
    net,card=fit(xs[ix[folds!=k]],xt[tsel],tgt.w0[tsel],config,config['seed']+100*direction+k,'cuda')
    logits=predict_logits(net,xs[selected],card,'cuda',config['batch_size'])
    out=Path('transformer-out');out.mkdir()
    card.update(features=BDT.FEATURES,tag=tag,fold=k,config=config,input_selection=manifest)
    torch.save(card,out/'model.pt')
    loaded=torch.load(out/'model.pt',map_location='cpu',weights_only=True)
    net.load_state_dict(loaded['state'])
    z=predict_logits(net,xs[selected[:128]],loaded,'cuda')
    np.testing.assert_allclose(logits[:128],z,rtol=1e-4,atol=2e-5)
    np.savez_compressed(out/'logits.npz',source_index=selected,logits=logits)
    record=dict(job=job,tag=tag,fold=k,config=config,features=BDT.FEATURES,input_selection=manifest,
                source_events=src.n,target_events=tgt.n,heldout_events=len(selected),
                validation_loss=card['validation_loss'],checkpoint_reload='PASS',
                package=Path('VERSION').read_text().strip(),torch=str(torch.__version__),gpu=torch.cuda.get_device_name(0))
    Path('transformer-part.json').write_text(C.dumps(record))
    with tarfile.open('transformer-part.tgz','w:gz') as t:
        for p in out.iterdir():t.add(p,arcname=p.name)
    print(C.dumps(record),flush=True)
elif task=='merge':
    out=Path('transformer-out');out.mkdir()
    records=[]
    weights=dict(run=src['run'],event=src['event'],in_window=src.in_window)
    result=dict(method='public-ParT-feature-token-C',config=config,features=BDT.FEATURES,input_selection=manifest,
                provenance=dict(package=Path('VERSION').read_text().strip(),source_events=src.n,target_events=tgt.n))
    models=out/'models';models.mkdir()
    for direction,tag in enumerate(tags):
        logits=np.full(len(ix),np.nan)
        for k in range(5):
            j=5*direction+k;directory=out/f'part-{j}';directory.mkdir()
            with tarfile.open(f'parts/fold-{j}.tgz') as t:
                assert all(m.name in ['model.pt','logits.npz'] for m in t.getmembers())
                t.extractall(directory)
            z=np.load(directory/'logits.npz',allow_pickle=False)
            assert np.array_equal(z['source_index'],ix[folds==k])
            card=torch.load(directory/'model.pt',map_location='cpu',weights_only=True)
            assert card['tag']==tag and card['fold']==k and card['config']==config and card['input_selection']==manifest
            net=network(17,[64]).cuda();net.load_state_dict(card['state'])
            check=predict_logits(net,xs[ix[folds==k][:128]],card,'cuda')
            np.testing.assert_allclose(z['logits'][:128],check,rtol=1e-4,atol=2e-5)
            logits[folds==k]=z['logits']
            (directory/'model.pt').rename(models/f'{tag}-fold-{k}.pt')
            del net,z,card
        assert np.isfinite(logits).all()
        w=np.ones(src.n);w[ix],norm=normalized_ratios(logits);weights['w_'+tag]=w
        test=tgt.in_window.copy()
        if direction:test &= np.arange(tgt.n)%2 == (1 if direction==1 else 0)
        result[tag]=BDT.evaluate(src,w,src.in_window,tgt,test,tag,print)
        result[tag]['log_normalization']=norm
    r=result['full']['reweighted']
    result['summary']={'1_stage13':r['pass_1'],'2_variables':r['pass_2'],'2b_features':r['pass_2b'],
        '3_stability':result['A_to_B']['pass'] and result['B_to_A']['pass']}
    result['summary']['all_pass']=all(result['summary'].values())
    result['diagnostics']=diagnostics(src,weights['w_full'],tgt)
    result['claim_reweighting_works']=bool(result['summary']['all_pass'] and result['diagnostics']['hardjet']['pth450']['within_band'])
    result['target_neff']={str(c):float(tgt.w0[tgt.in_window&(tgt['ptz']>=c)].sum()**2/(tgt.w0[tgt.in_window&(tgt['ptz']>=c)]**2).sum()) for c in [0,400,600]}
    (out/'result.json').write_text(C.dumps(result));np.savez_compressed(out/'weights.npz',**weights)
    with tarfile.open('transformer-models.tgz','w:gz') as t:t.add(models,arcname='models')
else:raise ValueError(task)
