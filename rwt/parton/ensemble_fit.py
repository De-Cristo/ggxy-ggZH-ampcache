"""One ensemble member: train on a whole frozen 1M/class subset, predict the calibration and evaluation sets.

Same network, optimizer, loss and early stopping as fit.py; no source folds or target halves.
Inputs: parton-data/ (the member subset), calib/parton-data/ and eval/parton-data/ (disjoint sets).
"""
import json,os,tarfile,time
from pathlib import Path
from common import worker,configuration,sha
worker()
import numpy as np
import torch
from torch.nn.functional import binary_cross_entropy_with_logits
from model import ZHParT

c=configuration();seed=c['initialization_seed'];assert isinstance(seed,int) and c.get('rapidity_convention')=='transverse-mass-rapidity-v1'
assert torch.cuda.is_available() and torch.cuda.device_count()==1,'Exactly one assigned CUDA GPU required'
torch.set_num_threads(int(os.environ.get('RWT_THREADS','4')));torch.manual_seed(seed);torch.cuda.manual_seed_all(seed)
rng=np.random.default_rng(seed)

def arrays(root,order):
    return {p.stem:np.load(p,mmap_mode='r',allow_pickle=False) for p in (Path(root)/'parton-data'/order).glob('*.npy')}
source,target=arrays('.','lo-cms'),arrays('.','nlo')
assert len(source['weight'])==len(target['weight'])==1000000
held={(name,order):arrays(name,order) for name in ('calib','eval') for order in ('lo-cms','nlo')}
for name in ('calib','eval'):
    other=json.loads(Path(name,'parton-data','selection.json').read_text())
    for order in ('lo-cms','nlo'):
        mine={x['seed'] for x in json.loads(Path('parton-data/selection.json').read_text())['samples'][order]}
        assert not mine & {x['seed'] for x in other['samples'][order]},f'{name} {order} overlaps the member subset'
s_pool=np.flatnonzero(source['in_window']);t_pool=np.flatnonzero(target['in_window'])
def split(pool):
    order=rng.permutation(pool);nv=max(1,int(len(pool)*c['validation_fraction']));assert len(pool)-nv>1
    return order[nv:],order[:nv]
si,sv=split(s_pool);ti,tv=split(t_pool)
smoke=int(os.environ.get('RWT_ENSEMBLE_SMOKE_ROWS','0'))
if smoke:
    # Technical smoke test of the complete chain only: tiny pools, one epoch, never a result.
    si,sv,ti,tv=si[:smoke],sv[:smoke//10],ti[:smoke],tv[:smoke//10];c['epochs']=1

def moments():
    ps=np.zeros(7);ps2=np.zeros(7);pn=0;gs=np.zeros(c['model']['global_features']);gs2=np.zeros(c['model']['global_features']);gn=0
    for data,indices in ((source,si),(target,ti)):
        for start in range(0,len(indices),8192):
            ix=indices[start:start+8192];x=np.asarray(data['pf_features'][ix],dtype=np.float64);mask=data['pf_mask'][ix]
            x=x*mask;ps+=x.sum(axis=(0,2));ps2+=(x*x).sum(axis=(0,2));pn+=int(mask.sum())
            g=np.asarray(data['globals'][ix],dtype=np.float64);gs+=g.sum(0);gs2+=(g*g).sum(0);gn+=len(ix)
    pm=ps/pn;gm=gs/gn;pstd=np.sqrt(np.maximum(ps2/pn-pm*pm,0));gstd=np.sqrt(np.maximum(gs2/gn-gm*gm,0))
    return dict(particle_mean=pm.tolist(),particle_std=np.where(pstd>1e-8,pstd,1).tolist(),global_mean=gm.tolist(),global_std=np.where(gstd>1e-8,gstd,1).tolist(),particle_tokens=pn,events=gn)
norm=moments()

def batch(data,indices):
    return [torch.as_tensor(np.asarray(data[k][indices]),device='cuda') for k in ('pf_features','pf_vectors','pf_mask','globals')]

def pools(six,tix):
    sw=np.full(len(six),.5/len(six));tw=np.maximum(np.asarray(target['weight'][tix]),0);assert tw.sum()>0
    return np.concatenate([six,tix]),np.concatenate([np.zeros(len(six),np.int8),np.ones(len(tix),np.int8)]),np.concatenate([sw,.5*tw/tw.sum()])
train_ids,labels,lossweights=pools(si,ti);val_ids,vlabels,vweights=pools(sv,tv)

def mixed(ids,y,w):
    # Preserve the shuffled batch order while drawing from separate memory maps.
    values=[]
    for key in ('pf_features','pf_vectors','pf_mask','globals'):
        shape=(len(ids),*source[key].shape[1:]);a=np.empty(shape,dtype=source[key].dtype)
        s=y==0;a[s]=source[key][ids[s]];a[~s]=target[key][ids[~s]]
        values.append(torch.as_tensor(a,device='cuda'))
    return values,torch.as_tensor(y,dtype=torch.float32,device='cuda'),torch.as_tensor(w,dtype=torch.float32,device='cuda')

net=ZHParT(c,norm).cuda();opt=torch.optim.AdamW(net.parameters(),lr=c['learning_rate'],weight_decay=c['weight_decay'])
best=float('inf');best_state=None;bad=0;history=[];start_time=time.time();bs=c['batch_size']
for epoch in range(c['epochs']):
    net.train();order=rng.permutation(len(train_ids));trainloss=0.
    for start in range(0,len(order),bs):
        ix=order[start:start+bs];x,y,w=mixed(train_ids[ix],labels[ix],lossweights[ix]);opt.zero_grad(set_to_none=True)
        loss=(binary_cross_entropy_with_logits(net(*x),y,reduction='none')*w).sum()*len(train_ids)/len(ix)
        assert torch.isfinite(loss),'Nonfinite training loss'
        loss.backward();torch.nn.utils.clip_grad_norm_(net.parameters(),c['gradient_clip'],error_if_nonfinite=True);opt.step()
        trainloss+=float(loss.detach())*len(ix)/len(train_ids)
    net.eval();score=0.
    with torch.no_grad():
        for start in range(0,len(val_ids),bs):
            sl=slice(start,start+bs);x,y,w=mixed(val_ids[sl],vlabels[sl],vweights[sl])
            score+=float((binary_cross_entropy_with_logits(net(*x),y,reduction='none')*w).sum())
    assert np.isfinite(score)
    history.append(dict(epoch=epoch+1,train_loss=trainloss,validation_loss=score,elapsed_seconds=time.time()-start_time))
    print(json.dumps(history[-1]),flush=True)
    if score<best-c['min_delta']:
        best=score;bad=0;best_state={k:v.detach().cpu().clone() for k,v in net.state_dict().items()}
    else:
        bad+=1
        if bad>=c['patience']:break
assert best_state is not None
net.load_state_dict(best_state);net.eval()

def predict(model,data):
    # Every row, in stored order, 256-event batches (the certified inference protocol of fit.py).
    n=len(data['weight']);out=np.empty(n,np.float64)
    with torch.no_grad():
        for start in range(0,n,bs):
            sl=np.arange(start,min(start+bs,n));out[sl]=model(*batch(data,sl)).cpu().numpy()
    assert np.isfinite(out).all()
    return out
out=Path('parton-fit');assert not out.exists();out.mkdir()
predictions={};inference_start=time.time()
for (name,order),data in held.items():
    predictions[name,order]=predict(net,data)
    np.savez_compressed(out/f'logits-{name}-{order}.npz',seed=np.asarray(data['seed']),lhe_ordinal=np.asarray(data['lhe_ordinal']),
        table_row=np.asarray(data['table_row']),in_window=np.asarray(data['in_window']),logits=predictions[name,order])
inference_seconds=time.time()-inference_start
card=dict(state=best_state,config=c,normalization=norm,member=c['study'],initialization_seed=seed,selection_sha256=c['selection_sha256'],package=Path('VERSION').read_text().strip())
torch.save(card,out/'model.pt')
loaded=torch.load(out/'model.pt',map_location='cpu',weights_only=True);check=ZHParT(loaded['config'],loaded['normalization']).cuda();check.load_state_dict(loaded['state']);check.eval()
reload_audit=dict(protocol='identical-batch-calibration-lo-v1',batch_size=bs,rtol=1e-4,atol=2e-5,events=len(held['calib','lo-cms']['weight']),state_tensors_equal=False,same_batch_pass=False)
try:
    assert set(loaded['state'])==set(best_state) and all(torch.equal(best_state[k],loaded['state'][k]) for k in best_state)
    reload_audit['state_tensors_equal']=True
    fresh=predict(check,held['calib','lo-cms']);reference=predictions['calib','lo-cms']
    reload_audit['max_absolute_difference']=float(np.max(np.abs(fresh-reference)))
    np.testing.assert_allclose(fresh,reference,rtol=1e-4,atol=2e-5);reload_audit['same_batch_pass']=True
except Exception as exc:
    reload_audit['failure']=repr(exc)
    raise
finally:
    (out/'reload-audit.json').write_text(json.dumps(reload_audit,indent=2)+'\n')

def selection_sha(name):
    return json.loads(Path(name,'parton-data','audit.json').read_text())['selection_sha256']
def summary(name):
    # Unit-normalized in-window weights of this member alone: a quick health check, not a result.
    lo,hi=held[name,'lo-cms'],held[name,'nlo'];w=np.where(np.asarray(lo['in_window']),np.exp(predictions[name,'lo-cms']-predictions[name,'lo-cms'].max()),0)
    w*=np.asarray(lo['weight'])
    return dict(in_window_lo=int(np.asarray(lo['in_window']).sum()),effective_fraction=float(w.sum()**2/np.sum(w*w)/max(1,int(np.asarray(lo['in_window']).sum()))),
                logit_lo_mean=float(predictions[name,'lo-cms'][np.asarray(lo['in_window'])].mean()),logit_nlo_mean=float(predictions[name,'nlo'][np.asarray(hi['in_window'])].mean()))
record=dict(member=c['study'],channel=c['channel'],selection_sha256=c['selection_sha256'],calibration_selection_sha256=selection_sha('calib'),
    evaluation_selection_sha256=selection_sha('eval'),package=card['package'],initialization_seed=seed,rapidity_convention=c['rapidity_convention'],
    source_train=len(si),source_validation=len(sv),target_train=len(ti),target_validation=len(tv),history=history,best_validation_loss=best,
    epochs_run=len(history),training_seconds=inference_start-start_time,inference_seconds=inference_seconds,checkpoint_reload='PASS',reload_audit=reload_audit,
    model_sha256=sha(out/'model.pt'),torch=str(torch.__version__),gpu=torch.cuda.get_device_name(0),parameters=sum(p.numel() for p in net.parameters()),
    smoke_rows=smoke,global_features=c['model']['global_features'],health={name:summary(name) for name in ('calib','eval')})
(out/'record.json').write_text(json.dumps(record,indent=2)+'\n');Path('ensemble-fit.json').write_text(json.dumps(record,indent=2)+'\n')
with tarfile.open('ensemble-fit.tgz','w:gz') as t:t.add(out,arcname='parton-fit')
print('ENSEMBLE_MEMBER_COMPLETED',c['study'],flush=True)
