"""One complete real-subset GPU fold, with training-only normalization and reload validation."""
import copy,json,os,sys,tarfile,time
from pathlib import Path
from common import worker,configuration,arrays,sha
worker()
import numpy as np
import torch
from torch.nn.functional import binary_cross_entropy_with_logits
from model import ZHParT

c=configuration();task=int(sys.argv[1]);assert 0<=task<15
fold=task%5;direction=task//5;tag=c['tags'][direction];seed=c['seed_base']+task
assert torch.cuda.is_available() and torch.cuda.device_count()==1,'Exactly one assigned CUDA GPU required'
torch.set_num_threads(int(os.environ.get('RWT_THREADS','4')));torch.manual_seed(seed);torch.cuda.manual_seed_all(seed)
rng=np.random.default_rng(seed);source,target=arrays('lo-cms'),arrays('nlo')
assert len(source['weight'])==len(target['weight'])==1000000
s_pool=np.flatnonzero(source['in_window']&(source['source_fold']!=fold));test=np.flatnonzero(source['in_window']&(source['source_fold']==fold))
tmask=target['in_window'].copy()
if direction:tmask &= target['target_half']==direction-1
t_pool=np.flatnonzero(tmask)
def split(pool):
    order=rng.permutation(pool);nv=max(1,int(len(pool)*c['validation_fraction']));assert len(pool)-nv>1
    return order[nv:],order[:nv]
si,sv=split(s_pool);ti,tv=split(t_pool)
assert not np.intersect1d(si,test).size and not np.intersect1d(sv,test).size
if direction:assert (target['target_half'][ti]==direction-1).all() and (target['target_half'][tv]==direction-1).all()

def moments():
    ps=np.zeros(7);ps2=np.zeros(7);pn=0;gs=np.zeros(24);gs2=np.zeros(24);gn=0
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
logits=np.empty(len(test),np.float64)
with torch.no_grad():
    for start in range(0,len(test),bs):
        sl=slice(start,start+bs);logits[sl]=net(*batch(source,test[sl])).cpu().numpy()
assert np.isfinite(logits).all()
out=Path('parton-fit');assert not out.exists();out.mkdir()
card=dict(state=best_state,config=c,normalization=norm,task=task,tag=tag,fold=fold,selection_sha256=c['selection_sha256'],package=Path('VERSION').read_text().strip())
torch.save(card,out/'model.pt')
loaded=torch.load(out/'model.pt',map_location='cpu',weights_only=True);check=ZHParT(loaded['config'],loaded['normalization']).cuda();check.load_state_dict(loaded['state']);check.eval()
reload_audit=dict(protocol='identical-batch-full-heldout-v2',batch_size=bs,rtol=1e-4,atol=2e-5,
                  events=len(test),state_tensors_equal=False,same_batch_pass=False)
try:
    assert set(loaded['state'])==set(best_state)
    assert all(torch.equal(best_state[k],loaded['state'][k]) for k in best_state)
    reload_audit['state_tensors_equal']=True
    maximum=0.
    with torch.no_grad():
        for start in range(0,len(test),bs):
            sl=slice(start,start+bs);fresh=check(*batch(source,test[sl])).cpu().numpy()
            maximum=max(maximum,float(np.max(np.abs(fresh-logits[sl]))))
            np.testing.assert_allclose(fresh,logits[sl],rtol=1e-4,atol=2e-5)
        short=check(*batch(source,test[:128])).cpu().numpy()
    reload_audit.update(same_batch_pass=True,max_absolute_difference=maximum,
        different_batch_max_absolute_difference=float(np.max(np.abs(short-logits[:128]))),
        different_batch_pass=bool(np.allclose(short,logits[:128],rtol=1e-4,atol=2e-5)),
        supported_inference='256-event batches in original order; final remainder as stored; arbitrary batch portability is not certified')
except Exception as exc:
    reload_audit['failure']=repr(exc)
    raise
finally:
    (out/'reload-audit.json').write_text(json.dumps(reload_audit,indent=2)+'\n')
    Path('parton-reload-audit.json').write_text(json.dumps(reload_audit,indent=2)+'\n')
np.savez_compressed(out/'logits.npz',source_index=test,seed=np.asarray(source['seed'][test]),lhe_ordinal=np.asarray(source['lhe_ordinal'][test]),table_row=np.asarray(source['table_row'][test]),logits=logits)
record=dict(task=task,tag=tag,fold=fold,selection_sha256=c['selection_sha256'],package=card['package'],source_train=len(si),source_validation=len(sv),source_test=len(test),target_train=len(ti),target_validation=len(tv),target_half=None if not direction else direction-1,history=history,best_validation_loss=best,checkpoint_reload='PASS',model_sha256=sha(out/'model.pt'),torch=str(torch.__version__),gpu=torch.cuda.get_device_name(0),training_seconds=time.time()-start_time,parameters=sum(p.numel() for p in net.parameters()),particle_attention_layers=8,scalar_token_substitute=False)
record['reload_audit']=reload_audit
record['channel']=c['channel']
(out/'record.json').write_text(json.dumps(record,indent=2)+'\n');Path('parton-fit.json').write_text(json.dumps(record,indent=2)+'\n')
with tarfile.open('parton-fit.tgz','w:gz') as t:t.add(out,arcname='parton-fit')
print('REAL_PARTICLE_FOLD_COMPLETED',task,flush=True)
