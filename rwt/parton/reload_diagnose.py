"""Real trained-checkpoint GPU diagnosis, without fitting, tolerance changes or fake events."""
import json,sys,tarfile
from pathlib import Path
from common import worker,sha
worker()
import numpy as np
import torch
from model import ZHParT
assert torch.cuda.is_available() and torch.cuda.device_count()==1
torch.set_num_threads(2)
with tarfile.open('reference-fit.tgz') as t:
    assert all((m.name=='parton-fit' or m.name.startswith('parton-fit/')) and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers());t.extractall('reference')
p=Path('reference/parton-fit');r=json.loads((p/'record.json').read_text());assert r['task']==0 and r['checkpoint_reload']=='PASS'
assert sha(p/'model.pt')==r['model_sha256']
card=torch.load(p/'model.pt',map_location='cpu',weights_only=True);config=card['config']
assert config['selection_sha256']=='e656b2bef79777d978111a0e0fa106320d5bd584cb3f7d8600d6ab89674f4160'
with np.load(p/'logits.npz',allow_pickle=False) as z:
    ids={k:np.array(z[k][:256]) for k in ('seed','table_row','lhe_ordinal')};saved=np.array(z['logits'][:256])
assert (ids['seed']==110011).all()
with np.load('reference-input.npz',allow_pickle=False) as z:
    ix=ids['table_row'].astype(int);assert np.array_equal(z['lhe_ordinal'][ix],ids['lhe_ordinal'])
    assert z['in_window'][ix].all() and (z['source_fold'][ix]==0).all()
    assert not z['overflow'][ix].any() and not z['pf_mask'][ix,:,44:].any()
    values=[np.array(z[k][ix,:,:44]) for k in ('pf_features','pf_vectors')]+[np.array(z['pf_mask'][ix,:,:44],dtype=bool),np.concatenate([z['event_features'][ix],z['ps_encoding_features'][ix]],axis=1)]
inputs=[torch.as_tensor(a,device='cuda') for a in values]
def net():
    n=ZHParT(config,card['normalization']).cuda();n.load_state_dict(card['state']);n.eval();return n

def compare(a,b):
    assert np.isfinite(a).all() and np.isfinite(b).all()
    return dict(pass_=bool(np.allclose(a,b,rtol=1e-4,atol=2e-5)),max_absolute_difference=float(np.max(np.abs(a-b))),max_relative_difference=float(np.max(np.abs(a-b)/np.maximum(np.abs(b),1e-30))),max_relative_odds_shift=float(np.max(np.abs(np.expm1(a-b)))))

def evaluate(model,values):
    with torch.no_grad():return model(*values).cpu().numpy().astype(np.float64)

a=net();b=net();state_equal=all(torch.equal(v,b.state_dict()[k]) for k,v in a.state_dict().items())
full=evaluate(a,inputs);repeat=evaluate(a,inputs);reload_full=evaluate(b,inputs);half=evaluate(a,[v[:128] for v in inputs])
report=dict(study='parton-zll-1m-v1',selection_sha256=config['selection_sha256'],diagnostic_package=Path('VERSION').read_text().strip(),trained_package=r['package'],actual_real_records=256,source_seed=110011,reference_model_sha256=r['model_sha256'],torch=str(torch.__version__),gpu=torch.cuda.get_device_name(0),rtol=1e-4,atol=2e-5,state_tensors_equal=state_equal,default_flags=dict(matmul_tf32=torch.backends.cuda.matmul.allow_tf32,cudnn_tf32=torch.backends.cudnn.allow_tf32,cudnn_benchmark=torch.backends.cudnn.benchmark),same_batch_repeat=compare(repeat,full),same_batch_reload=compare(reload_full,full),different_batch_same_model=compare(half,full[:128]),saved_reference_cross_device=compare(full,saved),interpretation='Diagnostic only. No new fit or accepted physics weights; failed tasks1/2 require their own valid checkpoint evidence.')
torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
precise=evaluate(a,inputs);precise_half=evaluate(a,[v[:128] for v in inputs]);precise_reload=evaluate(b,inputs)
report['full_fp32_same_batch_reload']=compare(precise_reload,precise)
report['full_fp32_different_batch']=compare(precise_half,precise[:128])
report['default_vs_full_fp32']=compare(full,precise)
Path('parton-reload-diagnosis.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
