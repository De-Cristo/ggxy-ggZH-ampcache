"""One BDT-17 ensemble member (feature set C, R6 hyperparameters) on a frozen 1M/class subset, CPU only.

Inputs: parton-data/ (the member subset), calib/parton-data/ and eval/parton-data/ (disjoint sets); only
their scalar tables (lo-cms.npz, nlo-cp5.npz) and selection/audit files are read.  The classifier is the
one of ggzh_reweight_bdt.py with --retuned (R6): balanced classes, every LO event weight 1, NLO events with
their nominal weight, negative weights zero in the loss.  It predicts logit = log(p/(1-p)) for every
calibration and evaluation row; the output layout is the one ensemble_assess.py reads.
"""
import json,os,pickle,sys,tarfile,time
from pathlib import Path
from common import worker,configuration,sha
worker()
import numpy as np
sys.path.insert(0,str(Path('rwt').resolve()))
import ggzh_reweight_common as C
import ggzh_reweight_bdt as BDT

c=configuration();seed=c['initialization_seed'];assert isinstance(seed,int)
BDT.configure(C.FEATURE_SETS['C']);BDT.CLF.update(max_iter=1000,max_leaf_nodes=63,min_samples_leaf=500,random_state=seed)
FEATURES=list(BDT.FEATURES);assert len(FEATURES)==17

def scalar(root,label):
    return C.Sample(label,Path(root)/'parton-data'/f'{label}.npz')
src,tgt=scalar('.','lo-cms'),scalar('.','nlo-cp5')
assert src.n==tgt.n==1000000
mine={o:{x['seed'] for x in json.loads(Path('parton-data/selection.json').read_text())['samples'][o]} for o in ('lo-cms','nlo')}
held={}
for name in ('calib','eval'):
    other=json.loads(Path(name,'parton-data','selection.json').read_text())
    for o in ('lo-cms','nlo'):assert not mine[o]&{x['seed'] for x in other['samples'][o]},f'{name} {o} overlaps the member subset'
    held[name,'lo-cms']=scalar(name,'lo-cms');held[name,'nlo']=scalar(name,'nlo-cp5')

start=time.time()
Xs=BDT.prepare(src)[src.in_window];Xt=BDT.prepare(tgt)[tgt.in_window];wt=tgt.w0[tgt.in_window]
assert np.isfinite(Xs).all() and np.isfinite(Xt).all()
clf=BDT.fit(Xs,Xt,wt);training_seconds=time.time()-start

def logits(model,s):
    p=model.predict_proba(BDT.prepare(s))[:,1];p=np.clip(p,1e-15,1-1e-15)
    out=np.log(p)-np.log1p(-p);assert np.isfinite(out).all();return out
out=Path('member-fit');assert not out.exists();out.mkdir()
predictions={}
for (name,order),s in held.items():
    predictions[name,order]=logits(clf,s)
    np.savez_compressed(out/f'logits-{name}-{order}.npz',seed=np.asarray(s['run']),lhe_ordinal=np.asarray(s['lhe_ordinal']) if 'lhe_ordinal' in s.cols else np.full(s.n,-1),
        table_row=np.asarray(s['event']),in_window=s.in_window,logits=predictions[name,order])
with (out/'model.pkl').open('wb') as f:pickle.dump(dict(model=clf,features=FEATURES,classifier=BDT.CLF,member=c['study'],initialization_seed=seed),f)
with (out/'model.pkl').open('rb') as f:reloaded=pickle.load(f)['model']
fresh=logits(reloaded,held['calib','lo-cms']);difference=float(np.max(np.abs(fresh-predictions['calib','lo-cms'])))
reload_audit=dict(protocol='pickle-reload-calibration-lo-v1',events=held['calib','lo-cms'].n,max_absolute_difference=difference,same_batch_pass=difference==0.)
assert reload_audit['same_batch_pass'],reload_audit
(out/'reload-audit.json').write_text(json.dumps(reload_audit,indent=2)+'\n')
def selection_sha(name):
    return json.loads(Path(name,'parton-data','audit.json').read_text())['selection_sha256']
import sklearn
record=dict(member=c['study'],model='BDT-17 (feature set C), R6 hyperparameters',features=FEATURES,classifier=BDT.CLF,channel=c['channel'],
    selection_sha256=c['selection_sha256'],calibration_selection_sha256=selection_sha('calib'),evaluation_selection_sha256=selection_sha('eval'),
    package=Path('VERSION').read_text().strip(),initialization_seed=seed,source_train=int(len(Xs)),target_train=int(len(Xt)),
    n_iter=int(clf.n_iter_),capacity_limited=bool(clf.n_iter_>=BDT.CLF['max_iter']),epochs_run=int(clf.n_iter_),
    best_validation_loss=float(-clf.validation_score_[-1]) if getattr(clf,'validation_score_',None) is not None and len(clf.validation_score_) else None,
    training_seconds=training_seconds,checkpoint_reload='PASS',reload_audit=reload_audit,model_sha256=sha(out/'model.pkl'),
    sklearn=sklearn.__version__,threads=os.environ.get('OMP_NUM_THREADS'))
(out/'record.json').write_text(json.dumps(record,indent=2)+'\n');Path('bdt-member.json').write_text(json.dumps(record,indent=2)+'\n')
with tarfile.open('bdt-member.tgz','w:gz') as t:t.add(out,arcname='member-fit')
print('BDT_MEMBER_COMPLETED',c['study'],'n_iter',clf.n_iter_,flush=True)
