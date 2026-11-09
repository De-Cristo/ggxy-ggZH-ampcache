"""Matched C BDT: original classifier/tests, with the saved source folds retained."""
import json,tarfile
from pathlib import Path
from common import worker,configuration,arrays,scalar_samples,result
worker()
import numpy as np
import joblib
c=configuration();src,tgt=scalar_samples()
import ggzh_reweight_bdt as BDT
import ggzh_reweight_common as C
BDT.configure(C.FEATURE_SETS['C']);xs,xt=BDT.prepare(src),BDT.prepare(tgt)
weights={};models=Path('parton-bdt-models');models.mkdir()
for direction,tag in enumerate(c['tags']):
    tsel=tgt.in_window.copy()
    if direction:tsel &= tgt['target_half']==direction-1
    w=np.ones(src.n)
    for fold in range(5):
        train=src.in_window&(src['source_fold']!=fold);test=src.in_window&(src['source_fold']==fold)
        model=BDT.fit(xs[train],xt[tsel],tgt.w0[tsel]);p=model.predict_proba(xs[test])[:,1]
        assert ((p>0)&(p<1)).all(),'No silent probability clipping'
        w[test]=p/(1-p)
        joblib.dump(dict(model=model,fold=fold,tag=tag,selection_sha256=c['selection_sha256'],features=BDT.FEATURES),models/f'{tag}-{fold}.joblib')
        print(tag,fold,'real BDT fit completed',flush=True)
    w[src.in_window]/=w[src.in_window].mean();assert np.isfinite(w).all();weights[tag]=w
out=result(src,tgt,weights);out.update(method='matched17C-BDT',config=BDT.CLF,selection_sha256=c['selection_sha256'],package=Path('VERSION').read_text().strip(),saved_folds_retained=True)
Path('parton-bdt.json').write_text(C.dumps(out));np.savez_compressed('parton-bdt-weights.npz',run=src['run'],event=src['event'],**weights)
with tarfile.open('parton-bdt.tgz','w:gz') as t:
    for name in ('parton-bdt.json','parton-bdt-weights.npz','parton-bdt-models'):t.add(name,arcname=name)
