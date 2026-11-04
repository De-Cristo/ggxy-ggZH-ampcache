"""Real-record forward validation only; no fit and no synthetic input."""
import json, sys
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0,str(Path('rwt/transformer/vendor').resolve()))
from ParticleTransformer import ParticleTransformer
# Random initialization is only an interface check, never a training or physics result.
torch.set_num_threads(1)
report={}
for view in ('jets-v1','partons-v1'):
    with np.load(f'particle-out/{view}.npz') as d:
        x=torch.from_numpy(d['pf_features'][:8]);v=torch.from_numpy(d['pf_vectors'][:8]);m=torch.from_numpy(d['pf_mask'][:8]).bool()
    model=ParticleTransformer(input_dim=x.shape[1],num_classes=1,embed_dims=[64],pair_input_dim=4,
        num_heads=4,num_layers=2,num_cls_layers=2,trim=False)
    model.eval()
    with torch.no_grad():y=model(x,v=v,mask=m)
    assert y.shape==(len(x),1) and torch.isfinite(y).all()
    report[view]={'real_events':len(x),'output_shape':list(y.shape),'finite':True}
p=Path('particle-out/validation.json');s=json.loads(p.read_text());s['forward_validation']=report;p.write_text(json.dumps(s,indent=2)+'\n')
