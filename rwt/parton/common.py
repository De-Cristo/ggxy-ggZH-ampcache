"""Public real-particle study interfaces. All data work is worker-only."""
import hashlib,json,os,platform,sys
from pathlib import Path

def worker():
    assert platform.node().split('.')[0]!='cmscon' and os.environ.get('_CONDOR_SCRATCH_DIR'),'Condor worker required'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()

def configuration():
    c=json.loads(Path('rwt/parton/config-1m.json').read_text())
    assert sha('production/parton-zll-1m-v1.json')==c['selection_sha256']
    return c

def arrays(order):
    import numpy as np
    d=Path('parton-data')/order
    return {p.stem:np.load(p,mmap_mode='r',allow_pickle=False) for p in d.glob('*.npy')}

def scalar_samples():
    sys.path.insert(0,str(Path('rwt').resolve()))
    import ggzh_reweight_common as C
    return C.Sample('source',Path('parton-data/lo-cms.npz')),C.Sample('target',Path('parton-data/nlo-cp5.npz'))

def result(src,tgt,weights):
    import numpy as np
    import ggzh_reweight_bdt as BDT
    import ggzh_reweight_common as C
    from ggzh_rwt10m import diagnostics
    BDT.configure(C.FEATURE_SETS['C']);BDT.prepare(src);BDT.prepare(tgt)
    out={}
    for tag in ('full','A_to_B','B_to_A'):
        test=tgt.in_window.copy()
        if tag!='full':test &= np.asarray(tgt['target_half'])==(1 if tag=='A_to_B' else 0)
        out[tag]=BDT.evaluate(src,weights[tag],src.in_window,tgt,test,tag,print)
    r=out['full']['reweighted']
    out['summary']=dict(criterion1=r['pass_1'],criterion2=r['pass_2'],scalar_17C_projections=r['pass_2b'],criterion3=out['A_to_B']['pass'] and out['B_to_A']['pass'])
    out['summary']['all_scalar_tests_pass']=all(out['summary'].values())
    out['diagnostics']=diagnostics(src,weights['full'],tgt)
    out['interpretation']='Exploratory 1M per sample; scalar/STXS/half tests do not certify full10M acceptance or the full joint particle-input distribution'
    return out
