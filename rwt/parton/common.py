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
    channel=os.environ.get('RWT_PARTON_CHANNEL','zll')
    assert channel in ('zll','zvv')
    name='config-1m.json' if channel=='zll' else 'config-zvv-1m.json'
    c=json.loads((Path('rwt/parton')/name).read_text())
    assert sha(f'production/parton-{channel}-1m-v1.json')==c['selection_sha256']
    c['channel']=channel
    variant=os.environ.get('RWT_PARTON_CONTEXT','baseline24')
    assert variant in ('baseline24','analysis-jets-v1')
    if variant=='analysis-jets-v1':
        from jet_context import SCHEMA
        c['study']=f'parton-{channel}-jets-1m-v1'
        c['context_variant']=variant
        c['jet_context_schema']=SCHEMA
        c['model']['global_features']=35
        c['normalization']='training partitions only; masked particles; global35; physical four-vectors unchanged'
    return c

def accepted_training_package(record, config):
    current=Path('VERSION').read_text().split()[0]
    pin=record['package'].split()[0]
    legacy='4205b4d0af70628c0d7fbd053911ff3c34aa09b9'
    return pin==current or (config['channel']=='zll' and pin==legacy)

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
    out['channel']=configuration()['channel']
    out['decay_columns']='Neutrino-pair mass/neutrino kinematics' if out['channel']=='zvv' else 'Dilepton mass/lepton kinematics'
    out['interpretation']='Exploratory 1M per sample; scalar/STXS/half tests do not certify full10M acceptance or the full joint particle-input distribution'
    return out
