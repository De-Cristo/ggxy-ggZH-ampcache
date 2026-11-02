#!/usr/bin/env python3
"""Worker steps for the pre-registered rwt10m campaign. Shared region definitions come from ggzh_rwt_summary.
Physics failures are recorded as failures; an execution error is a separate failure of the job.
"""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from ggzh_reweight_common import Sample, SCALE7, FEATURE_SETS, dumps
from ggzh_rwt_summary import REGIONS, HARDJET_EDGES, frac_region
from ggzh_reweight_closure import test_pois
import ggzh_reweight_bdt as BDT


def variant():
    return 'R6' if os.environ.get('RWT_RETUNED') == '1' else 'default'


def aligned(path, src):
    with np.load(path) as z:
        assert np.array_equal(z['run'], src['run']) and np.array_equal(z['event'], src['event']), str(path)
        return z['w_full']


def ratio(src, w, tgt, var, a, b):
    fs, es = frac_region(src, w, src.in_window, var, a, b)
    ft, et = frac_region(tgt, tgt.w0, tgt.in_window, var, a, b)
    if ft == 0:
        return {'ratio': None, 'error': None, 'target_share': float(ft)}
    r = fs/ft
    e = np.sqrt((es/ft)**2 + (fs*et/ft**2)**2)
    return {'ratio': float(r), 'error': float(e), 'target_share': float(ft)}


def hardjet(s, w, sel):
    hard = s['ptj1'] > s['mzh']
    den = w[sel].sum()
    if den == 0:
        return {'share': None, 'error': None, 'raw_hard': int((sel & hard).sum()), 'raw_total': int(sel.sum())}
    return {'share': float(w[sel & hard].sum()/den),
            'error': float(np.sqrt(np.sum(w[sel & hard]**2))/abs(den)),
            'raw_hard': int((sel & hard).sum()), 'raw_total': int(sel.sum())}


def diagnostics(src, w, tgt):
    ss, ts = src.in_window, tgt.in_window
    out = {'regions': {f'{v}_{a:g}_{b:g}': ratio(src,w,tgt,v,a,b) for v,a,b in REGIONS},
           'stage13': test_pois(src,w,ss,tgt,tgt.w0,ts),
           'neff': float(w[ss].sum()**2/np.sum(w[ss]**2)), 'hardjet': {}}
    for label, sel_s, sel_t in [('all',ss,ts), ('ptz400',ss & (src['ptz']>400),ts & (tgt['ptz']>400)),
                              ('pth450',ss & (src['pth']>450),ts & (tgt['pth']>450))]:
        target = hardjet(tgt,tgt.w0,sel_t)
        varied = [hardjet(tgt,tgt.weight(wid),sel_t)['share'] for wid in SCALE7]
        shares = [x for x in [target['share']]+varied if x is not None]
        band = [min(shares),max(shares)] if shares else [None,None]
        source = hardjet(src,w,sel_s)
        out['hardjet'][label] = {'source':source,'target':target,'scale_band':band,
            'within_band': bool(source['share'] is not None and shares and band[0] <= source['share'] <= band[1])}
    for var in ('ptz','pth'):
        out['hardjet']['vs_'+var] = [
            {'range':[a,b], 'source':hardjet(src,w,ss & (src[var]>=a) & (src[var]<b)),
             'target':hardjet(tgt,tgt.w0,ts & (tgt[var]>=a) & (tgt[var]<b))}
            for a,b in zip(HARDJET_EDGES[:-1],HARDJET_EDGES[1:])]
    # Fixed mass bins of the original closure and STXS pT(Z) boundaries, including overflow.
    mass_edges = np.arange(15.,150.01,5.)
    pt_edges = np.array([0.,75.,150.,250.,400.,600.,np.inf])
    def h2(s,ww):
        sel=s.in_window
        return [np.histogram2d(s['lhe_mll'][sel],s['ptz'][sel],bins=(mass_edges,pt_edges),weights=z[sel])[0]
                for z in (ww,ww**2,np.ones(s.n))]
    hs,ht=h2(src,w),h2(tgt,tgt.w0)
    norm_s,norm_t=w[ss].sum(),tgt.w0[ts].sum()
    fs,ft=hs[0]/norm_s,ht[0]/norm_t
    with np.errstate(divide='ignore',invalid='ignore'):
        rr=fs/ft
        er=np.sqrt(hs[1]/norm_s**2/ft**2 + fs**2*(ht[1]/norm_t**2)/ft**4)
    out['mll_x_ptz']={'mass_edges':mass_edges,'ptz_edges':pt_edges,'ratio':rr,'error':er,
                      'raw_source':hs[2],'raw_target':ht[2]}
    return out


def summary(D,I):
    src,tgt=Sample('source',D/'lo-cms.npz'),Sample('target',D/'nlo-cp5.npz')
    out={'provenance':{'package':Path('VERSION').read_text().strip(),'study':'rwt10m',
                       'source_events':src.n,'target_events':tgt.n,'variant':variant()},
         'unweighted':diagnostics(src,np.ones(src.n),tgt),'sets':{}}
    for letter in 'ABC':
        label='clf-'+letter.lower()
        c=json.loads((I/f'bdt-{label}.json').read_text())
        assert c['features']==FEATURE_SETS[letter]
        assert c['classifier']==BDT.CLF, 'closure classifier settings do not match this variant'
        d=diagnostics(src,aligned(I/f'weights-{label}.npz',src),tgt)
        d['criteria']=c['summary']; d['features']=c['features']; d['classifier']=c['classifier']
        d['claim_reweighting_works']=bool(c['summary']['all_pass'] and d['hardjet']['pth450']['within_band'])
        out['sets'][letter]=d
    out['target_neff']={}
    for cut in (0,400,600):
        sel=tgt.in_window & (tgt['ptz']>=cut)
        out['target_neff'][str(cut)]=float(tgt.w0[sel].sum()**2/np.sum(tgt.w0[sel]**2))
    return out


def transfer(D,I):
    import joblib
    src,tgt=Sample('new',D/'lo-cms.npz'),Sample('new-target',D/'nlo-cp5.npz')
    old,oldt=Sample('old',D/'old-lo-cms.npz'),Sample('old-target',D/'old-nlo-cp5.npz')
    import tarfile
    with tarfile.open(I/'ensemble-clf-a.tgz') as f:
        if any(Path(x.name).is_absolute() or '..' in Path(x.name).parts for x in f.getmembers()):
            raise ValueError('unsafe ensemble archive')
        f.extractall('ensemble-transfer')
    models=sorted(Path('ensemble-transfer').rglob('fold-*.joblib'))
    assert len(models)==5
    xo=BDT.prepare(old,FEATURE_SETS['A'])
    # Match out-of-fold source weights: assign an independent old event to one of five trained models.
    wnew=np.ones(old.n); idx=np.flatnonzero(old.in_window)
    for k,path in enumerate(models):
        saved=joblib.load(path); assert saved['features']==FEATURE_SETS['A']
        ix=idx[np.arange(len(idx))%5==k]
        p=saved['model'].predict_proba(xo[ix])[:,1]; wnew[ix]=p/(1-p)
    wnew[old.in_window]/=wnew[old.in_window].mean()
    # The actual existing six-input all-event model (not a fresh retraining).
    from ggzh_rwt_export import onnx_p
    with tarfile.open(I/'model-baseline.tgz') as f:
        member = next(x for x in f.getmembers() if x.name.endswith('/clf6-nominal.onnx'))
        model = I/'baseline-clf6-nominal.onnx'
        model.write_bytes(f.extractfile(member).read())
    p=onnx_p(model.read_bytes(),BDT.prepare(src,FEATURE_SETS['A']))
    wold=np.where(src.in_window,p/(1-p),1.)
    wold[src.in_window]/=wold[src.in_window].mean()
    return {'old_model_on_new_source':diagnostics(src,wold,tgt),
            'new_model_on_old_source_vs_new_target':diagnostics(old,wnew,tgt),
            'new_model_on_old_source_vs_old_target':diagnostics(old,wnew,oldt),
            'provenance':{'old_model_sha256':__import__('hashlib').sha256(model.read_bytes()).hexdigest(),
                          'variant':variant(),'new_model_classifier':saved['classifier'],
                          'new_model':'full-target five-fold ensemble, assigned by independent source index'}}


def export(D,I):
    study=json.loads((I/'summary.json').read_text())
    passing=[k for k,v in study['sets'].items() if v['claim_reweighting_works']]
    assert study.get('provenance',{}).get('variant','default')==variant()
    outdir=Path('model'); outdir.mkdir()
    if passing:
        import ggzh_rwt_export as E
        E.MODELS={f'clf{k}':FEATURE_SETS[k] for k in passing}
        for k in passing:
            assert study['sets'][k]['classifier']==E.CLF, 'export settings differ from the closure test'
        extra={**__import__('ggzh_reweight_common').EXTRA_VARS,
               'lhe_mll':(None,'dilepton invariant mass','GeV')}
        for key,value in extra.items():
            E.INPUTS[key]=(value[1], 'GeV' if key.startswith('pt') or key=='lhe_mll' else '')
        rc=E.main(['--source',str(D/'lo-cms.npz'),'--target',str(D/'nlo-cp5.npz'),
                   '--out-dir',str(outdir)])
        if rc: raise RuntimeError('export checks failed')
        card=json.loads((outdir/'model-card.json').read_text())
        card['closure_tests']={k:study['sets'][k]['criteria'] for k in passing}
        card['validity']='Shape reweighting for the rwt10m CMS-style source; scope and measured closure attached.'
        (outdir/'model-card.json').write_text(dumps(card))
    else:
        (outdir/'model-card.json').write_text(dumps({'exported':False,'reason':'No set passed all pre-registered conditions',
                                                   'closure_tests':{k:v['criteria'] for k,v in study['sets'].items()}}))
    import tarfile
    with tarfile.open('model.tgz','w:gz') as f: f.add(outdir,arcname='model')
    return {'passing_sets':passing,'exported':bool(passing),'variant':variant(),'classifier':BDT.CLF}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('step',choices=['summary-10m','transfer-10m','model-10m'])
    ap.add_argument('--dir',type=Path,required=True); ap.add_argument('--inputs',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if variant()=='R6':
        BDT.CLF.update(max_iter=1000,max_leaf_nodes=63,min_samples_leaf=500)
    result={'summary-10m':summary,'transfer-10m':transfer,'model-10m':export}[a.step](a.dir,a.inputs)
    a.out.write_text(dumps(result))
    print('Wrote',a.out,flush=True)
if __name__=='__main__': main()
