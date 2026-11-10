"""Exploratory performance on the accepted ParT fold0 held-out source, not full OOF/stability."""
import html,json,sys,tarfile
from pathlib import Path
from common import worker,configuration,scalar_samples
worker()
import numpy as np
from scipy.special import logsumexp
c=configuration();source,target=scalar_samples()
import ggzh_reweight_common as C
import ggzh_reweight_bdt as BDT
from ggzh_rwt10m import diagnostics,ratio
for archive,dest,allowed in [('reference-fit.tgz','fold0',{'record.json','logits.npz'}),('reference-bdt.tgz','control',{'parton-bdt.json','parton-bdt-weights.npz'})]:
    with tarfile.open(archive) as t:
        assert all(not m.name.startswith('/') and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers())
        t.extractall(dest,members=[m for m in t.getmembers() if Path(m.name).name in allowed])
r=json.loads(Path('fold0/parton-fit/record.json').read_text())
assert r['task']==0 and r['tag']=='full' and r['fold']==0 and r['checkpoint_reload']=='PASS'
assert r['selection_sha256']==c['selection_sha256']
with np.load('fold0/parton-fit/logits.npz',allow_pickle=False) as z:
    ix=np.array(z['source_index']);expected=np.flatnonzero(source.in_window&(source['source_fold']==0))
    assert np.array_equal(ix,expected) and len(np.unique(ix))==len(ix)
    assert np.array_equal(z['seed'],source['run'][ix]) and np.array_equal(z['table_row'],source['event'][ix])
    logits=np.array(z['logits']);assert np.isfinite(logits).all()
normalization=float(logsumexp(logits)-np.log(len(logits)));wp=np.exp(logits-normalization)
assert np.isfinite(wp).all() and (wp>0).all()
with np.load('control/parton-bdt-weights.npz',allow_pickle=False) as z:
    assert np.array_equal(z['run'],source['run']) and np.array_equal(z['event'],source['event'])
    wb=np.array(z['full'][ix]);assert np.isfinite(wb).all() and (wb>0).all();wb/=wb.mean()
control=json.loads(Path('control/parton-bdt.json').read_text());assert control['selection_sha256']==c['selection_sha256'] and control['saved_folds_retained']
heldout=source.resample(ix);assert (heldout['source_fold']==0).all()
BDT.configure(C.FEATURE_SETS['C']);BDT.prepare(heldout);BDT.prepare(target)
part=BDT.evaluate(heldout,wp,heldout.in_window,target,target.in_window,'ParT source fold0 only',print)
base=BDT.evaluate(heldout,wb,heldout.in_window,target,target.in_window,'matched BDT source fold0 only',print)
def summary(x):
    d=x['reweighted'];return dict(criterion1=d['pass_1'],criterion2=d['pass_2'],scalar_17C_projections=d['pass_2b'],nominal_scalar_tests_pass=d['pass_1'] and d['pass_2'] and d['pass_2b'],criterion3=None)
result=dict(study=c['study'],status='SINGLE_FOLD_EXPLORATORY',selection_sha256=c['selection_sha256'],package=Path('VERSION').read_text().strip(),training_package=r['package'],source_test_events=heldout.n,target_events=target.n,source_holdout_fold=0,target_scope='Full nominal NLO pool includes gradient-training and early-stopping events; independent target-half stability not yet evaluated',parton_summary=summary(part),matched_bdt_summary=summary(base),parton=part,matched_bdt=base,log_normalization=normalization,diagnostics=dict(parton=diagnostics(heldout,wp,target),matched_bdt=diagnostics(heldout,wb,target),unweighted=diagnostics(heldout,np.ones(heldout.n),target)),mll_15_40=dict(parton=ratio(heldout,wp,target,'lhe_mll',15,40),matched_bdt=ratio(heldout,wb,target,'lhe_mll',15,40)),interpretation='Preliminary performance on one excluded LO fold only. This does not replace complete five-fold OOF results or both target-half tests; a pass/fail here can change with full sample statistics.',checkpoint_reload='PASS')
out=Path('parton-fold0-report');out.mkdir();(out/'result.json').write_text(C.dumps(result));Path('parton-fold0-summary.json').write_text(C.dumps(result))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path('tools').resolve()));import prd_style
prd_style.apply(plt);plots=[]
labels={'pth':r'$p_T(H)$ (GeV)','ptj1':r'$p_T(j_1)$ (GeV)','ptz':r'$p_T(Z)$ (GeV)','lhe_mll':r'$m_{\ell\ell}$ (GeV)'}
for name,label in labels.items():
    section='features' if name in part['reweighted']['features'] else 'variables';v=part['reweighted'][section][name];b=base['reweighted'][section][name];u=part['unweighted'][section][name]
    edges=np.asarray(v['edges'],float)
    if np.isinf(edges[-1]):edges[-1]=edges[-2]+edges[-2]-edges[-3]
    x=(edges[:-1]+edges[1:])/2;ft=np.asarray(v['f_target'])
    fig,(ax,ra)=plt.subplots(2,1,sharex=True,figsize=(prd_style.DOUBLE,4.7),gridspec_kw={'height_ratios':[3,1]})
    ax.fill_between(x,v['band_lo'],v['band_hi'],color='#2a78d6',alpha=.2,label='NLO scale band');ax.errorbar(x,ft,yerr=v['e_target'],fmt='o',color='#2a78d6',label='NLO+PS')
    for text,d,color in [('LO raw',u,'#eb6834'),('LO x ParT fold0',v,'#1baf7a'),('LO x BDT same fold',b,'#9a51a0')]:
        de=np.asarray(d['edges'],float)
        if np.isinf(de[-1]):de[-1]=de[-2]+de[-2]-de[-3]
        dx=(de[:-1]+de[1:])/2;df=np.asarray(d['f_test']);dt=np.asarray(d['f_target']);good=dt!=0
        ax.errorbar(dx,df,yerr=d['e_test'],fmt='.',color=color,label=text);ra.errorbar(dx[good],df[good]/dt[good],yerr=np.asarray(d['e_test'])[good]/np.abs(dt[good]),fmt='.',color=color)
    bad=np.asarray(v['flags'],bool)&(ft!=0);ra.plot(x[bad],np.asarray(v['f_test'])[bad]/ft[bad],'x',color='#c22b36')
    if name!='lhe_mll':ax.set_yscale('log')
    ax.set_ylabel('Normalized bin fraction');ax.set_title(f'Single source fold: {heldout.n:,} LO vs {target.n:,} NLO; ParT p={v["pvalue"]:.3g}')
    ax.legend(ncol=2);ra.axhline(1,color='black',lw=.6);ra.set_ylim(0,2);ra.set_xlabel(label);ra.set_ylabel('/ NLO');fig.tight_layout()
    fig.savefig(out/(name+'.png'));fig.savefig(out/(name+'.pdf'));plt.close(fig);plots.append(name+'.png')
rows=[]
for section in ('features','variables'):
    for name,v in part['reweighted'][section].items():rows.append(f'<tr><td>{html.escape(name)}</td><td>{v["pvalue"]:.5g}</td><td>{v["bins_outside_band_2sigma"]}</td><td>{"PASS" if v["pass"] else "FAIL"}</td></tr>')
body=f'<h1>Single-fold ParT performance: {heldout.n:,} held-out LO events</h1><p>{html.escape(result["target_scope"])}. Full OOF and criterion3 remain pending.</p><h2>ParT</h2><pre>'+html.escape(json.dumps(result['parton_summary'],indent=2))+'</pre><h2>BDT on the same events</h2><pre>'+html.escape(json.dumps(result['matched_bdt_summary'],indent=2))+'</pre><p><a href="result.json">Complete results and diagnostics</a></p><table border="1"><tr><th>Variable</th><th>p</th><th>Bins outside band</th><th>Pass</th></tr>'+''.join(rows)+'</table>'+''.join('<p><img width="900" src="'+p+'"></p>' for p in plots)
(out/'index.html').write_text('<!doctype html><html><meta charset="utf-8"><title>Single-fold ParT exploratory performance</title><body>'+body+'</body></html>')
with tarfile.open('parton-fold0-report.tgz','w:gz') as t:t.add(out,arcname='fold0-report')
print('SINGLE_FOLD_PRELIMINARY_PERFORMANCE',json.dumps(result['parton_summary']),flush=True)
