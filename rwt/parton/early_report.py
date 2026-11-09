"""Early nominal closure from all five real OOF particle fits. Stability remains pending."""
import html,json,sys,tarfile
from pathlib import Path
from common import worker,configuration,arrays,scalar_samples
worker()
import numpy as np
from scipy.special import logsumexp
c=configuration();src,tgt=scalar_samples();data=arrays('lo-cms')
import ggzh_reweight_common as C
import ggzh_reweight_bdt as BDT
from ggzh_rwt10m import diagnostics
BDT.configure(C.FEATURE_SETS['C']);BDT.prepare(src);BDT.prepare(tgt)
logits=np.full(src.n,np.nan);records=[]
for fold in range(5):
    directory=Path('early-parts')/str(fold);directory.mkdir(parents=True)
    with tarfile.open(f'parts/fit-{fold}.tgz') as t:
        assert all((m.name=='parton-fit' or m.name.startswith('parton-fit/')) and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers())
        # Models remain in their original archives; only real fit records/logits are needed.
        t.extractall(directory,members=[m for m in t.getmembers() if Path(m.name).name in ('logits.npz','record.json')])
    p=directory/'parton-fit';record=json.loads((p/'record.json').read_text())
    assert record['task']==fold and record['fold']==fold and record['tag']=='full'
    assert record['selection_sha256']==c['selection_sha256'] and record['checkpoint_reload']=='PASS'
    assert record['package'].split()[0]=='4205b4d0af70628c0d7fbd053911ff3c34aa09b9'
    with np.load(p/'logits.npz',allow_pickle=False) as z:
        ix=np.flatnonzero(src.in_window&(data['source_fold']==fold));assert np.array_equal(ix,z['source_index'])
        for key in ('seed','table_row','lhe_ordinal'):assert np.array_equal(z[key],data[key][ix])
        assert np.isfinite(z['logits']).all() and np.isnan(logits[ix]).all();logits[ix]=z['logits']
    records.append(record)
assert np.isfinite(logits[src.in_window]).all()
normalization=float(logsumexp(logits[src.in_window])-np.log(src.in_window.sum()));weights=np.ones(src.n);weights[src.in_window]=np.exp(logits[src.in_window]-normalization)
assert np.isfinite(weights).all() and (weights>0).all()
closure=BDT.evaluate(src,weights,src.in_window,tgt,tgt.in_window,'early full nominal',print);r=closure['reweighted']
result=dict(study=c['study'],status='INTERMEDIATE_NOMINAL_CLOSURE',selection_sha256=c['selection_sha256'],package=Path('VERSION').read_text().strip(),training_package=records[0]['package'],summary=dict(criterion1=r['pass_1'],criterion2=r['pass_2'],scalar_17C_projections=r['pass_2b'],criterion3=None,criterion3_status='PENDING independent target-half fits; no final success claim'),full=closure,diagnostics=diagnostics(src,weights,tgt),fold_records=records,log_normalization=normalization,interpretation='All five nominal OOF source folds included. Nominal exploratory closure is assessable; target-half stability and full-inventory acceptance are not yet established.')
out=Path('parton-early-report');out.mkdir();(out/'result.json').write_text(C.dumps(result));Path('parton-early-summary.json').write_text(C.dumps(result))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path('tools').resolve()));import prd_style
prd_style.apply(plt)
plots=[]
for name in ('pth','ptj1','ptz','lhe_mll'):
    section='features' if name in r['features'] else 'variables';v=r[section][name];u=closure['unweighted'][section][name]
    edges=np.asarray(v['edges'],float)
    if np.isinf(edges[-1]):edges[-1]=edges[-2]+(edges[-2]-edges[-3])
    x=(edges[:-1]+edges[1:])/2;target=np.asarray(v['f_target'])
    fig,(ax,ratio)=plt.subplots(2,1,sharex=True,figsize=(prd_style.DOUBLE,4.5),gridspec_kw={'height_ratios':[3,1]})
    ax.fill_between(x,v['band_lo'],v['band_hi'],color='#2a78d6',alpha=.2,label='NLO scale band')
    ax.errorbar(x,target,yerr=v['e_target'],fmt='o',color='#2a78d6',label='NLO+PS')
    for label,values,errors,color in [('LO+PS',u['f_test'],u['e_test'],'#eb6834'),('LO x ParT',v['f_test'],v['e_test'],'#1baf7a')]:
        ax.errorbar(x,values,yerr=errors,fmt='.',color=color,label=label);good=target!=0
        ratio.errorbar(x[good],np.asarray(values)[good]/target[good],yerr=np.asarray(errors)[good]/np.abs(target[good]),fmt='.',color=color)
    ax.set_ylabel('Normalized bin fraction');ax.set_title(f'Early nominal 1M/class: {name}; p={v["pvalue"]:.3g}; outside band={v["bins_outside_band_2sigma"]}')
    if name!='lhe_mll':ax.set_yscale('log')
    ax.legend();ratio.axhline(1,color='black',lw=.6);ratio.set_xlabel(name);ratio.set_ylabel('/ NLO');ratio.set_ylim(0,2);fig.tight_layout()
    fig.savefig(out/(name+'.png'));fig.savefig(out/(name+'.pdf'));plt.close(fig);plots.append(name+'.png')
body='<h1>Early nominal particle ParT closure: 1M LO +1M NLO</h1><p>All five OOF source folds; independent target-half stability tests are PENDING. This is not the final physics conclusion.</p><pre>'+html.escape(json.dumps(result['summary'],indent=2))+'</pre><p><a href="result.json">Complete results and diagnostics</a></p>'+''.join('<p><img width="900" src="'+p+'"></p>' for p in plots)
(out/'index.html').write_text('<!doctype html><html><meta charset="utf-8"><title>Early nominal ZH ParT closure</title><body>'+body+'</body></html>')
with tarfile.open('parton-early-report.tgz','w:gz') as t:t.add(out,arcname='early-report')
print('EARLY_REAL_NOMINAL_CLOSURE',json.dumps(result['summary']),flush=True)
