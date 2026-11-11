"""Paired real-event OOF comparison of original ParT and analysis-jet context.

Worker-only. Preserves registered scalar tests; additional conditional plots are
exploratory and do not define a new acceptance test. NLO nominal reference is
signed and includes training/validation events. No independent-target claim.
"""
import base64,html,json,sys,tarfile
from pathlib import Path
from common import worker,configuration,arrays,scalar_samples
worker()
import numpy as np
from scipy.special import logsumexp
c=configuration();assert c['context_variant']=='analysis-jets-v1'
from jet_context import SCHEMA
src,tgt=scalar_samples();data=arrays('lo-cms')
import ggzh_reweight_common as C
import ggzh_reweight_bdt as BDT
from ggzh_reweight_2d import shapes,test_vars,values
from ggzh_rwt10m import diagnostics
BDT.configure(C.FEATURE_SETS['C']);BDT.prepare(src);BDT.prepare(tgt)
CURRENT=Path('VERSION').read_text().split()[0]
BASE='d8219891f20b25220882dc5479538aa67a1549e6'
LEGACY='4205b4d0af70628c0d7fbd053911ff3c34aa09b9'
weights={'raw':np.ones(src.n)};records={};normalizations={}
for variant in ('baseline','jets'):
    logits=np.full(src.n,np.nan);records[variant]=[]
    for fold in range(5):
        directory=Path('compare-parts')/variant/str(fold);directory.mkdir(parents=True)
        with tarfile.open(f'{variant}/fit-{fold}.tgz') as t:
            assert all((m.name=='parton-fit' or m.name.startswith('parton-fit/')) and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers())
            t.extractall(directory,members=[m for m in t.getmembers() if Path(m.name).name in ('logits.npz','record.json')])
        p=directory/'parton-fit';r=json.loads((p/'record.json').read_text())
        assert r['task']==r['fold']==fold and r['tag']=='full'
        assert r['selection_sha256']==c['selection_sha256'] and r['checkpoint_reload']=='PASS'
        assert r['particle_attention_layers']==8 and not r['scalar_token_substitute']
        if variant=='jets':
            assert r['package'].split()[0]==CURRENT and r['context_variant']=='analysis-jets-v1'
            assert r['global_features']==35 and r['jet_context_schema']==SCHEMA
            assert r['initialization_seed']==c['seed_base']+fold
        else:
            allowed={BASE,LEGACY} if c['channel']=='zll' else {BASE}
            if c.get('rapidity_convention')=='transverse-mass-rapidity-v1':
                assert r.get('rapidity_convention')==c['rapidity_convention']
                allowed={CURRENT}
            # External comparisons may use a baseline from this exact code pin.
            # Other iteration pins require an explicit provenance review.
            if r.get('execution_context')=='external-gpu':allowed.add(CURRENT)
            assert r['package'].split()[0] in allowed
            assert r.get('context_variant','baseline24')=='baseline24'
        if c.get('rapidity_convention')=='transverse-mass-rapidity-v1':
            assert r.get('rapidity_convention')==c['rapidity_convention']
        with np.load(p/'logits.npz',allow_pickle=False) as z:
            ix=np.flatnonzero(src.in_window&(data['source_fold']==fold));assert np.array_equal(ix,z['source_index'])
            for key in ('seed','table_row','lhe_ordinal'):assert np.array_equal(z[key],data[key][ix])
            assert np.isfinite(z['logits']).all() and np.isnan(logits[ix]).all();logits[ix]=z['logits']
        records[variant].append(r)
    assert np.isfinite(logits[src.in_window]).all()
    norm=float(logsumexp(logits[src.in_window])-np.log(src.in_window.sum()))
    w=np.ones(src.n);w[src.in_window]=np.exp(logits[src.in_window]-norm)
    assert np.isfinite(w).all() and (w>0).all()
    weights[variant]=w;normalizations[variant]=norm
for a,b in zip(records['baseline'],records['jets']):
    for key in ('source_train','source_validation','source_test','target_train','target_validation','target_half'):
        assert a[key]==b[key],(key,a[key],b[key])
closure={v:BDT.evaluate(src,weights[v],src.in_window,tgt,tgt.in_window,v,print) for v in ('baseline','jets')}
summary={}
for variant in closure:
    r=closure[variant]['reweighted']
    summary[variant]=dict(criterion1=r['pass_1'],criterion2=r['pass_2'],scalar_17C_projections=r['pass_2b'],criterion3=None,
        criterion3_status='Not tested for this nominal paired comparison',n_eff=closure[variant]['n_eff'])
summary['raw']={k:closure['baseline']['unweighted'][v] for k,v in [('criterion1','pass_1'),('criterion2','pass_2'),('scalar_17C_projections','pass_2b')]}
diag={v:diagnostics(src,w,tgt) for v,w in weights.items()}
band=[tgt.weight(wid) for wid in C.SCALE7]
jetvars={k:({**C.CLOSURE_VARS,**C.EXTRA_VARS})[k] for k in ('ptj1','etaj1','ptj2','etaj2')}
ptedges=[0,75,150,250,400,600,np.inf]

def fraction(w,sel,full):
    den=float(w[full].sum());assert den>0
    return dict(fraction=float(w[sel].sum()/den),error=float(np.sqrt(np.square(w[sel]).sum())/den),
                raw=int(sel.sum()),sumw=float(w[sel].sum()),sumw2=float(np.square(w[sel]).sum()))

def conditional(name,ss,tt):
    definition={name:jetvars[name]};edges=jetvars[name][0]
    counts=[C.hist(values(s,name)[mask],np.ones(int(mask.sum())),edges)[2] for s,mask in ((src,ss),(tgt,tt))]
    reason=None
    if any(x.sum()<20 for x in counts):reason='Fewer than20 in-range raw events in source or target'
    if reason is None:
        for s,w,mask in [(src,w,ss) for w in weights.values()]+[(tgt,w,tt) for w in [tgt.w0,*band]]:
            if C.hist(values(s,name)[mask],w[mask],edges)[0].sum()<=0:
                reason='Nonpositive conditional normalization';break
    if reason:return dict(status='UNUSABLE',reason=reason,raw_source=int(counts[0].sum()),raw_target=int(counts[1].sum()))
    target=shapes(tgt,tgt.w0,tt,definition,band);out={}
    for label,w in weights.items():
        result,_=test_vars(shapes(src,w,ss,definition),target,definition);r=result[name]
        # No new acceptance threshold for these post-registration diagnostics.
        r.pop('pass');out[label]=r
    status='OK' if all(r['ndof']>0 and np.isfinite(r['pvalue']) for r in out.values()) else 'UNUSABLE'
    return dict(status=status,reason=None if status=='OK' else 'Insufficient merged bins for shape comparison',comparisons=out,
                raw_source=int(counts[0].sum()),raw_target=int(counts[1].sum()),normalization='Separate unit area within category AND plotted variable range')

categories=[]
for ip,(lo,hi) in enumerate(zip(ptedges[:-1],ptedges[1:])):
    for nj in range(4):
        ss=src.in_window&(src['ptz']>=lo)&(src['ptz']<hi)&(np.minimum(src['njets_isr'],3)==nj)
        tt=tgt.in_window&(tgt['ptz']>=lo)&(tgt['ptz']<hi)&(np.minimum(tgt['njets_isr'],3)==nj)
        target=fraction(tgt.w0,tt,tgt.in_window);fracs={v:fraction(w,ss,src.in_window) for v,w in weights.items()}
        for f in fracs.values():
            ft=target['fraction'];f['ratio_to_nlo']=f['fraction']/ft if ft>0 else None
            f['ratio_error']=float(np.sqrt((f['error']/ft)**2+(f['fraction']*target['error']/ft**2)**2)) if ft>0 else None
        varied=[target['fraction']]+[fraction(w,tt,tgt.in_window)['fraction'] for w in band]
        cat=dict(vpt_index=ip,vpt_range=[lo,hi],njets=nj,njets_label='>=3' if nj==3 else str(nj),target=target,
                 source=fracs,target_scale_band=[min(varied),max(varied)],
                 status='OK' if ss.sum()>=20 and tt.sum()>=20 and target['sumw']>0 else 'SPARSE_OR_NONPOSITIVE',shapes={})
        for name in jetvars:
            if nj<int(name[-1]):cat['shapes'][name]=dict(status='NOT_APPLICABLE',reason='Jet absent in this category')
            else:cat['shapes'][name]=conditional(name,ss,tt)
        categories.append(cat)
result=dict(study=c['study'],status='EXPLORATORY_PAIRED_NOMINAL_COMPARISON',channel=c['channel'],package=Path('VERSION').read_text().strip(),
    selection_sha256=c['selection_sha256'],context_schema=SCHEMA,summary=summary,full=closure,diagnostics=diag,categories=categories,
    fold_records=records,log_normalization=normalizations,
    interpretation='Same1M/class source-OOF comparison. Nominal NLO includes training/validation events. No independent target-half claim; no bootstrap/PCA uncertainty estimated. Signed target and existing fixed-normalization MC error conventions retained. Baseline and jet results are correlated, not independent experiments.',
    conditional_statistics='Legacy diagonal chi-square with nbins-1 degrees of freedom and raw-count merging>=20. Descriptive p-values, no multiple-testing correction, not independent-target goodness-of-fit certification.')
out=Path('parton-jet-report');out.mkdir();(out/'result.json').write_text(C.dumps(result))
compact={k:result[k] for k in ('study','status','channel','package','selection_sha256','context_schema','summary','diagnostics','interpretation')}
compact['training_packages']={v:sorted({r['package'] for r in rec}) for v,rec in records.items()}
Path('parton-jet-summary.json').write_text(C.dumps(compact))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path('tools').resolve()));import prd_style
prd_style.apply(plt)
labels={'raw':'LO+PS unreweighted','baseline':'LO+PS x original ParT','jets':'LO+PS x jet-context ParT'}
colors={'raw':'#eb6834','baseline':'#2477b5','jets':'#139760'}
images=[]
def save(fig,name):
    fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=140);fig.savefig(out/(name+'.pdf'));plt.close(fig);images.append(name)
def ratio_panel(ax,comparisons,title):
    first=next(iter(comparisons.values()));edges=np.asarray(first['edges'],float);x=(edges[:-1]+edges[1:])/2;t=np.asarray(first['f_target']);good=t>0
    ax.fill_between(x,np.divide(first['band_lo'],t,out=np.full_like(t,np.nan),where=good),np.divide(first['band_hi'],t,out=np.full_like(t,np.nan),where=good),alpha=.17,color='gray',label='NLO scale band')
    for v,r in comparisons.items():
        assert np.array_equal(r['edges'],first['edges']);y=np.asarray(r['f_test']);e=np.asarray(r['e_test']);et=np.asarray(r['e_target'])
        error=np.sqrt((e[good]/t[good])**2+(y[good]*et[good]/t[good]**2)**2)
        ax.errorbar(x[good],y[good]/t[good],yerr=error,fmt='.',color=colors[v],label=labels[v])
    ax.axhline(1,color='black',lw=.6);ax.set_title(title);ax.set_ylabel('/ NLO+PS');ax.set_ylim(0,2.5)
for name in ('ptz','pth','ptj1','lhe_mll'):
    def get(v,kind='reweighted'):
        r=closure[v][kind];return r['features' if name in r['features'] else 'variables'][name]
    fig,ax=plt.subplots(figsize=(prd_style.DOUBLE,3.8));ratio_panel(ax,{'raw':get('baseline','unweighted'),'baseline':get('baseline'),'jets':get('jets')},f'{c["channel"]}: {name}, all five source folds')
    ax.set_xlabel('m(nu nu) (GeV)' if c['channel']=='zvv' and name=='lhe_mll' else name);ax.legend();save(fig,'inclusive-'+name)
for ip in range(6):
    for name in jetvars:
        fig,axes=plt.subplots(2,2,figsize=(prd_style.DOUBLE,6.3));lo,hi=ptedges[ip:ip+2]
        for nj,ax in enumerate(axes.flat):
            cat=categories[4*ip+nj];r=cat['shapes'][name];title=f'{lo:g} <= VpT < {hi:g}; {cat["njets_label"]} jets'
            if 'comparisons' in r:ratio_panel(ax,r['comparisons'],title)
            else:ax.text(.5,.5,r['status'],ha='center',transform=ax.transAxes);ax.set_title(title)
            ax.set_xlabel(jetvars[name][1])
        handles,ls=next(((ax.get_legend_handles_labels()) for ax in axes.flat if ax.get_legend_handles_labels()[0]),([],[]))
        if handles:fig.legend(handles,ls,loc='upper center',ncol=2)
        fig.subplots_adjust(top=.9);save(fig,f'conditional-{name}-vpt{ip}')
body='<h1>Analysis-jet context ParT comparison: '+c['channel']+'</h1><p>'+html.escape(result['interpretation'])+'</p><pre>'+html.escape(json.dumps(summary,indent=2))+'</pre>'
body+='<p>Conditional shapes have unit area inside each category and plotted variable range. Category yields below use the full selected sample normalization. Error bars use fixed-normalization MC statistics; scale band is separate. Sparse cells are flagged. Plot ratio axes are limited to0–2.5; complete values are in embedded JSON.</p>'
body+='<details><summary>Category yield fractions and ratios</summary><pre>'+html.escape(C.dumps([{k:v for k,v in cat.items() if k!='shapes'} for cat in categories]))+'</pre></details>'
for name in images:
    encoded=base64.b64encode((out/(name+'.png')).read_bytes()).decode()
    body+='<h2>'+html.escape(name)+'</h2><img style="max-width:100%" width="900" src="data:image/png;base64,'+encoded+'">'
encoded=base64.b64encode((out/'result.json').read_bytes()).decode()
body+='<p><a download="result.json" href="data:application/json;base64,'+encoded+'">Download full results JSON</a></p>'
(out/'index.html').write_text('<!doctype html><html><meta charset="utf-8"><title>ParT analysis-jet comparison</title><body>'+body+'</body></html>')
with tarfile.open('parton-jet-report.tgz','w:gz') as t:t.add(out,arcname='jet-report')
print('REAL_PAIRED_JET_CONTEXT_COMPARISON',json.dumps(summary),flush=True)
