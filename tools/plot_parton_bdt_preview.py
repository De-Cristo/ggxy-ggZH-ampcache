"""Small prebinned-result plots only. No event arrays, model loads or fits."""
import argparse,html,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import prd_style
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
assert a.input.suffix=='.json' and a.input.stat().st_size<8*1024**2
r=json.loads(a.input.read_text());assert r['method']=='matched17C-BDT' and r['saved_folds_retained']
assert not a.out.exists();a.out.mkdir(parents=True)
prd_style.apply(plt);plots=[]
labels={'pth':r'$p_T(H)$ (GeV)','ptj1':r'$p_T(j_1)$ (GeV)','ptz':r'$p_T(Z)$ (GeV)','lhe_mll':r'$m_{\ell\ell}$ (GeV)'}
for name,label in labels.items():
    section='features' if name in r['full']['reweighted']['features'] else 'variables'
    v=r['full']['reweighted'][section][name];u=r['full']['unweighted'][section][name]
    edges=np.asarray(v['edges'],float)
    if np.isinf(edges[-1]):edges[-1]=edges[-2]+edges[-2]-edges[-3]
    x=(edges[:-1]+edges[1:])/2;target=np.asarray(v['f_target']);values=np.asarray(v['f_test']);errors=np.asarray(v['e_test']);flags=np.asarray(v['flags'],bool)
    fig,(ax,ratio)=plt.subplots(2,1,sharex=True,figsize=(prd_style.DOUBLE,4.5),gridspec_kw={'height_ratios':[3,1]})
    ax.fill_between(x,v['band_lo'],v['band_hi'],color='#2a78d6',alpha=.2,label='NLO scale band')
    ax.errorbar(x,target,yerr=v['e_target'],fmt='o',color='#2a78d6',label='NLO+PS')
    ax.errorbar(x,u['f_test'],yerr=u['e_test'],fmt='.',color='#eb6834',label='LO+PS')
    ax.errorbar(x,values,yerr=errors,fmt='.',color='#1baf7a',label='LO x BDT C')
    good=target!=0
    ratio.errorbar(x[good],values[good]/target[good],yerr=errors[good]/np.abs(target[good]),fmt='.',color='#1baf7a')
    bad=good&flags;ratio.plot(x[bad],values[bad]/target[bad],'x',color='#c22b36',label='outside band by >2 MC sigma')
    ratio.axhline(1,color='black',lw=.6);ratio.set_ylim(0,2);ratio.set_xlabel(label);ratio.set_ylabel('/ NLO')
    ax.set_ylabel('Normalized bin fraction');ax.set_title(f'Matched BDT control: p={v["pvalue"]:.3g}; {v["bins_outside_band_2sigma"]} bins outside band')
    if name!='lhe_mll':ax.set_yscale('log')
    ax.legend(ncol=2);ratio.legend(fontsize=8,loc='upper right');fig.tight_layout()
    fig.savefig(a.out/(name+'.png'));fig.savefig(a.out/(name+'.pdf'));plt.close(fig);plots.append(name+'.png')
summary=dict(method=r['method'],selection_sha256=r['selection_sha256'],criteria=r['summary'],hardjet_pth450=r['diagnostics']['hardjet']['pth450'],interpretation='BDT control only. No ParT closure conclusion before complete accepted OOF predictions.')
(a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
body='<h1>Matched BDT control on real Zll1M/class</h1><p>This is the scalar BDT control. ParT closure is pending.</p><pre>'+html.escape(json.dumps(summary,indent=2))+'</pre>'+''.join('<p><img width="900" src="'+name+'"></p>' for name in plots)
(a.out/'index.html').write_text('<!doctype html><html><meta charset="utf-8"><title>Matched Zll BDT control</title><body>'+body+'</body></html>')
print('Prebinned BDT preview written',a.out)
