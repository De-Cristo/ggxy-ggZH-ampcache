"""Summarize existing small fold report; no event/model reads or new fits."""
import argparse,hashlib,json,math
from pathlib import Path
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
assert a.input.stat().st_size<4*1024**2
raw=a.input.read_bytes();s=json.loads(raw);assert s['status']=='SINGLE_FOLD_EXPLORATORY'
a.out.mkdir(parents=True,exist_ok=False)
rows=[]
for name in s['diagnostics']['parton']['regions']:
 rows.append(dict(region=name,**{k:s['diagnostics'][k]['regions'][name] for k in ('unweighted','matched_bdt','parton')}))
n=s['source_test_events'];ne=s['parton']['n_eff']
r=dict(source_sha256=hashlib.sha256(raw).hexdigest(),source_test_events=n,target_events=s['target_events'],training_package=s['training_package'],report_package=s['package'],regions=rows,n_eff=ne,efficiency=ne/n,fixed_weight_total_mc_error_factor=math.sqrt(n/ne),criteria=s['parton_summary'],hardjet=s['diagnostics']['parton']['hardjet'],interpretation='Single source fold; NLO includes training/validation events. Region ratios are normalized fractions, errors are stored MC-only errors. Methods share events; no independent-method significance calculated. Total Neff comparison is not a binwise or fitted-parameter precision claim.')
(a.out/'evidence.json').write_text(json.dumps(r,indent=2)+'\n')
lines=['| Region (GeV for momenta/mass) | LO+PS | LO+PS × BDT | LO+PS × ParT |','|---|---:|---:|---:|']
for x in rows:lines.append('| '+x['region']+' | '+' | '.join(f"{x[k]['ratio']:.3f} ± {x[k]['error']:.3f}" for k in ('unweighted','matched_bdt','parton'))+' |')
(a.out/'region-table.md').write_text('\n'.join(lines)+'\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import prd_style
prd_style.apply(plt)
fig,ax=plt.subplots(figsize=(prd_style.DOUBLE,5.2))
for key,label,color,offset in [('unweighted','LO+PS','#d46c30',-.19),('matched_bdt','LO+PS × BDT','#5667bb',0),('parton','LO+PS × ParT','#118860',.19)]:
 ax.errorbar([x[key]['ratio'] for x in rows],[i+offset for i in range(len(rows))],xerr=[x[key]['error'] for x in rows],fmt='o',color=color,label=label)
ax.axvline(1,color='black',lw=.8);ax.set_yticks(range(len(rows)),[x['region'] for x in rows]);ax.invert_yaxis();ax.set_xlabel('Region fraction / NLO+PS region fraction');ax.set_title('Single held-out source fold; MC errors only');ax.legend(loc='lower right');fig.tight_layout()
fig.savefig(a.out/'region-comparison.png');fig.savefig(a.out/'region-comparison.pdf');plt.close(fig)
print(json.dumps({k:r[k] for k in ['source_sha256','n_eff','efficiency','fixed_weight_total_mc_error_factor']}))
