"""Deliverable of a BDT-17 ensemble: ONNX members, calibration tables, apply/PCA scripts, test vectors, model card.

Inputs: tables.json (bdt-assess), pdf-transport.npz, members/member-XX.tgz (bdt-fit), eval/parton-data/ (test
vectors).  The ONNX thresholds are snapped to float32 first (ggzh_rwt_export.snap_thresholds), so scikit-learn
and ONNX take identical branches for float32 inputs; the test vectors check apply.py against scikit-learn.
"""
import json,pickle,shutil,sys,tarfile
from pathlib import Path
from common import worker,configuration,sha
worker()
import numpy as np
sys.path.insert(0,str(Path('rwt').resolve()))
import ggzh_reweight_common as C
import ggzh_reweight_bdt as BDT
from ggzh_rwt_export import to_onnx,snap_thresholds
from ensemble_assess_bins import MZH_EDGES,RHO_EDGES

c=configuration();ch=c['channel'];BDT.configure(C.FEATURE_SETS['C'])
tables=json.loads(Path('tables.json').read_text());names=tables['members'];T=tables['tables']
out=Path('deliver');assert not out.exists();(out/'members').mkdir(parents=True)
models={}
for n in names:
    with tarfile.open(f'members/{n}.tgz') as t:
        assert all('..' not in m.name and not m.issym() and not m.islnk() for m in t.getmembers());t.extractall(Path('members-x')/n)
    with (Path('members-x')/n/'member-fit'/'model.pkl').open('rb') as f:card=pickle.load(f)
    clf=card['model'];features=card['features'];snap_thresholds(clf);models[n]=clf
    (out/'members'/f'{n}.onnx').write_bytes(to_onnx(clf,len(features)).SerializeToString())
assert features==list(BDT.FEATURES),(features,BDT.FEATURES)
jet=T['jet_nonclosure'];ngroups_rho=len(jet['c']);nom=T['nominal']
cell_of=lambda m,r:np.asarray(nom['mzh_lookup'])[m]*ngroups_rho+np.asarray(jet['rho_lookup'])[r]
ncell=len(nom['k'])*ngroups_rho
# PDF4LHC21 transport, merged into the same (m_ZH group x rho group) cells as the scale transport.
d=np.load('pdf-transport.npz');sums=d['sums'];assert np.allclose(d['mzh_edges'],MZH_EDGES) and np.allclose(d['rho_edges'],RHO_EDGES)
mi,ri=np.meshgrid(np.arange(sums.shape[1]),np.arange(sums.shape[2]),indexing='ij');cells=cell_of(mi.ravel(),ri.ravel())
merged=np.stack([np.bincount(cells,weights=s.ravel(),minlength=ncell) for s in sums])
raw=np.bincount(cells,weights=d['raw'].ravel(),minlength=ncell)
# Cells with fewer than 50 events of the 5e5 sample take the ratio of their m_ZH group summed over rho, then
# the inclusive ratio; the number of such cells is recorded.
group=np.arange(ncell)//ngroups_rho;gsum=np.stack([np.bincount(group,weights=m_,minlength=len(nom['k'])) for m_ in merged])
graw=np.bincount(group,weights=raw,minlength=len(nom['k']));ratio=np.empty((42,ncell));fallback=dict(group=0,inclusive=0)
for i in range(ncell):
    if raw[i]>=50 and merged[0,i]>0:ratio[:,i]=merged[1:,i]/merged[0,i]
    elif graw[group[i]]>=50 and gsum[0,group[i]]>0:ratio[:,i]=gsum[1:,group[i]]/gsum[0,group[i]];fallback['group']+=1
    else:ratio[:,i]=merged[1:].sum(1)/merged[0].sum();fallback['inclusive']+=1
pdf={(f'pdf_{k-1:02d}' if k<=40 else ('alphas_dn' if k==41 else 'alphas_up')):ratio[k-1].tolist() for k in range(1,43)}
scale={'scale_{:g}_{:g}'.format(*C.MU_OF_ID[w]):v for w,v in T['scale_transport']['T'].items()}
member_tables=[dict(name=n,onnx=f'members/{n}.onnx',**{k:T[n][k] for k in ('offset','norm','k','mzh_lookup')}) for n in names]
calibration=dict(schema='ggzh-bdt17-ensemble-v1',channel=ch,features=features,window=[15.,150.],
    fine_mzh_edges=[float(x) for x in MZH_EDGES],fine_rho_edges=[float(x) for x in RHO_EDGES],
    nominal={k:nom[k] for k in ('offset','norm','k','mzh_lookup')},members=member_tables,
    jet=dict(rho_lookup=jet['rho_lookup'],c=jet['c'],rho_group_edges=jet['rho_edges']),
    scale=dict(T=scale,source=T['scale_transport']['note']),
    pdf=dict(T=pdf,fallback_cells=fallback,source='PDF4LHC21_40_pdfas members 1-40 (Hessian) and 41/42 (alpha_s), ratios to member 0 from the 5e5-event 13.6 TeV Z->ll NLO+PS production; same cells as scale'),
    cells='nominal m_ZH group (mzh_lookup over fine_mzh_edges) x rho group (jet.rho_lookup over fine rho cells; cell 0 = no jet)')
(out/'calibration.json').write_text(json.dumps(calibration,indent=1)+'\n')
for f in ('apply.py','pca.py'):shutil.copy(Path('rwt/parton/deliver')/f,out/f)
# Test vectors: 2000 in-window evaluation LO+PS events; expected weights from scikit-learn.
s=C.Sample('lo-cms',Path('eval/parton-data/lo-cms.npz'));X=BDT.prepare(s).astype(np.float32);ix=np.flatnonzero(s.in_window)[:2000]
ev={k:np.asarray(X[ix,i],np.float32) for i,k in enumerate(features)}
sys.path.insert(0,str(out));import apply as A
model=A.Model(out);got=model.weights(ev)
Lsk=[]
for n in names:
    p=np.clip(models[n].predict_proba(X[ix].astype(np.float64))[:,1],1e-15,1-1e-15);Lsk.append(np.log(p)-np.log1p(-p))
Lsk=np.asarray(Lsk);Lonnx=model.logits(model.matrix(ev))
check=dict(events=len(ix),max_abs_logit_difference=float(np.max(np.abs(Lsk-Lonnx))),pdf_fallback_cells=fallback,cells=ncell)
assert check['max_abs_logit_difference']<1e-4,check
np.savez_compressed(out/'test-vectors.npz',**ev,**{f'expected_{k}':v for k,v in got.items()})
card=dict(model='gg->ZH LO+PS -> NLO+PS shape weights, BDT-17 ensemble (20 members)',channel=ch,
    source='CMS-style LO+PS 13.6 TeV: POWHEG ggHZ (ggxy_ggZH, bornonly) NNPDF3.1 NNLO, mu_R = mu_F = mZ + mH, SCALUP = m_ZH, Pythia 8 CP5 + POWHEG veto at SCALUP, 15 < m_ll < 150 GeV',
    target='NLO+PS 13.6 TeV: ggxy_ggZH POWHEG NLO (Born suppression 250 GeV, hdamp 250 GeV), NNPDF3.1 central member as nominal, Pythia 8 CP5 + POWHEG veto',
    training='each member: HistGradientBoostingClassifier (feature set C, R6 hyperparameters) on 1e6 LO+PS and 1e6 NLO+PS events drawn from 1e7 per class; calibration and evaluation on two further disjoint 1e6-event sets',
    inputs={k:BDT.FEATURE_VARS.get(k,BDT.TEST_VARS.get(k,('', '', '')))[1] for k in features},
    input_conventions='generator level after the shower; Z, H = last copies; jets anti-kT 0.4 of all stable particles not from Z or H decays, pT > 30 GeV; ptj1/etaj1/drj1zh refer to the leading jet with |eta| < 4.7 and take the values of the training tables when absent (see test-vectors.npz); njets = min(N_jets(pT>30, any eta), 3); ayzh = |y(ZH)|; lhe_mll = Z mass (last copy); for Z->nunu the lepton variables are the neutrinos',
    validity='15 < m_ll < 150 GeV (NaN outside). Jet-inclusive observables. The LO+PS source has almost no events with a leading jet harder than m_ZH (POWHEG veto at SCALUP = m_ZH); the weights do not create them. The jet_up/jet_down variation covers this region; analyses that select or veto hard jets at high pT should use NLO+PS events.',
    uncertainties=dict(nn='fill each member_XX weight; pca.py -> one nuisance per retained mode',jet='jet_up / jet_down, two-sided',
        scale='scale_* weights, 7-point envelope with the nominal', pdf='PDF4LHC21 Hessian: sqrt(sum_k (w_pdf_k - w_nominal)^2); alpha_s: (alphas_up - alphas_dn)/2',
        combination='independent sources; do not multiply the widths'),
    normalization='weights have mean 1 on the training LO+PS sample inside the window; scale the total to sigma(NLO) of the model card note',
    checks=check,files={p.name:sha(p) for p in sorted(out.rglob('*')) if p.is_file()})
(out/'model-card.json').write_text(json.dumps(card,indent=1)+'\n')
with tarfile.open('bdt-deliver.tgz','w:gz') as t:t.add(out,arcname=f'ggzh-bdt17-{ch}')
Path('bdt-deliver.json').write_text(json.dumps(dict(channel=ch,members=names,checks=check,deliver_sha256=sha('bdt-deliver.tgz')),indent=1)+'\n')
print('DELIVERED',json.dumps(check),flush=True)
