"""Ensemble calibration (calibration set) and the registered coverage test (evaluation set), CPU only.

Inputs: calib/parton-data/, eval/parton-data/ (prepared bundles with their scalar tables) and
members/<name>.tgz (ensemble_fit.py products).  Plan: VH-TnP docs/handoff 2026-11-11 plan A.

Weights are shape weights on LO+PS: mean 1 over the in-window LO+PS events of the calibration set;
the inclusive normalization to the NLO+PS cross section is separate.  Every histogram below is a
fraction of the in-window total, as in the registered closure reports.

  nominal     w   = w0 * exp(mean_i logit_i) * k(m_ZH)
  members     w_i = w0 * exp(logit_i) * k_i(m_ZH)                 -> standardized histogram PCA (99 %)
  jet nuis.   w * c(rho)^-1 (up), w * c(rho) (down), rho = pT(j1)/m_ZH, events without a jet in their own bin
  scale       w * T_k(m_ZH, rho), T_k = sum(w0 sw_k)/sum(w0) of NLO+PS (calibration set)
"""
import json,sys,tarfile
from pathlib import Path
from common import worker,configuration,sha
worker()
import numpy as np
sys.path.insert(0,str(Path('rwt').resolve()))
import ggzh_reweight_common as C
import ggzh_reweight_bdt as BDT

c=configuration()
BDT.configure(C.FEATURE_SETS['C'])
VARS={**BDT.FEATURE_VARS,**BDT.TEST_VARS}
MZH_EDGES=np.concatenate([[0.],np.arange(150.,1250.+1e-9,25.),[np.inf]])
RHO_EDGES=np.array([0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.25,1.5,2.,np.inf])
VARS['rho']=(np.concatenate([[-1.5,-.5],RHO_EDGES[:-1],[3.]]),'pT(j1)/m(ZH)','-1 = no jet')
MIN_RAW=100;K_REL=.02;C_REL=.05;PCA_KEEP=.99

def scalar(name,label):
    return C.Sample(label,Path(name)/'parton-data'/f'{label}.npz')
def rho(s):
    j=np.asarray(s['ptj1'],float);return np.where(j>0,j/np.asarray(s['mzh'],float),-1.)
def rho_bin(s):
    r=rho(s);return np.where(r<0,0,1+np.clip(np.searchsorted(RHO_EDGES,r,side='right')-1,0,len(RHO_EDGES)-2))

def logits(path):
    with tarfile.open(path) as t:
        roots={m.name.split('/')[0] for m in t.getmembers()};assert len(roots)==1 and roots<={'parton-fit','member-fit'}
        assert all('..' not in m.name and not m.issym() and not m.islnk() for m in t.getmembers())
        name=path.stem;t.extractall(Path('members-x')/name)
    d=Path('members-x')/name/roots.pop();record=json.loads((d/'record.json').read_text())
    assert record['checkpoint_reload']=='PASS' and record['reload_audit']['same_batch_pass']
    return record,{(s,o):dict(np.load(d/f'logits-{s}-{o}.npz')) for s in ('calib','eval') for o in ('lo-cms','nlo')}

def merged_edges(edges,index_r,w_r,index_t,w_t,rel):
    """Merge adjacent bins from below until both samples have relative MC error < rel in every bin."""
    nb=len(edges)-1
    def sums(index,w):
        return np.bincount(index,weights=w,minlength=nb),np.bincount(index,weights=w*w,minlength=nb)
    (sr,sr2),(st,st2)=sums(index_r,w_r),sums(index_t,w_t)
    groups=[];cur=[];a=b=a2=b2=0.
    for i in range(nb):
        cur.append(i);a+=sr[i];a2+=sr2[i];b+=st[i];b2+=st2[i]
        if a>0 and b>0 and np.sqrt(a2)/a<rel and np.sqrt(b2)/b<rel:groups.append(cur);cur=[];a=b=a2=b2=0.
    if cur:
        if groups:groups[-1]+=cur
        else:groups.append(cur)
    lookup=np.empty(nb,int)
    for g,members in enumerate(groups):lookup[members]=g
    return [edges[g[0]] for g in groups]+[edges[-1]],lookup

def ratio_table(index,w_r,index_t,w_t,n):
    """Per-group NLO/reweighted fraction ratio (fractions of the in-window totals)."""
    fr=np.bincount(index,weights=w_r,minlength=n)/w_r.sum();ft=np.bincount(index_t,weights=w_t,minlength=n)/w_t.sum()
    assert (fr>0).all() and (ft>0).all()
    return ft/fr

def pca_sigma(replicas):
    """Standardized histogram PCA (check_reweighting upstream convention): per-bin sigma from the retained modes."""
    r=np.asarray(replicas,float);n=len(r)
    if n<2:return np.zeros(r.shape[1]),dict(members=n,modes=0)
    mean=r.mean(0);std=r.std(0,ddof=1);ok=std>0;z=np.zeros_like(r);z[:,ok]=(r[:,ok]-mean[ok])/std[ok]
    cov=np.cov(z,rowvar=False,ddof=1);val,vec=np.linalg.eigh(cov);val=np.clip(val[::-1],0,None);vec=vec[:,::-1]
    frac=np.cumsum(val)/val.sum() if val.sum()>0 else np.ones_like(val);keep=int(np.searchsorted(frac,PCA_KEEP)+1)
    shifts=(np.sqrt(val[:keep])[:,None]*vec[:,:keep].T)*std[None,:]
    return np.sqrt((shifts**2).sum(0)),dict(members=n,modes=keep,retained_variance=float(frac[keep-1]),plain_std_max=float(std.max()))

def histogram(s,w,sel,name):
    if name=='stxs':
        a,a2,raw=C.cell_sums(s,w,sel)
    else:
        x=rho(s) if name=='rho' else np.asarray(s[name],float);a,a2,_=C.hist(x,w*sel,VARS[name][0])
        raw=C.hist(x[sel],np.ones(int(sel.sum())),VARS[name][0])[2]
    total=w[sel].sum();return a/total,np.sqrt(a2)/total,raw

def main():
    members=sorted(Path('members').glob('*.tgz'));assert members
    records={};pred={}
    for path in members:records[path.stem],pred[path.stem]=logits(path)
    names=sorted(records);print('members',names,flush=True)
    samples={(s,o):scalar(s,lab) for s in ('calib','eval') for o,lab in (('lo-cms','lo-cms'),('nlo','nlo-cp5'))}
    for (s,o),smp in samples.items():
        for n in names:
            p=pred[n][s,o];assert np.array_equal(p['seed'],smp['run']) and np.array_equal(p['table_row'],smp['event']) and np.array_equal(p['in_window'],smp.in_window)
    L={key:np.stack([pred[n][key]['logits'] for n in names]) for key in samples}
    out=dict(members=names,member_records={n:dict(selection_sha256=r['selection_sha256'],initialization_seed=r['initialization_seed'],epochs_run=r['epochs_run'],best_validation_loss=r['best_validation_loss'],smoke_rows=r.get('smoke_rows',0)) for n,r in records.items()},
             calibration_selection_sha256={r['calibration_selection_sha256'] for r in records.values()}.pop(),evaluation_selection_sha256={r['evaluation_selection_sha256'] for r in records.values()}.pop())
    lo_c,nlo_c,lo_e,nlo_e=samples['calib','lo-cms'],samples['calib','nlo'],samples['eval','lo-cms'],samples['eval','nlo']
    sc,tc,se,te=lo_c.in_window,nlo_c.in_window,lo_e.in_window,nlo_e.in_window
    mzh_index=lambda s:np.clip(np.searchsorted(MZH_EDGES,np.asarray(s['mzh'],float),side='right')-1,0,len(MZH_EDGES)-2)

    def raw_weights(s,logit,sel,norm=None):
        # exp of (logit - offset) for stability; normalized to mean 1 over the calibration in-window LO.
        w=s.w0*np.exp(logit-np.max(logit[sel]));return w,(norm if norm is not None else sel.sum()/w[sel].sum())
    tables={}
    def calibrate(tag,logit_c):
        w,n=raw_weights(lo_c,logit_c,sc);w=w*n
        edges,lookup=merged_edges(MZH_EDGES,mzh_index(lo_c)[sc],w[sc],mzh_index(nlo_c)[tc],nlo_c.w0[tc],K_REL)
        g=lookup[mzh_index(lo_c)];k=ratio_table(g[sc],w[sc],lookup[mzh_index(nlo_c)][tc],nlo_c.w0[tc],len(edges)-1)
        n2=sc.sum()/(w*k[g])[sc].sum()
        tables[tag]=dict(offset=float(np.max(logit_c[sc])),norm=float(n*n2),mzh_edges=[float(x) for x in edges],k=k.tolist(),mzh_lookup=lookup.tolist())
        return tables[tag]
    def apply(tag,s,logit):
        t=tables[tag];lookup=np.asarray(t['mzh_lookup']);k=np.asarray(t['k'])
        return s.w0*np.exp(logit-t['offset'])*t['norm']*k[lookup[mzh_index(s)]]
    nominal_c=L['calib','lo-cms'].mean(0);calibrate('nominal',nominal_c)
    for i,n in enumerate(names):calibrate(n,L['calib','lo-cms'][i])
    w_nom_c=apply('nominal',lo_c,nominal_c)
    # Leading-jet non-closure c(rho) after k, on the calibration set.
    rb_c=rho_bin(lo_c);rb_t=rho_bin(nlo_c);nrho=len(RHO_EDGES)
    redges=np.concatenate([[-1.],RHO_EDGES]);redges,rlookup=merged_edges(redges,rb_c[sc],w_nom_c[sc],rb_t[tc],nlo_c.w0[tc],C_REL)
    cfac=1/ratio_table(rlookup[rb_c][sc],w_nom_c[sc],rlookup[rb_t][tc],nlo_c.w0[tc],len(redges)-1)
    tables['jet_nonclosure']=dict(rho_edges=[float(x) for x in redges],c=cfac.tolist(),rho_lookup=rlookup.tolist(),no_jet_bin='first edge -1: events without a jet')
    # Scale transport T_k(m_ZH group, rho group) from NLO+PS (stored ratios ids 1001-1009 at PDF4LHC21; see plan WP3).
    nom_lookup=np.asarray(tables['nominal']['mzh_lookup']);ng=len(tables['nominal']['k'])
    cell_t=nom_lookup[mzh_index(nlo_c)]*(len(redges)-1)+rlookup[rb_t];ncell=ng*(len(redges)-1)
    base=np.bincount(cell_t[tc],weights=nlo_c.w0[tc],minlength=ncell)
    T={}
    for wid in C.SCALE7:
        v=np.bincount(cell_t[tc],weights=nlo_c.weight(wid)[tc],minlength=ncell);T[wid]=np.where(base!=0,v/np.where(base!=0,base,1),1.)
    tables['scale_transport']=dict(ids=C.SCALE7,cells='mzh group (nominal k groups) x rho group',T={k:v.tolist() for k,v in T.items()},
                                   note='ids 1001-1009 at PDF4LHC21 divided by the stored nominal; target central = NNPDF3.1 id 2000')
    # Evaluation set: the registered coverage test.
    w_nom=apply('nominal',lo_e,L['eval','lo-cms'].mean(0))
    w_mem=[apply(n,lo_e,L['eval','lo-cms'][i]) for i,n in enumerate(names)]
    cidx=rlookup[rho_bin(lo_e)];w_up=w_nom/cfac[cidx];w_dn=w_nom*cfac[cidx]
    cell_e=nom_lookup[mzh_index(lo_e)]*(len(redges)-1)+cidx;cell_te=nom_lookup[mzh_index(nlo_e)]*(len(redges)-1)+rlookup[rho_bin(nlo_e)]
    coverage={};pulls=[]
    for name in ['stxs']+list(VARS):
        fr,er,raw=histogram(lo_e,w_nom,se,name);ft,et,_=histogram(nlo_e,nlo_e.w0,te,name)
        fu,_,_=histogram(lo_e,w_up,se,name);fd,_,_=histogram(lo_e,w_dn,se,name)
        reps=[histogram(lo_e,w,se,name)[0] for w in w_mem];snn,pinfo=pca_sigma(reps)
        snc=np.maximum(np.abs(fu-fr),np.abs(fd-fr))
        ok=raw>=MIN_RAW;den=np.sqrt(er**2+et**2+snn**2+snc**2);pull=np.where(ok&(den>0),(fr-ft)/np.where(den>0,den,1),np.nan)
        scale_rw=[histogram(lo_e,w_nom*T[wid][cell_e],se,name)[0] for wid in C.SCALE7]
        scale_t=[histogram(nlo_e,nlo_e.weight(wid),te,name)[0] for wid in C.SCALE7]
        scale_tt=[histogram(nlo_e,nlo_e.w0*T[wid][cell_te],te,name)[0] for wid in C.SCALE7]
        half=np.maximum(np.max(scale_t,0)-ft,ft-np.min(scale_t,0))
        # WP3 gate, NLO+PS side only: the binned transport T_k against the exact stored variation.
        transport=np.max(np.abs(np.asarray(scale_tt)-np.asarray(scale_t)),0)/np.where(half>0,half,np.inf)
        rel_half=half/np.where(ft>0,ft,1)
        coverage[name]=dict(edges=None if name=='stxs' else [float(x) for x in VARS[name][0]],raw_lo=raw.tolist(),reweighted=fr.tolist(),target=ft.tolist(),mc_reweighted=er.tolist(),mc_target=et.tolist(),
            nn=snn.tolist(),pca=pinfo,nonclosure=snc.tolist(),jet_up=fu.tolist(),jet_down=fd.tolist(),pull=[None if not np.isfinite(x) else float(x) for x in pull],
            ratio=(fr/np.where(ft>0,ft,np.nan)).tolist(),scale_reweighted=np.asarray(scale_rw).tolist(),scale_target=np.asarray(scale_t).tolist(),
            scale_transport_error_over_half_width=[float(x) for x in transport],scale_target_rel_half_width=[float(x) for x in rel_half])
        pulls+=[float(x) for x in pull if np.isfinite(x)]
    pulls=np.abs(np.asarray(pulls))
    out['coverage_test']=dict(bins=int(len(pulls)),max_abs_pull=float(pulls.max()),beyond_2=int((pulls>2).sum()),beyond_3=int((pulls>3).sum()),
        fraction_beyond_2=float((pulls>2).mean()),pass_=bool((pulls<=3).all() and (pulls>2).mean()<=.05),
        rule='bins with >=100 raw LO+PS events; pull=(R-T)/sqrt(MC_R^2+MC_T^2+NN^2+nonclosure^2); pass: none >3 and <=5% >2; scale not in the denominator')
    out['effective_fraction_eval']=float(w_nom[se].sum()**2/np.sum(w_nom[se]**2)/se.sum())
    out['tables']=tables;out['coverage']=coverage
    Path('ensemble-assess.json').write_text(json.dumps(out,indent=1)+'\n')
    Path('ensemble-tables.json').write_text(json.dumps(dict(members=names,tables=tables,vars=dict(rho=VARS['rho'][2])),indent=1)+'\n')
    print('COVERAGE',json.dumps(out['coverage_test']),flush=True)

if __name__=='__main__':main()
