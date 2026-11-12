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
import json,os,sys,tarfile
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
from ensemble_assess_bins import MZH_EDGES,RHO_EDGES,mzh_cell,rho_cell
VARS['rho']=(np.concatenate([[-1.5,-.5],RHO_EDGES[:-1],[3.]]),'pT(j1)/m(ZH)','-1 = no jet')
MIN_RAW=100;K_REL=.02;C_REL=.05;PCA_KEEP=.99

def scalar(name,label):
    return C.Sample(label,Path(name)/'parton-data'/f'{label}.npz')
def rho(s):
    j=np.asarray(s['ptj1'],float);return np.where(j>0,j/np.asarray(s['mzh'],float),-1.)
def rho_bin(s):
    return rho_cell(s['ptj1'],s['mzh'])

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
    if r.shape[1]==1:return std,dict(members=n,modes=1,retained_variance=1.,plain_std_max=float(std.max()))
    cov=np.cov(z,rowvar=False,ddof=1);val,vec=np.linalg.eigh(cov);val=np.clip(val[::-1],0,None);vec=vec[:,::-1]
    frac=np.cumsum(val)/val.sum() if val.sum()>0 else np.ones_like(val);keep=int(np.searchsorted(frac,PCA_KEEP)+1)
    shifts=(np.sqrt(val[:keep])[:,None]*vec[:,:keep].T)*std[None,:]
    return np.sqrt((shifts**2).sum(0)),dict(members=n,modes=keep,retained_variance=float(frac[keep-1]),plain_std_max=float(std.max()))

def histogram(s,w,sel,name,edges=None):
    if name=='stxs':
        a,a2,raw=C.cell_sums(s,w,sel)
    else:
        e=VARS[name][0] if edges is None else edges
        x=rho(s) if name=='rho' else np.asarray(s[name],float);a,a2,_=C.hist(x,w*sel,e)
        raw=C.hist(x[sel],np.ones(int(sel.sum())),e)[2]
    total=w[sel].sum();return a/total,np.sqrt(a2)/total,raw

def variants(lo_e,nlo_e,se,te,w_nom,w_mem,w_up,w_dn):
    """Post-hoc diagnostics on the already inspected evaluation set (not the registered test): wider bins
    (fine bins merged from below until each has >= 1000 raw LO+PS events; STXS cells unchanged) and/or the
    analysis domain m_ll > 40 GeV (charged leptons only).  Same pull definition and pass rule."""
    out={};zll=configuration()['channel']=='zll'
    def mll40(s):return np.asarray(s['lhe_mll'],float)>40
    def domain(s):
        # VH(bb)-like kinematic domain at generator level: Z and H back to back and balanced, Z mass window (ll).
        dphi=np.abs(np.arctan2(np.sin(np.asarray(s['phiz'],float)-np.asarray(s['phih'],float)),np.cos(np.asarray(s['phiz'],float)-np.asarray(s['phih'],float))))
        r=np.asarray(s['pth'],float)/np.maximum(np.asarray(s['ptz'],float),1e-9)
        d=(dphi>(2.5 if zll else 2.2))&(r>=.5)&(r<=2.)
        if zll:m=np.asarray(s['lhe_mll'],float);d&=(m>75)&(m<105)
        return d
    specs=(('wide',1000,None),('mll40',None,mll40),('wide+mll40',1000,mll40))
    if os.environ.get('RWT_ASSESS_VARIANT_SET')=='v2':specs=(('wide2k',2000,None),('domain',None,domain),('domain+wide2k',2000,domain))
    if os.environ.get('RWT_ASSESS_VARIANT_SET')=='tail':specs=(('tail500',-500,None),('tail1k',-1000,None),('tail2k',-2000,None))
    for label,minraw,cut in specs:
        if cut is mll40 and not zll:continue
        wide=minraw is not None;tail=wide and minraw<0;minraw=abs(minraw) if wide else None
        s2=se&cut(lo_e) if cut else se;t2=te&cut(nlo_e) if cut else te
        pulls=[];fails=[];hists={}
        for name in ['stxs']+list(VARS):
            edges=None
            if wide and name!='stxs':
                fine=np.asarray(VARS[name][0],float);raw=histogram(lo_e,lo_e.w0,s2,name)[2]
                if tail:
                    # Merge only sparse tail bins, from the high end down; populated bins stay as they are.
                    keep=[fine[-1]];acc=0
                    for i in range(len(raw)-1,-1,-1):
                        acc+=raw[i]
                        if acc>=minraw:keep.append(fine[i]);acc=0
                    if keep[-1]!=fine[0]:
                        if len(keep)>1:keep[-1]=fine[0]
                        else:keep.append(fine[0])
                    edges=np.asarray(keep[::-1])
                else:
                    keep=[fine[0]];acc=0
                    for i in range(len(raw)):
                        acc+=raw[i]
                        if acc>=minraw:keep.append(fine[i+1]);acc=0
                    if keep[-1]!=fine[-1]:
                        if len(keep)>1:keep[-1]=fine[-1]
                        else:keep.append(fine[-1])
                    edges=np.asarray(keep)
            fr,er,raw=histogram(lo_e,w_nom,s2,name,edges);ft,et,_=histogram(nlo_e,nlo_e.w0,t2,name,edges)
            fu=histogram(lo_e,w_up,s2,name,edges)[0];fd=histogram(lo_e,w_dn,s2,name,edges)[0]
            snn,_=pca_sigma([histogram(lo_e,w,s2,name,edges)[0] for w in w_mem]);snc=np.maximum(np.abs(fu-fr),np.abs(fd-fr))
            ok=raw>=MIN_RAW;den=np.sqrt(er**2+et**2+snn**2+snc**2);pull=np.where(ok&(den>0),(fr-ft)/np.where(den>0,den,1),np.nan)
            e=VARS[name][0] if (edges is None and name!='stxs') else edges
            hists[name]=dict(edges=None if e is None else [float(x) for x in e],raw_lo=raw.tolist(),reweighted=fr.tolist(),target=ft.tolist(),
                unc=np.sqrt(er**2+et**2+snn**2+snc**2).tolist(),pull=[None if not np.isfinite(x) else float(x) for x in pull],
                selected_lo=int(s2.sum()),selected_nlo=int(t2.sum()))
            for i,p in enumerate(pull):
                if np.isfinite(p):
                    pulls.append(abs(p))
                    if abs(p)>2:fails.append(dict(var=name,bin=[float(e[i]),float(e[i+1])] if e is not None else i,pull=float(p),ratio=float(fr[i]/ft[i]) if ft[i]>0 else None,raw_lo=int(raw[i])))
        pulls=np.asarray(pulls)
        out[label]=dict(bins=int(len(pulls)),beyond_2=int((pulls>2).sum()),beyond_3=int((pulls>3).sum()),fraction_beyond_2=float((pulls>2).mean()),
            max_abs_pull=float(pulls.max()),pass_=bool((pulls<=3).all() and (pulls>2).mean()<=.05),failing_bins=sorted(fails,key=lambda f:-abs(f['pull'])),histograms=hists)
        print('VARIANT',label,json.dumps({k:v for k,v in out[label].items() if k not in ('failing_bins','histograms')}),flush=True)
    return out

def scaleband(lo_c,nlo_c,Lc,names,tables,T,cell_of,apply,w0n):
    """Absolute scale bands on the calibration set (window only, domain convention of the deliverable): raw LO+PS
    with its Born-formula scale points (scale-inputs companion), the reweighted sample with the transported NLO scale
    points, members and jet variations, and NLO+PS with its scale points.  Sums are in units of the nominal totals."""
    with tarfile.open('scale-inputs.tgz') as t:
        assert all('..' not in m.name and not m.issym() and not m.islnk() for m in t.getmembers());t.extractall('.')
    comp=dict(np.load('scale-inputs/lo-cms.npz'))
    assert np.array_equal(comp['seed'],lo_c['run']) and np.array_equal(comp['table_row'],lo_c['event']),'companion rows do not align'
    ratio=comp['scale_nnpdf']/comp['nominal_lhe'][:,None];ids=list(comp['scale_ids']);assert ids[0]=='1001'
    sel_lo=lo_c.in_window;sel_nlo=nlo_c.in_window
    from ensemble_assess_bins import in_domain
    inside=in_domain(lo_c['ptj1'],lo_c['mzh']) if os.environ.get('RWT_DOMAIN')=='rho1' else np.ones(lo_c.n,bool)
    def rw(logit,tag):
        w=apply(tag,lo_c,logit);return np.where(inside,w,w0n(lo_c))
    w_nom=rw(Lc.mean(0),'nominal');cell=cell_of(lo_c)
    jet=tables.get('jet_nonclosure');cf=np.asarray(jet['c'])[np.asarray(jet['rho_lookup'])[rho_cell(lo_c['ptj1'],lo_c['mzh'])]]
    out={}
    # Outside the domain: weight 1 with the hard-jet share-ratio variation of the deliverable (up = r, down = 1/r).
    o=tables.get('domain',{}).get('outside')
    x=np.maximum(np.asarray(lo_c['ptz'],float),np.asarray(lo_c['pth'],float))
    r=np.asarray(o['share_ratio'])[np.clip(np.searchsorted(o['edges'],x,side='right')-1,0,len(o['share_ratio'])-1)] if o else np.ones(lo_c.n)
    sets=dict(lo_raw={'nominal':lo_c.w0,**{f'scale_{i}':lo_c.w0*ratio[:,k] for k,i in enumerate(ids) if i!='1001'}},
              rw={'nominal':w_nom,'jet_up':np.where(inside,w_nom/cf,w_nom*r),'jet_down':np.where(inside,w_nom*cf,w_nom/r),
                  **{f'member_{i:02d}':rw(Lc[i],n) for i,n in enumerate(names)},**{f'scale_{wid}':w_nom*T[wid][cell] for wid in T}},
              nlo={'nominal':nlo_c.w0,**{f'scale_{wid}':nlo_c.weight(wid) for wid in C.SCALE7}})
    for name in ['stxs']+list(VARS):
        entry=dict(edges=None if name=='stxs' else [float(x) for x in VARS[name][0]])
        for group,ws in sets.items():
            s,sel=(nlo_c,sel_nlo) if group=='nlo' else (lo_c,sel_lo);entry[group]={}
            for k,w in ws.items():
                if name=='stxs':a,a2,raw=C.cell_sums(s,w,sel)
                else:
                    x=rho(s) if name=='rho' else np.asarray(s[name],float);a,a2,_=C.hist(x,w*sel,VARS[name][0]);raw=C.hist(x[sel],np.ones(int(sel.sum())),VARS[name][0])[2]
                entry[group][k]=dict(sumw=a.tolist(),sumw2=a2.tolist(),nominal_total=float(ws['nominal'][sel].sum()))
                if group=='lo_raw' and k=='nominal':entry['raw_lo']=raw.tolist()
        out[name]=entry
    Path('scaleband-cache.json').write_text(C.dumps(dict(set='calibration',domain=os.environ.get('RWT_DOMAIN','none'),
        note='absolute: sum(weights)/nominal_total of the same group; LO scale points = Born-formula companion (NNPDF3.1, fixed mu0); rw scale = transported NLO points',
        lo_scale_ids=ids,histograms=out)))
    print('SCALEBAND_CACHE_WRITTEN',flush=True)

def main():
    members=sorted(Path('members').glob('*.tgz'));assert members
    records={};pred={}
    for path in members:records[path.stem],pred[path.stem]=logits(path)
    names=sorted(records);print('members',names,flush=True)
    samples={(s,o):scalar(s,lab) for s in ('calib','eval') for o,lab in (('lo-cms','lo-cms'),('nlo','nlo-cp5'))}
    for smp in samples.values():BDT.prepare(smp)   # adds the derived classifier columns njets and ayzh
    for (s,o),smp in samples.items():
        for n in names:
            p=pred[n][s,o];assert np.array_equal(p['seed'],smp['run']) and np.array_equal(p['table_row'],smp['event']) and np.array_equal(p['in_window'],smp.in_window)
    L={key:np.stack([pred[n][key]['logits'] for n in names]) for key in samples}
    out=dict(members=names,member_records={n:dict(selection_sha256=r['selection_sha256'],initialization_seed=r['initialization_seed'],epochs_run=r['epochs_run'],best_validation_loss=r['best_validation_loss'],smoke_rows=r.get('smoke_rows',0)) for n,r in records.items()},
             calibration_selection_sha256={r['calibration_selection_sha256'] for r in records.values()}.pop(),evaluation_selection_sha256={r['evaluation_selection_sha256'] for r in records.values()}.pop())
    lo_c,nlo_c,lo_e,nlo_e=samples['calib','lo-cms'],samples['calib','nlo'],samples['eval','lo-cms'],samples['eval','nlo']
    sc,tc,se,te=lo_c.in_window,nlo_c.in_window,lo_e.in_window,nlo_e.in_window
    domain=os.environ.get('RWT_DOMAIN','none');assert domain in ('none','rho1')
    if domain=='rho1':
        # Calibration, transport and coverage inside the training domain; the rest is outside the model's validity.
        from ensemble_assess_bins import in_domain
        dom={k:in_domain(v['ptj1'],v['mzh']) for k,v in samples.items()}
        out['domain']=dict(name='rho1',definition='no jet, or pT(j1) < m_ZH (jet: |eta| < 4.7, pT > 30 GeV)',
            nlo_fraction_outside_eval=float(nlo_e.w0[te&~dom['eval','nlo']].sum()/nlo_e.w0[te].sum()),
            lo_fraction_outside_eval=float(lo_e.w0[se&~dom['eval','lo-cms']].sum()/lo_e.w0[se].sum()))
        # Outside the domain the delivered weight is 1 (LO+PS x K).  Its two-sided variation is the ratio of the
        # NLO+PS and LO+PS shares outside the domain, in pT(Z) bins of the calibration set: up = r, down = 1/r.
        # Binned in max(pT(Z), pT(H)): events with a hard Higgs boson and a moderate Z need the larger ratio.
        edges=[0.,150.,250.,400.,600.,np.inf];r=[]
        ptmax=lambda s:np.maximum(np.asarray(s['ptz'],float),np.asarray(s['pth'],float))
        for a,b in zip(edges[:-1],edges[1:]):
            bl=lambda s,sel:sel&(ptmax(s)>=a)&(ptmax(s)<b)
            fl=lo_c.w0[bl(lo_c,sc)&~dom['calib','lo-cms']].sum()/lo_c.w0[bl(lo_c,sc)].sum()
            fn=nlo_c.w0[bl(nlo_c,tc)&~dom['calib','nlo']].sum()/nlo_c.w0[bl(nlo_c,tc)].sum()
            r.append(float(fn/fl) if fl>0 else 1.)
        outside=dict(variable='max(ptz,pth)',edges=edges[:-1]+[1e9],share_ratio=r,nominal=1.0,note='weight 1 outside the domain (LO+PS x K); up = share_ratio, down = 1/share_ratio, binned in max(pT(Z), pT(H))')
        out['domain']['outside']=outside
        sc=sc&dom['calib','lo-cms'];tc=tc&dom['calib','nlo'];se=se&dom['eval','lo-cms'];te=te&dom['eval','nlo']
    mzh_index=lambda s:mzh_cell(s['mzh'])

    # LO+PS event weights in units of their mean over the calibration window, so that a delivered column multiplies
    # an event weight of about 1 and the 'weight 1 outside the domain' convention has the same scale as inside.
    W=float(lo_c.w0[lo_c.in_window].mean())
    def w0n(s):return s.w0/W
    def raw_weights(s,logit,sel,norm=None):
        # exp of (logit - offset) for stability; normalized to mean 1 over the calibration in-window LO.
        w=w0n(s)*np.exp(logit-np.max(logit[sel]));return w,(norm if norm is not None else sel.sum()/w[sel].sum())
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
        return w0n(s)*np.exp(logit-t['offset'])*t['norm']*k[lookup[mzh_index(s)]]
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
    if 'domain' in out:tables['domain']=out['domain']
    tables['lo_weight_scale']=dict(mean_w0_calibration_window=W,note='columns multiply LO weights expressed in units of this mean (about 1 per event)')
    tables['scale_transport']=dict(ids=C.SCALE7,cells='mzh group (nominal k groups) x rho group',T={k:v.tolist() for k,v in T.items()},
                                   note='ids 1001-1009 at PDF4LHC21 divided by the stored nominal; target central = NNPDF3.1 id 2000')
    if os.environ.get('RWT_TABLES_ONLY')=='1':
        Path('ensemble-tables.json').write_text(C.dumps(dict(members=names,tables=tables,vars=dict(rho=VARS['rho'][2]))))
        print('TABLES_WRITTEN',json.dumps(tables.get('domain',{}).get('outside')),flush=True);return
    if os.environ.get('RWT_SCALEBAND')=='1':
        scaleband(lo_c,nlo_c,L['calib','lo-cms'],names,tables,T,cell_of=lambda s:nom_lookup[mzh_index(s)]*(len(redges)-1)+rlookup[rho_bin(s)],apply=apply,w0n=w0n)
        return
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
    # Histogram cache: fine-binned sums (and sums of squares) of every weight set on the evaluation set, so that
    # any later rebinning (merging adjacent bins) and its PCA can be redone from this file without the events.
    cache={}
    for name in ['stxs']+list(VARS):
        def sums(s,w,sel):
            if name=='stxs':
                a,a2,raw=C.cell_sums(s,w,sel)
            else:
                x=rho(s) if name=='rho' else np.asarray(s[name],float);a,a2,_=C.hist(x,w*sel,VARS[name][0]);raw=C.hist(x[sel],np.ones(int(sel.sum())),VARS[name][0])[2]
            return a.tolist(),a2.tolist(),raw.tolist(),float(w[sel].sum())
        lo_sets={'nominal':w_nom,'jet_up':w_up,'jet_down':w_dn,**{f'member_{i:02d}':w for i,w in enumerate(w_mem)},**{f'scale_{wid}':w_nom*T[wid][cell_e] for wid in C.SCALE7}}
        nlo_sets={'nominal':nlo_e.w0,**{f'scale_{wid}':nlo_e.weight(wid) for wid in C.SCALE7}}
        entry=dict(edges=None if name=='stxs' else [float(x) for x in VARS[name][0]],lo={},nlo={})
        for k,w in lo_sets.items():a,a2,raw,tot=sums(lo_e,w,se);entry['lo'][k]=dict(sumw=a,sumw2=a2,total=tot);entry['raw_lo']=raw
        for k,w in nlo_sets.items():a,a2,raw,tot=sums(nlo_e,w,te);entry['nlo'][k]=dict(sumw=a,sumw2=a2,total=tot)
        cache[name]=entry
    Path('histogram-cache.json').write_text(C.dumps(dict(set=os.environ.get('RWT_PARTON_SELECTION'),domain=os.environ.get('RWT_DOMAIN','none'),
        note='evaluation set; fractions = sumw/total; MC error = sqrt(sumw2)/total; NN = PCA over member_XX fractions',histograms=cache)))
    pulls=np.abs(np.asarray(pulls))
    out['coverage_test']=dict(bins=int(len(pulls)),max_abs_pull=float(pulls.max()),beyond_2=int((pulls>2).sum()),beyond_3=int((pulls>3).sum()),
        fraction_beyond_2=float((pulls>2).mean()),pass_=bool((pulls<=3).all() and (pulls>2).mean()<=.05),
        rule='bins with >=100 raw LO+PS events; pull=(R-T)/sqrt(MC_R^2+MC_T^2+NN^2+nonclosure^2); pass: none >3 and <=5% >2; scale not in the denominator')
    out['effective_fraction_eval']=float(w_nom[se].sum()**2/np.sum(w_nom[se]**2)/se.sum())
    if os.environ.get('RWT_ASSESS_VARIANTS')=='1':
        out['variants']=variants(lo_e,nlo_e,se,te,w_nom,w_mem,w_up,w_dn)
    # The predetermined criteria 1-3 of the paper (reported, not the pass rule) and the region/hard-jet
    # diagnostics, nominal weight on the evaluation set; 'unweighted' is the raw LO+PS sample.
    from ggzh_rwt10m import diagnostics
    out['closure_criteria']=BDT.evaluate(lo_e,w_nom,se,nlo_e,te,'full',print)
    out['diagnostics']=diagnostics(lo_e,w_nom,nlo_e)
    out['diagnostics_unweighted']=diagnostics(lo_e,np.ones(lo_e.n),nlo_e)
    for name in coverage:
        coverage[name]['raw']=histogram(lo_e,lo_e.w0,se,name)[0].tolist()
    out['tables']=tables;out['coverage']=coverage
    Path('ensemble-assess.json').write_text(C.dumps(out)+'\n')   # numpy scalars and NaN from the criteria code
    Path('ensemble-tables.json').write_text(C.dumps(dict(members=names,tables=tables,vars=dict(rho=VARS['rho'][2])))+'\n')
    print('COVERAGE',json.dumps(out['coverage_test']),flush=True)

if __name__=='__main__':main()
