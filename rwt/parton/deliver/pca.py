"""NN (ensemble) uncertainty in any binning: standardized histogram PCA of the member histograms.

Fill the same histogram once per member weight (member_00 .. member_19 from apply.py), then

    shifts, info = pca_shifts(hists)        # hists: (n_members, n_bins) bin contents

gives one shift vector per retained mode; use nominal +- shift as the up/down templates of independent
nuisance parameters.  Convention of arXiv:2510.16217 (Palmer, Kronheim): each bin is standardized across the
members, the modes explaining 99 % of the variance are kept, and the shifts are scaled back to bin contents.
Bin contents, not normalized shapes: normalization changes are part of the uncertainty.
"""
import numpy as np

def pca_shifts(hists,keep=0.99):
    r=np.asarray(hists,float);n=len(r)
    if n<2:return np.zeros((0,r.shape[1])),dict(members=n,modes=0)
    mean=r.mean(0);std=r.std(0,ddof=1);ok=std>0;z=np.zeros_like(r);z[:,ok]=(r[:,ok]-mean[ok])/std[ok]
    val,vec=np.linalg.eigh(np.cov(z,rowvar=False,ddof=1));val=np.clip(val[::-1],0,None);vec=vec[:,::-1]
    frac=np.cumsum(val)/val.sum() if val.sum()>0 else np.ones_like(val);k=int(np.searchsorted(frac,keep)+1)
    return (np.sqrt(val[:k])[:,None]*vec[:,:k].T)*std[None,:],dict(members=n,modes=k,retained_variance=float(frac[k-1]))
