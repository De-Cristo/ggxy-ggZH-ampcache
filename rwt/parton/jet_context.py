"""Registered analysis-jet context, shared by training and future model application.

Inputs are the exact scalar-table jet collection, pT>30 GeV and |eta|<4.7.
No jet masses or particle reclustering are inferred. Call extra_context with
aligned table columns, concatenate to the original24 context, and use the
normalization saved in the model checkpoint. All event-array work is worker-only.
"""
from common import worker

SCHEMA=['log1p_njets_isr'] + [f'j{j}_{name}' for j in (1,2)
    for name in ('present','log1p_pt_GeV','eta','sin_phi_minus_phiZ','cos_phi_minus_phiZ')]


def extra_context(cols):
    worker()
    import numpy as np
    nj=np.asarray(cols['njets_isr']);assert np.isfinite(nj).all() and (nj>=0).all() and (nj==np.floor(nj)).all()
    zphi=np.asarray(cols['phiz']);assert np.isfinite(zphi).all()
    entries=[np.log1p(nj.astype(float))]
    for j in (1,2):
        present=nj>=j;pt=np.asarray(cols[f'ptj{j}']);eta=np.asarray(cols[f'etaj{j}']);phi=np.asarray(cols[f'phij{j}'])
        assert np.isfinite(pt).all() and np.array_equal(pt>0,present),'Jet presence disagrees with multiplicity'
        assert (pt[present]>=30-1e-4).all() and (np.abs(eta[present])<=4.7+1e-5).all()
        assert np.isfinite(eta[present]).all() and np.isfinite(phi[present]).all()
        delta=phi[present]-zphi[present]
        transformed=[np.log1p(pt[present]),eta[present],np.sin(delta),np.cos(delta)]
        entries.append(present.astype(float))
        for x in transformed:
            a=np.zeros(len(nj));a[present]=x;entries.append(a)
    out=np.stack(entries,axis=1).astype(np.float32)
    assert out.shape==(len(nj),11) and np.isfinite(out).all()
    return out


def attach(data,scalar):
    worker()
    import numpy as np
    assert np.array_equal(scalar['run'],data['seed']) and np.array_equal(scalar['event'],data['table_row'])
    assert np.array_equal(scalar.in_window,data['in_window'])
    assert np.array_equal(scalar['w0'],data['weight'])
    assert data['globals'].shape==(scalar.n,24)
    extras=extra_context(scalar)
    data['globals']=np.concatenate([data['globals'],extras],axis=1)
    return dict(rows=scalar.n,event_identity='seed/table_row exact',weights_exact=True,
                original_context=24,added_context=11,global_features=35,schema=SCHEMA,
                j1_present=int((scalar['njets_isr']>=1).sum()),j2_present=int((scalar['njets_isr']>=2).sum()))
