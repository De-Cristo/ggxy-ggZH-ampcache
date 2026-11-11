"""Fine (m_ZH, rho) cells shared by ensemble_assess.py, pdf_transport.py and the delivered apply script.

rho = pT(j1)/m_ZH with the leading jet of |eta| < 4.7 and pT > 30 GeV; events without a jet are rho cell 0.
"""
import numpy as np
MZH_EDGES=np.concatenate([[0.],np.arange(150.,1250.+1e-9,25.),[np.inf]])
RHO_EDGES=np.array([0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.25,1.5,2.,np.inf])

def mzh_cell(mzh):
    return np.clip(np.searchsorted(MZH_EDGES,np.asarray(mzh,float),side='right')-1,0,len(MZH_EDGES)-2)

def rho_cell(ptj1,mzh):
    j=np.asarray(ptj1,float);r=np.where(j>0,j/np.asarray(mzh,float),-1.)
    return np.where(r<0,0,1+np.clip(np.searchsorted(RHO_EDGES,r,side='right')-1,0,len(RHO_EDGES)-2))

def in_domain(ptj1,mzh):
    """Training domain 'rho1': no jet, or leading jet softer than m_ZH (where the vetoed LO+PS source has events)."""
    j=np.asarray(ptj1,float);return (j<=0)|(j<np.asarray(mzh,float))

def fine_cells(mzh,ptj1):
    return mzh_cell(mzh),rho_cell(ptj1,mzh)
