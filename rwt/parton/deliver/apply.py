"""Per-event LO+PS -> NLO+PS weights and their variations for gg->ZH (BDT-17 ensemble).  Needs numpy and onnxruntime.

    python apply.py --model DIR --inputs events.npz --out weights.npz

events.npz holds one array per input named in DIR/calibration.json "features" (generator level, model card),
in GeV.  Events outside 15 < m_ll < 150 GeV (column lhe_mll) get NaN: the model is not defined there.
Output columns (multiply the event's own LO weight by them):
  nominal                       NLO+PS shape weight (mean 1 on the training LO+PS sample)
  member_00 .. member_19        ensemble members -> NN uncertainty with pca.py
  jet_up, jet_down              leading-jet non-closure (two-sided)
  scale_<muR>_<muF>             NLO+PS 7-point scale variations (6 varied points)
  pdf_00 .. pdf_39, alphas_dn, alphas_up   PDF4LHC21 Hessian members and alpha_s variations
Normalise the nominal total to the NLO+PS cross section of the model card and apply the same factor to every column:
the scale and PDF variations then also carry the change of the total cross section, the others are shape-only.
"""
import argparse,json
from pathlib import Path
import numpy as np

class Model:
    def __init__(self,directory):
        import onnxruntime as ort
        d=Path(directory);self.cal=json.loads((d/'calibration.json').read_text())
        so=ort.SessionOptions();so.intra_op_num_threads=so.inter_op_num_threads=1
        self.sessions=[ort.InferenceSession(str(d/m['onnx']),sess_options=so,providers=['CPUExecutionProvider']) for m in self.cal['members']]
        self.fine_mzh=np.asarray(self.cal['fine_mzh_edges'],float);self.fine_rho=np.asarray(self.cal['fine_rho_edges'],float)

    def matrix(self,ev):
        return np.stack([np.asarray(ev[k],np.float64) for k in self.cal['features']],axis=1).astype(np.float32)

    def logits(self,X):
        out=[]
        for s in self.sessions:
            names=[o.name for o in s.get_outputs()];p=s.run(None,{'X':X})[names.index('probabilities')][:,1].astype(np.float64)
            p=np.clip(p,1e-15,1-1e-15);out.append(np.log(p)-np.log1p(-p))
        return np.asarray(out)

    def cells(self,ev):
        mzh=np.asarray(ev['mzh'],float);j=np.asarray(ev['ptj1'],float)
        m=np.clip(np.searchsorted(self.fine_mzh,mzh,side='right')-1,0,len(self.fine_mzh)-2)
        r=np.where(j>0,j/mzh,-1.)
        rc=np.where(r<0,0,1+np.clip(np.searchsorted(self.fine_rho,r,side='right')-1,0,len(self.fine_rho)-2))
        return m,rc

    def weights(self,ev):
        c=self.cal;L=self.logits(self.matrix(ev));m,rc=self.cells(ev)
        def calibrated(t,logit):
            return np.exp(logit-t['offset'])*t['norm']*np.asarray(t['k'])[np.asarray(t['mzh_lookup'])[m]]
        out={'nominal':calibrated(c['nominal'],L.mean(0))}
        for i,t in enumerate(c['members']):out[f'member_{i:02d}']=calibrated(t,L[i])
        g=np.asarray(c['jet']['rho_lookup'])[rc];cf=np.asarray(c['jet']['c'])[g]
        out['jet_up']=out['nominal']/cf;out['jet_down']=out['nominal']*cf
        cell=np.asarray(c['nominal']['mzh_lookup'])[m]*len(c['jet']['c'])+g
        for name,T in c['scale']['T'].items():out[name]=out['nominal']*np.asarray(T)[cell]
        for name,T in c['pdf']['T'].items():out[name]=out['nominal']*np.asarray(T)[cell]
        window=(np.asarray(ev['lhe_mll'],float)>c['window'][0])&(np.asarray(ev['lhe_mll'],float)<c['window'][1])
        return {k:np.where(window,v,np.nan) for k,v in out.items()}

def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--model',required=True);p.add_argument('--inputs',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();ev=dict(np.load(a.inputs));np.savez_compressed(a.out,**Model(a.model).weights(ev))

if __name__=='__main__':main()
