"""Real stored-weight control of the dynamic HT Born-scale formula, worker only."""
import json, math
from pathlib import Path
import scale_inputs as base

def ht(event):
    outgoing=[p for p in event['particles'] if p['status']==1]
    higgs=[p for p in outgoing if p['id']==25]
    daughters=[p for p in outgoing if abs(p['id']) in (11,12,13,14,15,16)]
    assert len(outgoing)==3 and len(higgs)==1 and len(daughters)==2
    x,y,z,e=higgs[0]['v']
    m2=max(0.,e*e-x*x-y*y-z*z)
    return math.sqrt(m2+x*x+y*y)+sum(math.hypot(p['v'][0],p['v'][1]) for p in daughters)

def main():
    path=Path('atlas-scale-reference.lhe.gz')
    assert base.sha(path)==Path(str(path)+'.sha256.ref').read_text().split()[0]
    card,defs=base.header(path)
    required=dict(lhans1=93300,lhans2=93300,dynscale=3,bornonly=1,loevents=1,
                  ebeam1=6800,ebeam2=6800,vdecaymode=11,mtscheme=0,alphas_from_pdf=1,
                  renscfact=1,facscfact=1)
    for k,v in required.items():assert card[k]==v,(k,card.get(k))
    assert card.get('fixedscale',0)!=1
    pdf=base.lhapdf.mkPDF(base.PDFS[93300],0)
    audit=dict(accepted=False,package=Path('VERSION').read_text().strip(),
               lhe_sha256=base.sha(path),card=card,weight_definitions=defs,
               scale='HT=mT(H)+pT(daughter1)+pT(daughter2), no factor 1/2',
               rtol=1e-4,atol=1e-8,events=0,max_relative_error=0.)
    try:
        for event in base.events(path):
            expected=base.np.array([event['weights'][k]/event['weight'] for k in base.IDS])
            predicted=base.np.array(base.born_ratio(event,pdf,pdf,base.POINTS,False,mu0=ht(event)))
            audit['max_relative_error']=max(audit['max_relative_error'],float(base.np.max(abs(predicted-expected)/abs(expected))))
            base.np.testing.assert_allclose(predicted,expected,rtol=1e-4,atol=1e-8)
            audit['events']+=1
        assert audit['events']==5000
        audit['accepted']=True
    except Exception as exc:
        audit['error']=repr(exc)
        raise
    finally:
        Path('scale-audit.json').write_text(json.dumps(audit,indent=2)+'\n')

if __name__=='__main__':main()
