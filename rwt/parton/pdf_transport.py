"""PDF4LHC21 PDF+alpha_s ratios of NLO+PS in fine (m_ZH, rho) cells, CPU only.

The 10M production stores only the scale weights.  The 5e5-event 13.6 TeV NLO production of the paper
(set 13600, seeds 10001-10200, weight list pwg-rwl-pdf4lhc21.dat) carries PDF4LHC21 members 93300-93342 as LHE ids 3000-3042;
its CP5 shower tables (tables/nlo-cp5-13600) give m_ZH and the leading jet of the same events.  Rows are
joined by event order and checked on the stored scale weight 1001.  Output: per fine cell the sums of the
central weight (id 3000) and of each member, for ensemble_deliver.py to merge into the calibration groups.
"""
import gzip,json,os,re,subprocess,sys
from pathlib import Path
from common import worker,sha
worker()
import numpy as np
from ensemble_assess_bins import MZH_EDGES,RHO_EDGES,fine_cells

ROOT=os.environ['PAPER_EOS_ROOT'].rstrip('/')   # product storage, set by the job environment as for the other parton tasks
SEEDS=range(10001,10201);IDS=[str(3000+k) for k in range(43)]   # LHE ids 3000-3042 = LHAPDF 93300-93342
WGT=re.compile(r"<wgt id='(\d+)'>\s*([-+0-9.eEdD]+)\s*</wgt>")

def fetch(product,local):
    Path(local).parent.mkdir(parents=True,exist_ok=True)
    for url,target in ((f'{ROOT}/{product}.sha256',local+'.sha256'),(f'{ROOT}/{product}',local)):
        subprocess.run(['xrdcp','--retry','3','-s','-f',url,target],check=True)
    assert Path(local+'.sha256').read_text().split()[0]==sha(local),product

def lhe_weights(path):
    rows=[];event=None
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('<event'):event={}
            elif event is not None and '<wgt' in line:
                for k,v in WGT.findall(line):event[k]=float(v.replace('d','e').replace('D','e'))
            elif line.startswith('</event'):
                rows.append([event[k] for k in ['1001']+IDS]);event=None
    return np.asarray(rows)

def main():
    sums=np.zeros((len(IDS),len(MZH_EDGES)-1,len(RHO_EDGES)));count=np.zeros((len(MZH_EDGES)-1,len(RHO_EDGES)))
    audit=dict(seeds=[SEEDS.start,SEEDS.stop-1],files=0,events=0,max_rel_1001=0.)
    for s in SEEDS:
        fetch(f'lhe/nlo-13600/events-{s}.lhe.gz',f'in/events-{s}.lhe.gz');fetch(f'tables/nlo-cp5-13600/table-{s}.csv.gz',f'in/table-{s}.csv.gz')
        w=lhe_weights(f'in/events-{s}.lhe.gz')
        with gzip.open(f'in/table-{s}.csv.gz','rt') as f:
            head=f.readline().strip().split(',');t=np.loadtxt(f,delimiter=',',ndmin=2)
        col={k:i for i,k in enumerate(head)};ievt=t[:,col['ievt']].astype(int)
        assert len(t)==len(w) and np.array_equal(ievt,np.arange(len(t)))
        rel=np.max(np.abs(t[:,col['w:1001']]/w[:,0]-1));assert rel<1e-4,(s,rel)
        mi,ri=fine_cells(t[:,col['mzh']],t[:,col['ptj1eta47']])
        mll=t[:,col['mz']];keep=(mll>15)&(mll<150)
        for k in range(len(IDS)):np.add.at(sums[k],(mi[keep],ri[keep]),w[keep,k+1])
        np.add.at(count,(mi[keep],ri[keep]),1)
        audit['files']+=1;audit['events']+=int(len(t));audit['max_rel_1001']=max(audit['max_rel_1001'],float(rel))
        for p in Path('in').glob('*'):p.unlink()
    np.savez_compressed('pdf-transport.npz',ids=np.asarray(IDS),sums=sums,raw=count,mzh_edges=MZH_EDGES,rho_edges=RHO_EDGES)
    audit.update(central='3000 = LHAPDF 93300',members='3001-3040 Hessian (93301-93340), 3041/3042 alpha_s (93341/93342)',window='15<m_ll<150 (last-copy Z mass)',
                 output_sha256=sha('pdf-transport.npz'))
    Path('pdf-transport.json').write_text(json.dumps(audit,indent=2)+'\n');print('PDF_TRANSPORT',json.dumps(audit),flush=True)

if __name__=='__main__':main()
