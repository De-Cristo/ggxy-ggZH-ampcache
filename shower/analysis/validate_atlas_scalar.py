"""Worker-only full scalar validation; exact decompressed parity, no tolerance."""
import csv, gzip, hashlib, io, json, math, os, sys
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    assert os.environ.get('_CONDOR_SCRATCH_DIR') and os.uname().nodename.split('.')[0]!='cmscon'
    ch,order,seed,expected,ss=sys.argv[1:]
    seed,expected,ss=map(int,(seed,expected,ss))
    assert ch in ('zll','zvv') and order in ('lo-atlas','nlo-a14')
    first=(500001 if ch=='zll' else 800001) if order=='lo-atlas' else (100001 if ch=='zll' else 300001)
    shower_first=(600001 if ch=='zll' else 900001) if order=='lo-atlas' else (700001 if ch=='zll' else 1000001)
    assert ss==shower_first+seed-first and 0<=seed-first<(2000 if order=='lo-atlas' else 4000)
    assert expected==(5000 if order=='lo-atlas' else 2500)
    report=dict(schema='ggzh-atlas-scalar-validation-v1',channel='atlas-'+ch,order=order,
                seed=seed,shower_seed=ss,expected=expected,package=Path('VERSION').read_text().strip(),
                primary_included=True,accepted=False,lhe_sha256=sha('input.lhe.gz'))
    try:
        blob=gzip.decompress(Path('table.csv.gz').read_bytes())
        reader=csv.DictReader(io.StringIO(blob.decode()));cols=reader.fieldnames
        assert len(cols)==len(set(cols))
        required={'w:nominal','mz','mzh','ptz','pth','njets30','ptj1','etaj1','phij1','ptj2','etaj2','phij2','ptl1','etal1','phil1','ptl2','etal2','phil2'}
        if order=='nlo-a14':required.update('w:'+str(x) for x in range(1001,1010));required.add('w:2000')
        assert required<=set(cols),(required-set(cols))
        rows=list(reader);assert len(rows)==expected,(len(rows),expected)
        assert all(set(row)==set(cols) and all(v is not None and math.isfinite(float(v)) for v in row.values()) for row in rows)
        nominal='w:nominal' if order=='lo-atlas' else 'w:1001'
        ww=[float(r[nominal]) for r in rows]
        report.update(events=len(rows),columns=cols,table_sha256=sha('table.csv.gz'),
                      decompressed_sha256=hashlib.sha256(blob).hexdigest(),
                      sumw=sum(ww),sumw2=sum(w*w for w in ww),
                      in_window_events=sum(15<float(r['mz'])<150 for r in rows))
        if Path('parity-reference.csv.gz').exists():
            spec=json.loads(Path('production/atlas-tables-parity-v1.json').read_text())
            entry=next(x for x in spec[order] if x['seed']==seed)
            assert sha('parity-reference.csv.gz')==entry['sha256']
            old=gzip.decompress(Path('parity-reference.csv.gz').read_bytes())
            report['parity']=dict(reference_sha256=entry['sha256'],original_package=entry['original_package'],
                                 exact_decompressed_content=(old==blob),reference_decompressed_sha256=hashlib.sha256(old).hexdigest())
            if old!=blob:
                before=old.splitlines();after=blob.splitlines()
                report['parity']['first_different_row']=next((i for i,(a,b) in enumerate(zip(before,after)) if a!=b),min(len(before),len(after)))
            assert old==blob,'Exact table-content identity gate failed; no rounding relaxation'
        report['accepted']=True
        schema=dict(schema='ggzh-atlas-scalar-table-v1',columns=cols,channel='atlas-'+ch,
                    lhe_mll='mz: unchanged existing table reader mapping',nominal_weight=nominal,
                    daughter_columns='neutrinos' if ch=='zvv' else 'charged leptons',
                    weights='all stored legacy scalar weights retained',producer=report['package'])
        Path('table-schema.json').write_text(json.dumps(schema,indent=2)+'\n')
    except Exception as exc:
        report['error']=repr(exc)
        raise
    finally:
        Path('table-validation.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
