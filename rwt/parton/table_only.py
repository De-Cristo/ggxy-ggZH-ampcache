"""Table-only copy of a prepared bundle for the data record, CPU worker.

The BDT-17 tasks (bdt_member.py, ensemble_assess.py, ensemble_deliver.py) read only the scalar tables lo-cms.npz and
nlo-cp5.npz and the files selection.json and audit.json of a prepared bundle.  The CMS-style bundles also carry the
particle tensors of the earlier transformer study.  This task keeps the four files byte for byte, with the
corresponding lines of the bundle's own SHA256SUMS, and writes parton-tables.tgz (same layout, parton-data/) and
table-only.json.  The input bundle is unpacked and checked against its SHA256SUMS by run_parton_study.sh.
"""
import json,tarfile
from pathlib import Path
from common import worker,sha
worker()

KEEP=('lo-cms.npz','nlo-cp5.npz','selection.json','audit.json')
src=Path('parton-data');sums={}
for line in (src/'SHA256SUMS').read_text().splitlines():
    digest,name=line.split(maxsplit=1);sums[name.strip().lstrip('*')]=digest
out=Path('tables')/'parton-data';assert not out.parent.exists();out.mkdir(parents=True)
kept=[]
for name in KEEP:
    f=src/name;assert f.is_file(),name;assert sums[name]==sha(f),name
    (out/name).write_bytes(f.read_bytes());assert sha(out/name)==sums[name],name
    kept.append(dict(name=name,sha256=sums[name],bytes=f.stat().st_size))
(out/'SHA256SUMS').write_text(''.join(f"{k['sha256']}  {k['name']}\n" for k in kept))
with tarfile.open('parton-tables.tgz','w:gz',compresslevel=6) as t:t.add(out,arcname='parton-data')
dropped=sorted(str(p.relative_to(src)) for p in src.rglob('*') if p.is_file() and p.name not in KEEP+('SHA256SUMS',))
record=dict(input_sha256=sha('parton-prepared.tgz'),output_sha256=sha('parton-tables.tgz'),kept=kept,dropped_files=len(dropped),
            dropped_bytes=sum((src/d).stat().st_size for d in dropped),note='kept files are byte-identical to the input bundle (its SHA256SUMS)')
Path('table-only.json').write_text(json.dumps(record,indent=1)+'\n')
print('TABLE_ONLY',json.dumps({k:v for k,v in record.items() if k!='kept'}),flush=True)
