#!/bin/bash
set -eo pipefail
[[ -n ${_CONDOR_SCRATCH_DIR:-} && $(hostname -s) != cmscon ]] || exit 2
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
bash generator/condor/run_job.sh grids 8 3 || echo 'Grid process returned nonzero; inspect product checks'
python3 - <<'PY'
import json,tarfile,math
from pathlib import Path
with tarfile.open('grids.tgz') as t:
    card=t.extractfile('powheg.input').read().decode()
    xs=json.loads(t.extractfile('grids-xsec.json').read())
    names=t.getnames();assert sum(n.startswith('pwggrid-') and n.endswith('.dat') for n in names)==8
values={r[0].lower():r[1] for line in card.splitlines() if len(r:=line.split())>=2}
expected=dict(lhans1=93300,lhans2=93300,dynscale=3,bornonly=1,loevents=1,vdecaymode=12,ebeam1=6800,ebeam2=6800)
for k,v in expected.items():assert float(values[k].replace('d','e').replace('D','E'))==v,(k,values.get(k))
assert xs['sigma_fb']>0 and math.isfinite(xs['sigma_fb']) and xs['mc_err_fb']>=0
Path('grid-check.json').write_text(json.dumps(dict(accepted=True,expected=expected,card=card,xsec=xs,package=Path('VERSION').read_text().strip()),indent=2)+'\n')
PY
