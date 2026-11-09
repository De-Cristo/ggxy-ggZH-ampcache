#!/bin/bash
# The effective card comes from the validated ATLAS grid archive, not environment overrides.
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
SEED=${1:?}
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 - <<'CARD'
import json,re,tarfile
from pathlib import Path
with tarfile.open('grids.tgz') as t:
    card=t.extractfile('powheg.input').read().decode()
values={}
for line in card.splitlines():
    line=re.split(r'[!#]',line,1)[0].split()
    if len(line)>=2:values[line[0].lower()]=line[1]
expected={'lhans1':93300,'lhans2':93300,'dynscale':3,'bornonly':1,'loevents':1,'ebeam1':6800,'ebeam2':6800,'vdecaymode':11}
for key,want in expected.items():
    assert float(values[key].replace('d','e').replace('D','E'))==want,(key,values.get(key),want)
Path('atlas-grid-card-check.json').write_text(json.dumps({'pass':True,'expected':expected,'card':card},indent=2)+'\n')
CARD
# POWHEG can report an unreliable status: the wrapper and payload verify the LHE product.
bash generator/condor/run_job.sh events 5000 "$SEED" none
test -s "events-$SEED.lhe.gz"
test -s "out/xsec-raw-$SEED.json"
