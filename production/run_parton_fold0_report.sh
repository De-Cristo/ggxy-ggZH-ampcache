#!/bin/bash
# Real single-held-out-fold report only; no fitting or new events.
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh >/dev/null 2>&1
unset PYTHONPATH PYTHONHOME
mkdir -p parton-logs
python3 -m venv early-env
source early-env/bin/activate
python -m pip install --no-cache-dir --disable-pip-version-check numpy==1.26.4 scipy==1.13.1 matplotlib==3.8.4 scikit-learn==1.2.2 >parton-logs/pip-analysis.log 2>&1
python -m pip freeze >parton-logs/environment.txt
python - <<'PY'
from pathlib import Path
import tarfile
with tarfile.open('parton-prepared.tgz') as t:
    assert all((m.name=='parton-data' or m.name.startswith('parton-data/')) and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers())
    needed={'parton-data/lo-cms.npz','parton-data/nlo-cp5.npz','parton-data/audit.json','parton-data/SHA256SUMS'}
    t.extractall('.',members=[m for m in t.getmembers() if m.name in needed])
PY
python - <<'VERIFY'
import hashlib
from pathlib import Path
root=Path('parton-data');entries={line.split(maxsplit=1)[1]:line.split()[0] for line in (root/'SHA256SUMS').read_text().splitlines()}
for name in ('lo-cms.npz','nlo-cp5.npz','audit.json'):
    h=hashlib.sha256()
    with (root/name).open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    assert h.hexdigest()==entries[name]
VERIFY
/usr/bin/time -v -o parton-logs/resources.txt python rwt/parton/fold0_report.py 2>&1 | tee parton-logs/task.log
