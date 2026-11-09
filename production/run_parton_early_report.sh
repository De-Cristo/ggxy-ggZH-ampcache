#!/bin/bash
# Independent early nominal report; no existing training/production wrapper is changed.
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
    assert all((m.name=='parton-data' or m.name.startswith('parton-data/')) and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers());t.extractall('.')
PY
(cd parton-data && sha256sum -c --quiet SHA256SUMS)
/usr/bin/time -v -o parton-logs/resources.txt python rwt/parton/early_report.py 2>&1 | tee parton-logs/task.log
