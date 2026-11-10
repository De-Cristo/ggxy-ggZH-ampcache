#!/bin/bash
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh >/dev/null 2>&1
unset PYTHONPATH PYTHONHOME
mkdir -p parton-logs
python3 -m venv diagnose-env
source diagnose-env/bin/activate
python -m pip install --no-cache-dir --disable-pip-version-check numpy==1.26.4 >parton-logs/pip-analysis.log 2>&1
python -m pip install --no-cache-dir --disable-pip-version-check torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121 >parton-logs/pip-torch.log 2>&1
python -m pip freeze >parton-logs/environment.txt
/usr/bin/time -v -o parton-logs/resources.txt python rwt/parton/reload_diagnose.py 2>&1 | tee parton-logs/task.log
