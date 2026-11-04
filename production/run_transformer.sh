#!/bin/bash
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh >/dev/null 2>&1
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
unset PYTHONPATH PYTHONHOME
mkdir -p transformer-logs
python3 -m venv transformer-env
source transformer-env/bin/activate
python -m pip install --no-cache-dir --disable-pip-version-check torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121 > transformer-logs/pip-torch.log 2>&1
python -m pip install --no-cache-dir --disable-pip-version-check numpy==1.26.4 scipy==1.13.1 matplotlib==3.8.4 > transformer-logs/pip-analysis.log 2>&1
python -m pip freeze >transformer-logs/environment.txt
/usr/bin/time -v -o transformer-logs/resources.txt python rwt/transformer/task.py "$@" 2>&1 | tee transformer-logs/train.log
tar czf out-rwt-transformer.tgz transformer-logs transformer-part.json transformer-out/result.json 2>/dev/null || true
