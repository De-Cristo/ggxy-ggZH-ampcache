#!/bin/bash
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh >/dev/null 2>&1
unset PYTHONPATH PYTHONHOME
python3 -m venv audit-env
source audit-env/bin/activate
python -m pip install --disable-pip-version-check --no-cache-dir numpy==1.26.4 pyarrow==20.0.0 >audit-install.log 2>&1
python -m pip install --disable-pip-version-check --no-cache-dir torch==2.5.1 --index-url https://download.pytorch.org/whl/cpu >>audit-install.log 2>&1
# LCG also exports a different libtorch. Resolve this pinned wheel's shared libraries first.
TORCH_LIB=$(python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"]+"/torch/lib")')
export LD_LIBRARY_PATH="$TORCH_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
python rwt/parton/stable_audit.py "${1:?channel}"
