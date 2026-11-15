#!/bin/bash
set -eo pipefail
[[ -n ${_CONDOR_SCRATCH_DIR:-} && $(hostname -s) != cmscon ]] || exit 2
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source generator/env/lcg.sh >/dev/null 2>&1
: "${LHAPDF_DATA_PATH:?}"
python3 rwt/parton/scale_selection.py "${1:?selection}"
