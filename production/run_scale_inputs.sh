#!/bin/bash
# Complete frozen real subsets; no event generation or shower rerun.
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source generator/env/lcg.sh >/dev/null 2>&1
# Keep the LCG-provided data AND library configuration directories.
# Replacing this path with only PDF grids hides lhapdf.conf.
: "${LHAPDF_DATA_PATH:?LCG LHAPDF search path missing}"
python3 rwt/parton/scale_inputs.py "${1:?channel}"
