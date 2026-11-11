#!/bin/bash
# All original Zll NLO events; existing accepted shard reused without recalculation.
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir install
tar xzf ggzh-install.tgz -C install
export GGZH_PREFIX=$PWD/install
source "$GGZH_PREFIX/env.sh"
python3 generator/analysis/nlo_topup_all.py prepare
if [[ ${NLO_TOPUP_REUSE:-0} != 1 ]]; then
  (
    cd topup-run
    ulimit -t 3240
    /usr/bin/time -f '%U %S %e' -o ../topup-time.txt pwhg_main < seed-index.txt > ../topup-powheg.log 2>&1
  ) || echo "Generator returned nonzero; product validation determines outcome"
fi
python3 generator/analysis/nlo_topup_all.py validate
