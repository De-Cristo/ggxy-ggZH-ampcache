#!/bin/bash
# Same reference shower/cards/seeds as particle-atlas-zll; no particle writer.
set -eo pipefail
[[ -n ${_CONDOR_SCRATCH_DIR:-} && $(hostname -s) != cmscon ]] || exit 2
CHANNEL=${1:?}; ORDER=${2:?}; SEED=${3:?}; EXPECTED=${4:?}; SHOWER_SEED=${5:?}
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
MODE=nlo
case $ORDER in lo-atlas) MODE=atlas-lo-veto;; nlo-a14) ;; *) exit 2;; esac
mkdir reference-install
tar xzf reference-install.tgz -C reference-install
export GGZH_SHOWER_PREFIX=$PWD/reference-install GGZH_SHOWER_TABLE=1
bash shower/run/shower.sh input.lhe.gz scalar-out "$MODE" a14-athena "$SHOWER_SEED" atlas-common
cp scalar-out/*.table.csv.gz table.csv.gz
python3 shower/analysis/validate_atlas_scalar.py "$CHANNEL" "$ORDER" "$SEED" "$EXPECTED" "$SHOWER_SEED"
