#!/bin/bash
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
CHANNEL=${1:?}; ORDER=${2:?}; SEED=${3:?}; EXPECTED=${4:?}
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p particle-install particle-logs
tar xzf particle-install.tgz -C particle-install
source particle-install/env.sh >/dev/null 2>&1
source shower/particle/python_env.sh
# Preserve all current shower cards and order; only output and LHE ordinal instrumentation added.
gunzip -c input.lhe.gz > input.lhe
MODE=nlo; [[ $ORDER != lo-cms ]] || MODE=cms-lo-veto
if [[ $CHANNEL = zvv ]]; then
  mkdir -p reference-install
  tar xzf reference-install.tgz -C reference-install
  # A full same-seed reference replay establishes table alignment for new Zvv shards.
  export GGZH_SHOWER_PREFIX=$PWD/reference-install GGZH_SHOWER_TABLE=1
  bash shower/run/shower.sh input.lhe.gz reference-out "$MODE" cp5 "$((SEED+100000))" cms-common >particle-logs/reference.log 2>&1
  cp reference-out/*.table.csv.gz legacy.csv.gz
  source particle-install/env.sh >/dev/null 2>&1
source shower/particle/python_env.sh
fi
/usr/bin/time -v -o particle-logs/shower-resources.txt ggzh-shower --lhe input.lhe --out particle \
 --cmnd shower/cards/common.cmnd --cmnd shower/cards/tune-cp5.cmnd --cmnd shower/cards/tune-cp5-ecm13600.cmnd \
 --cmnd shower/cards/$MODE.cmnd --cmnd shower/cards/cms-common.cmnd --seed "$((SEED+100000))" \
 --analysis 'SimplifiedTemplateCrossSections:HiggsProdMode=GG2ZH' --analysis GGZH_NLOPS \
 --analysis GGZH_EVENTTABLE:OUT=particle.table.csv >particle-logs/shower.log 2>&1
gzip particle.table.csv
# Dependencies pinned by the production LCG view; avoid per-shard package downloads.
python3 - <<'VERSIONS' >particle-logs/environment.txt
from importlib.metadata import version
for name,want in [('numpy','2.1.3'),('pyarrow','20.0.0'),('torch','2.7.1')]:
 actual=version(name);print(name,actual);assert actual.split('+')[0]==want,(name,actual,want)
VERSIONS
/usr/bin/time -v -o particle-logs/convert-resources.txt python shower/particle/convert.py --seed "$SEED" --channel "$CHANNEL" --order "$ORDER" --expected "$EXPECTED" >particle-logs/convert.log 2>&1
python shower/particle/forward.py >particle-logs/forward.log 2>&1
cp shower/particle/schema.json particle-out/schema.json
cp particle.table.csv.gz particle-out/legacy.csv.gz
cp -p particle-out/validation.json particle-logs/
tar czf out-particle-logs.tgz particle-logs
