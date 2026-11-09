#!/bin/bash
# Full real-shard ATLAS inputs. The reference and passive writer use identical cards/seeds.
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
ORDER=${1:?}; SEED=${2:?}; EXPECTED=${3:?}; SHOWER_SEED=${4:?}
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
MODE=nlo
case $ORDER in lo-atlas) MODE=atlas-lo-veto ;; nlo-a14) ;; *) exit 2 ;; esac
mkdir -p particle-install reference-install particle-logs
tar xzf particle-install.tgz -C particle-install
tar xzf reference-install.tgz -C reference-install
export GGZH_SHOWER_PREFIX=$PWD/reference-install GGZH_SHOWER_TABLE=1
bash shower/run/shower.sh input.lhe.gz reference-out "$MODE" a14-athena "$SHOWER_SEED" atlas-common >particle-logs/reference.log 2>&1
cp reference-out/*.table.csv.gz legacy.csv.gz
source particle-install/env.sh >/dev/null 2>&1
source shower/particle/python_env.sh
gunzip -c input.lhe.gz > input.lhe
/usr/bin/time -v -o particle-logs/shower-resources.txt ggzh-shower --lhe input.lhe --out particle \
 --cmnd shower/cards/common.cmnd --cmnd shower/cards/tune-a14-athena.cmnd \
 --cmnd "shower/cards/$MODE.cmnd" --cmnd shower/cards/atlas-common.cmnd --seed "$SHOWER_SEED" \
 --analysis 'SimplifiedTemplateCrossSections:HiggsProdMode=GG2ZH' --analysis GGZH_NLOPS \
 --analysis GGZH_EVENTTABLE:OUT=particle.table.csv >particle-logs/shower.log 2>&1
gzip particle.table.csv
python3 - <<'VERSIONS' >particle-logs/environment.txt
from importlib.metadata import version
for name,want in [('numpy','2.1.3'),('pyarrow','20.0.0'),('torch','2.7.1')]:
    actual=version(name);print(name,actual);assert actual.split('+')[0]==want,(name,actual,want)
VERSIONS
/usr/bin/time -v -o particle-logs/convert-resources.txt python3 shower/particle/convert_atlas.py \
 --seed "$SEED" --shower-seed "$SHOWER_SEED" --channel zll --order "$ORDER" --expected "$EXPECTED" >particle-logs/convert.log 2>&1
python3 shower/particle/forward.py >particle-logs/forward.log 2>&1
cp shower/particle/schema_atlas.json particle-out/schema.json
cp particle.table.csv.gz particle-out/legacy.csv.gz
cp -p particle-out/validation.json particle-logs/
tar czf out-particle-logs.tgz particle-logs
