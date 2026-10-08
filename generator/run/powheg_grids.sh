#!/bin/bash
# POWHEG NLO integration grids (stages 1-3) for event generation, all seeds on this node.
#
#   GGZH_PREFIX=<prefix> run/powheg_grids.sh <workdir> [nseeds = 16] [x-grid iterations = 3]
#
# The settings of the CMS gridpack: ncall1 20000 x 3, ncall2 40000 x 2, nubound 20000 per seed,
# 16 seeds, 3 x-grid iterations; took 4 h 16 min on 16 cores.  Run it on a node with that
# many cores (one condor job with request_cpus = 16); do not run it on a login node.
# Output: <workdir> with the grids and upper bounds (pwg*.dat) and grids.tgz, everything
# run/powheg_events.sh needs.  The NLO cross section of the grids is printed at the end
# (error-weighted mean of the seeds' stage-2 results, and POWHEG's stage-3 combination).
set -o pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
WD=${1:?workdir}; NS=${2:-16}; NX=${3:-3}
load_env
mkdir -p "$WD"; WD=$(cd "$WD" && pwd); cd "$WD"
manifest manifest.txt "powheg_grids.sh seeds $NS xgrid $NX"
seq 1 "$NS" > pwgseeds.dat
card() { make_card "$GEN/cards/ggzh_13p6TeV.input" powheg.input numevts=1000 manyseeds=1 "$@"; }
for xg in $(seq 1 "$NX"); do
  card parallelstage=1 xgriditeration=$xg
  run_parallel "$NS" "run-st1-xg$xg" || die "stage 1 iteration $xg failed (run-st1-xg$xg-*.log)"
  echo "$(stamp) stage 1 iteration $xg done"
done
card parallelstage=2 xgriditeration=$NX
run_parallel "$NS" run-st2 || die "stage 2 failed (run-st2-*.log)"
echo "$(stamp) stage 2 done"
card parallelstage=3 xgriditeration=$NX
run_parallel "$NS" run-st3 || die "stage 3 failed (run-st3-*.log)"
echo "$(stamp) stage 3 done"
for f in pwggrid pwgubound pwgbtildeupb pwgremnupb; do
  [ "$(ls $f-*.dat 2>/dev/null | wc -l)" -eq "$NS" ] || die "expected $NS files $f-*.dat"
done
# single-run event generation reads pwgubound.dat and pwg-stat.dat: provide them exactly as the
# CMS gridpack tarball step does (createTarBall_template.sh: the first seed's files)
cp -p "$(ls pwgubound-*.dat | head -1)" pwgubound.dat
f=$(ls pwg-[0-9]*-stat.dat 2>/dev/null | head -1); [ -n "$f" ] && cp -p "$f" pwg-stat.dat
python3 "$GEN/analysis/powheg_stat.py" --points r1f1="$WD" --label grids --out "$WD/grids-xsec.json"
tar czf grids.tgz pwg*.dat powheg.input manifest.txt grids-xsec.json
echo "$(stamp) GRIDS-DONE $WD/grids.tgz"
