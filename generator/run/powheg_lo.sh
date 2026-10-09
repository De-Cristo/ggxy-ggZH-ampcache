#!/bin/bash
# POWHEG fixed-order LO (Born only) at the seven scale points: stage 1 (3 x-grid iterations)
# and stage 2 with <nseeds> seeds in parallel on this node, one directory per point.
#
#   GGZH_PREFIX=<prefix> [GRID_SETTINGS="key=value ..."] run/powheg_lo.sh <outdir> [nseeds = 8] [points ...]
#
# Points: r1f1 r0p5f0p5 r0p5f1 r1f0p5 r1f2 r2f1 r2f2 (renscfact/facscfact).  Integration
# statistics of the reference: ncall1 50000 x 3, ncall2 200000 x 3; about 20 min per point
# with 8 seeds on 8 cores.  GRID_SETTINGS as for run/powheg_grids.sh (e.g. "lhans1=93300
# lhans2=93300 ebeam1=6500d0 ebeam2=6500d0" for PDF4LHC21 at 13 TeV).  Combine with
#   analysis/powheg_stat.py --points r1f1=<outdir>/r1f1 r0p5f0p5=<outdir>/r0p5f0p5 ...
set -o pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
OUT=${1:?outdir}; NS=${2:-8}; shift 2 2>/dev/null || shift $#
POINTS=("$@"); [ ${#POINTS[@]} -gt 0 ] || POINTS=(r1f1 r0p5f0p5 r0p5f1 r1f0p5 r1f2 r2f1 r2f2)
load_env
read -r -a EXTRA <<< "${GRID_SETTINGS:-}"
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
manifest "$OUT/manifest.txt" "powheg_lo.sh points ${POINTS[*]} seeds $NS settings: ${GRID_SETTINGS:-none}"
fac() { local x=${1/p/.}; echo "${x}d0"; }
for P in "${POINTS[@]}"; do
  [[ $P =~ ^r([0-9p]+)f([0-9p]+)$ ]] || die "bad point name $P"
  R=$(fac "${BASH_REMATCH[1]}"); F=$(fac "${BASH_REMATCH[2]}")
  D=$OUT/$P; mkdir -p "$D"; cd "$D"
  seq 1 "$NS" | awk '{print $1*1000+17}' > pwgseeds.dat
  t=$SECONDS
  for xg in 1 2 3; do
    make_card "$GEN/cards/ggzh_13p6TeV.input" powheg.input "${EXTRA[@]}" numevts=1000 bornonly=1 renscfact=$R facscfact=$F \
      ncall1=50000 itmx1=3 ncall2=200000 itmx2=3 manyseeds=1 parallelstage=1 xgriditeration=$xg
    run_parallel "$NS" "run-st1-xg$xg" || die "$P stage 1 iteration $xg failed (see $D/run-st1-xg$xg-*.log)"
  done
  make_card "$GEN/cards/ggzh_13p6TeV.input" powheg.input "${EXTRA[@]}" numevts=1000 bornonly=1 renscfact=$R facscfact=$F \
    ncall1=50000 itmx1=3 ncall2=200000 itmx2=3 manyseeds=1 parallelstage=2 xgriditeration=3
  run_parallel "$NS" run-st2 || die "$P stage 2 failed (see $D/run-st2-*.log)"
  n=$(ls pwg-st2-*-stat.dat 2>/dev/null | wc -l)
  echo "$(stamp) $P (renscfact $R, facscfact $F): $n/$NS stage-2 results in $((SECONDS - t)) s"
  [ "$n" -eq "$NS" ] || die "$P: $n of $NS stage-2 stat files"
done
