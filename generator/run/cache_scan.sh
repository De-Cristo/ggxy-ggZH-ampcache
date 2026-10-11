#!/bin/bash
# Amplitude-cache timing scan: the same events with weight lists of increasing size, each with the
# cache on and off (run/cache_ab.sh, which also checks that both runs give identical LHE files), all
# in one job, so that the times of all runs compare on the same machine.
#
#   GGZH_PREFIX=<prefix> run/cache_scan.sh <grids dir | grids.tgz> <outdir> <nevents> <seed> <weight list|none> ...
#
# Output per list, with N its number of weights: <outdir>/w<N>/ (timing.txt, compare-<seed>.json,
# logs and manifests of the cached and uncached runs; the LHE files are removed after the comparison).
# Lists: cards/rwl-scan/ (stock POWHEG-BOX-V2 accepts at most 500 weights).
set -o pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
GRIDS=${1:?grids}; OUT=${2:?outdir}; NEV=${3:?nevents}; SEED=${4:?seed}; shift 4
[ $# -gt 0 ] || die "no weight lists"
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
rc=0
for L in "$@"; do
  n=$([ "$L" = none ] && echo 0 || grep -c '<weight ' "$L")
  echo "$(stamp) weight list $L: $n weights"
  "$GEN/run/cache_ab.sh" "$GRIDS" "$OUT/w$n" "$NEV" "$SEED" "$L" || rc=1
  rm -f "$OUT/w$n"/cached/events-*.lhe.gz "$OUT/w$n"/uncached/events-*.lhe.gz
done
exit $rc
