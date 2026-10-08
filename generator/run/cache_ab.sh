#!/bin/bash
# Amplitude-cache closure: the same events twice, same seed, with the cache on and off.
#
#   GGZH_PREFIX=<prefix> run/cache_ab.sh <grids dir | grids.tgz> <outdir> <nevents> <seed> [weight list]
#
# Both runs use the same weight list (cards/pwg-rwl.dat by default).  The LHE files must be
# identical, event by event, in the kinematics and in every weight as written (analysis/
# compare_lhe.py; only POWHEG's '#rwgt' comment line may differ, at 1e-15 relative), and the
# cached run should be faster.  Reference (CMS 884-weight list, 9 x 50 events): 450/450 events
# and 397 800/397 800 weights identical, 3.3 s vs 79 s per event.  With the 112-weight list the
# gain is smaller, roughly in proportion to the number of weights.
set -o pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
GRIDS=${1:?grids}; OUT=${2:?outdir}; NEV=${3:?nevents}; SEED=${4:?seed}; RWL=${5:-$GEN/cards/pwg-rwl.dat}
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
GGXY_AMPCACHE=1 "$GEN/run/powheg_events.sh" "$GRIDS" "$OUT/cached" "$NEV" "$SEED" "$RWL" || die "cached run failed"
GGXY_AMPCACHE=0 "$GEN/run/powheg_events.sh" "$GRIDS" "$OUT/uncached" "$NEV" "$SEED" "$RWL" || die "uncached run failed"
cat "$OUT/cached/timing.txt" "$OUT/uncached/timing.txt" | tee "$OUT/timing.txt"
python3 "$GEN/analysis/compare_lhe.py" "$OUT/cached/events-$SEED.lhe.gz" "$OUT/uncached/events-$SEED.lhe.gz" \
  --out "$OUT/compare-$SEED.json"
