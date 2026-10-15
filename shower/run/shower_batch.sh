#!/bin/bash
# Shower a set of LHE files with one mode and tune, several files in parallel, and sum the Rivet
# results exactly (analysis/yoda_summary.py: raw sums per file, then merged).
#
#   GGZH_SHOWER_PREFIX=<prefix> run/shower_batch.sh <outdir> <mode> <tune> [extra cards ...] -- <events-<seed>.lhe.gz> ...
#
# The shower of events-<seed>.lhe.gz uses the random seed <seed> + 100000, so the result does not
# depend on how the files are grouped into jobs or on the order in which they run.  NPAR files at a
# time (default NCPU, else nproc).  Output: <outdir>/<file>.raw.json per file, <outdir>/sum.raw.json
# (all files; combine several jobs with yoda_summary.py --raw --merge), logs in <outdir>/logs/.
# GGZH_SHOWER_TABLE=1 (run/shower.sh): also table-<seed>.csv.gz per file, one row per event, in
# GGZH_SHOWER_TABLE_DIR (default <outdir>; production/paper.sh keeps them out of the result tarball).
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
OUT=${1:?outdir}; MODE=${2:?mode}; TUNE=${3:?tune}; shift 3
EXTRA=(); while [ $# -gt 0 ] && [ "$1" != -- ]; do EXTRA+=("$1"); shift; done
[ "${1:-}" = -- ] && shift
FILES=("$@"); [ ${#FILES[@]} -gt 0 ] || { echo "FATAL: no LHE files" >&2; exit 2; }
[ -n "${GGZH_SHOWER_PREFIX:-}" ] && [ -f "$GGZH_SHOWER_PREFIX/env.sh" ] || { echo "FATAL: set GGZH_SHOWER_PREFIX" >&2; exit 1; }
NPAR=${NPAR:-${NCPU:-$(nproc)}}
mkdir -p "$OUT/logs"; OUT=$(cd "$OUT" && pwd)
TDIR=${GGZH_SHOWER_TABLE_DIR:-$OUT}; mkdir -p "$TDIR"; TDIR=$(cd "$TDIR" && pwd)
echo "$(date -u +%FT%TZ) ${#FILES[@]} files, mode $MODE, tune $TUNE, extra ${EXTRA[*]:-none}, $NPAR in parallel"
one() {  # one <lhe file>
  local f=$1 base seed
  base=$(basename "$f"); base=${base%.gz}; base=${base%.lhe}
  seed=${base#events-}; [[ $seed =~ ^[0-9]+$ ]] || { echo "FATAL: no seed in $f" >&2; return 1; }
  "$HERE/run/shower.sh" "$f" "$OUT/logs" "$MODE" "$TUNE" "$((seed + 100000))" "${EXTRA[@]}" >> "$OUT/logs/batch.txt" 2>&1 || return 1
  ( source "$GGZH_SHOWER_PREFIX/env.sh"
    python3 "$HERE/analysis/yoda_summary.py" "$OUT/logs/$base.$MODE.$TUNE.yoda" --raw --out "$OUT/$base.raw.json" > /dev/null ) \
    && rm -f "$OUT/logs/$base.$MODE.$TUNE.yoda" || return 1
  if [ -n "${GGZH_SHOWER_TABLE:-}" ]; then mv "$OUT/logs/$base.$MODE.$TUNE.table.csv.gz" "$TDIR/table-$seed.csv.gz" || return 1; fi
}
rc=0; running=0
for f in "${FILES[@]}"; do
  one "$f" & running=$((running + 1))
  if [ $running -ge "$NPAR" ]; then wait -n || rc=1; running=$((running - 1)); fi
done
while [ $running -gt 0 ]; do wait -n || rc=1; running=$((running - 1)); done
n=$(ls "$OUT"/*.raw.json 2>/dev/null | grep -vc '/sum.raw.json$')
echo "$(date -u +%FT%TZ) $n of ${#FILES[@]} files summarised"
[ "$n" = ${#FILES[@]} ] || rc=1
( source "$GGZH_SHOWER_PREFIX/env.sh"
  python3 "$HERE/analysis/yoda_summary.py" --raw --merge $(ls "$OUT"/*.raw.json | grep -v '/sum.raw.json$') --out "$OUT/sum.raw.json" ) || rc=1
exit $rc
