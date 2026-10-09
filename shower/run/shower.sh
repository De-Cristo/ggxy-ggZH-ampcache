#!/bin/bash
# Shower one LHE file with Pythia 8 and analyse it with Rivet (STXS + GGZH_NLOPS).
#
#   GGZH_SHOWER_PREFIX=<prefix> run/shower.sh <events.lhe[.gz]> <outdir> <mode> [tune=monash] [seed] [extra cards...]
#
# mode: nlo (POWHEG veto; for NLO events) or lo-default, lo-scalup, lo-scalup-half,
#       lo-scalup-quarter, lo-power (no veto; for Born-only LO events; cards/lo-*.cmnd), or
#       noshower (test: LHE kinematics only).
# extra cards: e.g. ps-variations (cards/ps-variations.cmnd).
# Output: <outdir>/<base>.<mode>.<tune>.yoda and the log; <base> is the LHE file name.
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LHE=${1:?lhe file}; OUT=${2:?outdir}; MODE=${3:?mode}; TUNE=${4:-monash}; SEED=${5:-}
shift 5 2>/dev/null || shift $#
[ -n "${GGZH_SHOWER_PREFIX:-}" ] && [ -f "$GGZH_SHOWER_PREFIX/env.sh" ] || { echo "FATAL: set GGZH_SHOWER_PREFIX (shower/build.sh)" >&2; exit 1; }
source "$GGZH_SHOWER_PREFIX/env.sh"
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
base=$(basename "$LHE"); base=${base%.gz}; base=${base%.lhe}
work=$(mktemp -d "${TMPDIR:-/tmp}/ggzh-shower.XXXXXX")
trap 'rm -rf "$work"' EXIT
case $LHE in *.gz) gunzip -c "$LHE" > "$work/events.lhe";; *) ln -s "$(readlink -f "$LHE")" "$work/events.lhe";; esac
cmnds=(--cmnd "$HERE/cards/common.cmnd" --cmnd "$HERE/cards/tune-$TUNE.cmnd")
case $MODE in
  noshower) cmnds+=(--cmnd "$HERE/cards/lo-default.cmnd" --cmnd "$HERE/cards/noshower.cmnd");;
  *) [ -f "$HERE/cards/$MODE.cmnd" ] || { echo "FATAL: no cards/$MODE.cmnd" >&2; exit 2; }
     cmnds+=(--cmnd "$HERE/cards/$MODE.cmnd");;
esac
for x in "$@"; do cmnds+=(--cmnd "$HERE/cards/$x.cmnd"); done
name=$base.$MODE.$TUNE
t=$SECONDS
( cd "$work" && ggzh-shower --lhe events.lhe --out "$OUT/$name" "${cmnds[@]}" ${SEED:+--seed "$SEED"} ) > "$OUT/$name.log" 2>&1
rc=$?
echo "$name rc $rc wall $((SECONDS - t)) s; $(grep -h '^ggzh-shower:' "$OUT/$name.log" | tail -1)"
[ $rc = 0 ] && [ -s "$OUT/$name.yoda" ] || { tail -20 "$OUT/$name.log" >&2; exit 3; }
