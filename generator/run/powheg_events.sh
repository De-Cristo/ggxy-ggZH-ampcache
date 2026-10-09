#!/bin/bash
# NLO+PS events (LHE) from the grids of run/powheg_grids.sh: one seed, one core.
#
#   GGZH_PREFIX=<prefix> run/powheg_events.sh <grids dir | grids.tgz> <outdir> <nevents> <seed> [weight list]
#
# Weight list: cards/pwg-rwl.dat by default (9 scale points + NNPDF3.1 Hessian and alpha_s,
# 112 weights; stock POWHEG-BOX-V2 allows up to 500), "none" for no reweighting.
# The amplitude cache is on unless GGXY_AMPCACHE=0 is set.
# Run as the CMS gridpack does: all pwg*.dat of the grids in a fresh directory (with
# pwgubound.dat = the first seed's upper bounds), the grids' own card without the
# parallel-stage lines, numevts and iseed set.  LO grids (GRID_SETTINGS="bornonly=1
# LOevents=1") give Born-only LO events; grids with bornsuppfact give weighted events.  Output:
# <outdir>/events-<seed>.lhe.gz, the POWHEG log and timing.txt.
# Cost per event with the 112-weight list and the cache: a few seconds (reference: 3.3 s with
# the 884-weight CMS list, 2-4 s without any list).
set -o pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
GRIDS=${1:?grids}; OUT=${2:?outdir}; NEV=${3:?nevents}; SEED=${4:?seed}; RWL=${5:-$GEN/cards/pwg-rwl.dat}
load_env
GRIDS=$(readlink -f "$GRIDS"); [ -e "$GRIDS" ] || die "no grids $GRIDS"   # absolute: the run happens in its own directory
RWL=$( [ "$RWL" = none ] && echo none || readlink -f "$RWL" )
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
RUN=$OUT/run-$SEED; rm -rf "$RUN"; mkdir -p "$RUN"; cd "$RUN"
if [ -f "$GRIDS" ]; then tar xzf "$GRIDS" --wildcards 'pwg*.dat'; else cp -p "$GRIDS"/pwg*.dat .; fi
rm -f pwg-rwl.dat pwgseeds.dat pwgevents*.lhe
[ "$(ls pwggrid-*.dat 2>/dev/null | wc -l)" -gt 0 ] || die "no pwggrid-*.dat in $GRIDS"
[ -f pwgubound.dat ] || die "no pwgubound.dat in $GRIDS (made by run/powheg_grids.sh)"
args=(numevts=$NEV iseed=$SEED storeinfo_rwgt=0 pdfreweight=0)
if [ "$RWL" != none ]; then
  cp -p "$RWL" pwg-rwl.dat
  args+=("rwl_file='pwg-rwl.dat'" rwl_format_rwgt=1 rwl_group_events=2000 lhapdf6maxsets=50)
fi
# the card of the grids (with any GRID_SETTINGS), without the parallel-stage lines: single run
if [ -f "$GRIDS" ]; then tar xzf "$GRIDS" powheg.input; else cp -p "$GRIDS/powheg.input" .; fi
mv powheg.input powheg.input-grids
sed -i -e '/^[[:space:]]*manyseeds[[:space:]]/d' -e '/^[[:space:]]*parallelstage[[:space:]]/d' \
       -e '/^[[:space:]]*xgriditeration[[:space:]]/d' powheg.input-grids
make_card powheg.input-grids powheg.input "${args[@]}"
manifest "$OUT/manifest.txt" "powheg_events.sh seed $SEED nevents $NEV weights $(grep -c '<weight ' pwg-rwl.dat 2>/dev/null || echo 0) GGXY_AMPCACHE=${GGXY_AMPCACHE:-1}"
t=$SECONDS
pwhg_main < /dev/null > pwhg.log 2>&1 || { tail -20 pwhg.log >&2; die "pwhg_main failed (seed $SEED)"; }
dt=$((SECONDS - t))
N=$(grep -c "<event>" pwgevents.lhe 2>/dev/null || echo 0)
[ "$N" = "$NEV" ] || die "seed $SEED: $N events, wanted $NEV"
gzip -c pwgevents.lhe > "$OUT/events-$SEED.lhe.gz"
cp -p pwhg.log "$OUT/pwhg-$SEED.log"
echo "seed $SEED cache=${GGXY_AMPCACHE:-1} events $N wall ${dt}s per_event $(awk -v d=$dt -v n=$N 'BEGIN{printf "%.2f", d/n}')s $(grep -h 'ggxy ampcache' pwhg.log | tail -1)" \
  | tee -a "$OUT/timing.txt"
cd "$OUT" && rm -rf "$RUN"
