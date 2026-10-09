#!/bin/bash
# Generic HTCondor job wrapper for the generator-level steps (EL9 worker with /cvmfs/sft.cern.ch).
#
#   run_job.sh <step> [arguments of the step]
#
# Transferred inputs: generator.tgz (tar czf generator.tgz -C <repository> generator) always;
# ggzh-install.tgz (the product of the build step) for every other step; sources.tgz
# (fetch_sources.sh) for the build step if the worker cannot reach gitlab.com; grids.tgz for
# the events and cache steps.  Output: one tarball out-<step>[-<tag>].tgz in the job directory.
#
#   build                                  -> ggzh-install.tgz
#   grids  [nseeds] [x-grid iterations]    -> grids.tgz (+ out-grids.tgz with logs)
#   events <nevents> <seed> [none]         -> out-events-<seed>.tgz (LHE, log, xsec-raw-<seed>.json)
#   cache  <nevents> <seed> [weight list]  -> out-cache-<seed>.tgz
#   fo     <lo|nlo> <first seed> <nseeds> [np] [np_opt] [pdf] [program] -> out-fo-<mode>-<first seed>.tgz
#   lo     [nseeds] [points...]            -> out-lo.tgz
set -o pipefail
STEP=${1:?step}; shift
export NCPU=${NCPU:-${OMP_NUM_THREADS:-1}}
echo "run_job.sh $STEP $* on $(hostname) $(date -u +%FT%TZ), NCPU=$NCPU"
tar xzf generator.tgz || { echo "FATAL: generator.tgz"; exit 2; }
G=$PWD/generator
if [ "$STEP" != build ]; then
  [ -f ggzh-install.tgz ] || { echo "FATAL: ggzh-install.tgz not transferred"; exit 2; }
  mkdir -p install && tar xzf ggzh-install.tgz -C install || exit 2
  export GGZH_PREFIX=$PWD/install
  source "$GGZH_PREFIX/env.sh" || exit 3
fi
case $STEP in
  build)
    if [ -f sources.tgz ]; then tar xzf sources.tgz && export SRC_DIR=$PWD/$(tar tzf sources.tgz | head -1 | cut -d/ -f1); fi
    source "$G/env/lcg.sh" || exit 3
    "$G/build.sh" install "$NCPU" || { tar czf out-build.tgz install/build-logs 2>/dev/null; exit 4; }
    tar czf ggzh-install.tgz -C install ggxy bin env.sh lcg.sh BUILD-MANIFEST.txt
    tar czf out-build.tgz install/build-logs install/BUILD-MANIFEST.txt ;;
  grids)
    "$G/run/powheg_grids.sh" grids "${1:-$NCPU}" "${2:-3}"; rc=$?
    cp -p grids/grids.tgz . 2>/dev/null
    tar czf out-grids.tgz grids/*.log grids/*.json grids/manifest.txt 2>/dev/null; exit $rc ;;
  events)
    "$G/run/powheg_events.sh" grids.tgz out "$1" "$2" "${3:-$G/cards/pwg-rwl.dat}"; rc=$?
    # per-job sums for analysis/lhe_xsec.py --merge
    [ $rc = 0 ] && python3 "$G/analysis/lhe_xsec.py" --files "out/events-$2.lhe.gz" --raw --out "out/xsec-raw-$2.json"
    tar czf "out-events-$2.tgz" out; exit $rc ;;
  cache)
    "$G/run/cache_ab.sh" grids.tgz out "$1" "$2" "${3:-$G/cards/pwg-rwl.dat}"; rc=$?
    tar czf "out-cache-$2.tgz" out; exit $rc ;;
  fo)
    "$G/run/fo.sh" "$1" out "$2" "$3" "${@:4}"; rc=$?
    tar czf "out-fo-$1-$2.tgz" out; exit $rc ;;
  lo)
    "$G/run/powheg_lo.sh" out "${1:-$NCPU}" "${@:2}"; rc=$?
    tar czf out-lo.tgz out/*/pwg-st2-*-stat.dat out/*/powheg.input out/*/run-st2-1.log out/manifest.txt 2>/dev/null; exit $rc ;;
  *) echo "unknown step $STEP"; exit 2 ;;
esac
