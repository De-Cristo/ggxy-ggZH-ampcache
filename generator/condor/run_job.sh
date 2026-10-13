#!/bin/bash
# Generic HTCondor job wrapper for the generator-level steps (EL9 worker with /cvmfs/sft.cern.ch).
#
#   run_job.sh <step> [arguments of the step]
#
# Inputs in the job directory: the package, either unpacked (generator/ present, e.g. by
# production/paper.sh) or as generator.tgz (git archive HEAD VERSION generator | gzip); for every
# step except build and mcfm, ggzh-install.tgz (the product of the build step); sources.tgz
# (fetch_sources.sh) for the build step if the worker cannot reach gitlab.com; grids.tgz for the
# events and cache steps; MCFM-10.3.1.tar.gz for mcfm if the worker cannot reach mcfm.fnal.gov.
# Output: a tarball out-<step>[-<tag>].tgz with results and logs in the job directory.
#
#   build                                  -> ggzh-install.tgz
#   grids  [nseeds] [x-grid iterations]    -> grids.tgz (+ out-grids.tgz with logs)
#   events <nevents> <seed> [weight list|none] -> events-<seed>.lhe.gz + out-events-<seed>.tgz (log, xsec-raw-<seed>.json)
#   cache  <nevents> <seed> [weight list]  -> out-cache-<seed>.tgz
#   cache-scan <nevents> <seed> <weight list|none> ...   -> out-cache-scan-<seed>.tgz (run/cache_scan.sh)
#   fo-members <lo|nlo> <seed> <np> <np_opt> <PDF set> <first member> <n members> [program] [sqrt(s)]
#                                          -> out-fo-members-<first member>.tgz (run/fo_members.sh)
#   fo     <lo|nlo> <first seed> <nseeds> [np] [np_opt] [pdf] [program] [sqrt(s)] [mtscheme] [mu_t/m_ZH]
#                                          -> out-fo-<mode>-<first seed>.tgz
#   lo     [nseeds] [points...]            -> out-lo.tgz
#   mcfm                                   -> out-mcfm.tgz (MCFM_SQRTS from the environment)
set -o pipefail
STEP=${1:?step}; shift
export NCPU=${NCPU:-${OMP_NUM_THREADS:-1}}
echo "run_job.sh $STEP $* on $(hostname) $(date -u +%FT%TZ), NCPU=$NCPU"
[ -d generator ] || tar xzf generator.tgz || { echo "FATAL: generator.tgz"; exit 2; }
G=$PWD/generator
echo "package $(head -1 VERSION 2>/dev/null || echo unknown)"
if [ "$STEP" != build ] && [ "$STEP" != mcfm ]; then
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
    # the LHE file separately (it is the large product); results and logs in the tarball
    [ -f "out/events-$2.lhe.gz" ] && mv "out/events-$2.lhe.gz" .
    tar czf "out-events-$2.tgz" out; exit $rc ;;
  cache)
    "$G/run/cache_ab.sh" grids.tgz out "$1" "$2" "${3:-$G/cards/pwg-rwl.dat}"; rc=$?
    tar czf "out-cache-$2.tgz" out; exit $rc ;;
  fo-members)
    "$G/run/fo_members.sh" "$1" out "$2" "$3" "$4" "$5" "$6" "$7" "${8:-fo-ggllh}" "${9:-13600}"; rc=$?
    tar czf "out-fo-members-$6.tgz" out; exit $rc ;;
  cache-scan)
    "$G/run/cache_scan.sh" grids.tgz out "$1" "$2" "${@:3}"; rc=$?
    tar czf "out-cache-scan-$2.tgz" out; exit $rc ;;
  fo)
    "$G/run/fo.sh" "$1" out "$2" "$3" "${@:4}"; rc=$?
    tar czf "out-fo-$1-$2.tgz" out; exit $rc ;;
  lo)
    "$G/run/powheg_lo.sh" out "${1:-$NCPU}" "${@:2}"; rc=$?
    tar czf out-lo.tgz out/*/pwg-st2-*-stat.dat out/*/powheg.input out/*/run-st2-1.log out/manifest.txt 2>/dev/null; exit $rc ;;
  mcfm)
    mkdir -p mcfm; [ -f MCFM-10.3.1.tar.gz ] && mv MCFM-10.3.1.tar.gz mcfm/
    "$G/mcfm/run_mcfm.sh" mcfm "$NCPU"; rc=$?
    mkdir -p out; for d in mcfm/out/*/; do [ -f "$d/stdout.log" ] && python3 "$G/analysis/mcfm_xsec.py" "$d" --out "out/mcfm-$(basename "$d").json"; done
    tar czf out-mcfm.tgz out mcfm/out; exit $rc ;;
  *) echo "unknown step $STEP"; exit 2 ;;
esac
