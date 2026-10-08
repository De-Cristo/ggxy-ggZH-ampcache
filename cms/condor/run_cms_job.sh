#!/bin/bash
# Generic HTCondor job wrapper for the CMS steps.  Runs the step inside the EL8 container
# (cmssw-el8) unless the worker already is EL8.
#
#   run_cms_job.sh gridpack <card dir name> <name> [nseeds = 16] [nxgrid = 3]
#   run_cms_job.sh nanogen  <nevents> <seed> [scale|full] [nlo|lo|lo-gridpack]
#
# Transferred inputs: cms.tgz (tar czf cms.tgz -C <repository> cms) always; for nanogen the
# gridpack (*.tgz other than cms.tgz).  For gridpack, the genproductions checkout is cloned from
# GENPROD_URL at GENPROD_REF (environment; defaults below) unless genproductions.tgz is transferred.
# Output: the gridpack and out-gridpack.tgz (logs), or out-nanogen-<seed>.tgz.
set -o pipefail
STEP=${1:?step}; shift
GENPROD_URL=${GENPROD_URL:-https://gitlab.cern.ch/lichengz/genproductions_scripts.git}
GENPROD_REF=${GENPROD_REF:-powheg_ggxy_ggzh}
tar xzf cms.tgz || { echo "FATAL: cms.tgz"; exit 2; }
C=$PWD/cms
inner() {
  case $STEP in
    gridpack)
      if [ -f genproductions.tgz ]; then mkdir -p genprod && tar xzf genproductions.tgz -C genprod --strip-components=1
      else git clone -q --branch "$GENPROD_REF" "$GENPROD_URL" genprod || { echo "FATAL: clone $GENPROD_URL"; exit 3; }; fi
      "$C/run/gridpack_single_node.sh" genprod "$C/card/$1" "$2" "${3:-16}" "${4:-3}"; rc=$?
      tar czf out-gridpack.tgz work/*.log work/"$2"/run-st*-1.log work/"$2"/pwg-st*-0001-stat.dat 2>/dev/null
      exit $rc ;;
    nanogen)
      GPK=$(ls *.tgz | grep -v -e '^cms.tgz$' -e '^out-' | head -1)
      [ -n "$GPK" ] || { echo "FATAL: no gridpack transferred"; exit 2; }
      "$C/run/gridpack_to_nanogen.sh" "$GPK" "$1" "$2" "${3:-scale}" "${4:-nlo}" out; rc=$?
      tar czf "out-nanogen-$2.tgz" out; exit $rc ;;
    *) echo "unknown step $STEP"; exit 2 ;;
  esac
}
if grep -q 'release 8' /etc/redhat-release 2>/dev/null || [ -n "${INSIDE_EL8:-}" ]; then
  inner "$@"
else
  export INSIDE_EL8=1
  /cvmfs/cms.cern.ch/common/cmssw-el8 --command-to-run "$(readlink -f "$0")" "$STEP" "$@"
fi
