#!/bin/bash
# Validation chain for one job: gridpack -> LHE (runcmsgrid.sh) -> Pythia8 -> NanoGEN with HTXS.
#
#   gridpack_to_nanogen.sh <gridpack.tgz> <nevents> <seed> [weights: scale|full] [order: nlo|lo] [outdir = out]
#
# weights=scale keeps only the nine scale weights of the gridpack's list (fast; enough for the
# validation), full keeps the gridpack's list (884 weights for the CMS list: 3.3 s/event with
# the cache).  order=lo makes Born-only events from the NLO gridpack (as for the reference LO
# sample: the stored grids are dropped, bornonly 1 and LOevents 1, the Born is integrated in
# the job, a few minutes) and showers them without the POWHEG veto; with a dedicated LO
# gridpack (card/ggxy_ggZH_ZToLL_LO_13p6TeV) use order=lo-gridpack.
# Environment: CMSSW_12_4_11, el8_amd64_gcc10 (inside cmssw-el8); the fragments of ../fragment.
# HTXS: rivetProducerHTXS.ProductionMode must be GG2ZH; the default AUTO reads the LHE process
# id, does not know ggxy_ggZH and puts every event in category 0.
# Output: <outdir>/nanogen_<seed>.root and the logs.
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GPK=$(readlink -f "${1:?gridpack}"); NEV=${2:?nevents}; SEED=${3:?seed}; RWL=${4:-scale}; ORDER=${5:-nlo}; OUT=${6:-out}
export SCRAM_ARCH=${SCRAM_ARCH:-el8_amd64_gcc10}
CMSSW_RELEASE=${CMSSW_RELEASE:-CMSSW_12_4_11}
die() { echo "FATAL: $*" >&2; exit 1; }
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
T0=$SECONDS

# --- 1. LHE --------------------------------------------------------------------------------
rm -rf gp && mkdir gp && tar xzf "$GPK" -C gp || die "cannot unpack $GPK"
if [ "$RWL" = scale ]; then
  awk 'BEGIN{keep=1} /<weightgroup/ && !/scale_variation/{keep=0} keep{print} /<\/weightgroup>/{if(!keep){keep=1}}' gp/pwg-rwl.dat > gp/pwg-rwl.scale
  grep -q "</initrwgt>" gp/pwg-rwl.scale || echo "</initrwgt>" >> gp/pwg-rwl.scale
  mv gp/pwg-rwl.scale gp/pwg-rwl.dat
fi
echo "weights in the list: $(grep -c '<weight ' gp/pwg-rwl.dat)"
if [ "$ORDER" = lo ]; then
  find gp -maxdepth 1 -name "pwg*" ! -name "pwg-rwl.dat" -delete; rm -f gp/mint_upb_*
  sed -i "/^bornonly/d; /^LOevents/d" gp/powheg.input
  printf "bornonly 1\nLOevents 1\n" >> gp/powheg.input
fi
( cd gp && ./runcmsgrid.sh "$NEV" "$SEED" 1 > lhe.log 2>&1 ) || { tail -20 gp/lhe.log; cp gp/lhe.log "$OUT/lhe-$SEED.log"; die "runcmsgrid"; }
N=$(grep -c "<event>" gp/cmsgrid_final.lhe)
[ "$N" = "$NEV" ] || die "LHE has $N events, wanted $NEV"
echo "LHE: $N events in $((SECONDS - T0)) s; $(grep -h 'ggxy ampcache' gp/lhe.log | tail -1)"
cp gp/lhe.log "$OUT/lhe-$SEED.log"; mv gp/cmsgrid_final.lhe cmsgrid_final.lhe; rm -rf gp

# --- 2. GEN + NanoGEN --------------------------------------------------------------------
source /cvmfs/cms.cern.ch/cmsset_default.sh
[ -d "$CMSSW_RELEASE" ] || scram project CMSSW "$CMSSW_RELEASE" > /dev/null || die "scram project"
FRAG=ggZH_ZToLL_HToBB_ggxy-NLO_pythia8_CP5_hadronizer_cff.py
case $ORDER in lo|lo-gridpack) FRAG=ggZH_ZToLL_HToBB_ggxy-LO_pythia8_CP5_hadronizer_cff.py;; esac
mkdir -p "$CMSSW_RELEASE/src/Configuration/GenProduction/python"
cp -p "$HERE/../fragment/$FRAG" "$CMSSW_RELEASE/src/Configuration/GenProduction/python/"
( cd "$CMSSW_RELEASE/src" && eval "$(scram runtime -sh)" && scram b -j 1 > /dev/null 2>&1 ) || die "scram build"
cd "$CMSSW_RELEASE/src" && eval "$(scram runtime -sh)" && cd - > /dev/null
cmsDriver.py "Configuration/GenProduction/python/$FRAG" \
  --filein file:cmsgrid_final.lhe --filetype LHE --fileout "file:nanogen_$SEED.root" \
  --mc --eventcontent NANOAODGEN --datatier NANOAOD --conditions auto:mc --step GEN,NANOGEN \
  --beamspot Realistic25ns13p6TeVEarly2022Collision --era Run3 --nThreads 1 \
  --customise_commands "process.RandomNumberGeneratorService.generator.initialSeed = ${SEED}\nprocess.rivetProducerHTXS.ProductionMode = 'GG2ZH'" \
  --python_filename gen_cfg.py --no_exec -n -1 > "$OUT/cmsdriver-$SEED.log" 2>&1 || die "cmsDriver ($OUT/cmsdriver-$SEED.log)"
cmsRun gen_cfg.py > "$OUT/cmsrun-$SEED.log" 2>&1 || { tail -30 "$OUT/cmsrun-$SEED.log"; die "cmsRun"; }
[ -s "nanogen_$SEED.root" ] || die "no NanoGEN output"
mv "nanogen_$SEED.root" "$OUT/" && rm -f cmsgrid_final.lhe
echo "nanogen_$SEED.root: $(stat -c %s "$OUT/nanogen_$SEED.root") bytes; total $((SECONDS - T0)) s"
