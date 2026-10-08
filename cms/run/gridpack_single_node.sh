#!/bin/bash
# Build the CMS gridpack for ggxy_ggZH with CMS genproductions on ONE multicore node: compile
# (run_pwg_condor.py -p 0, which clones ggxy_ggZH and builds ggxy), stages 1-3 with <nseeds>
# seeds in parallel on this node, tarball (run_pwg_condor.py -p 9).  This is how the reference
# gridpack was integrated (16 seeds, 3 x-grid iterations, 4 h 16 min on 16 cores).  Run it as a
# batch job with request_cpus = nseeds, inside an EL8 container (cmssw-el8); never on a login node.
# The standard alternative on lxplus is the condor route of genproductions (README.md).
#
#   gridpack_single_node.sh <genproductions_scripts checkout> <card dir> <name> [nseeds = 16] [nxgrid = 3]
#
#   <genproductions_scripts checkout>: with ggxy_ggZH support (until merged: the fork branch in README.md)
#   <card dir>: e.g. ../card/ggxy_ggZH_ZToLL_13p6TeV (powheg.input, process.dat)
# Environment: CMSSW_12_4_11, el8_amd64_gcc10 (CMSSW_RELEASE / SCRAM_ARCH to change).
# Output: <CMSSW release>/work/<name>/ (the work area, inside the release as genproductions
# requires) and the gridpack ggxy_ggZH_<SCRAM_ARCH>_<CMSSW>_<name>.tgz in the current directory.
set -o pipefail
GP=$(cd "${1:?genproductions_scripts checkout}" && pwd); CARD=$(cd "${2:?card dir}" && pwd); NAME=${3:?name}
NS=${4:-16}; NX=${5:-3}
export SCRAM_ARCH=${SCRAM_ARCH:-el8_amd64_gcc10}
CMSSW_RELEASE=${CMSSW_RELEASE:-CMSSW_12_4_11}
PROCESS=$(cat "$CARD/process.dat")
TOP=$PWD
stamp() { date -u +%FT%TZ; }
die() { echo "FATAL: $*" >&2; exit 1; }
[ -d "$GP/bin/Powheg" ] || die "no bin/Powheg in $GP"
grep -q '"ggxy_ggZH"' "$GP/bin/Powheg/Utilities/helpers.py" || die "this genproductions checkout has no ggxy_ggZH recipe"

echo "[$(stamp)] 1. CMSSW area and genproductions copy"
source /cvmfs/cms.cern.ch/cmsset_default.sh
[ -d "$CMSSW_RELEASE" ] || scram project -n "$CMSSW_RELEASE" CMSSW "$CMSSW_RELEASE" > /dev/null || die "scram project"
cd "$CMSSW_RELEASE/src" && eval "$(scram runtime -sh)" || die "scram runtime"
cd "$TOP"
# the work area must be inside the release: genproductions calls "scram tool info" (LHAPDF,
# Boost, Eigen) from there, which fails anywhere else
W=$TOP/$CMSSW_RELEASE/work; rm -rf "$W"; mkdir -p "$W"
cp -rp "$GP/bin/Powheg/." "$W/"
mkdir -p "$W/cards"; cp -rp "$CARD" "$W/cards/"
CARDIN=$W/cards/$(basename "$CARD")/powheg.input
cd "$W"

echo "[$(stamp)] 2. compile (run_pwg_condor.py -p 0)"
python3 ./run_pwg_condor.py -d 1 -p 0 -i "$CARDIN" -m "$PROCESS" -f "$NAME" > "compile-$NAME.log" 2>&1
[ -x "$W/$NAME/pwhg_main" ] || { grep -n -i error "run_src_$NAME.log" 2>/dev/null | head -20; die "compilation failed (compile-$NAME.log, run_src_$NAME.log)"; }
[ -d "$W/$NAME/lib" ] || die "the ggxy libraries were not shipped into $NAME/lib"
echo "pwhg_main sha256 $(sha256sum "$W/$NAME/pwhg_main" | cut -c1-64)"

echo "[$(stamp)] 3. stages 1-3: $NS seeds, $NX x-grid iterations, on this node"
cd "$W/$NAME"
export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH      # the ggxy libraries first (CMSSW has its own CutTools)
LHAPDF_BASE=$(scram tool info lhapdf | grep LHAPDF_BASE= | cut -d= -f2)
export LHAPDF_DATA_PATH=$("$LHAPDF_BASE/bin/lhapdf-config" --datadir)
sed -e "s/NEVENTS/2000/; s/SEED/1/" "$CARDIN" > powheg.input-save
seq 1 "$NS" > pwgseeds.dat
runthem() {  # runthem <log prefix>
  local i rc=0 pids=()
  for i in $(seq 1 "$NS"); do ( echo "$i" | ./pwhg_main > "$1-$i.log" 2>&1 ) & pids+=($!); done
  for i in "${pids[@]}"; do wait "$i" || rc=1; done
  return $rc
}
for xg in $(seq 1 "$NX"); do
  sed "s/^parallelstage.*/parallelstage 1/; s/^xgriditeration.*/xgriditeration $xg/" powheg.input-save > powheg.input
  runthem "run-st1-$xg" || die "stage 1 iteration $xg"
  echo "[$(stamp)] stage 1 iteration $xg done"
done
sed "s/^parallelstage.*/parallelstage 2/; s/^xgriditeration.*/xgriditeration $NX/" powheg.input-save > powheg.input
runthem run-st2 || die "stage 2"; echo "[$(stamp)] stage 2 done"
sed "s/^parallelstage.*/parallelstage 3/; s/^xgriditeration.*/xgriditeration $NX/" powheg.input-save > powheg.input
runthem run-st3 || die "stage 3"; echo "[$(stamp)] stage 3 done"
grep -h "total (btilde" pwg-st3-0001-stat.dat
cp -p "$CARDIN" powheg.input                         # the tarball step sets the NEVENTS/SEED placeholders itself

echo "[$(stamp)] 4. tarball (run_pwg_condor.py -p 9, keep the test plots)"
cd "$W"
python3 ./run_pwg_condor.py -d 1 -p 9 -i "$CARDIN" -m "$PROCESS" -f "$NAME" -k 1 > "tarball-$NAME.log" 2>&1
T=$(ls "$W"/"${PROCESS}"_*_"$NAME".tgz 2>/dev/null | head -1)
[ -f "$T" ] || die "no tarball (tarball-$NAME.log)"
mv "$T" "$TOP/" && echo "$(sha256sum "$TOP/$(basename "$T")")"
echo "[$(stamp)] GRIDPACK-DONE $TOP/$(basename "$T")"
