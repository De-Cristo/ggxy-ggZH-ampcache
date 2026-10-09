#!/bin/bash
# Optional LO cross-check with MCFM 10.3.1: loop-induced gg -> ZH with MCFM's own one-loop
# amplitudes (gg_zh), independent of the Recola amplitudes of ggxy_ggZH.
#
#   [MCFM_SQRTS=<GeV>] run_mcfm.sh <workdir> [threads = 8] [cards = cards/*.ini]
#
# Downloads MCFM-10.3.1.tar.gz (or uses one in <workdir>), checks its sha256, applies
# mcfm-10.3.1-only-ggzh.patch (input switch extra%only_ggzh: the LO ZH matrix element returns
# only gg -> ZH), builds with the LCG_106 view (gcc 13; MCFM's internal LHAPDF) and runs every
# card at LO with the seven-point scale variation.  Build about 20 min, about 100 s per card
# on 8 threads.  Results: analysis/mcfm_xsec.py <workdir>/out/<card>.
# The cards use removebr: the stable-Z results are sigma(gg -> ZH); for the off-shell cards
# multiply by 3 x BR(Z -> ee) (printed as 'only_ggzh: brzee') for Z -> l+l-, l = e, mu, tau.
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WD=${1:?workdir}; NT=${2:-8}; shift 2 2>/dev/null || shift $#
CARDS=("$@"); [ ${#CARDS[@]} -gt 0 ] || CARDS=("$HERE"/cards/*.ini)
mkdir -p "$WD"; WD=$(cd "$WD" && pwd); cd "$WD"; mkdir -p out
URL=https://mcfm.fnal.gov/downloads/MCFM-10.3.1.tar.gz
SHA=db52a8d824309b3573eca3bcfc46016a076092bca31068ceabd63f026b4913fd
[ -f MCFM-10.3.1.tar.gz ] || curl -sSfL -o MCFM-10.3.1.tar.gz "$URL" || { echo "FATAL: download $URL"; exit 2; }
echo "$SHA  MCFM-10.3.1.tar.gz" | sha256sum -c - || { echo "FATAL: MCFM tarball checksum"; exit 2; }

LCGVIEW=${MCFM_LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt}
source "$LCGVIEW/setup.sh" || { echo "FATAL: LCG view $LCGVIEW"; exit 3; }
if [ ! -x MCFM-10.3.1/build/mcfm ]; then
  rm -rf MCFM-10.3.1 && tar xzf MCFM-10.3.1.tar.gz
  ( cd MCFM-10.3.1 && patch -p1 --batch --forward -i "$HERE/mcfm-10.3.1-only-ggzh.patch" ) > out/patch.log 2>&1 \
    || { echo "FATAL: patch"; cat out/patch.log; exit 4; }
  ( cd MCFM-10.3.1 && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -Duse_internal_lhapdf=ON \
      -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran \
    && cmake --build build -j"$NT" ) > out/build.log 2>&1 || { echo "FATAL: MCFM build, see out/build.log"; tail -20 out/build.log; exit 5; }
fi
# LHAPDF sets from the LCG area, plus a global lhapdf.conf (the sets area has none)
mkdir -p lhapdf-cfg
cat > lhapdf-cfg/lhapdf.conf <<'LHCONF'
Verbosity: 1
Interpolator: logcubic
Extrapolator: continuation
ForcePositive: 0
AlphaS_Type: analytic
MZ: 91.1876
MUp: 0.002
MDown: 0.005
MStrange: 0.10
MCharm: 1.29
MBottom: 4.19
MTop: 172.9
Pythia6LambdaV5Compat: true
LHCONF
export LHAPDF_DATA_PATH=$WD/lhapdf-cfg:/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current
export OMP_NUM_THREADS=$NT OMP_STACKSIZE=256M
RC=0
for C in "${CARDS[@]}"; do
  N=$(basename "$C" .ini); R=run-$N; mkdir -p "$R/runs" "out/$N"
  cp -p MCFM-10.3.1/Bin/br.sm1 MCFM-10.3.1/Bin/br.sm2 MCFM-10.3.1/Bin/hto_output.dat MCFM-10.3.1/Bin/process.DAT "$R/"
  cp -p "$C" "$R/input.ini"
  [ -n "${MCFM_SQRTS:-}" ] && sed -i "s/^\([[:space:]]*sqrts[[:space:]]*=\).*/\1 ${MCFM_SQRTS}/" "$R/input.ini"
  grep -m1 sqrts "$R/input.ini"
  ( cd "$R" && "$WD/MCFM-10.3.1/build/mcfm" input.ini > "$WD/out/$N/stdout.log" 2> "$WD/out/$N/stderr.log" ) || RC=1
  find "$R/runs" -name "*total_cross*.txt" -exec cp -p {} "out/$N/" \;
  grep -a "only_ggzh\|Value of integral" "out/$N/stdout.log" | tail -2
done
exit $RC
