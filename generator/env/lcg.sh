# Compiler, CMake, Boost, Eigen and LHAPDF from an LCG view on /cvmfs/sft.cern.ch.
# Source it before build.sh, or set LCG_VIEW to another view.
#
# Default: LCG_107, gcc 11. ggxy has been compiled with gcc 10 (CMSSW_12_4_11, used for
# the reference numbers) and gcc 11 (CMSSW_13_0_13); later compilers are untested.
#
# The view ships its own librecola, libcollier, libcuttools and liboneloop, which are
# NOT the versions bundled with ggxy (ggxy modifies Recola and CutTools).  build.sh
# links every program with an RPATH to the ggxy libraries and checks with ldd that
# they resolve there; the env.sh it writes also puts them first in LD_LIBRARY_PATH.

: "${LCG_VIEW:=/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt}"
if [ ! -f "$LCG_VIEW/setup.sh" ]; then
  echo "env/lcg.sh: no LCG view at $LCG_VIEW (is /cvmfs/sft.cern.ch mounted?)" >&2
  return 1 2>/dev/null || exit 1
fi
source "$LCG_VIEW/setup.sh"
export LCG_VIEW
