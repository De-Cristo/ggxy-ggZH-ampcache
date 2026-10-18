#!/bin/bash
# Build the pinned ggxy library, POWHEG-BOX-V2 with the ggxy_ggZH process (pwhg_main) and the
# three ggxy fixed-order programs (bin/fo-ggllh, bin/fo-ggzh, and bin/fo-ggvvh = fo-ggllh for Z -> nu nubar) into <prefix>.
#
#   source env/lcg.sh            # or provide your own toolchain, see README.md
#   ./build.sh <prefix> [make jobs = 8]
#
# Requirements: gcc/g++/gfortran, cmake >= 3.16, Boost, Eigen 3, LHAPDF 6 (lhapdf-config in
# PATH), git, and access to gitlab.com.  Without network, clone the three repositories of
# versions.sh into $SRC_DIR beforehand (ggxy-release, POWHEG-BOX-V2, POWHEG-BOX-V2/ggxy_ggZH);
# the commits are checked out from there.  The ggxy build needs about 2.5 GB of memory per
# make job.  Takes 20-40 min on 8 cores.
#
# Output: <prefix>/ggxy (library install), <prefix>/bin/{pwhg_main,fo-ggllh,fo-ggzh,fo-ggvvh},
# <prefix>/env.sh (source it before running anything) and <prefix>/BUILD-MANIFEST.txt.
set -Eeo pipefail
trap 'echo "FATAL: build.sh line $LINENO: $BASH_COMMAND" >&2' ERR
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/versions.sh"
[ -n "$1" ] || { echo "usage: $0 <prefix> [make jobs]" >&2; exit 2; }
mkdir -p "$1"; PREFIX=$(cd "$1" && pwd)
JOBS=${2:-8}
SRC=${SRC_DIR:-$PREFIX/src}
mkdir -p "$SRC" "$PREFIX/bin"
LOG=$PREFIX/build-logs; mkdir -p "$LOG"
stamp() { date -u +%FT%TZ; }

for t in gcc g++ gfortran cmake git lhapdf-config; do
  command -v $t >/dev/null || { echo "FATAL: $t not found (source env/lcg.sh?)" >&2; exit 2; }
done
LHAPDF_PREFIX=$(lhapdf-config --prefix)
# Boost and Eigen: from the LCG view by default, otherwise BOOST_ROOT / EIGEN_ROOT
BOOST_ROOT=${BOOST_ROOT:-${LCG_VIEW:-}}
EIGEN_ROOT=${EIGEN_ROOT:-${LCG_VIEW:-}}
[ -f "$BOOST_ROOT/include/boost/version.hpp" ] || { echo "FATAL: Boost headers not found; set BOOST_ROOT" >&2; exit 2; }
[ -d "$EIGEN_ROOT/include/eigen3" ] || { echo "FATAL: Eigen not found; set EIGEN_ROOT" >&2; exit 2; }
EIGEN_CMAKE=$(ls -d "$EIGEN_ROOT"/share/eigen3/cmake "$EIGEN_ROOT"/lib*/cmake/eigen3 2>/dev/null | head -1 || true)
BOOST_CMAKE=$(ls -d "$BOOST_ROOT"/lib*/cmake/Boost-* 2>/dev/null | sort -V | tail -1 || true)
# The LCG views also contain Recola and Collier (with CMake configs).  Recola's CMake looks
# for Collier with find_package(collier HINTS <ggxy's>), and the CMAKE_PREFIX_PATH environment
# variable takes precedence over HINTS: Recola would be compiled against the wrong Collier.
# So ggxy is configured with CMAKE_PREFIX_PATH unset and Boost/Eigen given explicitly.
NOPREFIX=(env -u CMAKE_PREFIX_PATH -u collier_DIR -u recola_DIR -u COLLIER_PATH)

fetch() {  # fetch <dir> <url> <commit>
  [ -d "$1/.git" ] || git clone -q "$2" "$1"
  git -C "$1" checkout -q "$3" || { echo "FATAL: cannot check out $3 in $1" >&2; exit 3; }
  echo "$(basename "$1") $(git -C "$1" rev-parse HEAD)"
}

echo "[$(stamp)] 1. ggxy $GGXY_COMMIT"
fetch "$SRC/ggxy-release" "$GGXY_URL" "$GGXY_COMMIT"
CMAKE_EXTRA=()
# Recola and Collier need this with CMake >= 4 (ggxy CMakeLists.txt says so)
[ "$(cmake --version | head -1 | awk '{print $3}' | cut -d. -f1)" -ge 4 ] && CMAKE_EXTRA+=(-DCMAKE_POLICY_VERSION_MINIMUM=3.5)
rm -rf "$SRC/ggxy-release/build"; mkdir -p "$SRC/ggxy-release/build"
(
  cd "$SRC/ggxy-release/build"
  # ggxy's CMake finds Boost but does not add its include directory to the target
  export CXXFLAGS="-I$BOOST_ROOT/include ${CXXFLAGS:-}"
  "${NOPREFIX[@]}" cmake .. -DCMAKE_INSTALL_PREFIX="$PREFIX/ggxy" -DCMAKE_BUILD_TYPE=Release \
    -DLHAPDFPATH="$LHAPDF_PREFIX" -DggHH=Off -DggZH=On -DonlyFF=Off \
    -DBOOST_ROOT="$BOOST_ROOT" ${BOOST_CMAKE:+-DBoost_DIR="$BOOST_CMAKE"} ${EIGEN_CMAKE:+-DEigen3_DIR="$EIGEN_CMAKE"} \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran "${CMAKE_EXTRA[@]}" \
    > "$LOG/ggxy-cmake.log" 2>&1
  "${NOPREFIX[@]}" make -j "$JOBS" > "$LOG/ggxy-make.log" 2>&1
  "${NOPREFIX[@]}" make install > "$LOG/ggxy-install.log" 2>&1
) || { echo "FATAL: ggxy build failed, see $LOG/ggxy-*.log" >&2; tail -20 "$LOG"/ggxy-make.log >&2; exit 4; }
GGLIB=$PREFIX/ggxy/lib
for l in ggxy crundec recola collier avhlib cuttools oneloop; do
  [ -f "$GGLIB/lib$l.so" ] || { echo "FATAL: $GGLIB/lib$l.so missing" >&2; exit 4; }
done
# Recola prints the Collier it was configured against; it must be ggxy's own
CLL=$(grep -h "Collier library path" "$LOG/ggxy-make.log" | head -1 || true)
case "$CLL" in
  *"$PREFIX"*|*"$SRC/ggxy-release"*) echo "Recola uses ggxy's Collier: ${CLL##*: }";;
  "") echo "WARNING: could not find Recola's Collier path in the build log; check $LOG/ggxy-make.log" >&2;;
  *) echo "FATAL: Recola was configured against another Collier: $CLL" >&2; exit 4;;
esac

echo "[$(stamp)] 2. POWHEG-BOX-V2 $POWHEG_COMMIT + ggxy_ggZH $GGZH_COMMIT"
fetch "$SRC/POWHEG-BOX-V2" "$POWHEG_URL" "$POWHEG_COMMIT"
fetch "$SRC/POWHEG-BOX-V2/ggxy_ggZH" "$GGZH_URL" "$GGZH_COMMIT"
P=$SRC/POWHEG-BOX-V2/ggxy_ggZH
git -C "$P" checkout -q -- .          # pristine process sources (re-runs)
# local patches (patches/*.patch); each is inert unless enabled in the card:
#   ggxy_ggZH-born-suppression.patch  Born suppression in pT(H) for weighted events (bornsuppfact)
for pf in "$HERE"/patches/*.patch; do
  [ -f "$pf" ] || continue
  patch -p1 -d "$P" --forward --batch < "$pf" > "$LOG/patch-$(basename "$pf").log" 2>&1 \
    || { echo "FATAL: patch $(basename "$pf") does not apply" >&2; cat "$LOG/patch-$(basename "$pf").log" >&2; exit 5; }
  echo "applied $(basename "$pf")"
done
grep -q '^ANALYSIS = default' "$P/Makefile" && grep -q '^LIBS+=-L$(GGXYLIB) -lggxy' "$P/Makefile" && grep -q '^CFLAGS=' "$P/Makefile" \
  && grep -q '^LIBHEPMC=' "$P/Makefile" \
  || { echo "FATAL: unexpected ggxy_ggZH Makefile; the edits below need revisiting" >&2; exit 5; }
# - no analysis (no FastJet); the results come from the LHE files and the stat files
# - RPATH to the ggxy libraries, as DT_RPATH, so that other Recola/CutTools/Collier/OneLOop
#   libraries in LD_LIBRARY_PATH (LCG views, CMSSW) are not picked up
# - LHAPDF include directory (ggxy headers include LHAPDF/LHAPDF.h)
# - no HepMC flags: they are only for the Pythia8 interface, and "pythia8-config --hepmc" of
#   some installations (LCG views) prints a warning that would end up on the compiler line
sed -i -e 's|^ANALYSIS = default|ANALYSIS = none|' \
       -e 's|^LIBHEPMC=.*|LIBHEPMC=|' \
       -e 's|^LIBS+=-L$(GGXYLIB) |LIBS+=-L$(GGXYLIB) -Wl,-rpath,$(GGXYLIB) -Wl,--disable-new-dtags |' \
       -e 's|^CFLAGS=\(.*\)$|CFLAGS=\1 -I$(shell $(LHAPDF_CONFIG) --incdir)|' "$P/Makefile"
(
  cd "$P"
  export GGXYPATH=$PREFIX/ggxy LD_LIBRARY_PATH=$GGLIB:${LD_LIBRARY_PATH:-}
  make clean > /dev/null 2>&1 || true
  mkdir -p obj-gfortran mod-gfortran              # the .f90 rule does not create the module directory
  make ggxy.o > "$LOG/pwhg-make.log" 2>&1        # Fortran module first, other objects use it
  make -j "$JOBS" pwhg_main >> "$LOG/pwhg-make.log" 2>&1
) || { echo "FATAL: pwhg_main build failed, see $LOG/pwhg-make.log" >&2; grep -n -i error "$LOG/pwhg-make.log" | head >&2; exit 6; }
cp -p "$P/pwhg_main" "$PREFIX/bin/pwhg_main"

echo "[$(stamp)] 3. ggxy fixed-order programs"
# Compiled directly (the upstream example CMakeLists uses find_package(recola), which an LCG
# view on CMAKE_PREFIX_PATH would satisfy with its own Recola).
INC=$PREFIX/ggxy/include
FOFLAGS=(-std=c++14 -O3 -ffast-math -fno-finite-math-only -fpic -D GNU_NAME_MANGLING)
for d in ggxy ggxy/modules crundec avhlib oneloop cuttools collier recola; do FOFLAGS+=(-I"$INC/$d"); done
FOFLAGS+=(-I"$(lhapdf-config --incdir)" -I"$BOOST_ROOT/include" -I"$EIGEN_ROOT/include" -I"$EIGEN_ROOT/include/eigen3")
for prog in fo-ggllh fo-ggzh fo-ggvvh; do
  src=$prog; def=()
  # fo-ggvvh: gg -> nu nubar H, the l-l+ program with the decay mode switched at compile time
  [ "$prog" = fo-ggvvh ] && { src=fo-ggllh; def=(-DGGZH_DECAYMODE='"vv"'); }
  g++ "${FOFLAGS[@]}" "${def[@]}" "$HERE/fo/$src.cpp" -o "$PREFIX/bin/$prog" \
      -L"$GGLIB" -lggxy -lcrundec -lavhlib -loneloop -lcuttools -lcollier -lrecola \
      -L"$(lhapdf-config --libdir)" -lLHAPDF -lgfortran \
      -Wl,-rpath,"$GGLIB" -Wl,--disable-new-dtags > "$LOG/$prog.log" 2>&1 \
    || { echo "FATAL: $prog failed, see $LOG/$prog.log" >&2; tail -20 "$LOG/$prog.log" >&2; exit 7; }
done

echo "[$(stamp)] 4. environment and checks"
# env.sh is relocatable: the installation (ggxy/, bin/, env.sh, lcg.sh) can be packed and
# unpacked elsewhere, e.g. on batch workers with the same LCG view (condor/run_job.sh).
cp -p "$HERE/env/lcg.sh" "$PREFIX/lcg.sh"
{
  echo "# generated by build.sh $(stamp); source before running"
  if [ -n "${LCG_VIEW:-}" ]; then echo "export LCG_VIEW=$LCG_VIEW"; echo "source \"\$(dirname \"\${BASH_SOURCE[0]}\")/lcg.sh\""; fi
  echo "export GGZH_PREFIX=\$(cd \"\$(dirname \"\${BASH_SOURCE[0]}\")\" && pwd)"
  echo "export PATH=\$GGZH_PREFIX/bin:\$PATH"
  echo "# the ggxy libraries first: other Recola/CutTools/Collier/OneLOop builds may be in the path"
  echo "export LD_LIBRARY_PATH=\$GGZH_PREFIX/ggxy/lib:\${LD_LIBRARY_PATH:-}"
} > "$PREFIX/env.sh"
source "$PREFIX/env.sh"
for b in pwhg_main fo-ggllh fo-ggzh fo-ggvvh; do
  [ -x "$PREFIX/bin/$b" ] || { echo "FATAL: $b not built" >&2; exit 8; }
  if ldd "$PREFIX/bin/$b" | grep -q "not found"; then echo "FATAL: $b has unresolved libraries" >&2; ldd "$PREFIX/bin/$b" | grep "not found" >&2; exit 8; fi
  for l in recola collier cuttools oneloop avhlib crundec ggxy; do
    where=$(ldd "$PREFIX/bin/$b" | awk -v l="lib$l.so" '$1 ~ "^"l {print $3}')
    [ -z "$where" ] || [ "$(dirname "$where")" = "$GGLIB" ] \
      || { echo "FATAL: $b loads lib$l.so from $where, not from $GGLIB" >&2; exit 8; }
  done
done
{
  echo "built $(stamp) on $(hostname)"
  v=$(head -1 "$HERE/../VERSION" 2>/dev/null)
  case $v in ''|*Format*) v="$(git -C "$HERE" rev-parse HEAD 2>/dev/null || echo unknown)";; esac
  echo "package $v"
  echo "ggxy $(git -C "$SRC/ggxy-release" rev-parse HEAD)"
  echo "POWHEG-BOX-V2 $(git -C "$SRC/POWHEG-BOX-V2" rev-parse HEAD)"
  echo "ggxy_ggZH $(git -C "$P" rev-parse HEAD)"
  for pf in "$HERE"/patches/*.patch; do [ -f "$pf" ] && echo "patch $(basename "$pf") sha256 $(sha256sum "$pf" | cut -c1-16)"; done
  echo "LCG_VIEW ${LCG_VIEW:-none}"
  echo "gcc $(gcc --version | head -1)"
  echo "gfortran $(gfortran --version | head -1)"
  echo "cmake $(cmake --version | head -1)"
  echo "LHAPDF $(lhapdf-config --version) at $LHAPDF_PREFIX"
  echo "Boost $(grep -m1 '#define BOOST_LIB_VERSION' "$BOOST_ROOT/include/boost/version.hpp" | awk '{print $3}')"
  for b in pwhg_main fo-ggllh fo-ggzh fo-ggvvh; do echo "sha256 $(sha256sum "$PREFIX/bin/$b")"; done
  for l in "$GGLIB"/lib*.so; do echo "sha256 $(sha256sum "$l")"; done
} > "$PREFIX/BUILD-MANIFEST.txt"
cat "$PREFIX/BUILD-MANIFEST.txt"
echo "[$(stamp)] BUILD-DONE  (source $PREFIX/env.sh)"
