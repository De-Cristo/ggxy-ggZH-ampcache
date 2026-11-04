#!/bin/bash
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
source shower/versions.sh
source "$SHOWER_LCG_VIEW/setup.sh" >/dev/null 2>&1
mkdir -p particle-install/bin particle-install/lib
# Retain exactly the already validated production Rivet plugin.
tar xzf shower-install.tgz -C particle-install lib
# Compile only the passive-writer driver against the identical Pythia/Rivet stack.
g++ -O2 shower/src/ggzh_particle_shower.cc -o particle-install/bin/ggzh-shower \
 $(pythia8-config --cxxflags) $(rivet-config --cppflags) $(HepMC3-config --cflags) $(fastjet-config --cxxflags) -std=c++17 \
 $(pythia8-config --libs) $(rivet-config --libs) $(HepMC3-config --libs) $(fastjet-config --libs) -ldl
cat > particle-install/env.sh <<'ENV'
source /cvmfs/sft.cern.ch/lcg/views/LCG_108a/x86_64-el9-gcc13-opt/setup.sh
export GGZH_SHOWER_PREFIX=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export PATH=$GGZH_SHOWER_PREFIX/bin:$PATH
export RIVET_ANALYSIS_PATH=$GGZH_SHOWER_PREFIX/lib${RIVET_ANALYSIS_PATH:+:$RIVET_ANALYSIS_PATH}
ENV
{ cat VERSION; pythia8-config --version; rivet-config --version; fastjet-config --version; sha256sum particle-install/bin/ggzh-shower particle-install/lib/RivetGGZH.so; } > particle-install/BUILD-MANIFEST.txt
tar czf particle-install.tgz -C particle-install bin lib env.sh BUILD-MANIFEST.txt
