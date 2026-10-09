# Pinned inputs of the shower and analysis step (sourced by the other scripts).

# LCG view with Rivet >= 4.1.2 (needed by the LHC Higgs WG STXS routine), Pythia 8.315,
# Herwig 7.3.0, HepMC3, FastJet, LHAPDF 6.5.5.  The generator step may use another view: the
# LHE files are the interface.
: "${SHOWER_LCG_VIEW:=/cvmfs/sft.cern.ch/lcg/views/LCG_108a/x86_64-el9-gcc13-opt}"

# LHC Higgs WG STXS classification (Rivet routine SimplifiedTemplateCrossSections; stages 0,
# 1.2 and 1.3, coarse and fine).
STXS_URL=https://gitlab.cern.ch/LHCHIGGSXS/LHCHXSWG2/STXS/Classification.git
STXS_COMMIT=38eea57811211fd5895842033653489dc357fb71        # 2026-06-24, Stage 1.3
