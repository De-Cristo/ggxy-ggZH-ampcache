# Pinned inputs of the generator-level reproduction (sourced by the other scripts).
#
# The reference numbers in reference/generator.json were produced with these commits.
# ggxy_ggZH 0fbaccf is upstream main after the merge of the amplitude cache (MR !1,
# 2026-10-07): commit a22cb4f of the original release plus the cache, nothing else in
# the physics code.

GGXY_URL=https://gitlab.com/ggxy/ggxy-release.git
GGXY_COMMIT=9edd1d5658ebd9f56951f5f48276bf36c41fa5d7          # 2026-03-18

POWHEG_URL=https://gitlab.com/POWHEG-BOX/V2/POWHEG-BOX-V2.git
POWHEG_COMMIT=6654ea3e8cc364b3a628f56aa0ab5efd3c3b114f        # 2026-08-17

GGZH_URL=https://gitlab.com/POWHEG-BOX/V2/User-Processes/ggxy_ggZH.git
GGZH_COMMIT=0fbaccfca6afd56b41c1f4cdcc1876076a3cd27a          # 2026-10-07, includes the cache

# PDF sets (LHAPDF names); both are on /cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current
PDF_MAIN=NNPDF31_nnlo_as_0118_mc_hessian_pdfas                # LHAPDF id 325300
PDF_PAPER=PDF4LHC21_40_pdfas                                  # arXiv:2603.15762, Table 3
