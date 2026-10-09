# gg → ZH at NLO: generator-level reproduction (no CMS software)

This directory rebuilds the public codes at pinned commits and reproduces the fixed-order and
event-level results for gg → ℓ⁺ℓ⁻H (ℓ = e, μ, τ) at √s = 13.6 TeV, 15 < m_ℓℓ < 150 GeV,
μ_R = μ_F = m_ZH/2, on-shell m_t = 172.5 GeV, NNPDF3.1 NNLO (LHAPDF 325300).  It needs only an
LCG view from `/cvmfs/sft.cern.ch`, or any compiler + CMake + Boost + Eigen + LHAPDF stack.

| Code | Commit | |
|---|---|---|
| [ggxy](https://gitlab.com/ggxy/ggxy-release) | `9edd1d5` | two-loop amplitudes; bundles Recola, Collier, CutTools, OneLOop, avhlib, CRunDec |
| [POWHEG-BOX-V2](https://gitlab.com/POWHEG-BOX/V2/POWHEG-BOX-V2) | `6654ea3` | |
| [ggxy_ggZH](https://gitlab.com/POWHEG-BOX/V2/User-Processes/ggxy_ggZH) | `0fbaccf` | includes the amplitude cache (merged 2026-10-07) |

**Status: verified end to end** (2026-10-08/09). A fresh build of the pinned commits with the
LCG_107 view (gcc 11.3, LHAPDF 6.5.5) on HTCondor workers re-ran every step with the reference
seeds, and `check.py` passes all 48 checks (`reference/verification-2026-10-09.json`):

- LO (ggxy program, stable Z, POWHEG LO at 7 points): identical to 10⁻¹⁵ with the same seeds,
  although the reference was built with gcc 10.3 / LHAPDF 6.4.0.
- NLO: statistically equal. ggxy 13.0432(38) vs 13.0461(38) fb; half of the 40 seeds agree to
  ~10⁻⁵, the others take a different path in the adaptive integration after floating-point
  differences (up to 0.5 % per seed). Grids 13.0312(70) fb, as the reference. 50 000 events:
  13.054(59) fb, 7-point band +16.7/−14.0 %, PDF 0.96 %, α_s 2.8 %, 0.14 % negative weights.
- Cache closure (cache on/off, same seed, 112 weights): 150/150 events and 16 800/16 800 weights
  identical; 5.5–6× faster with 112 weights (23.9× with the 884-weight CMS list).
- MCFM LO cross-check: reproduced with the standalone patch.

The verification found and fixed five problems of the extracted scripts (silent `set -e` exit,
module directory, HepMC flags of the LCG `pythia8-config`, relative grids path, and the CMS
work-area location in `../cms`).

## What is reproduced

| Result | Reference (fb) | Steps |
|---|---|---|
| ggxy fixed order, LO | 7.0324(10), +26.8/−20.0 % | `run/fo.sh lo`, 8 seeds × 2·10⁶ points |
| ggxy fixed order, NLO | 13.0461(38), +16.6/−14.0 % | `run/fo.sh nlo`, 40 seeds × 10⁶ points |
| K = NLO/LO; K(m_ZH), K(p_T(Z)) | 1.8552(6); K(m_ZH) from 2.94 (200–225 GeV) to 1.52 (700–800 GeV) | `analysis/fo_combine.py` |
| stable Z, LO (`fo-ggzh`) | 70.705(7) NNPDF3.1; 70.601(7) PDF4LHC21 (paper: 70.59(2)) | `run/fo.sh lo … fo-ggzh` |
| POWHEG fixed-order LO, 7 points | 7.0306(19), +26.8/−20.0 % | `run/powheg_lo.sh` |
| POWHEG NLO grids (integral) | 13.031(7) | `run/powheg_grids.sh` |
| NLO events: σ from the weights, 9 scale points | 13.060(30), +16.7/−14.0 % (200 000 events) | `run/powheg_events.sh` |
| NLO events: PDF (Hessian), α_s ± 0.002 | 0.95 %, 2.8 % | `run/powheg_events.sh` (112-weight list) |
| Amplitude cache closure | identical events and weights, cache on vs off | `run/cache_ab.sh` |
| MCFM LO, stable Z (optional) | 71.111(18) NNPDF3.1; 71.003(17) PDF4LHC21 | `mcfm/run_mcfm.sh` |
| arXiv:2603.15762 (PDF4LHC21) | LO 7.026(4), NLO 13.02(1), K 1.853 | comparison only |

All values with their errors are in `reference/generator.json`.  Two corrections to numbers
quoted earlier: the POWHEG LO cross section 7.037(5) fb and the gridpack cross section
13.04 fb were the first seed's results alone; combined over the seeds they are 7.0306(19) and
13.031(7) fb.  The +0.6 % between MCFM and ggxy at LO is not explained (it is not from m_b,
BR(Z→ℓℓ) or α_s).

## Quick start (one machine, EL9, `/cvmfs/sft.cern.ch` mounted)

Do not run any of this on a shared login node; use a batch job (below) or your own machine.

```sh
source env/lcg.sh                       # LCG_107, gcc 11 (LCG_VIEW=... for another view)
./build.sh $PWD/install 8               # 20-40 min, ~20 GB memory at 8 jobs
export GGZH_PREFIX=$PWD/install

run/fo.sh lo  out/fo-lo  2001 8         # ggxy LO, 8 seeds in parallel
run/fo.sh nlo out/fo-nlo 3001 40        # ggxy NLO (use 40 cores, or several machines with other first seeds)
run/powheg_lo.sh out/powheg-lo 8        # POWHEG LO, 7 points x 8 seeds
run/powheg_grids.sh out/grids 16        # NLO grids, 16 seeds (4-5 h on 16 cores)
for s in $(seq 1001 1050); do run/powheg_events.sh out/grids/grids.tgz out/events 1000 $s; done
run/cache_ab.sh out/grids/grids.tgz out/cache 50 51000

python3 analysis/fo_combine.py --lo out/fo-lo/histo_lo_*.txt --nlo out/fo-nlo/histo_nlo_*.txt --out res/fo.json
python3 analysis/powheg_stat.py --points $(for p in r1f1 r0p5f0p5 r0p5f1 r1f0p5 r1f2 r2f1 r2f2; do echo $p=out/powheg-lo/$p; done) --out res/powheg-lo.json
python3 analysis/lhe_xsec.py --files out/events/events-*.lhe.gz --out res/events.json
python3 check.py --fo res/fo.json --powheg-lo res/powheg-lo.json --grids out/grids/grids-xsec.json \
                 --events res/events.json --cache out/cache/compare-51000.json
```

`check.py` passes a value if it agrees with the reference within three standard deviations
(combined MC errors), at every scale point; it also compares the codes with each other
(POWHEG LO vs ggxy LO, grids and events vs ggxy NLO).

## On HTCondor

`condor/run_job.sh` runs one step on an EL9 worker with `/cvmfs/sft.cern.ch`; the submit files
are examples (add your site's runtime and accounting settings).

```sh
tar czf generator.tgz -C <this repository> generator
mkdir -p logs
condor_submit generator/condor/build.sub            # -> ggzh-install.tgz (needs gitlab.com, or sources.tgz)
condor_submit generator/condor/grids.sub            # -> grids.tgz
condor_submit generator/condor/fo.sub               # ggxy LO (2 jobs) and NLO (10 jobs)
condor_submit generator/condor/lo.sub               # POWHEG LO
condor_submit generator/condor/events.sub           # 50 x 1000 NLO events
condor_submit generator/condor/cache.sub            # 3 x 50 events, cache on/off
```

If the workers cannot reach gitlab.com: `./fetch_sources.sh sources` on a machine that can,
and add `sources.tgz` to the build job's input files.

**Cost** (core-hours, measured in the verification): build 2.3 (17 min on 8 cores for ggxy);
ggxy LO 2 and NLO 60; POWHEG LO 20; grids 70–155 (16 seeds; 4.3 h on fast, 9.6 h on slow
cores); 50 000 NLO events with 112 weights 18 (1.3 s per event with the cache); cache closure 3;
MCFM 2. About 200 core-hours in total.

## Options for production samples

- **LO events** (Born only): grids made with `GRID_SETTINGS="bornonly=1 LOevents=1"
  run/powheg_grids.sh <dir>`; `run/powheg_events.sh` then gives LO events (it always uses the
  card stored with the grids).
- **Born suppression** (weighted events that populate the high-p_T tail):
  `GRID_SETTINGS="bornsuppfact=250"`, optionally `bornsuppexp` (default 1) and
  `bornsuppfloor` (default 0.1): the generation density is multiplied by
  ((p_T² + (rL)²)/(p_T² + L²))ⁿ in p_T(H) of the underlying Born, and every event carries the
  weight 1/fact. Implemented by `patches/ggxy_ggZH-born-suppression.patch` (the process ships
  a stub that always returns 1), which `build.sh` applies; without `bornsuppfact` it changes
  nothing. With suppression, the integrals printed by POWHEG are those of the suppressed
  function; the cross section is Σw/N of the events. *Not yet validated (paper plan WP1).*
- **Another PDF**: e.g. `GRID_SETTINGS="lhans1=93300 lhans2=93300"` (PDF4LHC21_40_pdfas) and a
  weight list with that set: `cards/pwg-rwl-pdf4lhc21.dat` (9 scale points at PDF4LHC21,
  PDF4LHC21_40_pdfas 43 members, NNPDF3.1 103 members: 155 weights); fixed order with the
  set name as the 7th argument of `run/fo.sh`.
- **Another energy**: POWHEG with `GRID_SETTINGS="ebeam1=6500d0 ebeam2=6500d0"` (13 TeV; also
  for `run/powheg_lo.sh`), the ggxy programs with √s in GeV as the 9th argument of
  `run/fo.sh`, MCFM with `MCFM_SQRTS=13000`. `check.py --paper-energy 13600|13000` compares a
  PDF4LHC21 reproduction with Tables 1 and 3 of arXiv:2603.15762 at that energy. The shower
  picks the energy-dependent part of the CP5 tune from the beams of the LHE file.

## Notes

- **Other Recola/Collier/CutTools/OneLOop builds.** ggxy bundles modified versions. LCG views
  (and CMSSW) ship their own, which must not be used: `build.sh` configures ggxy without
  `CMAKE_PREFIX_PATH` (otherwise Recola is compiled against the view's Collier), links every
  program with an RPATH to the ggxy libraries, and checks with `ldd` where they resolve;
  `env.sh` puts the ggxy libraries first in `LD_LIBRARY_PATH`.
- **Compilers.** ggxy has been built with gcc 10 and 11. Fortran module files are compiler
  specific: build ggxy and POWHEG with the same compiler.
- **Event normalisation.** POWHEG events carry the weight σ_abs × sign; the cross section is the
  sum of the weights over the number of events, for the nominal and every alternative weight
  (`analysis/lhe_xsec.py`). Event generation reads the first seed's upper bounds
  (`pwgubound.dat`), as CMS gridpacks do.
- **Weight list.** `cards/pwg-rwl.dat`: 9 scale points and the 103 NNPDF3.1 members (Hessian +
  α_s), 112 weights; stock POWHEG-BOX-V2 allows 500. The 884-weight CMS list needs the CMS
  changes to POWHEG (see `../cms`).
- **Cache.** On by default; `GGXY_AMPCACHE=0` disables it. Valid for reweighting PDFs, μ_R, μ_F
  and the m_t scheme only. "Identical" means the same text in the LHE file (weights are written
  with 5 significant digits); POWHEG's internal `#rwgt` line differs at the 10⁻¹⁵ level.
- **Kinematic limits.** The two-loop amplitudes are expanded in the external masses: keep m_ℓℓ
  and m_H below 150 GeV.
- **ggxy output units.** `fo-ggllh` writes nb per unit of the observable for one lepton flavour;
  `fo_combine.py --nlep 3` converts to fb for ℓ = e, μ, τ.
