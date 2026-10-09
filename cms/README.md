# gg → ZH at NLO+PS in CMS: gridpack, generator chain and validation

This directory holds what a CMS user needs to produce and validate gg → Z(ℓ⁺ℓ⁻)H samples at
NLO+PS with the POWHEG-BOX-V2 process `ggxy_ggZH`:

| | |
|---|---|
| `card/ggxy_ggZH_ZToLL_13p6TeV/` | the gridpack card (NLO); `card/ggxy_ggZH_ZToLL_LO_13p6TeV/` its LO companion |
| `fragment/` | Pythia8 fragments: full (ExternalLHEProducer + Pythia8, for requests) and hadronizer-only (validation), NLO and LO |
| `run/gridpack_single_node.sh` | gridpack with CMS genproductions, integration on one multicore node |
| `run/gridpack_to_nanogen.sh` | one validation job: gridpack → LHE → Pythia8 → NanoGEN with HTXS |
| `validation/nanogen_summary.py`, `check_cms.py` | summary of a NanoGEN sample and comparison with `reference/cms.json` |
| `condor/` | HTCondor examples |

The generator itself is upstream: `ggxy_ggZH` contains the amplitude cache since commit `0fbaccf`
(2026-10-07). The gridpack recipe for CMS genproductions (`bin/Powheg/Utilities/helpers.py`:
clone `ggxy_ggZH` at `0fbaccf`, build ggxy `9edd1d5` with the release's compiler, ship its
libraries in `lib/`) is on the branch `powheg_ggxy_ggzh` of
`gitlab.cern.ch/lichengz/genproductions_scripts` until it is merged into CMS genproductions.

**Status.** The reference (`reference/cms.json`) comes from a gridpack built with the same card
and the same process code (`ggxy_ggZH` `a22cb4f` + the cache; `0fbaccf` differs only by the
authors' Born-key change, which has no effect for the on-shell top-mass scheme), but with a
POWHEG-BOX-V2 core repacked from git `6654ea3`. The genproductions recipe (fork branch
`powheg_ggxy_ggzh`: `0fbaccf`, ggxy built in the release, CMS's own POWHEG-BOX-V2 tarball
`powhegboxV2_rev4122_date20250213`) was verified on 2026-10-09 with `run/gridpack_single_node.sh`
(CMSSW_12_4_11, 16 seeds): it compiles, integrates and packs a gridpack whose grid integral is
13.0315(70) fb, equal to the reference gridpack's to 10⁻⁷. The NanoGEN validation of that
gridpack (`check_cms.py`) has not been run yet. The scripts here are extracted from the validated
workflow; `nanogen_summary.py` and `check_cms.py` were run on the original NanoGEN files (a
subset passes 21/21 checks; the LO subset checked against the NLO reference fails 19 checks, as
it should).

## 1. Gridpack

Environment of the reference: CMSSW_12_4_11, `el8_amd64_gcc10` (inside `cmssw-el8`).

**Standard genproductions route (lxplus condor):**

```sh
cmsrel CMSSW_12_4_11 && cd CMSSW_12_4_11/src && cmsenv && cd -
git clone -b powheg_ggxy_ggzh https://gitlab.cern.ch/lichengz/genproductions_scripts.git
cd genproductions_scripts/bin/Powheg
cp -r <this repository>/cms/card/ggxy_ggZH_ZToLL_13p6TeV .
python3 ./run_pwg_parallel_condor.py -p 0123 -i ggxy_ggZH_ZToLL_13p6TeV/powheg.input -m ggxy_ggZH \
        -f ggxy_ggZH_ZToLL_13p6TeV -q longlunch -j 16 -x 3
python3 ./run_pwg_condor.py -p 9 -i ggxy_ggZH_ZToLL_13p6TeV/powheg.input -m ggxy_ggZH \
        -f ggxy_ggZH_ZToLL_13p6TeV -k 1
```

**One multicore node** (how the reference gridpack was integrated: 16 seeds, 3 x-grid
iterations, 4 h 16 min on 16 cores): `run/gridpack_single_node.sh <genproductions checkout>
card/ggxy_ggZH_ZToLL_13p6TeV <name> 16 3` in a 16-core batch job (`condor/gridpack.sub`).
Never on a login node.

The grid integral is in the stage-3 statistics of the work directory
(`pwg-st3-0001-stat.dat`: POWHEG's combination of all seeds). The gridpack's `pwg-stat.dat` is
the **first seed's** result only (the tarball step copies `pwg-0001-stat.dat`); reference:
13.031(7) fb combined, 13.038(27) fb from seed 1.

## 2. Using the gridpack

- **Process.** gg → Z(→ℓ⁺ℓ⁻)H at NLO QCD with the full top-mass dependence (two-loop virtual from
  ggxy, one-loop Born and reals from Recola); ℓ = e, μ, τ (`vdecaymode 11`), 15 < m_ℓℓ < 150 GeV,
  on-shell undecayed Higgs boson (decayed by Pythia8), μ_R = μ_F = m_ZH/2, on-shell m_t = 172.5 GeV,
  `hdamp 250`, NNPDF3.1 NNLO. Cross section 13.03 fb (scale band +16.7/−14.0 %, PDF 0.95 %,
  α_s ± 0.002: 2.8 %); multiply by the Higgs branching ratio forced in Pythia8.
- **Limits.** The two-loop amplitudes are expanded in the external masses: m_ℓℓ and m_H must stay
  below 150 GeV. Z → νν is available in the process (`vdecaymode 12`) but has not been validated
  here; hadronic Z decays are not available.
- **Weights.** The CMS list has 884 weights (9 scale, 570 Hessian, 305 replica members). With the
  amplitude cache an event costs about 3.3 s with the full list (79 s without the cache, 2–4 s
  without any list). The cache is on by default (`GGXY_AMPCACHE=0` turns it off); it is exact for
  PDF, μ_R, μ_F and m_t-scheme reweighting, and must not be used when reweighting other
  parameters. The CMS POWHEG source allows up to 2000 weights (stock POWHEG-BOX-V2: 500).
- **Negative weights:** 0.12 % of the events.
- **Pythia8.** CP5, POWHEG emission veto with `POWHEG:nFinal = 3` (the validated setting),
  `vetoCount 100` (CMS default). The LO companion sample is showered without the veto.
- **HTXS.** Set `rivetProducerHTXS.ProductionMode = 'GG2ZH'`: the default `AUTO` reads the LHE
  process and does not recognise `ggxy_ggZH`, and every event ends up in category 0.
- **LO companion.** `card/ggxy_ggZH_ZToLL_LO_13p6TeV` (`bornonly 1`, `LOevents 1`) for a LO
  gridpack. The reference LO sample was made differently, from the NLO gridpack with the grids
  removed and the Born integrated in every job (`run/gridpack_to_nanogen.sh … lo`); both give
  Born-level events with the LO normalisation (7.04 fb).

## 3. Validation chain

```sh
# one job (inside cmssw-el8): 1000 NLO events with the nine scale weights
run/gridpack_to_nanogen.sh <gridpack>.tgz 1000 <seed> scale nlo out
# summary of all jobs, and comparison with the reference
python3 validation/nanogen_summary.py --files 'out/nanogen_*.root' --label nlo --out summary-nlo.json
python3 validation/check_cms.py --nlo summary-nlo.json [--lo summary-lo.json]
```

Steps of a job: `runcmsgrid.sh` (weight list cut to the nine scale points), then
`cmsDriver.py … --step GEN,NANOGEN --era Run3 --conditions auto:mc --beamspot
Realistic25ns13p6TeVEarly2022Collision` with the hadronizer fragment and HTXS `GG2ZH`.
`nanogen_summary.py` needs uproot and awkward (about 1 s per file); for many files run it in
chunks with `--raw` and combine with `--merge`.

Reference (`reference/cms.json`): 200 000 NLO+PS events (200 × 1000) and 200 000 LO+PS events
(100 × 2000).

| | NLO+PS | LO+PS |
|---|---|---|
| σ from the event weights | 13.060(30) fb, +16.7/−14.0 % | 7.042(16) fb, +26.8/−20.0 % |
| fixed-order value (`../generator`) | 13.046(4) fb | 7.032(1) fb |
| negative weights | 0.12 % | 0 |
| \|y_H\| > 2.5 (HTXS FWDH) | 3.4 % | |
| HTXS category 0 | 0 | 0 |

`check_cms.py` compares the cross section at the nine scale points (pull < 3), the
negative-weight and FWDH fractions, requires no event in HTXS category 0, and compares the 18
STXS Stage 1.3 cells and seven distributions (p_T(Z), p_T(H), m_ZH, p_T(ZH), y_ZH, N_jets,
leading-jet p_T) with χ² tests (p > 0.001).

**Cost** (reference): gridpack 4.3 h on 16 cores (~70 core-hours); NanoGEN 200 × 1000 NLO events
with the nine scale weights ~100 core-hours, the LO sample ~30.

## 4. Fragments

- `ggZH_ZToLL_HToBB_ggxy-NLO_pythia8_CP5_cff.py`: ExternalLHEProducer (`GRIDPACK` placeholder) +
  Pythia8 CP5 with the POWHEG veto and PS weights, H → bb. Change `25:onIfMatch` for other decays.
- `ggZH_ZToLL_HToBB_ggxy-NLO_pythia8_CP5_hadronizer_cff.py`: the same Pythia8 settings without
  the LHE producer and without PS weights (the validation configuration, LHE input).
- `ggZH_ZToLL_HToBB_ggxy-LO_pythia8_CP5_hadronizer_cff.py`: LO companion, no POWHEG veto.
