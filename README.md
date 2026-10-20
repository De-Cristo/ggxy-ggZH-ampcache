# ggxy-ggZH-ampcache

An exact amplitude cache for the POWHEG-BOX-V2 process
[`ggxy_ggZH`](https://gitlab.com/POWHEG-BOX/V2/User-Processes/ggxy_ggZH)
(gg → ZH at NLO QCD with full top-mass dependence, matched to parton showers;
J. Davies, K. Schönwald, M. Steinhauser, D. Stremmer, [arXiv:2603.15762](https://arxiv.org/abs/2603.15762)),
and a reproduction package for gg → ZH with Z → ℓ⁻ℓ⁺ or Z → νν̄ at NLO and NLO+PS
at √s = 13.6 and 13 TeV. Both channels include all three lepton flavours.
The `vv13600` and `vv13000` sets in [`production/`](production/README.md) reproduce the neutrino
channel; its independently trained CMS-style map and 16 ONNX models are included in
[`models/vv13600/`](models/vv13600/), with channel-specific metadata and original closure outcomes.
The `atlas13600` set adds a public ATLAS-style charged-lepton source (PDF4LHC21, H_T scale,
A14 and the ATLAS shower veto), its A14 NLO+PS target, and separate correctionlib/ONNX exports.
The validated export files are included in [`models/atlas13600/`](models/atlas13600/):
the correctionlib map, five-/six-input ONNX classifiers, nominal and seven scale variations,
model card and checksums. The emulation scope and failed physics closure criteria are documented
in [`rwt/README.md`](rwt/README.md); successful export checks do not imply physics closure.
CMS-style reweighting covers both ℓℓ and νν̄ in the paper and package. ATLAS-style reweighting
covers ℓℓ only; no ATLAS-style νν̄ samples or corrections have been produced.

**The cache is part of `ggxy_ggZH` since 2026-10-07** (merge request !1, commit `1a6dbc6`; the
process README documents it). Use `ggxy_ggZH` at `0fbaccf` or later; the patch in `patch/` is
only needed for the original release `a22cb4f`.

| Directory | Contents |
|---|---|
| `src/`, `patch/`, `tools/`, `validation/` | the cache: source, patch for `a22cb4f`, patch generator, closure test of the cache |
| [`generator/`](generator/README.md) | generator-level reproduction without CMS software: build of ggxy, POWHEG-BOX-V2 and `ggxy_ggZH` at pinned commits; fixed-order LO/NLO with the ggxy program, POWHEG LO, NLO grids and events, cache closure, optional MCFM LO cross-check; reference results and a check script |
| [`production/`](production/README.md) | the production of the gg → ZH NLO+PS paper: `paper.sh` defines every run (jobs, statistics, seeds, settings) at 13.6 and 13 TeV and the verification; `run_local.sh` executes a set |
| [`shower/`](shower/README.md) | Pythia showers, Rivet distributions and STXS analysis for both decay channels |
| [`rwt/`](rwt/README.md) | LO+PS reweighting, closure checks, and reusable ONNX/correctionlib models for each channel |
| [`cms/`](cms/README.md) | CMS usage: gridpack card and recipe, Pythia8 fragments, LHE → NanoGEN validation chain with HTXS, reference results and a check script |

## Why the cache

POWHEG computes every alternative event weight (`rwl_compute_new_weight` → `setrandom` →
`gen_btilderw`) by recomputing B̃ at exactly the same kinematics; only the PDF member and the
scale factors change. In `ggxy_ggZH` the Born and the real emission are one-loop amplitudes
from Recola, and the virtual is the ggxy two-loop amplitude. None of them depends on the PDF or on
μ_F, and each is an exact power of α_s times a function of the momenta. With the 884-weight
list used by CMS, the same amplitudes were evaluated 884 times per event: 69–125 s per event
instead of 2–4 s.

## What the cache does

- A process-wide hash map (`ampcache.cpp`) from (flavours, raw bit pattern of the momenta, top
  mass, and μ_R for the two-loop virtual) to the squared amplitude with the α_s power divided
  out; on a hit, the Fortran caller multiplies back the current α_s.
- The Born is stored together with its spin-correlated tensor B_μν (`Born.f`), because Recola's
  `compute_bmunu_rec` reads the amplitudes of its last evaluation; the one-loop routine caches
  only the six-leg real amplitudes (`pwhg_ggxy.f`).
- Every key contains the top mass that is actually passed to the amplitude: the on-shell mass
  for `mtscheme 0`, the running mass m_t(μ_t) for `mtscheme 1` (MS-bar). The Born key originally
  used the on-shell input mass; the ggxy_ggZH authors pointed this out and supplied the change
  to `Born.f`. With the parameters POWHEG-BOX-V2 can reweight (μ_R, μ_F, PDF) both versions give
  the same results, because m_t(μ_t) is then fixed by the momenta, but the key now follows the
  same rule as the other amplitudes.
- Controls (environment variables): `GGXY_AMPCACHE=0` disables the cache (stock behaviour);
  `GGXY_AMPCACHE_MAX` caps the number of entries (default 400000, about 100 MB). A hit/miss
  summary is printed at exit. The cache is exact for reweighting PDFs, μ_R, μ_F and the m_t
  scheme; it must not be used when reweighting other input parameters.

## Patch for the original release

The patch is made against `ggxy_ggZH` commit `a22cb4f` ("add arXiv number"):

```sh
cd POWHEG-BOX-V2
git -C ggxy_ggZH checkout a22cb4f
patch -p1 < /path/to/patch/ggxy_ggZH-ampcache-a22cb4f.patch
cd ggxy_ggZH && make pwhg_main
```

`tools/make_patched_source.py <POWHEG-BOX-V2>/ggxy_ggZH <output dir> new.patch` regenerates it;
every edit is checked to match exactly once, so a changed upstream source fails loudly.

## Validation of the cache

`validation/closure_summary.py` compares, job by job, events generated with the same seed
with the cache (R3) and without it (R4), both with the full 884-weight list, in the on-shell
top-mass scheme (`mtscheme 0`).
`validation/results/closure.json` holds the result of 9 jobs on 9 worker nodes at four sites:

| | |
|---|---|
| events | 450, all records identical |
| weights | 397 800 / 397 800 identical as written in the LHE file (5 significant digits) |
| time per event | 3.3 s cached vs 79.0 s uncached (23.9× in aggregate) |
| cache hit rate | 97 % |

The only difference is POWHEG's internal `#rwgt` comment, at 10⁻¹⁵ relative, from dividing
out and multiplying back α_s; it matters only for a later `rwl_add` re-reweighting of the LHE.
`generator/run/cache_ab.sh` and `generator/analysis/compare_lhe.py` repeat this test with a
fresh build.

## Cross sections

gg → ℓ⁻ℓ⁺H (ℓ = e, μ, τ) at √s = 13.6 TeV, 15 < m_ℓℓ < 150 GeV, μ₀ = m_ZH/2, on-shell
m_t = 172.5 GeV, seven-point scale bands (NNPDF3.1 here, PDF4LHC21 in the reference):

| | this setup | arXiv:2603.15762, Table 3 |
|---|---|---|
| LO, ggxy fixed-order program | 7.0324(10) fb, +26.8/−20.0 % | 7.026(4) fb, +26.2/−19.7 % |
| LO, POWHEG fixed-order integration (8 seeds) | 7.0306(19) fb, +26.8/−20.0 % | |
| NLO, ggxy fixed-order program | 13.046(4) fb, +16.6/−14.0 % | 13.02(1) fb, +16.4/−13.8 % |
| NLO, POWHEG grids (16 seeds) | 13.031(7) fb | |
| NLO, weights of 200 000 cached events | 13.060(30) fb, +16.7/−14.0 % | |

With PDF4LHC21 the ggxy program reproduces the reference's stable-Z LO cross section, 70.601(7)
vs 70.59(2) fb. Summed over 200 000 LO events, the event weights reproduce the fixed-order LO
integration to 0.1 % at every scale point. An earlier version of this README quoted the POWHEG
LO integration as 7.037(5) fb: that was the first of the eight seeds alone.
Details, inputs and scripts: [`generator/`](generator/README.md).

## Licence and credits

The cache code, the patch and the scripts are released under the GNU GPL v3 (see `LICENSE`),
the licence of the ggxy library. The patch modifies files of the `ggxy_ggZH` process, whose
authors and licence terms apply to those files; this repository does not redistribute the
process, POWHEG-BOX, ggxy, Recola or MCFM. `generator/fo/` contains programs derived from the
ggxy examples (`examples/ggzh-nlo`). Please cite arXiv:2603.15762 and the ggxy paper (Comput.
Phys. Commun. 320 (2026) 109933) when using the process, and this repository (see `CITATION.cff`)
when using the cache or the reproduction package.
