# ggxy-ggZH-ampcache

An exact amplitude cache for the POWHEG-BOX-V2 process
[`ggxy_ggZH`](https://gitlab.com/POWHEG-BOX/V2/User-Processes/ggxy_ggZH)
(gg → ZH at NLO QCD with full top-mass dependence, matched to parton showers;
J. Davies, K. Schönwald, M. Steinhauser, D. Stremmer, [arXiv:2603.15762](https://arxiv.org/abs/2603.15762)).

With a list of event weights (scale and PDF variations), the cache makes event generation
**about 24 times faster** and gives **bit-identical events and weights**.

## Why

POWHEG computes every alternative event weight (`rwl_compute_new_weight` → `setrandom` →
`gen_btilderw`) by recomputing B̃ at exactly the same kinematics; only the PDF member and the
scale factors change. In `ggxy_ggZH` the Born and the real emission are one-loop amplitudes
from Recola, and the virtual is the ggxy two-loop amplitude. None of them depends on the PDF or on
μ_F, and each is an exact power of α_s times a function of the momenta. With the 884-weight
list used by CMS, the same amplitudes were evaluated 884 times per event: 69–125 s per event
instead of 2–4 s.

## What the patch does

- A process-wide hash map (`ampcache.cpp`) from (flavours, raw bit pattern of the momenta, top
  mass, and μ_R for the two-loop virtual) to the squared amplitude with the α_s power divided
  out; on a hit, the Fortran caller multiplies back the current α_s.
- The Born is stored together with its spin-correlated tensor B_μν (`Born.f`), because Recola's
  `compute_bmunu_rec` reads the amplitudes of its last evaluation; the one-loop routine caches
  only the six-leg real amplitudes (`pwhg_ggxy.f`).
- Controls (environment variables): `GGXY_AMPCACHE=0` disables the cache (stock behaviour);
  `GGXY_AMPCACHE_MAX` caps the number of entries (default 400000, about 100 MB). A hit/miss
  summary is printed at exit.

## Apply

The patch is made against `ggxy_ggZH` commit `a22cb4f` ("add arXiv number"):

```sh
cd POWHEG-BOX-V2
git -C ggxy_ggZH checkout a22cb4f
patch -p1 < /path/to/patch/ggxy_ggZH-ampcache-a22cb4f.patch
cd ggxy_ggZH && make pwhg_main
```

For a newer upstream version, regenerate the patch; every edit is checked to match exactly once,
so a changed upstream source fails loudly instead of producing a wrong patch:

```sh
python3 tools/make_patched_source.py <POWHEG-BOX-V2>/ggxy_ggZH <output dir> new.patch
```

**CMS gridpacks.** Apply the patch to the process directory that `genproductions/bin/Powheg`
prepares (that tree also carries the CMS changes to POWHEG: `rwl_maxweights` 2000 and
`pdfweights_new.patch`), then build as usual. An existing gridpack can also be updated by
replacing only its `pwhg_main`; grids and cards are unchanged.

## Validation

`validation/closure_summary.py` compares, job by job, events generated with the same seed
with the cache (R3) and without it (R4), both with the full 884-weight list.
`validation/results/closure.json` holds the result of 9 jobs on 9 worker nodes at four sites:

| | |
|---|---|
| events | 450, all records identical |
| weights | 397 800 / 397 800 identical (string comparison) |
| time per event | 3.3 s cached vs 79.0 s uncached (23.9× in aggregate) |
| cache hit rate | 97 % |

The only difference is POWHEG's internal `#rwgt` comment, at 10⁻¹⁵ relative, from dividing
out and multiplying back α_s; it matters only for a later `rwl_add` re-reweighting of the LHE.

Physics validation (cross sections, scale bands, STXS Stage 1.3 predictions) is documented in
the accompanying CMS analysis note. At √s = 13.6 TeV, 15 < m_ℓℓ < 150 GeV, μ₀ = m_ZH/2, the
cached NLO events give 13.060 ± 0.030 fb with a seven-point band of +16.7/−14.0 %, against
13.02 fb and +16.4/−13.8 % in arXiv:2603.15762 (Table 3).

## Licence and credits

The cache code and the patch are released under the GNU GPL v3 (see `LICENSE`), the licence of
the ggxy library. The patch modifies files of the `ggxy_ggZH` process, whose authors and licence
terms apply to those files; this repository does not redistribute the process, POWHEG-BOX,
ggxy or Recola. Please cite arXiv:2603.15762 and the ggxy paper (Comput. Phys. Commun. 320
(2026) 109933) when using the process, and this repository (see `CITATION.cff`) when using the
cache.
