# Shower and analysis: Pythia 8 + Rivet (no experiment software)

Showers the LHE events of `../generator` with Pythia 8 and analyses them in the same process with
Rivet: the LHC Higgs WG STXS classification (`SimplifiedTemplateCrossSections`, Stage 0, 1.2 and
1.3, coarse and fine) and `GGZH_NLOPS` (p_T(Z), p_T(H), m_ZH, p_T(ZH), y_ZH, jets, and our own
18-cell Stage 1.3 fine classification as a cross-check). The LHE weights (scale, PDF, α_s) and
Pythia's shower-variation weights are carried to every histogram.

**Status:** built and run on NLO events (smoke test 2026-10-09: 14/14 checks of `check_shower.py`); LO+PS modes and the full samples in progress (paper plan WP2).

| | |
|---|---|
| `versions.sh` | LCG view (default LCG_108a gcc 13: Pythia 8.315, Rivet 4.1.2, HepMC3, FastJet), STXS routine commit |
| `build.sh <prefix>` | Rivet plugin `RivetGGZH.so`, driver `ggzh-shower`, relocatable `env.sh` |
| `src/ggzh_shower.cc` | Pythia 8 driver: LHE in, HepMC3 event handed to `Rivet::AnalysisHandler` in-process, POWHEG hooks when `POWHEG:veto > 0` |
| `rivet/GGZH_NLOPS.cc` | distributions; normalisation Σw/N per weight (pb per bin) |
| `cards/` | `common`, `tune-monash`, `nlo` (POWHEG veto), `lo-*` (shower-start variants for LO+PS), `ps-variations`, `noshower` (test) |
| `run/shower.sh` | one LHE file → `<name>.<mode>.<tune>.yoda` |
| `analysis/yoda_summary.py` | YODA (text) → JSON from Rivet's RAW sums: all histograms and weight variations with under/overflow, gg → ZH(ℓℓ) STXS bins; per job (`--raw`) and combined (`--merge`); plain Python |

```sh
./build.sh $PWD/shower-install
export GGZH_SHOWER_PREFIX=$PWD/shower-install
run/shower.sh events-1001.lhe.gz out nlo monash 1001                   # NLO+PS
run/shower.sh lo-events-2001.lhe.gz out lo-default monash 2001         # LO+PS, default shower start
python3 analysis/yoda_summary.py out/*.nlo.monash.yoda --out nlo.json                # or per job with --raw, then --merge
```

## Modes

- `nlo`: POWHEG matching for NLO events: `POWHEG:veto 1`, `nFinal 3` (ℓ⁺, ℓ⁻, H), `pTdef 1`,
  `pThard 0`, `pTemt 0`, `emitted 0`, `vetoCount 100`, `pTmaxMatch 2` (ISR and FSR).
- LO events (Born only, `bornonly 1`, `LOevents 1` in the generator) without veto, with the
  shower start as the variable under study:
  `lo-default` (`pTmaxMatch 0`: for a final state without coloured partons Pythia showers up to
  the kinematic limit), `lo-scalup` (`pTmaxMatch 1`: up to SCALUP of the event),
  `lo-scalup-half`, `lo-scalup-quarter` (`pTmaxFudge` 0.5, 0.25), `lo-power` (`pTmaxMatch 2`).
- `noshower`: no shower, MPI or hadronisation; for the tests.

## Tests (to pass before results are used)

1. `noshower` on NLO events: the GGZH_NLOPS `sigma` histogram equals `generator/analysis/lhe_xsec.py`
   of the same file for the nominal and every scale weight (weights and normalisation reach Rivet).
2. No event in the STXS routine's UNKNOWN category; the STXS fine bins agree with GGZH_NLOPS's own
   18-cell classification.
3. Same seed, same output; the POWHEG veto acts (vetoed emissions are counted in the log).
4. LO+PS with the shower off reproduces the fixed-order LO distributions of `generator/` (ggxy).

Tunes other than Monash (CP5, A14) are added as cards once their parameters are taken from the
publications (paper plan D2).
