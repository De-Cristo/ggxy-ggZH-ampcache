# Data of "gg → ZH at NLO matched to parton showers" (event samples, grids, analysis outputs, reweighting models)

These files were made with the workflows distributed in the reproduction package `ggxy-ggZH-ampcache`
(branch `reproduction`; the release is named in the record's metadata), and every number of the paper comes from
the analysis outputs below (`production/analyze.sh`).  The package can regenerate all of it; these files save the
computing and let the results be checked without rerunning.  Licence: CC BY 4.0.

Processes: gg → Z(→ℓ⁻ℓ⁺)H and gg → Z(→νν̄)H at NLO QCD with the full top-quark mass dependence (ggxy_ggZH in POWHEG-BOX-V2, with the
amplitude cache), Z decay and 15 < m_ℓℓ (or m_νν) < 150 GeV as in `generator/cards/ggzh_13p6TeV.input`, m_H = 125 GeV,
m_t = 172.5 GeV (on-shell), scale m_ZH/2 unless stated.

## Archives

| Archive | Content | Run (`production/paper.sh`) |
|---|---|---|
| `lhe-nlo-13600.tar`, `lhe-nlo-13000.tar` | NLO events, 200 files × 2 500, weighted (Born suppression, L = 250 GeV), PDF4LHC21, √s = 13.6 / 13 TeV | `events-nlo` |
| `lhe-lo-13600.tar`, `lhe-lo-13000.tar` | LO events (Born only), 100 files × 5 000, unweighted, PDF4LHC21 | `events-lo` |
| `lhe-nlo-hdampinf-13600.tar` | NLO events with h_damp = ∞ (no damping), 200 × 2 500, 13.6 TeV | `events-nlo-hdampinf` |
| `lhe-nlo-nnpdf31-13600.tar`, `lhe-lo-nnpdf31-13600.tar` | validation set: NLO unweighted (50 × 1 000), LO (25 × 2 000), NNPDF3.1, 13.6 TeV | set `verify`: `events-nlo`, `events-lo` |
| `lhe-nlo-bs250-nnpdf31-13600.tar` | validation set: NLO with Born suppression (50 × 1 000), NNPDF3.1 | `events-nlo-bs250` |
| `lhe-lo-cms-13600.tar` | CMS-style LO events of the reweighting study: NNPDF3.1, fixed scale m_Z + m_H, 100 × 5 000, no weight list | `events-lo-cms` |
| `grids.tar` | the POWHEG grids of every sample (`grids/<sample>.tgz`): more events of a sample need only `generator/run/powheg_events.sh` with these grids, not the integration; and the LO cross section of the CMS-style grids; Z → νν̄ grids with `-vv-` in the name | `grids-*` |
| `tables-lo-cms-13600.tar`, `tables-lo-cms-default-13600.tar`, `tables-nlo-cp5-13600.tar` | reweighting study: one row per showered event (Rivet `GGZH_EVENTTABLE`, columns in `shower/rivet/GGZH_EVENTTABLE.cc`): CMS-style LO+PS, the same events with Pythia's default start, NLO+PS target (CP5) | `table-*` |
| `rwt-weights.tar` | per-event classifier weights of the paper (`weights-<name>.npz`: run, event, in_window, w_full, w_A_to_B, w_B_to_A; out of fold) and the classifier test outputs (`bdt-<name>.json`); `-vv-13600` in the name: Z → νν̄; `clf-atlas`, `clf-jet-atlas`: the ATLAS-style source | `rwt-clf*` |
| `rwt-model.tgz` | reusable models: classifier with five and six inputs as ONNX and the p_T(Z) × N_jets map as a correctionlib table, nominal and seven scale variations, and `model-card.json` (inputs, weight formula, normalisations, checks, validity); `model-vv-13600.tgz`: the same for Z → νν̄; `model-atlas-13600.tgz`: for the ATLAS-style source | `rwt-model`, `rwt-model-atlas` |
| `lhe-nlo-vv-13600.tar`, `lhe-nlo-vv-13000.tar`, `lhe-lo-vv-13600.tar`, `lhe-lo-vv-13000.tar`, `lhe-nlo-hdampinf-vv-13600.tar`, `lhe-lo-cms-vv-13600.tar` | Z → νν̄ (ν_e, ν_μ, ν_τ): the samples of the rows above with `vdecaymode 12` (the m window applies to m_νν) and seeds independent of the ℓℓ ones | sets `vv13600`, `vv13000` |
| `tables-lo-cms-vv-13600.tar`, `tables-lo-cms-default-vv-13600.tar`, `tables-nlo-cp5-vv-13600.tar` | the reweighting tables for Z → νν̄; their rows carry 18 columns more than this package's `GGZH_EVENTTABLE` writes (the four-vectors of the Z-decay fermion and antifermion, from the production build), the other columns as written here | `table-*` (set `vv13600`) |
| `lhe-lo-atlas-13600.tar` | ATLAS-style LO events of the reweighting study: PDF4LHC21, dynamic scale H_T = m_T(H) + p_T(ℓ₁) + p_T(ℓ₂) (`dynscale 3`), 100 × 5 000, seeds 70001–70100, the PDF4LHC21 weight list | `events-lo` (set `atlas13600`) |
| `tables-lo-atlas-13600.tar`, `tables-nlo-a14-13600.tar` | reweighting tables of the ATLAS-style source (A14 as the ATLAS job options set it, the ATLAS base settings and the ATLAS POWHEG veto) and of its target (the production NLO events, seeds 10001–10200, with the same A14 settings; PDF4LHC21 nominal) | `table-lo-atlas`, `table-nlo-a14` |
| `analysis-13600.tar`, `analysis-13000.tar`, `analysis-verify.tar`, `analysis-vv13600.tar`, `analysis-vv13000.tar`, `analysis-atlas13600.tar` | outputs of `production/analyze.sh <set>`: fixed-order, event and shower summaries (Rivet histograms with every weight variation, `shower-*.json`, `parton-*.json`), the comparisons and checks behind every number of the paper, `STATUS.txt` (package commit of every run), `SHA256SUMS` | analysis of each set |

`README.md` is this file; `FILES.sha256` has the checksum of every file inside the archives and `MANIFEST.sha256` that
of every archive.

## Production provenance

The neutrino samples were generated on the development branch before its production code was
incorporated into `reproduction`. `STATUS.txt` preserves the original commit of each run; the
analysis header records the original commit used to combine them. The numerical analysis and
reweighting sources match those in the integrated package byte for byte; the original result
headers and checksums are retained. The paper's `results/vv-provenance.json` records the source
checksums and explains the superseded 13 TeV grid attempt in the original status report. The curated neutrino
job definitions were checked against those used for production, while the charged-lepton job
definitions were unchanged. The generator installation adds `fo-ggvvh`; the production shower
build also recorded decay-fermion information. The additional columns are documented above,
and the shared event-table columns and shower summaries passed the build comparisons.

## Event weights

The production LHE files carry 155 weights per event (`generator/cards/pwg-rwl-pdf4lhc21.dat`): ids 1001–1009 the
nine (μ_R, μ_F) ∈ {½, 1, 2}² combinations, ids 3000–3042 the PDF4LHC21 members 93300–93342 (Hessian and α_s), and
ids 2000–2102 the NNPDF3.1 members 325300–325402.  The validation samples carry 112 (`pwg-rwl.dat`: scale and
NNPDF3.1).  Weighted
samples: the cross section is the sum of the weights divided by the number of events (pb, summed over the three charged-lepton flavours, or the three neutrino flavours for Z → νν̄).

## Using the reweighting models

See `rwt/README.md` of the package and `model-card.json`: inputs in the card's order as float32, p = `probabilities`[:, 1],
w = p/(1−p), normalised to mean 1 over the sample's events inside 15 < m_ℓℓ (or m_νν) < 150 GeV.  The weights reshape a sample
defined like the CMS-style source; they do not create events where it has none (the paper's hard-jet test).
For Z → νν̄ samples use the separately trained `model-vv-13600` (m_νν in place of m_ℓℓ).
The two decay channels have different NLO corrections near threshold, where γ* → ℓℓ contributes
only to the charged-lepton channel.
For samples defined like the ATLAS-style source (μ = H_T, PDF4LHC21, A14 and the ATLAS POWHEG veto) use
`model-atlas-13600`; its target is NLO+PS with the same A14 settings.  The same caveat applies: the vetoed LO+PS
sample lacks events with a jet harder than m_ZH.

## Citation

The paper, and the package release named in this record's metadata.
