# Reweighting LO+PS to NLO+PS (the paper's section "Reweighting LO+PS to NLO+PS")

Can generator-level event weights turn an LO+PS sample of gg → ZH into NLO+PS?  The source is an
LO+PS sample made with public CMS settings ("CMS-style"); the target is the production NLO+PS sample of
this package.  The weights come from a classifier (the paper's method) or, for comparison, from a binned
map in p_T(Z) × N_jets.  Settings and pass criteria were fixed before each test was run.

## Samples (production/paper.sh, set 13600)

The same runs exist in set `vv13600` for Z → νν̄ (`vdecaymode 12`; products named `…-vv-13600`, seeds independent
of the ℓℓ ones: CMS-style LO events 45001+, NLO events 30001+).

Set `atlas13600` repeats the study with an ATLAS-style source, the settings of the public ATLAS Run-3 gg → ZH job
options (LHE 602411, shower 603335) with public tools: `events-lo` with the card's `dynscale 3` (μ_R = μ_F = H_T =
m_T(H) + p_T(ℓ₁) + p_T(ℓ₂), patch `generator/patches/ggxy_ggZH-scale-ht-atlas.patch`) and PDF4LHC21, seeds 70001+;
`table-lo-atlas`: shower with A14 as athena sets it (`shower/cards/tune-a14-athena.cmnd`), the ATLAS base settings
(`atlas-common.cmnd`) and the ATLAS POWHEG veto at SCALUP = m_ZH (`atlas-lo-veto.cmnd`); `table-nlo-a14`: the
production NLO events with the same A14 settings (target, PDF4LHC21 nominal); `rwt-clf-atlas`, `rwt-clf-jet-atlas`
(five and six inputs), `rwt-summary-atlas` (`ggzh_rwt_summary_atlas.py`), `rwt-model-atlas` (`model-atlas-13600.tgz`);
`powheg-lo` with the same card checks the scale option (the event weights at the seven scale points against
separate integrations, `generator/check.py --same-card`).  Differences to the ATLAS samples: ggxy's amplitudes
(m_b = 0) for ggHZ's top and bottom loops, 15 < m_ℓℓ < 150 GeV for 10–1000, no EvtGen, κ or shower weights. The LO veto uses `POWHEG:nFinal=-1` to accommodate the ggxy
event record, instead of the official configuration’s final-state count of 2. The NLO target retains
the package’s validated NLO veto (`nFinal=3`, `vetoCount=100`); the LO source uses `vetoCount=3`.
This is a public emulation, not an official ATLAS sample or a validation against its gridpack.

The ATLAS closure study counts every LO event once, including events carrying upper-bound-excess
weights. The summary records their number. NLO closure distributions retain signed event weights;
negative target weights are set to zero only for classifier training. Exported classifiers are
trained on all events; the reported closure tests use the separate out-of-fold and split-sample
weights. Their pass/fail outcomes accompany the ATLAS model card, independently of the numerical
ONNX and correctionlib export checks.

| Run | Jobs | What |
|---|---|---|
| `grids-lo-cms` | 1 | LO grids, NNPDF3.1 (the card's PDF), fixed scale μ_R = μ_F = m_Z + m_H (`fixedscale=1`), 8 seeds |
| `events-lo-cms` | 100 | 5 000 unweighted LO events each, seeds 40001–40100, no weight list |
| `table-lo-cms` | 100 | the LO events showered with CP5, the CMS common settings (`shower/cards/cms-common.cmnd`) and the CMS POWHEG veto at SCALUP = m_ZH (`cms-lo-veto.cmnd`): the source |
| `table-lo-cms-default` | 100 | the same events with Pythia's default start (`lo-default.cmnd`, no veto) |
| `table-nlo-cp5` | 200 | the production NLO events (`events-nlo`, seeds 10001–10200) showered with CP5 and the CMS common settings: the target, with the NNPDF3.1 central member (weight id 2000) as nominal weight |

The `table-*` runs write one row per event (Rivet routine `shower/rivet/GGZH_EVENTTABLE.cc`, objects as in
GGZH_NLOPS), staged as products `tables/<sample>-13600/table-<seed>.csv.gz`.  The showers use the same seeds as
every other shower of the package (LHE seed + 100000).

## Reweighting and tests

| Run | What |
|---|---|
| `rwt-clf` | classifier weights of the CMS-style source, five inputs: m_ZH, p_T(ZH), p_T(Z), N_jets (0, 1, 2, ≥3), \|y_ZH\| |
| `rwt-clf-default` | the same for the default-start source |
| `rwt-clf-jet`, `rwt-clf-jet-default` | six inputs: the leading-jet p_T added |
| `rwt-summary` | closure of the binned map, and the summary: every number of the paper's section and the data of its figures (`rwt-summary.json`) |
| `rwt-model` | the reusable models of the data release: the classifier (five and six inputs) as ONNX and the binned map as a correctionlib table, each for the nominal target and its seven scale variations, trained on all events; `model-card.json` gives inputs, weight formula, normalisations, training samples and the checks (ONNX against scikit-learn, correctionlib against the map code) |

Classifier (`ggzh_reweight_bdt.py`): scikit-learn `HistGradientBoostingClassifier` (learning rate 0.05, at most
400 iterations, 31 leaves, at least 100 events per leaf, L2 1.0, early stopping on 20 % of the events, fixed seed),
balanced classes, NLO events with their weights (negative weights set to zero); an LO event with classifier output
p gets the weight p/(1−p), out of fold (five folds by event index), normalised to mean 1 inside 15 < m_ℓℓ < 150 GeV.

Binned map (`ggzh_reweight_map.py`): ratio of the normalised NLO+PS and LO+PS fractions in p_T(Z) (25 GeV bins to
400 GeV, then 400–600, > 600) × N_jets (0, 1, ≥2); above 250 GeV a smooth fit per jet bin.

Pass criteria (`ggzh_reweight_closure.py`, `ggzh_reweight_bdt.py`):
1. the eight STXS Stage 1.3 POI fractions agree with NLO+PS, every pull below 2 (MC statistics);
2. every distribution not used as an input (p_T(H), the leading-jet p_T, ΔR(Z,H), y_H, m_ℓℓ; for the map also
   p_T(ZH), m_ZH), and for the classifier also every input, has a χ² p-value above 5 % against NLO+PS (bins merged to
   at least 20 events in both samples), and no merged bin lies outside the NLO+PS scale band by more than two
   standard deviations of the MC statistics.  The band is the minimum and maximum over the nominal and seven
   variations: all points of the (μ_R, μ_F) ∈ {½, 1, 2}² grid except (2, ½) (`SCALE7` in `ggzh_reweight_common.py`);
3. criteria 1 and 2 hold with weights from one half of the NLO+PS sample (even or odd event index) tested against
   the other half, in both directions.

## Running

On the pool, like every run of the package (`production/paper.sh show 13600 <run> <p>`), in the order of the
tables above.  The steps of the `rwt-*` runs are `rwt/run.sh clf …` and `rwt/run.sh summary`; they use the LCG_106
view (scikit-learn 1.2.2, the version the published weights were made with), single-threaded.
`production/analyze.sh 13600` copies the summary and the closure outputs to `analysis/rwt/`.

Locally, from per-event tables already made:

    rwt/ggzh_table.py --label lo-cms --tables 'tables/lo-cms/table-*.csv.gz' --out lo-cms.npz
    rwt/ggzh_table.py --label nlo-cp5 --tables 'tables/nlo-cp5/table-*.csv.gz' --out nlo-cp5.npz --nominal 2000
    rwt/ggzh_reweight_bdt.py --source lo-cms.npz --target nlo-cp5.npz --out bdt-clf.json --weights-out weights-clf.npz

Applying a released model to another LO+PS sample: build the inputs as in `model-card.json` (float32, in its order),
take p = `probabilities`[:, 1] of the ONNX model, w = p/(1−p), and normalise w to mean 1 over the sample's events with
15 < m_ℓℓ < 150 GeV (weight 1 outside); the map: `correctionlib` evaluate(systematic, p_T(Z), N_jets), already
normalised for the source sample.  The weights reshape a sample like the source; they do not create events where it
has none (the paper's hard-jet test).

`ggzh_reweight_2d.py` (a map in m_ZH × p_T(ZH)) is not used in the paper; the classifier uses its helpers.


## Included ATLAS model files

[`../models/atlas13600/`](../models/atlas13600/) contains the actual trained delivery files:
`map-ptz-njets.json`, `clf5-*.onnx`, `clf6-*.onnx`, and `model-card.json`. Each classifier has
one nominal model and seven target-scale variations. `SHA256SUMS` covers all 18 files;
`provenance.json` links them to the package commits and the original data-release archive
`rwt/model-atlas-13600.tgz`. These are byte-identical to that archive's contents.

The export was checked against scikit-learn and the map implementation on the source events:
maximum classifier probability difference 3.25e-7, maximum map-weight difference 1.19e-14.
The closure study does **not** pass: both classifiers fail the STXS, differential (including
inputs), and split-sample criteria. The map passes the STXS criterion but fails the other two.
Use the files as documented shape corrections within their stated scope, with these limitations.

For ONNX inference use float32 inputs in the card's order and `probabilities[:, 1]` as p;
apply p/(1-p)/c with the normalization c for that model and variation in the card. For a new
sample, normalize to unit mean inside the mass window. For correctionlib, evaluate
`ggzh_lo_to_nlo_ptz_njets(systematic, ptz, njets)` with real-valued ptz and njets; the map is
already normalized for the production source. Apply weight 1 outside the mass window for
both methods. Neither file format automatically tests the event's dilepton mass.

The exported classifiers use all training events. The paper's closure numbers instead use
out-of-fold weights; `closure_tests` in the card preserves those tests, while the export's
normalizations and numerical comparisons describe the all-event models. Production job
archives omit the root `models/` directory; regenerating the models uses `rwt-model-atlas`.
