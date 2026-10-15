# Reweighting LO+PS to NLO+PS (the paper's section "Reweighting LO+PS to NLO+PS")

Can generator-level event weights turn an LO+PS sample of gg → ZH into NLO+PS?  The source is an
LO+PS sample made with public CMS settings ("CMS-style"); the target is the production NLO+PS sample of
this package.  The weights come from a classifier (the paper's method) or, for comparison, from a binned
map in p_T(Z) × N_jets.  Settings and pass criteria were fixed before each test was run.

## Samples (production/paper.sh, set 13600)

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

`ggzh_reweight_2d.py` (a map in m_ZH × p_T(ZH)) is not used in the paper; the classifier uses its helpers.
