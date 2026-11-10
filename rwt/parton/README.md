# Real particle-input ZH ParT study

Set `parton-zll-1m-v1` in `production/paper.sh` defines the complete exploratory
1M LO +1M NLO chain: real CPU preparation, fifteen CUDA source-fold/target-half
fits, matched C BDT, and CPU closure/report merge. This is particle-input ParT,
not the scalar-token adaptation in `rwt/transformer/`. Its public eight-layer,
eight-head backbone uses particle features, four-vectors, padding masks and
physical pair attention, with event/shower context after class attention.
No private extractor or pretrained model is used.

Selection is frozen in `production/parton-zll-1m-v1.json`; configuration contains
its SHA. Preparation verifies every selected artifact and signed weight/ID/split,
retains every row and removes only verified empty padding. Bulk input locations
are supplied by `PAPER_EOS_ROOT` on a compute worker. Existing production pins,
raw records, frozen full-sample results and the original A/B/C criteria are untouched.

Checkpoint reload and event-aligned OOF prediction are mandatory. Negative target
weights are zeroed only in the classification loss; closure is signed. Normalization
is fitted only on the allowed training partitions. Use saved source folds and
opposite target halves. The matched BDT keeps its original classifier configuration
but reads those same saved folds, rather than recomputing ranks after sampling.

Scalar/STXS/half closure and all failures are reported, with effective counts,
region and hard-jet diagnostics, low-mass ratio and mass-by-pT(Z) check. These are
an exploratory subset result, not full-inventory acceptance or proof of the full
joint particle density. The report includes JSON, HTML and individual PNG/PDF
figures through the paper's exact `tools/prd_style.py` source. A trained checkpoint
can be reloaded through `model.py`; saved logits/weights preserve event identity.
No success is claimed before real worker fits and held-out closure finish.

### Analysis-jet context ablation

`paper.sh` sets `parton-zll-jets-1m-v1` and `parton-zvv-jets-1m-v1` reuse the
corresponding frozen prepared1M/class bundle. Each defines five nominal GPU fits
and one CPU `parton-jet-compare-1m` job. No new generation or shower is needed.
The comparison consumes the five original nominal fits and five variant fits;
it verifies package provenance, selection, reload status and exact OOF identities.
Original2L fits may have mixed accepted legacy/recovery provenance, recorded in
full. Every output uses a new variant namespace.

`jet_context.py` appends11 entries to the original24: log1p(analysis jet count),
and presence/log1p(pT/GeV)/eta/sin(delta phi to Z)/cos(delta phi to Z) for each of
j1,j2. Absent jets are zero before training-only normalization. These are the
existing scalar-table pT>30GeV, |eta|<4.7 jets, not reclustered partons. No jet mass
is available. For application, use `extra_context` on aligned table columns,
concatenate to the original24 context and load normalization from the checkpoint.
The model's saved configuration specifies35 globals. Particle inputs, pair
attention, folds, seeds and all other hyperparameters remain unchanged; the wider
head has additional parameters. Only nominal fits are registered for this ablation.

The comparison publishes unchanged scalar/STXS tests plus conditional j1/j2
pT/eta shapes in fixed VpT bins0/75/150/250/400/600/infinity and analysis jet
counts0/1/2/>=3. Category fractions use full-window normalization; conditional
shapes use unit area within each category and plotted variable range. Sparse or
nonpositive cells are flagged. Signed target weights, scale bands and legacy
fixed-normalization MC errors are retained. Extra conditional p-values are
descriptive, not new acceptance criteria. The NLO reference includes training
and validation events; no independent-target or target-half-stability success is
claimed. Reports include N_eff, region/hard-jet/low-mass diagnostics, JSON, a
self-contained HTML and PNG/PDF plots using `tools/prd_style.py`. A baseline/variant
comparison is correlated because both use the same events; it is not two
independent experiments. Bootstrap/PCA uncertainty production is a separate task.

### External GPU workers

The user-authorized external continuation uses the same public trainer on an
allocated compute worker. `tools/run_parton_external.py` snapshots a clean
committed checkout (or immutable exported source), checks prepared-data checksums
and audit, selects an explicitly assigned CUDA GPU, and runs a complete real
fold in a new work directory. Dependencies must already be installed. Declare
`RWT_EXTERNAL_COMPUTE_NODE=1`; pass channel/context/task/data/work/run-id and set
`CUDA_VISIBLE_DEVICES` to one GPU. Original Condor behavior is the default; cmscon
and lxplus login hosts remain forbidden. No Condor environment is impersonated.
`RWT_EXECUTION_CONTEXT=external-gpu` and the unique external campaign ID are saved
in checkpoints and fit records. Local iterations require new commits/campaigns;
publish only in a separate external namespace. Registered grid pins/results are
not replaced automatically. Worker runs retain failure artifacts, exact code
hashes, dependencies, GPU details, training logs and CPU/allocated GPU wall time.
The receiving agent must check `external-run.json` AND fit/reload products, rather
than treating exit status as sufficient. Actual real fits validate numerics; no
synthetic or disposable-pilot test is part of this adapter.
