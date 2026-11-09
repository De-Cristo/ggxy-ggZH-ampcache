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
