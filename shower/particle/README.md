# ggZH particle inputs v1

Defined by `production/paper.sh` sets `particle-zll` / `particle-zvv`, runs
`particle-build`, `particle-nlo`, `particle-lo-cms`. Run only on Condor workers.

The new driver uses exactly Pythia's public LHAupLHEF reader through frameType=5,
with a derived setEvent() that counts successfully consumed original LHE entries.
It records accepted-event and original-LHE ordinals separately. No random draw or
physics mutation is made by the writer. Production reference comparisons verify
that changing the reader route preserves the existing replay. The production
Rivet plugin is reused byte for byte from shower-install-rwt10m-0.tgz.

Records preserve native Pythia statuses, mother/daughter lists and all four-vectors
at float64 precision. Raw records include H->bb decay particles, Z descendants,
beam remnants and MPI. Radiation objects exclude all Z/H descendants; the parton
view selects terminal QCD copies with abs(native status)71..79, pT>=10, MPI included,
status63 remnants excluded. The pre-terminal selected index list is also saved.
This is an explicit public selection, not a reconstruction of opaque CMS statusFlags.

`schema.json` specifies input tensor axes, units, transformations and masks.
Caps64/128 are storage choices for the first real shards, subject to occupancy
review before freezing training. The complete raw records have no such cap.
No model is trained in this production. The forward check only loads real event
inputs into the MIT public ParT module with physical pair embeddings enabled.

Training-ready datasets live in `/scratch/lzhang/ggzh-reweighting/particle-v1/`.
The private transport stages them through new EOS particle-v1 namespaces and
publishes a per-shard ready manifest only after bounded checksum-verified copying.
Do not use a file merely because it exists without its ready manifest.

Example after reading a shard's manifest (on a worker):

```python
import numpy as np
from pathlib import Path
root = Path('/scratch/lzhang/ggzh-reweighting/particle-v1')
with np.load(root/'zll/lo-cms/inputs/jets-v1/110001.npz', allow_pickle=False) as d:
    x = {k: d[k] for k in ('pf_features','pf_vectors','pf_mask',
                           'event_features','ps_encoding_features')}
    ids = d['lhe_ordinal']
    weights = d['weight']  # never a classifier feature
    fold = d['source_fold']
```

Stage EOS shards to the worker when scratch is not mounted; no worker mount is
assumed. Exact sample inclusion/exclusions, package commits, cluster/index maps
and cost live in the planning repo ledger and campaign manifests. Physical
signed weights are preserved; future training loss and feature normalization
must be prospectively fixed separately.

ATLAS-style charged-lepton inputs use set `particle-atlas-zll` and the prospective
contract `production/atlas-parton10m.json`. The independent LO source uses
PDF4LHC21 and the validated ATLAS dynamic scale; both source and target use A14
with `atlas-common`. The target reuses the rwt10m NLO hard events and selects
PDF4LHC21 weight1001, while the source uses its nominal weight. Shower seeds are
independent of the CMS campaign. Run definitions, input/output names and resources
are in `production/paper.sh`; the `@atlas/` product prefix selects the separately
stored validated ATLAS prerequisites. A staging provider must resolve that prefix.

`run_particle_atlas.sh` compares the full real shard with an original same-seed,
same-card public shower replay before publication. `convert_atlas.py` retains the
same native records, feature functions,128-parton capacity and overflow counters;
the full raw record has no tensor cap. ATLAS splits are assigned by original row,
without a CMS split map. Inspect actual multiplicities and overflow before freezing
a training capacity. Production of these inputs does not establish closure or a
trained model. Existing CMS scripts, definitions and frozen split maps are unchanged.
