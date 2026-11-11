# User-approved NLO control gate revision, 2026-11-11

`paper.sh` set `nlo-topup-zll-10m-v2` selects
`NLO_TOPUP_CONTROL_GATE=user-20261111-rtol1e-3`. It changes only control RTOL
from1e-4 to1e-3. ATOL remains1e-8; all ten existing controls per event, exact
original event/weight preservation, replay identity, finite values, complete
2500-event count, checksums and0.90CPUh per-shard replay gate remain mandatory.
The user's stated reason is three comparisons out of25000 with maximum relative
difference1.42e-4 at stored-weight rounding. Original strict-gate failures remain
failures under that gate. Existing sets retain their original tolerance.

New products use rwt/nlo-topup-zll-10m-v2/. Already accepted v1 companion products
are reused by reference, never recalculated or overwritten. A provenance catalog
must record the actual package, tolerance and namespace for every accepted shard.
The native first shard101654 can reuse its existing accepted calculation exactly
as in v1. Seeds, generator, physics settings and scale definitions are unchanged.

Prioritize the NLO shards in production/parton-zll-1m-v1.json, then the separately
frozen evaluation selection supplied by the study owner. Do not invent evaluation
membership. This is weight-input production only, not a new model or fit.
