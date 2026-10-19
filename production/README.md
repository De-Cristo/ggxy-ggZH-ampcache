# Production of the gg → ZH NLO+PS paper

`paper.sh` is the single definition of every run behind the paper: which runs exist, how many
jobs each has, the statistics, random seeds and settings of every job.  Nothing else (no
submission script, no command history) defines them.  A run is a set of identical jobs that
differ only in the job index p; seeds are fixed functions of p.

| Set | Contents |
|---|---|
| `common` | the generator build (`generator/build.sh`) and the shower build (`shower/build.sh`) |
| `13600`, `13000` | PDF4LHC21 production at √s = 13.6 and 13 TeV, identical runs and seeds: ggxy fixed order (ℓℓH and stable Z, LO and NLO), POWHEG LO at 7 scale points, LO grids, NLO grids with Born suppression (L = 250 GeV), MCFM, 5·10⁵ NLO events (weighted) and 5·10⁵ LO events with 155 weights, cache closure, showers: NLO+PS and five LO+PS shower starts for the tunes Monash, CP5 and A14; the NLO and LO events at parton level (no shower) |
| `verify` | NNPDF3.1 at 13.6 TeV: the verification runs behind `generator/reference` (48 checks), the unweighted 5·10⁴ NLO and LO validation samples and their showers, and the Born-suppression test (same seeds with `bornsuppfact 250`) |
| `vv13600`, `vv13000` | Z → νν̄ (`vdecaymode 12`, three neutrino flavours): the runs of `13600` / `13000` except the decay-independent or code checks (stable Z, MCFM, POWHEG LO, cache, PDF members) and the m_t scheme scan (one MS-bar check, μ_t = m_ZH, at 13.6 TeV; the ℓℓ relative band applies); seeds independent of the ℓℓ sets; `vv13600` builds `ggzh-install-v3` (v2 plus `fo-ggvvh`); `production/compare_nunu.py` compares with the ℓℓ sets |
| `atlas13600` | the reweighting study with an ATLAS-style LO+PS source (`rwt/README.md`): LO events with the ATLAS scale H_T (`dynscale 3`) and PDF4LHC21, POWHEG LO with the same card (check of the scale option), the A14 tables of the source and of the target (the `13600` production NLO events), classifier weights, summary and model `rwt/model-atlas-13600.tgz`; builds `ggzh-install-v4` (v3 plus the scale option, inert unless a card sets `dynscale 3`) |

```sh
production/paper.sh list 13600                         # runs: <run> <jobs> <cpus> <memory GB> <hours>
production/paper.sh show 13600 events-nlo 7            # environment and command of job 7
production/paper.sh inputs 13600 shower-nlo-cp5 3      # files it needs: <file in job dir> <product>
production/paper.sh outputs 13600 grids-nlo 0          # products it makes
production/paper.sh run 13600 events-nlo 7             # execute it in the current directory
```

**Job directory.** The unpacked package (`git archive <commit> VERSION production generator
shower`) plus the inputs.  Products have storage-independent names (`install/ggzh-install.tgz`,
`grids/grids-nlo-13600.tgz`, `lhe/nlo-13600/events-10001.lhe.gz`, ...); the storage is up to the
executor.  Every other file a job writes (out-*.tgz, *.json, logs) is its result.  `VERSION`
holds the package commit (written by `git archive`); every job prints it and every run manifest
and build manifest records it.

**Executing.** `run_local.sh <set> <workdir>` runs all jobs of a set one after another with a
local product store (the reference implementation; hundreds of CPU hours per energy).  On a batch
system, submit each run as `<jobs>` jobs with the listed resources, each executing
`production/paper.sh run <set> <run> $(Process)` after copying its inputs in and before copying
its outputs to the store; runs depend on each other only through the products, in list order
(builds → grids → events → showers).  Seeds do not depend on where or in which order jobs run,
and the shower seed of a file is fixed by the file, so any distribution of the work gives the
same results (NLO integrations: statistically the same across compilers and CPUs, see
`generator/README.md`).

**Results.**  `analyze.sh <set> <results> <outdir>` combines the job results of a set into the
files the paper quotes and runs every check: `<results>/<run>/<p>/` holds what job p left
(`run_local.sh`'s `<workdir>/jobs/<set>` has this layout).  It calls the analysis scripts
(`generator/analysis/fo_combine.py`, `powheg_stat.py`, `lhe_xsec.py --merge`;
`shower/analysis/yoda_summary.py --merge` with exact sums over jobs, `compare_orders.py`,
`compare_samples.py`, `shower_effect.py`: fixed order vs parton level vs shower at each order and
the high-p_T test of the NLO events) and the checks (`generator/check.py`, with `--paper-energy` for the
production sets, and `shower/check_shower.py` for every shower).  `STATUS.txt` lists every run with
its jobs and the package commits they ran with; a run is used only when all its jobs are present.
Minutes on one core.

```sh
production/run_local.sh 13600 work && production/analyze.sh 13600 work/jobs/13600 results/13600
```
