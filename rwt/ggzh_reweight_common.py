#!/usr/bin/env python3
"""Shared definitions of the LO+PS -> NLO+PS reweighting study (rwt/README.md): sample loading, the m_ll
window, the map binning (pT(V) x N_jets30), the closure-variable binning, weighted histograms and the
Stage 1.3 POI layout.

Samples are the per-event .npz tables written by rwt/ggzh_table.py from the Rivet GGZH_EVENTTABLE output;
the column names (and the NanoAOD order of the nine scale weights) follow the per-event tables this code
was first written for.  The source sample is called "official" in the code (keys, labels).
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

NTUPLE_DIR = Path(os.environ.get("GGZH_RWT_DIR", "."))      # where Sample(name) looks for <name>.npz
MLL_WINDOW = (15.0, 150.0)                      # the NLO+PS generation cut on the LHE m_ll

# map binning: 25 GeV to 400 GeV (Stage 1.3 boundaries are edges), 400-600, > 600; x 0J / 1J / >=2J
PT_EDGES = np.array([25.0 * i for i in range(17)] + [600.0, np.inf])     # 18 pT bins
N_PT, N_JET = len(PT_EDGES) - 1, 3
FIT_PT_MIN = 250.0                             # above this the binned ratio is replaced by a smooth fit
STAGE13_EDGES = [0.0, 75.0, 150.0, 250.0, 400.0, 600.0, np.inf]

# NanoAOD LHEScaleWeight index -> LHE weight id, and the seven varied points of the target scale band:
# all except (mu_R, mu_F) = (2, 0.5); with the nominal, the eight-point band of the pass criteria
NANO_TO_ID = {0: "1009", 1: "1007", 2: "1008", 3: "1003", 4: "1001", 5: "1002", 6: "1006", 7: "1004", 8: "1005"}
ID_TO_NANO = {v: k for k, v in NANO_TO_ID.items()}
SCALE7 = ["1002", "1004", "1005", "1003", "1007", "1009", "1008"]
MU_OF_ID = {"1001": (1, 1), "1002": (1, 2), "1003": (1, 0.5), "1004": (2, 1), "1005": (2, 2), "1006": (2, 0.5),
            "1007": (0.5, 1), "1008": (0.5, 2), "1009": (0.5, 0.5)}

# closure / characterisation variables: name -> (edges, axis label, note)
CLOSURE_VARS = {
    "pth": (np.arange(0.0, 1000.0 + 1e-9, 25.0), "pT(H) [GeV]", "last-copy generator Higgs"),
    "ptzh": (np.arange(0.0, 300.0 + 1e-9, 10.0), "pT(ZH) [GeV]", "Z + H system"),
    "mzh": (np.arange(150.0, 1250.0 + 1e-9, 25.0), "m(ZH) [GeV]", "Z + H invariant mass"),
    "ptj1": (np.arange(25.0, 1000.0 + 1e-9, 25.0), "leading ISR-jet pT [GeV]", "GenJet pT>30, |eta|<4.7, away from H b quarks and Z leptons; events with a jet only"),
    "yh": (np.arange(-4.0, 4.0 + 1e-9, 0.4), "y(H)", "Higgs rapidity"),
    "drzh": (np.arange(0.0, 6.4 + 1e-9, 0.2), "dR(Z,H)", ""),
    "lhe_mll": (np.arange(15.0, 150.0 + 1e-9, 5.0), "m(ll) at LHE [GeV]", "the Z line shape, unchanged by the shower"),
}
OTHER_VARS = {
    "ptz": (np.arange(0.0, 1000.0 + 1e-9, 25.0), "pT(Z) [GeV]", "last-copy generator Z (the map variable)"),
    "lhe_vpt": (np.arange(0.0, 1000.0 + 1e-9, 25.0), "pT(V) at LHE [GeV]", "before the shower"),
    "yzh": (np.arange(-3.0, 3.0 + 1e-9, 0.3), "y(ZH)", ""),
    "njets": (np.array([-0.5, 0.5, 1.5, 2.5, 3.5]), "HTXS N_jets30", ">=3 merged"),
}
POIS = {"PTV_0_75": (0, [0, 1, 2]), "PTV_75_150": (1, [0, 1, 2]), "PTV_150_250_0J": (2, [0]), "PTV_150_250_GE1J": (2, [1, 2]),
        "PTV_250_400_0J": (3, [0]), "PTV_250_400_GE1J": (3, [1, 2]), "PTV_400_600": (4, [0, 1, 2]), "PTV_GE600": (5, [0, 1, 2])}
CELL_LABELS = [f"{p}_{j}" for p in ["0_75", "75_150", "150_250", "250_400", "400_600", "GE600"] for j in ["0J", "1J", "GE2J"]]

# fixed categorical colours (dataviz reference palette, light mode), assigned by series identity
COLORS = {"nlo": "#2a78d6", "official": "#eb6834", "reweighted": "#1baf7a", "lo": "#eda100", "band": "#2a78d6"}
LABELS = {"nlo": "NLO+PS (target)", "official": "LO+PS (source)", "reweighted": "LO+PS, reweighted",
          "lo": "LO+PS (reference)"}


class Sample:
    """One sample; weights: nominal w0 (1 for an unweighted sample), sw = scale weights / nominal."""

    def __init__(self, name: str, path: Path | None = None, cols: dict | None = None):
        self.name = name
        self.path = path or (NTUPLE_DIR / f"{name}.npz")
        if cols is None:
            d = np.load(self.path)
            cols = {k: d[k] for k in d.files}
            self.meta = json.loads(self.path.with_suffix(".json").read_text()) if self.path.with_suffix(".json").exists() else {}
        else:
            self.meta = {}
        self.cols = cols
        self.n = len(self.cols["ptz"])
        w0 = self.cols["w0"].astype(np.float64)
        self.unweighted = bool(np.all(w0 == w0[0]))
        self.w0 = np.ones(self.n) if self.unweighted else w0
        self.sw = self.cols["sw"].astype(np.float64)
        mll = self.cols["lhe_mll"]
        self.in_window = (mll > MLL_WINDOW[0]) & (mll < MLL_WINDOW[1])
        self.njet = np.minimum(self.cols["htxs_njets30"].astype(int), 2)
        self.ptz = self.cols["ptz"].astype(np.float64)
        self.ipt = np.clip(np.searchsorted(PT_EDGES, self.ptz, side="right") - 1, 0, N_PT - 1)
        self.mapbin = self.ipt * N_JET + self.njet              # 0..53
        self.cell = self.cols["cell"].astype(int)                # -1 = FWDH

    def __getitem__(self, k):
        return self.cols[k]

    def resample(self, idx: np.ndarray) -> "Sample":
        """Bootstrap replica: the events at the given indices (with repetition)."""
        return Sample(self.name, self.path, cols={k: v[idx] for k, v in self.cols.items()})

    def weight(self, wid: str = "1001") -> np.ndarray:
        """Nominal POWHEG weight times the scale-variation ratio for LHE id wid."""
        if wid == "1001":
            return self.w0
        return self.w0 * self.sw[:, ID_TO_NANO[wid]]

    def sigma_fb(self) -> float:
        """sum of weights / number of events, in fb (POWHEG weights in pb; meaningless for an unweighted sample)."""
        return 1e3 * float(np.sum(self.w0)) / self.n


def hist(x, w, edges, w2=None):
    """Weighted histogram inside [edges[0], edges[-1]]: (sum w, sum w^2, raw counts) per bin."""
    x = np.asarray(x, dtype=np.float64)
    ok = np.isfinite(x) & (x >= edges[0]) & (x < edges[-1])
    idx = np.searchsorted(edges, x[ok], side="right") - 1
    nb = len(edges) - 1
    s = np.bincount(idx, weights=w[ok], minlength=nb)
    s2 = np.bincount(idx, weights=(w[ok] ** 2 if w2 is None else w2[ok]), minlength=nb)
    c = np.bincount(idx, minlength=nb).astype(float)
    return s, s2, c


def fractions(s, s2):
    """Normalised shape and its MC statistical error (bin-wise, normalisation treated as fixed)."""
    tot = s.sum()
    return s / tot, np.sqrt(s2) / tot


def map_fractions(sample: Sample, w: np.ndarray, sel: np.ndarray):
    """Fraction of the selected events (weighted) in each of the 54 map bins, with sum w^2 and raw counts."""
    nb = N_PT * N_JET
    s = np.bincount(sample.mapbin[sel], weights=w[sel], minlength=nb)
    s2 = np.bincount(sample.mapbin[sel], weights=w[sel] ** 2, minlength=nb)
    c = np.bincount(sample.mapbin[sel], minlength=nb).astype(float)
    return s / s.sum(), np.sqrt(s2) / s.sum(), c, s


def cell_sums(sample: Sample, w: np.ndarray, sel: np.ndarray):
    """Sum of weights, sum w^2 and raw counts in the 18 Stage 1.3 cells (FWDH excluded)."""
    keep = sel & (sample.cell >= 0)
    s = np.bincount(sample.cell[keep], weights=w[keep], minlength=18)
    s2 = np.bincount(sample.cell[keep], weights=w[keep] ** 2, minlength=18)
    c = np.bincount(sample.cell[keep], minlength=18).astype(float)
    return s, s2, c


def poi_sums(cells: np.ndarray) -> dict:
    """Aggregate an 18-cell array into the 8 POIs (cells are pT-major: index = q * 3 + jet)."""
    return {p: float(sum(cells[q * 3 + j] for j in jets)) for p, (q, jets) in POIS.items()}


def merge_sparse_bins(edges, counts_list, min_count=20):
    """Merge adjacent bins so that every merged bin has >= min_count raw events in EVERY sample in
    counts_list; scanning left to right, a sparse bin absorbs its right neighbour; a sparse last bin
    merges into the previous one. Returns the list of (lo, hi) index slices into the original bins."""
    nb = len(edges) - 1
    groups, i = [], 0
    while i < nb:
        j = i + 1
        while j < nb and any(c[i:j].sum() < min_count for c in counts_list):
            j += 1
        groups.append((i, j))
        i = j
    if len(groups) > 1 and any(c[groups[-1][0]:groups[-1][1]].sum() < min_count for c in counts_list):
        lo, _ = groups.pop()
        groups[-1] = (groups[-1][0], nb)
    return groups


def regroup(arr, groups):
    return np.array([arr[lo:hi].sum() for lo, hi in groups])


def chi2_shape(f1, e1, f2, e2):
    """chi2 between two normalised shapes with independent bin errors; ndof = nbins - 1."""
    from scipy.stats import chi2 as chi2dist
    var = e1 ** 2 + e2 ** 2
    ok = var > 0
    c2 = float(np.sum((f1[ok] - f2[ok]) ** 2 / var[ok]))
    ndof = int(ok.sum()) - 1
    return c2, ndof, float(chi2dist.sf(c2, ndof)) if ndof > 0 else float("nan")


def dumps(obj) -> str:
    def conv(o):
        if isinstance(o, np.ndarray):
            return [conv(x) for x in o.tolist()]
        if isinstance(o, (np.floating, float)):
            return None if (isinstance(o, float) and (np.isnan(o) or np.isinf(o))) else float(o)
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, dict):
            return {str(k): conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [conv(x) for x in o]
        return o
    return json.dumps(conv(obj), indent=1) + "\n"
