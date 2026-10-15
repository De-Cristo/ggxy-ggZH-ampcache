#!/usr/bin/env python3
"""The binned LO+PS -> NLO+PS reweighting map in generator-level pT(V) x N_jets30 (rwt/README.md).

Method (fixed before the run):
  w(pT, Nj) = f_NLO+PS(pT, Nj) / f_source(pT, Nj),  f = fraction of the sample in the bin,
  bins 25 GeV up to 400 GeV, 400-600, > 600, times 0J / 1J / >=2J (HTXS_njets30);
  above 250 GeV the binned ratio is replaced per jet bin by a smooth fit: a polynomial in
  log pT of log(ratio), the lowest degree (1..3) whose chi2 probability exceeds 5 %, fitted to
  the bins at and above 250 GeV at the source sample's mean pT in each bin (if no degree
  passes, the best one is used and flagged);
  shape only: the weights are scaled to mean 1 over the source events inside
  15 < m_ll < 150 GeV (LHE); events outside the window get weight 1;
  seven more maps from the NLO+PS scale variations (source sample at nominal);
  statistical uncertainty per event from the bin errors (below 250 GeV) or the fit covariance.

    ggzh_reweight_map.py --official <source>.npz --nlo <target>.npz --out map.json [--weights-out weights.npz]
The build_map / apply_map functions are reused by ggzh_reweight_closure.py (half-sample test).
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2dist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ggzh_reweight_common import (SCALE7, FIT_PT_MIN, N_JET, N_PT, NTUPLE_DIR, PT_EDGES, Sample, dumps,  # noqa: E402
                                  map_fractions)

MAX_DEGREE = 3
P_MIN = 0.05
FIRST_FIT_BIN = int(np.searchsorted(PT_EDGES, FIT_PT_MIN))          # index of the 250-275 bin


@dataclass
class JetFit:
    degree: int
    coef: np.ndarray            # highest power first (np.polyfit convention), for log r vs log pT
    cov: np.ndarray
    chi2: float
    ndof: int
    pvalue: float
    passed: bool
    x: np.ndarray               # log of the mean official pT per fitted bin
    y: np.ndarray               # log r
    sy: np.ndarray
    bins: list
    tried: list = field(default_factory=list)


@dataclass
class Map:
    wid: str
    r: np.ndarray               # (54,) binned ratio, nan where undefined
    dr: np.ndarray
    f_nlo: np.ndarray
    df_nlo: np.ndarray
    c_nlo: np.ndarray
    f_off: np.ndarray
    df_off: np.ndarray
    c_off: np.ndarray
    fits: list                  # per jet bin
    norm: float = 1.0           # divisor that makes the mean weight 1 over the in-window official events

    def raw_weight(self, ptz, ipt, njet):
        """Un-normalised weight and its statistical uncertainty for events with (ptz, ipt, njet)."""
        w = np.ones(len(ptz)); dw = np.zeros(len(ptz))
        low = ipt < FIRST_FIT_BIN
        b = ipt * N_JET + njet
        w[low] = np.nan_to_num(self.r[b[low]], nan=1.0); dw[low] = np.nan_to_num(self.dr[b[low]], nan=0.0)
        for j in range(N_JET):
            m = (~low) & (njet == j)
            if not m.any():
                continue
            fit = self.fits[j]
            x = np.log(ptz[m])
            X = np.vander(x, fit.degree + 1)            # highest power first
            lr = X @ fit.coef
            var = np.einsum("ij,jk,ik->i", X, fit.cov, X)
            w[m] = np.exp(lr); dw[m] = w[m] * np.sqrt(np.maximum(var, 0.0))
        return w, dw

    def to_json(self):
        return {"wid": self.wid, "norm": self.norm, "r": self.r, "dr": self.dr, "f_nlo": self.f_nlo, "df_nlo": self.df_nlo, "c_nlo": self.c_nlo,
                "f_off": self.f_off, "df_off": self.df_off, "c_off": self.c_off,
                "fits": [{"jet_bin": j, "degree": f.degree, "coef_highest_first": f.coef, "cov": f.cov, "chi2": f.chi2, "ndof": f.ndof,
                          "pvalue": f.pvalue, "passed_p_gt_0p05": f.passed, "log_pt_points": f.x, "log_r_points": f.y, "log_r_err": f.sy,
                          "bins_used": f.bins, "tried": f.tried} for j, f in enumerate(self.fits)]}


def fit_jet_bin(x, y, sy):
    tried = []
    best = None
    if len(x) < 3:                      # too few tail points (tiny test samples): flagged weighted mean
        if len(x) == 0:
            return None, tried, False
        wgt = 1.0 / sy ** 2
        c0 = float(np.sum(wgt * y) / np.sum(wgt)); var = 1.0 / np.sum(wgt)
        c2 = float(np.sum(((y - c0) / sy) ** 2)); ndof = len(x) - 1
        tried.append({"degree": 0, "chi2": c2, "ndof": ndof, "pvalue": float(chi2dist.sf(c2, ndof)) if ndof > 0 else None})
        return (0, np.array([c0]), np.array([[var]]), c2, ndof, tried[-1]["pvalue"]), tried, False
    for d in range(1, MAX_DEGREE + 1):
        if len(x) < d + 2:
            break
        coef, cov = np.polyfit(x, y, d, w=1.0 / sy, cov="unscaled")
        res = (y - np.polyval(coef, x)) / sy
        c2 = float(np.sum(res ** 2)); ndof = len(x) - (d + 1); p = float(chi2dist.sf(c2, ndof))
        tried.append({"degree": d, "chi2": c2, "ndof": ndof, "pvalue": p})
        cand = (d, coef, cov, c2, ndof, p)
        if p > P_MIN:
            return cand, tried, True
        if best is None or p > best[5]:
            best = cand
    return best, tried, False


def build_map(official: Sample, nlo: Sample, wid: str = "1001", off_sel=None, nlo_sel=None) -> Map:
    off_sel = official.in_window if off_sel is None else off_sel
    nlo_sel = nlo.in_window if nlo_sel is None else nlo_sel
    f_nlo, df_nlo, c_nlo, _ = map_fractions(nlo, nlo.weight(wid), nlo_sel)
    f_off, df_off, c_off, _ = map_fractions(official, official.weight("1001"), off_sel)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where((c_off > 0) & (f_off > 0), f_nlo / f_off, np.nan)
        dr = np.abs(r) * np.sqrt((df_nlo / f_nlo) ** 2 + (df_off / f_off) ** 2)
    # mean official pT per bin (for the fit abscissa)
    sel = off_sel
    mean_pt = np.bincount(official.mapbin[sel], weights=official.ptz[sel], minlength=N_PT * N_JET)
    mean_pt = np.where(c_off > 0, mean_pt / np.maximum(c_off, 1), np.nan)
    fits = []
    for j in range(N_JET):
        bins = [i * N_JET + j for i in range(FIRST_FIT_BIN, N_PT)]
        use = [b for b in bins if c_off[b] > 0 and c_nlo[b] > 0 and f_nlo[b] > 0 and np.isfinite(r[b]) and r[b] > 0 and np.isfinite(dr[b]) and dr[b] > 0]
        x = np.log(mean_pt[use]); y = np.log(r[use]); sy = dr[use] / r[use]
        cand, tried, passed = fit_jet_bin(x, y, sy)
        if cand is None:
            raise RuntimeError(f"jet bin {j}: not enough points to fit ({len(use)})")
        d, coef, cov, c2, ndof, p = cand
        fits.append(JetFit(d, coef, cov, c2, ndof, p, passed, x, y, sy, use, tried))
    m = Map(wid, r, dr, f_nlo, df_nlo, c_nlo, f_off, df_off, c_off, fits)
    w, _ = m.raw_weight(official.ptz[sel], official.ipt[sel], official.njet[sel])
    m.norm = float(np.mean(w))
    return m


def apply_map(m: Map, official: Sample, sel=None):
    """Per-event weight (mean 1 over the selected events, 1 outside) and its statistical uncertainty."""
    sel = official.in_window if sel is None else sel
    w = np.ones(official.n); dw = np.zeros(official.n)
    wr, dwr = m.raw_weight(official.ptz[sel], official.ipt[sel], official.njet[sel])
    w[sel] = wr / m.norm; dw[sel] = dwr / m.norm
    return w, dw


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--official", type=Path, required=True, help="the source LO+PS sample (.npz)")
    ap.add_argument("--nlo", type=Path, required=True, help="the NLO+PS target (.npz)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--weights-out", type=Path, default=None)
    args = ap.parse_args(argv)
    official, nlo = Sample("official", args.official), Sample("nlo", args.nlo)
    maps = {wid: build_map(official, nlo, wid) for wid in ["1001"] + SCALE7}
    weights = {}
    out = {"pt_edges": PT_EDGES.tolist(), "fit_pt_min": FIT_PT_MIN, "max_degree": MAX_DEGREE, "p_min": P_MIN,
           "official_events": official.n, "official_in_window": int(official.in_window.sum()),
           "nlo_events": nlo.n, "nlo_in_window": int(nlo.in_window.sum()), "maps": {}}
    sel = official.in_window
    for wid, m in maps.items():
        w, dw = apply_map(m, official)
        weights[f"w_{wid}"] = w; weights[f"dw_{wid}"] = dw
        ws = w[sel]
        j = m.to_json()
        j.update({"n_eff_in_window": float(ws.sum() ** 2 / np.sum(ws ** 2)), "n_eff_all": float(w.sum() ** 2 / np.sum(w ** 2)),
                  "weight_quantiles": dict(zip(["0", "1", "5", "50", "95", "99", "100"], np.percentile(ws, [0, 1, 5, 50, 95, 99, 100]).tolist())),
                  "mean_weight_in_window": float(ws.mean()), "frac_events_weight_below_0p1": float(np.mean(ws < 0.1)),
                  "mean_rel_stat_err": float(np.mean(dw[sel] / w[sel]))})
        out["maps"][wid] = j
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(dumps(out))
    if args.weights_out:
        args.weights_out.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.weights_out, run=official["run"], lumi=official["lumi"], event=official["event"], ifile=official["ifile"],
                            in_window=sel, **weights)
    m = maps["1001"]
    print(f"official {official.n} events ({sel.sum()} in window), NLO {nlo.n}; nominal map: norm {m.norm:.4f}, "
          f"N_eff {out['maps']['1001']['n_eff_in_window']:.0f} of {sel.sum()} ({100*out['maps']['1001']['n_eff_in_window']/sel.sum():.1f} %)")
    for j, f in enumerate(m.fits):
        print(f"  jet bin {j}: degree {f.degree}, chi2/ndof {f.chi2:.1f}/{f.ndof}, p = {f.pvalue:.3f} {'ok' if f.passed else 'NO DEGREE PASSES'}; "
              f"points {len(f.x)}; r at 250-275 / 400-600 / >600: " + " ".join(f"{m.r[b]:.3f}" for b in (FIRST_FIT_BIN * N_JET + j, (N_PT - 2) * N_JET + j, (N_PT - 1) * N_JET + j)))
    print(f"  weight quantiles (1,50,99 %): " + " ".join(f"{out['maps']['1001']['weight_quantiles'][q]:.3f}" for q in ("1", "50", "99")))
    print(f"wrote {args.out}" + (f" and {args.weights_out}" if args.weights_out else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
