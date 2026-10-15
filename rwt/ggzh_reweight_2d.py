#!/usr/bin/env python3
"""LO+PS -> NLO+PS reweighting of gg -> ZH in generator-level m(ZH) x pT(ZH), with the pass criteria fixed
before this script was run on the full samples (rwt/README.md).  Not used in the paper; the classifier
(ggzh_reweight_bdt.py) uses its helpers shapes and test_vars:
  1. the 8 Stage 1.3 POI pulls < 2;  2. pT(Z), N_jets, pT(H), leading-jet pT, y(H), dR(Z,H), m_ll:
  chi2 p > 5 % and no merged bin outside the NLO 8-point band by > 2 sigma;  3. both half-sample
  directions;  D. optionally the map applied to another LO+PS sample (--apply-to), criteria 1-2.
Map: cells of m(ZH) edges M_EDGES x pT(ZH) edges P_EDGES; within each m(ZH) slice pT(ZH) cells are merged
from the high end down until each has >= MIN_CELL raw events in both samples and a positive NLO sum;
a slice that cannot be filled is merged with the previous slice.  Weight = f_NLO / f_LO per cell.

    ggzh_reweight_2d.py --source <source>.npz --target <target>.npz [--apply-to <other>.npz] --out 2d.json [--fig-dir figs]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ggzh_reweight_common as C  # noqa: E402
from ggzh_reweight_common import SCALE7, CLOSURE_VARS, OTHER_VARS, POIS, Sample, chi2_shape, dumps, fractions, hist, merge_sparse_bins, regroup  # noqa: E402
from ggzh_reweight_closure import test_pois  # noqa: E402

M_EDGES = np.array([0.0, 250, 300, 350, 400, 450, 500, 600, 800, np.inf])
P_EDGES = np.array([0.0, 10, 20, 30, 45, 60, 80, 100, 130, 170, 220, 300, np.inf])
MIN_CELL = 50
MIN_RAW = 20
TEST_VARS = {"ptz": OTHER_VARS["ptz"], "njets": OTHER_VARS["njets"], "pth": CLOSURE_VARS["pth"], "ptj1": CLOSURE_VARS["ptj1"],
             "yh": CLOSURE_VARS["yh"], "drzh": CLOSURE_VARS["drzh"], "lhe_mll": CLOSURE_VARS["lhe_mll"]}
CHECK_VARS = {"mzh": CLOSURE_VARS["mzh"], "ptzh": CLOSURE_VARS["ptzh"]}


def values(s: Sample, name):
    return np.minimum(s["htxs_njets30"].astype(float), 3) if name == "njets" else s[name].astype(float)


def fine_index(s: Sample):
    im = np.clip(np.searchsorted(M_EDGES, s["mzh"].astype(float), side="right") - 1, 0, len(M_EDGES) - 2)
    ip = np.clip(np.searchsorted(P_EDGES, s["ptzh"].astype(float), side="right") - 1, 0, len(P_EDGES) - 2)
    return im, ip


class Map2D:
    def __init__(self, src: Sample, tgt: Sample, src_sel, tgt_sel, wids=("1001",)):
        nm, npt = len(M_EDGES) - 1, len(P_EDGES) - 1
        ims, ips = fine_index(src); imt, ipt = fine_index(tgt)
        cs = np.zeros((nm, npt)); ct = np.zeros((nm, npt)); wt = np.zeros((nm, npt))
        np.add.at(cs, (ims[src_sel], ips[src_sel]), 1.0)
        np.add.at(ct, (imt[tgt_sel], ipt[tgt_sel]), 1.0)
        np.add.at(wt, (imt[tgt_sel], ipt[tgt_sel]), tgt.w0[tgt_sel])
        # slices of m(ZH): merge with the previous slice when a slice cannot form one valid cell
        slices, cur = [], []
        for i in range(nm):
            cur.append(i)
            if cs[cur].sum() >= MIN_CELL and ct[cur].sum() >= MIN_CELL and wt[cur].sum() > 0:
                slices.append(cur); cur = []
        if cur:
            slices[-1] = slices[-1] + cur
        self.cell = np.full((nm, npt), -1, dtype=int); self.cells = []
        for sl in slices:
            s_cs, s_ct, s_wt = cs[sl].sum(0), ct[sl].sum(0), wt[sl].sum(0)
            groups, g = [], []
            for j in reversed(range(npt)):
                g.append(j)
                if s_cs[g].sum() >= MIN_CELL and s_ct[g].sum() >= MIN_CELL and s_wt[g].sum() > 0:
                    groups.append(g); g = []
            if g:
                groups[-1] = groups[-1] + g
            for g in groups:
                k = len(self.cells)
                self.cells.append({"m_slices": [int(i) for i in sl], "pt_bins": sorted(int(j) for j in g),
                                   "m_range": [float(M_EDGES[min(sl)]), float(M_EDGES[max(sl) + 1])],
                                   "ptzh_range": [float(P_EDGES[min(g)]), float(P_EDGES[max(g) + 1])],
                                   "raw_source": float(cs[np.ix_(sl, g)].sum()), "raw_target": float(ct[np.ix_(sl, g)].sum())})
                for i in sl:
                    self.cell[i, g] = k
        ncell = len(self.cells)
        ks = self.cell[ims[src_sel], ips[src_sel]]; kt = self.cell[imt[tgt_sel], ipt[tgt_sel]]
        fs = np.bincount(ks, weights=src.w0[src_sel], minlength=ncell); fs /= fs.sum()
        fs2 = np.bincount(ks, weights=src.w0[src_sel] ** 2, minlength=ncell) / src.w0[src_sel].sum() ** 2
        self.r, self.dr = {}, None
        for wid in wids:
            w = tgt.weight(wid)[tgt_sel]
            ft = np.bincount(kt, weights=w, minlength=ncell); tot = ft.sum(); ft = ft / tot
            self.r[wid] = ft / fs
            if wid == "1001":
                ft2 = np.bincount(kt, weights=w ** 2, minlength=ncell) / tot ** 2
                self.dr = self.r[wid] * np.sqrt(ft2 / ft ** 2 + fs2 / fs ** 2)
        for k, c in enumerate(self.cells):
            c["ratio"] = float(self.r["1001"][k]); c["ratio_err"] = float(self.dr[k])

    def apply(self, s: Sample, sel, wid="1001"):
        im, ip = fine_index(s)
        w = np.ones(s.n)
        w[sel] = self.r[wid][self.cell[im[sel], ip[sel]]]
        w[sel] /= w[sel].mean()
        return w


def shapes(s: Sample, w, sel, vars_, band_weights=None):
    out = {}
    for name, (edges, _, _) in vars_.items():
        x = values(s, name)
        f, e = fractions(*hist(x[sel], w[sel], edges)[:2])
        d = {"f": f, "e": e, "raw": hist(x[sel], w[sel], edges)[2]}
        if band_weights is not None:
            fv = np.array([f] + [fractions(*hist(x[sel], bw[sel], edges)[:2])[0] for bw in band_weights])
            d["variations"] = fv
        out[name] = d
    return out


def test_vars(test, target, vars_):
    res, ok_all = {}, True
    for name, (edges, _, _) in vars_.items():
        t, g = test[name], target[name]
        groups = merge_sparse_bins(edges, [t["raw"], g["raw"]], MIN_RAW)
        ft, fg = regroup(t["f"], groups), regroup(g["f"], groups)
        et, eg = np.sqrt(regroup(t["e"] ** 2, groups)), np.sqrt(regroup(g["e"] ** 2, groups))
        fv = np.array([regroup(v, groups) for v in g["variations"]]); lo, hi = fv.min(0), fv.max(0)
        c2, ndof, p = chi2_shape(ft, et, fg, eg)
        sig = np.sqrt(et ** 2 + eg ** 2)
        flags = (ft > hi + 2 * sig) | (ft < lo - 2 * sig)
        ok = bool(p > 0.05 and not flags.any()); ok_all &= ok
        res[name] = {"chi2": c2, "ndof": ndof, "pvalue": p, "bins_outside_band_2sigma": int(flags.sum()), "pass": ok,
                     "max_abs_pull": float(np.max(np.abs((ft - fg) / sig))), "edges": [edges[a] for a, _ in groups] + [edges[groups[-1][1]]],
                     "f_test": ft, "e_test": et, "f_target": fg, "e_target": eg, "band_lo": lo, "band_hi": hi, "flags": flags}
    return res, bool(ok_all)


def run(src, tgt, src_sel, tgt_build, tgt_test, other=None, tag="", log=print):
    m = Map2D(src, tgt, src_sel, tgt_build, wids=["1001"] + SCALE7)
    w = m.apply(src, src_sel)
    band = [tgt.weight(wid) for wid in SCALE7]
    target = shapes(tgt, tgt.w0, tgt_test, {**TEST_VARS, **CHECK_VARS}, band)
    res = {"tag": tag, "cells": m.cells, "n_cells": len(m.cells),
           "n_eff": float(w[src_sel].sum() ** 2 / np.sum(w[src_sel] ** 2)),
           "weight_quantiles": dict(zip(["0", "1", "50", "99", "100"], np.percentile(w[src_sel], [0, 1, 50, 99, 100]).tolist()))}
    samples = {"reweighted": (src, w, src_sel), "unweighted": (src, np.ones(src.n), src_sel)}
    if other is not None:
        osel = other.in_window
        samples["applied_reweighted"] = (other, m.apply(other, osel), osel)
        samples["applied_unweighted"] = (other, np.ones(other.n), osel)
    for name, (s, ws, sel) in samples.items():
        sh = shapes(s, ws, sel, {**TEST_VARS, **CHECK_VARS})
        tv, ok = test_vars({k: sh[k] for k in TEST_VARS}, {k: target[k] for k in TEST_VARS}, TEST_VARS)
        cv, _ = test_vars({k: sh[k] for k in CHECK_VARS}, {k: target[k] for k in CHECK_VARS}, CHECK_VARS)
        st = test_pois(s, ws, sel, tgt, tgt.w0, tgt_test)
        res[name] = {"stage13": st, "variables": tv, "map_check": cv, "pass_1": st["pass"], "pass_2": ok}
    log(f"[{tag}] {len(m.cells)} cells; N_eff {res['n_eff']:.0f}")
    for name in samples:
        r = res[name]
        log(f"  {name:20s} criterion 1 {'PASS' if r['pass_1'] else 'FAIL'} (max |pull| {r['stage13']['max_abs_pull']:.2f}); criterion 2 {'PASS' if r['pass_2'] else 'FAIL'}"
            + "; " + ", ".join(f"{v} p={r['variables'][v]['pvalue']:.2g}/out {r['variables'][v]['bins_outside_band_2sigma']}" for v in TEST_VARS)
            + " | map check: " + ", ".join(f"{v} p={r['map_check'][v]['pvalue']:.2g}" for v in CHECK_VARS))
    return res, m


def figure_vars(path, res, keys, title, labels):
    from ggzh_reweight_plots import figure, plt, shape_panel
    saved = dict(C.LABELS); C.LABELS.update(labels)
    names = list(TEST_VARS) + list(CHECK_VARS)
    fig, pairs = figure(3, 3)
    for (ax, rax), name in zip(pairs, names):
        sec = "map_check" if name in CHECK_VARS else "variables"
        v = res[keys[0]][sec][name]; u = res[keys[1]][sec][name]
        ser = {"nlo": (v["f_target"], v["e_target"]), "official": (u["f_test"], u["e_test"]), "reweighted": (v["f_test"], v["e_test"])}
        edges = (TEST_VARS.get(name) or CHECK_VARS[name])
        shape_panel(ax, rax, v["edges"], ser, "nlo", band=(v["band_lo"], v["band_hi"]), xlabel=edges[1], logy=name not in ("yh", "drzh", "njets"),
                    title=f"{name}{' (map variable)' if name in CHECK_VARS else ''}: p = {v['pvalue']:.2f}", flags=v["flags"], legend=(name == names[0]))
    fig.suptitle(title, x=0.01, y=0.93, ha="left")
    fig.savefig(path, bbox_inches="tight"); plt.close(fig)
    C.LABELS.clear(); C.LABELS.update(saved)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--target", type=Path, required=True)
    ap.add_argument("--apply-to", type=Path, default=None)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--fig-dir", type=Path, default=None)
    args = ap.parse_args(argv)
    lines = []
    def log(s):
        print(s, flush=True); lines.append(s)
    src, tgt = Sample("source", args.source), Sample("target", args.target)
    other = Sample("applied", args.apply_to) if args.apply_to else None
    ss, ts = src.in_window, tgt.in_window
    idx = np.arange(tgt.n); A, B = ts & (idx % 2 == 0), ts & (idx % 2 == 1)
    full, mfull = run(src, tgt, ss, ts, ts, other, "full", log)
    ab, _ = run(src, tgt, ss, A, B, None, "map A -> test B", log)
    ba, _ = run(src, tgt, ss, B, A, None, "map B -> test A", log)
    summary = {"1_stage13": full["reweighted"]["pass_1"], "2_variables": full["reweighted"]["pass_2"],
               "3_stability": bool(ab["reweighted"]["pass_1"] and ab["reweighted"]["pass_2"] and ba["reweighted"]["pass_1"] and ba["reweighted"]["pass_2"]),
               "D_applied": (bool(full["applied_reweighted"]["pass_1"] and full["applied_reweighted"]["pass_2"]) if other else None)}
    summary["variables_failing"] = [v for v in TEST_VARS if not full["reweighted"]["variables"][v]["pass"]]
    log(f"summary: {summary}")
    out = {"plan": "rwt/README.md", "m_edges": M_EDGES, "ptzh_edges": P_EDGES, "min_cell": MIN_CELL,
           "inputs": {"source": str(args.source), "target": str(args.target), "applied": str(args.apply_to) if other else None},
           "full": full, "A_to_B": ab, "B_to_A": ba, "summary": summary, "log": lines}
    args.out.parent.mkdir(parents=True, exist_ok=True); args.out.write_text(dumps(out))
    if args.fig_dir:
        args.fig_dir.mkdir(parents=True, exist_ok=True)
        lab = {"nlo": "NLO+PS", "official": "CMS-style LO+PS", "reweighted": "CMS-style LO x rwt + PS"}
        figure_vars(args.fig_dir / "closure-2d.png", full, ("reweighted", "unweighted"), "m(ZH) x pT(ZH) map: CMS-style LO+PS reweighted to NLO+PS", lab)
        figure_vars(args.fig_dir / "closure-2d-halfA-to-B.png", ab, ("reweighted", "unweighted"), "m(ZH) x pT(ZH) map from target half A, tested on half B", lab)
        if other:
            figure_vars(args.fig_dir / "closure-2d-applied.png", full, ("applied_reweighted", "applied_unweighted"),
                        "the m(ZH) x pT(ZH) map applied to another LO+PS sample",
                        {"nlo": "NLO+PS (CP5)", "official": "other LO+PS", "reweighted": "other LO+PS, reweighted"})
        print(f"figures in {args.fig_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
