#!/usr/bin/env python3
"""Rescaled LO+PS against NLO+PS in the STXS Stage 1.3 bins and in p_T(Z).

    compare_orders.py --nlo nlo.json --lo lo-default.json [--lo lo-scalup-half.json ...] --out cmp.json

Inputs: finished summaries of analysis/yoda_summary.py (NLO+PS from NLO events with the POWHEG
veto; LO+PS from Born-only events with a given shower start).  K = sigma_NLO / sigma_LO from the
inclusive cross sections of the two samples (nominal weights); "rescaled LO+PS" = K x LO+PS.
Seven-point scale envelopes from the MUR/MUF weight variations of each sample.  The scale band of
rescaled LO+PS is defined as in the paper (B, "uncorrelated"): the LO+PS variations as shapes only
(each normalised to its own inclusive cross section), combined in quadrature with the inclusive
NLO scale band, which sets the normalisation; for comparison also (A) K times the full LO+PS
envelope and (C, "correlated") the envelope of sigma_NLO(point) x LO+PS shape(point).  Output per
LO+PS variant: for each of the 18 Stage 1.3 fine bins (FWDH excluded) and for p_T(Z) fractions
above 75, 150, 250, 400, 600 GeV: NLO+PS and rescaled LO+PS cross sections with bands and MC
errors, their ratio, and whether each band contains the NLO+PS central value ("coverage", counted
over the bins whose NLO+PS MC error is below 10 %).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

SEVEN = ["MUR1.0_MUF1.0", "MUR0.5_MUF0.5", "MUR2.0_MUF2.0", "MUR2.0_MUF1.0", "MUR0.5_MUF1.0", "MUR1.0_MUF2.0", "MUR1.0_MUF0.5"]
PTV = ["0-75", "75-150", "150-250", "250-400", "400-600", ">600"]
JETS = ["0J", "1J", ">=2J"]
CELLS = [f"{p} {j}" for j in JETS for p in PTV]          # the routine's (jet-major) order after FWDH
BOUNDS = [75.0, 150.0, 250.0, 400.0, 600.0]


def sigma(s, var=""):
    return sum(s["hists"]["GGZH_NLOPS/sigma"]["var"][var]["xs"])


def cells(s, var=""):
    d = s["stxs13_fine"][var]
    return d["xs"][1:], d["err"][1:]                    # drop FWDH


def fractions(s, var=""):
    """sigma above each boundary of p_T(Z) / sigma (incl. FWDH, overflow counted)."""
    h = s["hists"]["GGZH_NLOPS/ptz"]
    edges, xs = h["edges"], h["var"][var]["xs"]        # xs: underflow, bins..., overflow
    tot = sum(xs)
    out = {}
    for b in BOUNDS:
        above = sum(x for lo, x in zip(edges[:-1], xs[1:-1]) if lo >= b) + xs[-1]
        out[str(int(b))] = above / tot
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nlo", type=Path, required=True)
    ap.add_argument("--lo", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    nlo = json.loads(a.nlo.read_text())
    res = {"nlo": str(a.nlo), "sigma_nlo_fb": 1e3 * sigma(nlo), "variants": {}}
    nc, ne = cells(nlo)
    nenv = [cells(nlo, v)[0] for v in SEVEN]
    snlo = [sigma(nlo, v) for v in SEVEN]
    nrm_up, nrm_dn = max(snlo) / snlo[0] - 1, min(snlo) / snlo[0] - 1      # inclusive NLO band
    res["nlo_inclusive_band_rel"] = [nrm_up, nrm_dn]
    for lo_path in a.lo:
        lo = json.loads(lo_path.read_text())
        k = sigma(nlo) / sigma(lo)
        lc, le = cells(lo)
        lenv = [cells(lo, v)[0] for v in SEVEN]
        slo = [sigma(lo, v) for v in SEVEN]
        rows = []
        for i, name in enumerate(CELLS):
            n_up, n_dn = max(e[i] for e in nenv), min(e[i] for e in nenv)
            l_up, l_dn = k * max(e[i] for e in lenv), k * min(e[i] for e in lenv)
            ratio = k * lc[i] / nc[i] if nc[i] > 0 else float("nan")
            rerr = ratio * math.hypot(le[i] / lc[i] if lc[i] else 0, ne[i] / nc[i] if nc[i] else 0) if nc[i] > 0 else float("nan")
            c = k * lc[i]
            shape = [(e[i] / s) / (lc[i] / sigma(lo)) - 1 if lc[i] else 0.0 for e, s in zip(lenv, slo)]
            b_up = math.hypot(max(max(shape), 0.0), nrm_up); b_dn = math.hypot(min(min(shape), 0.0), nrm_dn)
            band_b = [c * (1 - b_dn), c * (1 + b_up)]
            corr = [s_n * e[i] / s for s_n, e, s in zip(snlo, lenv, slo)]
            band_c = [min(corr), max(corr)]
            inside = lambda b: b[0] <= nc[i] <= b[1]
            rows.append({"cell": name, "nlo_ps_fb": 1e3 * nc[i], "nlo_ps_err_fb": 1e3 * ne[i],
                         "nlo_ps_band_fb": [1e3 * n_dn, 1e3 * n_up],
                         "k_lo_ps_fb": 1e3 * k * lc[i], "k_lo_ps_err_fb": 1e3 * k * le[i],
                         "k_lo_ps_band_fb": [1e3 * x for x in band_b],
                         "k_lo_ps_band_A_fb": [1e3 * l_dn, 1e3 * l_up], "k_lo_ps_band_C_fb": [1e3 * x for x in band_c],
                         "covered": {"A": inside([l_dn, l_up]), "B": inside(band_b), "C": inside(band_c)},
                         "precise": bool(nc[i] > 0 and ne[i] < 0.1 * nc[i]),
                         "ratio_klo_over_nlo": ratio, "ratio_err": rerr})
        cov = {d: [r["cell"] for r in rows if r["precise"] and r["covered"][d]] for d in "ABC"}
        coverage = {"bins": sum(r["precise"] for r in rows), **{d: len(v) for d, v in cov.items()},
                    "not_covered_B": [r["cell"] for r in rows if r["precise"] and not r["covered"]["B"]]}
        fn, fl = fractions(nlo), fractions(lo)
        res["variants"][lo_path.stem] = {"lo": str(lo_path), "sigma_lo_fb": 1e3 * sigma(lo), "K": k, "cells": rows,
                                         "coverage": coverage,
                                         "fraction_above_nlo_ps": fn, "fraction_above_lo_ps": fl,
                                         "fraction_ratio_lo_over_nlo": {b: fl[b] / fn[b] if fn[b] else float("nan") for b in fn}}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    print(f"NLO+PS: {res['sigma_nlo_fb']:.3f} fb")
    for name, v in res["variants"].items():
        print(f"\n{name}: LO+PS {v['sigma_lo_fb']:.3f} fb, K = {v['K']:.3f}")
        print("  fraction above pT(Z) b, LO+PS / NLO+PS: " +
              ", ".join(f"{b}: {r:.2f}" for b, r in v["fraction_ratio_lo_over_nlo"].items()))
        for r in v["cells"]:
            print(f"  {r['cell']:16s} NLO+PS {r['nlo_ps_fb']:8.4f} fb   K x LO+PS {r['k_lo_ps_fb']:8.4f} fb   "
                  f"ratio {r['ratio_klo_over_nlo']:6.2f} +- {r['ratio_err']:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
