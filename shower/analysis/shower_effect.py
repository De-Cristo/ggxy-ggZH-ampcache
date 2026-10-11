#!/usr/bin/env python3
"""Where the excess of rescaled LO+PS comes from: fixed order, parton level and shower at each order.

    shower_effect.py --fo fo-llh.json --parton-nlo parton-nlo.json --parton-lo parton-lo.json
                     --nlo-ps shower-nlo-<tune>.json --lo-ps shower-lo-<start>-<tune>.json [...] --out effect.json

Inputs: generator/analysis/fo_combine.py of the ggxy fixed-order runs (l-l+H, central scale), and
finished summaries of analysis/yoda_summary.py: the NLO and LO events at parton level (run/shower.sh
mode noshower: the LHE kinematics), the same events after the shower (NLO+PS, LO+PS for each start).
Nominal weights; cross sections in fb.

  parton_vs_fo   p_T(Z) per bin of the fixed-order histogram (25-100 GeV wide, up to 1 TeV): events
                 at parton level against fixed order, at LO and at NLO, with pulls and chi2 (all bins,
                 and p_T(Z) > 250 GeV).  At LO the events carry the Born kinematics and must agree
                 within the MC errors.  At NLO the POWHEG emission moves p_T(Z) away from the fixed-
                 order value by terms beyond NLO, so agreement is expected to within a few per cent,
                 not exactly; with Born-suppressed events this is the test of their high-p_T tail.
  fractions      fraction of sigma above each STXS p_T(Z) boundary (75, 150, 250, 400, 600 GeV) for
                 fixed-order LO and NLO, parton-level LO and NLO, NLO+PS and every LO+PS start, and
                 the ratios: shower effect at each order (LO+PS / parton LO, NLO+PS / parton NLO),
                 parton level / fixed order, LO+PS / NLO+PS.
  mzh            m_ZH after / before the shower (normalised shapes, per bin) for NLO+PS and every
                 LO+PS start: the shower does not change m_ZH (Z and H keep their invariant mass up
                 to recoil), so a correction in m_ZH alone cannot remove an excess the shower makes.
Errors of ratios between different samples are combined in quadrature; LO+PS / parton LO and NLO+PS /
parton NLO use the same events, so their quoted errors are upper bounds.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

BOUNDS = [75.0, 150.0, 250.0, 400.0, 600.0]


def chi2_prob(chi2, ndf):
    """Upper tail of the chi2 distribution (series of the regularised gamma function)."""
    if ndf <= 0:
        return 1.0
    a, x = ndf / 2.0, chi2 / 2.0
    if x <= 0:
        return 1.0
    term = s = 1.0 / a
    k = 0
    while term > 1e-12 * s and k < 1000:
        k += 1
        term *= x / (a + k)
        s += term
    lower = s * math.exp(-x + a * math.log(x) - math.lgamma(a))
    return max(0.0, 1.0 - lower)


def hist(s, name):
    """-> edges, xs[fb], err[fb] with underflow first and overflow last (nominal weight)."""
    h = s["hists"][f"GGZH_NLOPS/{name}"]
    v = h["var"][""]
    return h["edges"], [1e3 * x for x in v["xs"]], [1e3 * e for e in v["err"]]


def rebin(edges, xs, err, new):
    """sum the (25 GeV) bins of a summary into the bins of new (all edges of new must be edges of the summary)"""
    out, oute = [], []
    for a, b in zip(new[:-1], new[1:]):
        idx = [i + 1 for i in range(len(edges) - 1) if edges[i] >= a - 1e-9 and edges[i + 1] <= b + 1e-9]
        if not idx or abs(edges[idx[0] - 1] - a) > 1e-9 or abs(edges[idx[-1]] - b) > 1e-9:
            raise SystemExit(f"bin {a}-{b} is not a union of the summary bins")
        out.append(sum(xs[i] for i in idx)); oute.append(math.sqrt(sum(err[i] ** 2 for i in idx)))
    return out, oute


def compare(name, new, ev, eve, fo, foe, tail_from=250.0):
    rows, chi2, ndf, chi2t, ndft = [], 0.0, 0, 0.0, 0
    for (a, b), x, ex, y, ey in zip(zip(new[:-1], new[1:]), ev, eve, fo, foe):
        e = math.hypot(ex, ey)
        p = (x - y) / e if e > 0 else 0.0
        rows.append({"bin": [a, b], "events_fb": x, "events_err_fb": ex, "fo_fb": y, "fo_err_fb": ey,
                     "ratio": x / y if y else float("nan"), "ratio_err": math.hypot(ex / y, x * ey / y ** 2) if y else float("nan"),
                     "pull": p})
        if e > 0:
            chi2 += p * p; ndf += 1
            if a >= tail_from:
                chi2t += p * p; ndft += 1
    return {"what": name, "bins": rows, "chi2": chi2, "ndf": ndf, "pvalue": chi2_prob(chi2, ndf),
            "tail_from_gev": tail_from, "chi2_tail": chi2t, "ndf_tail": ndft, "pvalue_tail": chi2_prob(chi2t, ndft)}


def frac_summary(s):
    """fraction of sigma above each boundary (incl. overflow), and its error (uncorrelated estimate)"""
    edges, xs, err = hist(s, "ptz")
    tot = sum(xs)
    out = {}
    for b in BOUNDS:
        idx = [i + 1 for i in range(len(edges) - 1) if edges[i] >= b] + [len(xs) - 1]
        above = sum(xs[i] for i in idx)
        out[str(int(b))] = [above / tot, math.sqrt(sum(err[i] ** 2 for i in idx)) / tot]
    return out


def frac_fo(fo, order):
    d, sig = fo["differential"]["ptz25"], fo[order]["sigma_fb"]
    edges, xs, err = d["edges"], d[f"{order}_fb"], d[f"{order}_err_fb"]
    out = {}
    for b in BOUNDS:
        idx = [i for i in range(len(xs)) if edges[i] >= b]
        out[str(int(b))] = [sum(xs[i] for i in idx) / sig, math.sqrt(sum(err[i] ** 2 for i in idx)) / sig]
    return out


def ratio(f, g):
    return {b: [f[b][0] / g[b][0], (f[b][0] / g[b][0]) * math.hypot(f[b][1] / f[b][0], g[b][1] / g[b][0])]
            if f[b][0] > 0 and g[b][0] > 0 else [float("nan"), float("nan")] for b in f}


def mzh_effect(after, before):
    ea, xa, _ = hist(after, "mzh")
    eb, xb, _ = hist(before, "mzh")
    if ea != eb:
        raise SystemExit("different m_ZH binning")
    ta, tb = sum(xa), sum(xb)
    rows = [{"bin": [ea[i], ea[i + 1]], "ratio": (xa[i + 1] / ta) / (xb[i + 1] / tb) if xb[i + 1] > 0 else float("nan"),
             "fraction_before": xb[i + 1] / tb} for i in range(len(ea) - 1)]
    big = [r for r in rows if r["fraction_before"] > 0.01]
    return {"bins": rows, "max_abs_deviation_bins_above_1pc": max(abs(r["ratio"] - 1) for r in big) if big else float("nan")}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fo", type=Path, required=True)
    ap.add_argument("--parton-nlo", type=Path, required=True)
    ap.add_argument("--parton-lo", type=Path, required=True)
    ap.add_argument("--nlo-ps", type=Path, required=True)
    ap.add_argument("--lo-ps", type=Path, nargs="*", default=[])
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    L = lambda p: json.loads(p.read_text())
    fo, pn, pl, nps = L(a.fo), L(a.parton_nlo), L(a.parton_lo), L(a.nlo_ps)
    d = fo["differential"]["ptz25"]
    new = d["edges"]
    res = {"fo": a.fo.name, "parton_nlo": a.parton_nlo.name, "parton_lo": a.parton_lo.name, "nlo_ps": a.nlo_ps.name,
           "lo_ps": [p.name for p in a.lo_ps], "parton_vs_fo": {}}
    for order, s in (("lo", pl), ("nlo", pn)):
        ev, eve = rebin(*hist(s, "ptz"), new)
        res["parton_vs_fo"][order] = compare(f"parton-level {order.upper()} events vs fixed-order {order.upper()}",
                                             new, ev, eve, d[f"{order}_fb"], d[f"{order}_err_fb"])
    f = {"fo_lo": frac_fo(fo, "lo"), "fo_nlo": frac_fo(fo, "nlo"), "parton_lo": frac_summary(pl),
         "parton_nlo": frac_summary(pn), "nlo_ps": frac_summary(nps)}
    lops = {p.stem: frac_summary(L(p)) for p in a.lo_ps}
    res["fractions"] = {**f, "lo_ps": lops}
    res["ratios"] = {"parton_lo_over_fo_lo": ratio(f["parton_lo"], f["fo_lo"]),
                     "parton_nlo_over_fo_nlo": ratio(f["parton_nlo"], f["fo_nlo"]),
                     "fo_lo_over_fo_nlo": ratio(f["fo_lo"], f["fo_nlo"]),
                     "nlo_ps_over_parton_nlo": ratio(f["nlo_ps"], f["parton_nlo"]),
                     "lo_ps_over_parton_lo": {k: ratio(v, f["parton_lo"]) for k, v in lops.items()},
                     "lo_ps_over_nlo_ps": {k: ratio(v, f["nlo_ps"]) for k, v in lops.items()}}
    res["mzh"] = {"nlo_ps_over_parton_nlo": mzh_effect(nps, pn),
                  "lo_ps_over_parton_lo": {p.stem: mzh_effect(L(p), pl) for p in a.lo_ps}}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1) + "\n")

    for order in ("lo", "nlo"):
        r = res["parton_vs_fo"][order]
        print(f"{r['what']}: chi2 {r['chi2']:.1f}/{r['ndf']} (p = {r['pvalue']:.3g}); "
              f"pT(Z) > {r['tail_from_gev']:.0f} GeV: {r['chi2_tail']:.1f}/{r['ndf_tail']} (p = {r['pvalue_tail']:.3g})")
        print("   ratio per bin above 250 GeV: " + ", ".join(
            f"{b['bin'][0]:.0f}-{b['bin'][1]:.0f}: {b['ratio']:.3f}+-{b['ratio_err']:.3f}" for b in r["bins"] if b["bin"][0] >= 250))
    print("fraction above pT(Z) boundary      " + "  ".join(f"{int(b):>13d}" for b in BOUNDS))
    show = lambda n, v: print(f"{n:34s} " + "  ".join(f"{v[str(int(b))][0]:.5f}({1e5 * v[str(int(b))][1]:3.0f})" for b in BOUNDS))
    for k in ("fo_lo", "parton_lo", "fo_nlo", "parton_nlo", "nlo_ps"):
        show(k, f[k])
    for k, v in lops.items():
        show(k, v)
    print("shower effect                      " + "  ".join(f"{int(b):>13d}" for b in BOUNDS))
    show("NLO+PS / parton NLO", res["ratios"]["nlo_ps_over_parton_nlo"])
    for k, v in res["ratios"]["lo_ps_over_parton_lo"].items():
        show(f"{k} / parton LO", v)
    print("m_ZH after/before the shower, largest deviation in bins with > 1 % of sigma: NLO+PS "
          f"{res['mzh']['nlo_ps_over_parton_nlo']['max_abs_deviation_bins_above_1pc']:.3f}; " + ", ".join(
              f"{k} {v['max_abs_deviation_bins_above_1pc']:.3f}" for k, v in res["mzh"]["lo_ps_over_parton_lo"].items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
