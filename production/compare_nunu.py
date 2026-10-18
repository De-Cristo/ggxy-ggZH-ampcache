#!/usr/bin/env python3
"""Z -> nu nubar against Z -> l+l-: the results of the vv sets of production/paper.sh next to those of the
13600 / 13000 sets, and the published values (arXiv:2603.15762, Table 3).  JSON in, JSON out; seconds.

    production/compare_nunu.py --ll-13600 <analysis 13600> --vv-13600 <analysis vv13600>
                               [--ll-13000 <analysis 13000> --vv-13000 <analysis vv13000>] --out nunu-vs-ll.json

<analysis X> is the output directory of production/analyze.sh X.  Per energy (fixed order, ggxy):
sigma_LO, sigma_NLO with the seven-point band for nu nubar H and l+l-H, the deviation from Table 3, the
ratio nu nubar / l+l- against the published ratio, K for both and the pull of the difference, and the
differential K (m_ZH, pT(Z), pT(H), pT(ZH)) as nu nubar / l+l- with a chi2.  At 13.6 TeV (NLO+PS): the
event-weight cross sections against the fixed order, and the STXS Stage 1.3 fine bins (FWDH included):
the fraction of sigma per bin for both channels, nu nubar / l+l-, the pull per bin and a chi2 over the bins
filled in both (n - 1 degrees of freedom: the fractions sum to one), plus the relative scale band per bin.
Expected (checked with an independent CMSSW production of the same process): equal shapes and STXS
fractions within the MC errors; K for nu nubar about 0.15 % below l+l-, because the l+l- NLO includes
gamma* -> l+l- (enhanced at small m_ll), absent for neutrinos.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

# arXiv:2603.15762, Table 3 (PDF4LHC21, mu0 = m_ZH/2, l = e, mu, tau, 15 < m_ll < 150 GeV): sigma [fb] (error), band [%]
TABLE3 = {
    "13000": {"ll": {"lo": (6.366, 0.004, 26.7, -20.0), "nlo": (11.80, 0.01, 16.6, -14.0)},
              "vv": {"lo": (12.581, 0.008, 26.7, -20.0), "nlo": (23.30, 0.02, 16.5, -14.0)}},
    "13600": {"ll": {"lo": (7.026, 0.004, 26.2, -19.7), "nlo": (13.02, 0.01, 16.4, -13.8)},
              "vv": {"lo": (13.883, 0.008, 26.2, -19.7), "nlo": (25.70, 0.02, 16.4, -13.8)}},
}


def chi2_sf(x: float, k: int) -> float:
    """Upper tail of the chi2 distribution."""
    a, z = k / 2.0, x / 2.0
    if k <= 0:
        return float("nan")
    if z <= 0:
        return 1.0
    if z < a + 1:
        s = t = 1.0 / a; n = a
        while t > 1e-15 * s:
            n += 1; t *= z / n; s += t
        return max(0.0, 1.0 - s * math.exp(-z + a * math.log(z) - math.lgamma(a)))
    b = z + 1 - a; c = 1e300; d = 1.0 / b; h = d
    for i in range(1, 500):
        an = -i * (i - a); b += 2; d = an * d + b; d = 1.0 / d if d else 1e300
        c = b + an / c; h *= d * c
        if abs(d * c - 1) < 1e-15:
            break
    return math.exp(-z + a * math.log(z) - math.lgamma(a)) * h


def load(d: Path, name: str):
    p = d / name
    return json.loads(p.read_text()) if p.exists() else None


def band(order: dict) -> list:
    c = order["sigma_fb"]; v = [p["sigma_fb"] for p in order["points"]]
    return [100 * (max(v) / c - 1), 100 * (min(v) / c - 1)]


def fixed_order(ll: dict, vv: dict, energy: str) -> dict:
    out = {}
    for ch, d in (("ll", ll), ("vv", vv)):
        out[ch] = {}
        for m in ("lo", "nlo"):
            ref = TABLE3[energy][ch][m]; o = d[m]
            out[ch][m] = {"sigma_fb": o["sigma_fb"], "mc_err_fb": o["mc_err_fb"], "band_pct": band(o),
                          "table3_fb": ref[0], "table3_err_fb": ref[1], "table3_band_pct": [ref[2], ref[3]],
                          "deviation": o["sigma_fb"] / ref[0] - 1,
                          "band_diff_points": [band(o)[0] - ref[2], band(o)[1] - ref[3]]}
        out[ch]["K"] = [d["K"], d["K_err"]]
    for m in ("lo", "nlo"):
        r = vv[m]["sigma_fb"] / ll[m]["sigma_fb"]
        er = r * math.hypot(vv[m]["mc_err_fb"] / vv[m]["sigma_fb"], ll[m]["mc_err_fb"] / ll[m]["sigma_fb"])
        rr = TABLE3[energy]["vv"][m][0] / TABLE3[energy]["ll"][m][0]
        out[f"ratio_vv_ll_{m}"] = {"ours": [r, er], "table3": rr}
    dk = vv["K"] - ll["K"]; sk = math.hypot(vv["K_err"], ll["K_err"])
    out["K_vv_minus_K_ll"] = {"difference": dk, "sigma": sk, "pull": dk / sk, "relative": vv["K"] / ll["K"] - 1,
                              "table3_relative": (TABLE3[energy]["vv"]["nlo"][0] / TABLE3[energy]["vv"]["lo"][0])
                              / (TABLE3[energy]["ll"]["nlo"][0] / TABLE3[energy]["ll"]["lo"][0]) - 1}
    diff = {}
    for obs in ll.get("differential", {}):
        a, b = vv["differential"].get(obs), ll["differential"][obs]
        if not a or a["edges"] != b["edges"]:
            continue
        rows, chi2, n = [], 0.0, 0
        for i in range(len(b["K"])):
            ka, kb, ea, eb = a["K"][i], b["K"][i], a["K_err"][i], b["K_err"][i]
            if ka is None or kb is None or not (ea and eb and math.isfinite(ka) and math.isfinite(kb)):
                rows.append(None); continue
            pull = (ka - kb) / math.hypot(ea, eb); chi2 += pull * pull; n += 1
            rows.append({"K_vv": ka, "K_ll": kb, "ratio": ka / kb if kb else None, "pull": pull})
        diff[obs] = {"edges": b["edges"], "bins": rows, "chi2": chi2, "ndf": n, "p_value": chi2_sf(chi2, n)}
    out["differential_K"] = diff
    return out


def stxs(ll: dict, vv: dict) -> dict:
    cl = {c["cell"]: c for c in ll["cells"]}; cv = {c["cell"]: c for c in vv["cells"]}
    tl = sum(c["xs_fb"] for c in ll["cells"]); tv = sum(c["xs_fb"] for c in vv["cells"])
    rows, chi2, n = {}, 0.0, 0
    for cell in [c["cell"] for c in ll["cells"]]:
        a, b = cv[cell], cl[cell]
        fa, fb = a["xs_fb"] / tv, b["xs_fb"] / tl
        ea, eb = a["mc_err_fb"] / tv, b["mc_err_fb"] / tl
        row = {"f_vv": fa, "f_ll": fb, "ratio": fa / fb if fb else None,
               "scale_band_vv": a["scale"], "scale_band_ll": b["scale"],
               "pdf_alphas_vv": a.get("pdf_alphas"), "pdf_alphas_ll": b.get("pdf_alphas")}
        if ea > 0 and eb > 0:
            row["pull"] = (fa - fb) / math.hypot(ea, eb); chi2 += row["pull"] ** 2; n += 1
        rows[cell] = row
    ndf = n - 1
    return {"sigma_fb": {"vv": vv["sigma_fb"], "ll": ll["sigma_fb"]}, "cells": rows, "chi2": chi2, "ndf": ndf,
            "p_value": chi2_sf(chi2, ndf), "max_abs_pull": max(abs(r["pull"]) for r in rows.values() if "pull" in r)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for e in ("13600", "13000"):
        ap.add_argument(f"--ll-{e}", type=Path); ap.add_argument(f"--vv-{e}", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    res = {"reference": "arXiv:2603.15762, Table 3", "energies": {}}
    for e in ("13600", "13000"):
        ld, vd = getattr(a, f"ll_{e}"), getattr(a, f"vv_{e}")
        if not (ld and vd):
            continue
        r = {}
        fl, fv = load(ld, "fo-llh.json"), load(vd, "fo-vvh.json")
        if fl and fv:
            r["fixed_order"] = fixed_order(fl, fv, e)
        el, ev = load(ld, "events-nlo.json"), load(vd, "events-nlo.json")
        if ev:
            r["events_nlo"] = {"vv": {k: ev[k] for k in ("sigma_fb", "mc_err_fb", "scale7_band_rel", "pdf_hessian_rel",
                                                       "alphas_rel_half_difference", "negative_weight_fraction", "events")}}
            if el:
                r["events_nlo"]["ll"] = {k: el[k] for k in ("sigma_fb", "mc_err_fb", "scale7_band_rel", "pdf_hessian_rel",
                                                           "alphas_rel_half_difference")}
            if fv:
                f = fv["nlo"]["sigma_fb"]; fe = fv["nlo"]["mc_err_fb"]
                r["events_nlo"]["vv_vs_fixed_order"] = {"deviation": ev["sigma_fb"] / f - 1,
                                                        "pull": (ev["sigma_fb"] - f) / math.hypot(ev["mc_err_fb"], fe)}
        g = load(vd, "grids-nlo.json")
        if g:
            r["grids_nlo_vv"] = g
        sl, sv = load(ld, "stxs-nlo.json"), load(vd, "stxs-nlo.json")
        if sl and sv:
            r["stxs_nlo_ps"] = stxs(sl, sv)
        cs = load(vd, "check-shower-nlo-monash.json")
        if cs:
            r["check_shower_vv"] = cs
        res["energies"][e] = r
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    for e, r in res["energies"].items():
        fo = r.get("fixed_order")
        if fo:
            for m in ("lo", "nlo"):
                v = fo["vv"][m]
                print(f"{e} {m:3s} nunu {v['sigma_fb']:.4f}({1e4 * v['mc_err_fb']:.0f}) fb  band {v['band_pct'][0]:+.2f}/{v['band_pct'][1]:+.2f} %"
                      f"  Table 3 {v['table3_fb']}  dev {100 * v['deviation']:+.2f} %  nunu/ll {fo[f'ratio_vv_ll_{m}']['ours'][0]:.4f} (Table 3 {fo[f'ratio_vv_ll_{m}']['table3']:.4f})")
            k = fo["K_vv_minus_K_ll"]
            print(f"{e} K nunu {fo['vv']['K'][0]:.4f} ll {fo['ll']['K'][0]:.4f}: {100 * k['relative']:+.3f} % (pull {k['pull']:+.2f}; Table 3 {100 * k['table3_relative']:+.3f} %)")
            for obs, d in fo["differential_K"].items():
                print(f"{e} differential K {obs}: chi2 {d['chi2']:.1f}/{d['ndf']} p {d['p_value']:.2f}")
        if "events_nlo" in r and "vv_vs_fixed_order" in r["events_nlo"]:
            x = r["events_nlo"]
            print(f"{e} NLO+PS events nunu {x['vv']['sigma_fb']:.3f} fb vs fixed order {100 * x['vv_vs_fixed_order']['deviation']:+.2f} % (pull {x['vv_vs_fixed_order']['pull']:+.2f})")
        if "stxs_nlo_ps" in r:
            s = r["stxs_nlo_ps"]
            print(f"{e} STXS fine fractions nunu vs ll: chi2 {s['chi2']:.1f}/{s['ndf']} p {s['p_value']:.2f}, max |pull| {s['max_abs_pull']:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
