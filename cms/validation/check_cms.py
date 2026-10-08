#!/usr/bin/env python3
"""Compare NanoGEN validation summaries (nanogen_summary.py) with reference/cms.json.

    check_cms.py [--nlo summary-nlo.json] [--lo summary-lo.json] [--out report.json]

Criteria (statistical; a new gridpack and new events agree within errors, not bit by bit):
  * inclusive cross section at each of the nine scale points: pull below 3;
  * negative-weight and |y_H| > 2.5 (FWDH) fractions: binomial pull below 3;
  * HTXS category 0 (unclassified) fraction: exactly 0 (ProductionMode must be GG2ZH);
  * the 18 Stage 1.3 cells and every distribution (incl. under/overflow): chi2 against the
    reference, p-value above 0.001; the bin variance is the reference's, scaled to the sample
    size (err^2 (1 + N_ref/N)), so that empty bins of a small sample are handled.
Exit status 1 if anything fails.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

REF = json.loads((Path(__file__).resolve().parents[1] / "reference" / "cms.json").read_text())
ROWS: list[dict] = []


def chi2_prob(chi2, ndf):
    if ndf <= 0 or chi2 <= 0:
        return 1.0
    a, x = ndf / 2.0, chi2 / 2.0
    term = s = 1.0 / a
    k = 0
    while term > 1e-12 * s and k < 2000:
        k += 1
        term *= x / (a + k)
        s += term
    return max(0.0, 1.0 - s * math.exp(-x + a * math.log(x) - math.lgamma(a)))


def pull(name, x, ex, r, er):
    p = (x - r) / math.hypot(ex, er) if math.hypot(ex, er) > 0 else (0.0 if x == r else float("inf"))
    ROWS.append({"check": name, "value": x, "err": ex, "ref": r, "ref_err": er, "pull": p, "pass": abs(p) < 3})


def binom(name, f, n, fr, nr):
    e = math.sqrt(max(f * (1 - f), 1e-12) / n); er = math.sqrt(max(fr * (1 - fr), 1e-12) / nr)
    pull(name, f, e, fr, er)


def chi2(name, got, gerr, ref, rerr, n, nref):
    """Bin variances from the reference, scaled to the sample size (err^2 ~ 1/N for the same
    distribution): sparse bins of a small sample often have no entries, so their own error is 0."""
    c = 0.0; ndf = 0
    for g, ge, r, re_ in zip(got, gerr, ref, rerr):
        v = re_ ** 2 * (1.0 + nref / n) if re_ > 0 else ge ** 2
        if v > 0:
            c += (g - r) ** 2 / v; ndf += 1
    p = chi2_prob(c, ndf)
    ROWS.append({"check": name, "chi2": c, "ndf": ndf, "pvalue": p, "pass": p > 1e-3})


def compare(label, s):
    r = REF[label]
    t, rt = s["totals"], r["totals"]
    key = lambda p: (float(p["muR"]), float(p["muF"]))
    got = {key(p): p for p in t["scale_points"]}
    for p in rt["scale_points"]:
        q = got[key(p)]
        pull(f"{label} sigma muR={p['muR']:g} muF={p['muF']:g}", q["sigma_fb"], q["mc_err_fb"], p["sigma_fb"], p["mc_err_fb"])
    binom(f"{label} negative-weight fraction", t["negative_weight_fraction"], s["events"], rt["negative_weight_fraction"], r["events"])
    st, rs = s["stxs"], r["stxs"]
    ROWS.append({"check": f"{label} HTXS category 0 fraction", "value": st["category0_fraction"], "ref": 0.0,
                 "pass": st["category0_fraction"] == 0.0})
    binom(f"{label} FWDH (|y_H| > 2.5) fraction", st["fwdh_fraction"], s["events"], rs["fwdh_fraction"], r["events"])
    chi2(f"{label} STXS 1.3 cells", st["sigma_fb"], st["mc_err_fb"], rs["sigma_fb"], rs["mc_err_fb"], s["events"], r["events"])
    for h, rh in r["hists"].items():
        gh = s["hists"][h]
        assert gh["edges"] == rh["edges"], h
        chi2(f"{label} {h}", gh["sigma_fb"], gh["mc_err_fb"], rh["sigma_fb"], rh["mc_err_fb"], s["events"], r["events"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nlo", type=Path)
    ap.add_argument("--lo", type=Path)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    for label, p in (("nlo", a.nlo), ("lo", a.lo)):
        if p:
            compare(label, json.loads(p.read_text()))
    if a.nlo:
        t = json.loads(a.nlo.read_text())["totals"]; g = REF["gridpack"]
        pull("nlo sigma vs reference grid integral", t["sigma_fb"], t["mc_err_fb"], g["sigma_fb"], g["mc_err_fb"])
    if not ROWS:
        raise SystemExit("give --nlo and/or --lo")
    for r in ROWS:
        s = "PASS" if r["pass"] else "FAIL"
        if "pull" in r:
            print(f"{s}  {r['check']:44s} {r['value']:10.5g} +- {r['err']:.3g}  ref {r['ref']:10.5g} +- {r['ref_err']:.3g}  pull {r['pull']:+5.2f}")
        elif "chi2" in r:
            print(f"{s}  {r['check']:44s} chi2 {r['chi2']:7.1f} / {r['ndf']:3d}  p = {r['pvalue']:.3g}")
        else:
            print(f"{s}  {r['check']:44s} {r['value']}  (required {r['ref']})")
    nfail = sum(not r["pass"] for r in ROWS)
    print(f"{len(ROWS) - nfail}/{len(ROWS)} checks pass")
    if a.out:
        a.out.write_text(json.dumps(ROWS, indent=1) + "\n")
    return 1 if nfail else 0


if __name__ == "__main__":
    raise SystemExit(main())
