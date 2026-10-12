#!/usr/bin/env python3
"""Two independent samples of the same prediction after the shower, bin by bin.

    compare_samples.py --a a.json --b b.json [--out cmp.json] [--report-only]

Inputs: finished summaries of analysis/yoda_summary.py (nominal weights).  Used for the
Born-suppression test: NLO+PS from weighted events (bornsuppfact) against NLO+PS from unweighted
events of the same generator settings.  Compared: the inclusive cross section (pull), the 19
STXS Stage 1.3 fine bins (FWDH + 18) and the p_T(Z) histogram including its overflow (chi2 over
the bins with an entry in either sample).  The samples agree if the pull is below 3 and both
chi2 p-values are above 0.001.  Exit status 1 otherwise, unless --report-only (for two samples that
are expected to differ, e.g. two matching settings: the per-bin values and ratios are the result).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


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


def bins(xa, ea, xb, eb):
    rows, chi2, ndf = [], 0.0, 0
    for i, (x, ex, y, ey) in enumerate(zip(xa, ea, xb, eb)):
        e = math.hypot(ex, ey)
        p = (x - y) / e if e > 0 else 0.0
        rows.append({"bin": i, "a_fb": 1e3 * x, "a_err_fb": 1e3 * ex, "b_fb": 1e3 * y, "b_err_fb": 1e3 * ey, "pull": p})
        if e > 0:
            chi2 += p * p; ndf += 1
    return {"bins": rows, "chi2": chi2, "ndf": ndf, "pvalue": chi2_prob(chi2, ndf)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", type=Path, required=True)
    ap.add_argument("--b", type=Path, required=True)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--report-only", action="store_true")
    a = ap.parse_args(argv)
    sa, sb = json.loads(a.a.read_text()), json.loads(a.b.read_text())
    res = {"a": str(a.a), "b": str(a.b)}

    sig = lambda s, k: sum(s["hists"]["GGZH_NLOPS/sigma"]["var"][""][k])
    sig_err = lambda s: math.sqrt(sum(e * e for e in s["hists"]["GGZH_NLOPS/sigma"]["var"][""]["err"]))
    x, ex, y, ey = sig(sa, "xs"), sig_err(sa), sig(sb, "xs"), sig_err(sb)
    pull = (x - y) / math.hypot(ex, ey)
    res["sigma"] = {"a_fb": 1e3 * x, "a_err_fb": 1e3 * ex, "b_fb": 1e3 * y, "b_err_fb": 1e3 * ey, "pull": pull}
    fa, fb = sa["stxs13_fine"][""], sb["stxs13_fine"][""]
    res["stxs13_fine"] = bins(fa["xs"], fa["err"], fb["xs"], fb["err"])
    ha, hb = sa["hists"]["GGZH_NLOPS/ptz"], sb["hists"]["GGZH_NLOPS/ptz"]
    res["ptz"] = bins(ha["var"][""]["xs"], ha["var"][""]["err"], hb["var"][""]["xs"], hb["var"][""]["err"])
    res["ptz"]["edges"] = ha["edges"]
    res["pass"] = abs(pull) < 3 and res["stxs13_fine"]["pvalue"] > 1e-3 and res["ptz"]["pvalue"] > 1e-3

    print(f"sigma  a {1e3 * x:.4f} +- {1e3 * ex:.4f} fb   b {1e3 * y:.4f} +- {1e3 * ey:.4f} fb   pull {pull:+.2f}")
    for k in ("stxs13_fine", "ptz"):
        r = res[k]
        print(f"{k:12s} chi2 {r['chi2']:.1f} / {r['ndf']}   p = {r['pvalue']:.3g}")
    print("PASS" if res["pass"] else "FAIL")
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(res, indent=1) + "\n")
    return 0 if res["pass"] or a.report_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
