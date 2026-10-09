#!/usr/bin/env python3
"""Compare a reproduction with reference/generator.json (and its parts with each other).

    check.py [--fo fo.json] [--stable-nnpdf31 st.json] [--stable-pdf4lhc21 st.json]
             [--powheg-lo powheg-lo.json] [--grids grids-xsec.json] [--events events.json]
             [--events-lo events-lo.json]
             [--cache compare-*.json ...] [--mcfm-stable-nnpdf31 m.json] [--mcfm-stable-pdf4lhc21 m.json]
             [--out report.json]

Inputs are the outputs of analysis/fo_combine.py, powheg_stat.py, lhe_xsec.py and
compare_lhe.py (and mcfm_xsec.py); give whichever exist.  Criteria:
  * a cross section (value and MC error) agrees with its reference if the pull
    |x - ref| / sqrt(err^2 + err_ref^2) is below 3, at every scale point;
  * differential K factors: chi2 over the bins, p-value above 0.001;
  * PDF and alpha_s uncertainties (from the event weights): within 0.3 and 0.5 percentage points;
  * cache closure: every compared file identical.
A rebuild with another compiler or library versions is expected to agree within these
statistical criteria, not bit by bit.  Exit status 1 if anything fails.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

REF = json.loads((Path(__file__).resolve().parent / "reference" / "generator.json").read_text())
ROWS: list[dict] = []


def pull(name, x, ex, r, er, what="vs reference"):
    p = (x - r) / math.hypot(ex, er) if math.hypot(ex, er) > 0 else float("inf")
    ROWS.append({"check": name, "what": what, "value": x, "err": ex, "ref": r, "ref_err": er, "pull": p, "pass": abs(p) < 3})


def within(name, x, r, tol, what="vs reference"):
    ROWS.append({"check": name, "what": what, "value": x, "ref": r, "tol": tol, "pass": abs(x - r) <= tol})


def flag(name, ok, what):
    ROWS.append({"check": name, "what": what, "pass": bool(ok)})


def points(name, got, ref):
    key = lambda p: (float(p["muR"]), float(p["muF"]))
    g = {key(p): p for p in got}
    for p in ref:
        if key(p) in g:
            q = g[key(p)]
            pull(f"{name} muR={p['muR']:g} muF={p['muF']:g}", q["sigma_fb"], q["mc_err_fb"], p["sigma_fb"], p["mc_err_fb"])
        else:
            flag(f"{name} muR={p['muR']:g} muF={p['muF']:g}", False, "scale point missing")


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


def kfactor(name, got, ref, skip_first=True):
    chi2, ndf = 0.0, 0
    for i, (k, e, rk, re_) in enumerate(zip(got["K"], got["K_err"], ref["K"], ref["K_err"])):
        if skip_first and i == 0:      # 150-200 GeV, below threshold: tiny LO, K undetermined
            continue
        chi2 += (k - rk) ** 2 / (e ** 2 + re_ ** 2); ndf += 1
    p = chi2_prob(chi2, ndf)
    ROWS.append({"check": name, "what": f"chi2 {chi2:.1f} / {ndf}", "pvalue": p, "pass": p > 1e-3})


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for opt in ("fo", "stable-nnpdf31", "stable-pdf4lhc21", "powheg-lo", "grids", "events", "events-lo",
                "mcfm-stable-nnpdf31", "mcfm-stable-pdf4lhc21"):
        ap.add_argument(f"--{opt}", type=Path)
    ap.add_argument("--cache", nargs="*", type=Path, default=[])
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    load = lambda p: json.loads(p.read_text())

    if a.fo:
        fo, r = load(a.fo), REF["fo_ggllh"]
        for o in ("lo", "nlo"):
            if o in fo:
                points(f"ggxy {o.upper()}", fo[o]["points"], r[o]["points"])
        if "K" in fo:
            pull("ggxy K = NLO/LO", fo["K"], fo["K_err"], r["K"], r["K_err"])
            within("ggxy K vs arXiv:2603.15762", fo["K"], REF["paper"]["K"], 0.01, "PDF4LHC21 in the paper")
        d = fo.get("differential", {})
        if "mzh25" in d:
            kfactor("ggxy K(m_ZH)", d["mzh25"], r["K_mzh"])
        if "ptz25" in d:
            kfactor("ggxy K(pT(Z))", d["ptz25"], r["K_ptz"], skip_first=False)
    for opt, key in (("stable_nnpdf31", "NNPDF31"), ("stable_pdf4lhc21", "PDF4LHC21")):
        p = getattr(a, opt)
        if p:
            st, r = load(p)["lo"], REF["fo_ggzh_stable_lo"][key]
            pull(f"ggxy stable-Z LO {key}", st["sigma_fb"], st["mc_err_fb"], r["sigma_fb"], r["mc_err_fb"])
            if key == "PDF4LHC21":
                pull("ggxy stable-Z LO PDF4LHC21", st["sigma_fb"], st["mc_err_fb"], REF["paper"]["stable_lo"]["sigma_fb"],
                     REF["paper"]["stable_lo"]["mc_err_fb"], "vs arXiv:2603.15762 Table 1")
    if a.powheg_lo:
        pl = load(a.powheg_lo)
        got = [{"muR": float(n.split("f")[0][1:].replace("p", ".")), "muF": float(n.split("f")[1].replace("p", ".")),
                "sigma_fb": p["sigma_fb"], "mc_err_fb": p["mc_err_fb"]} for n, p in pl["points"].items()]
        points("POWHEG LO", got, REF["powheg_lo"]["points"])
        if a.fo and "lo" in load(a.fo):
            lo = load(a.fo)["lo"]
            pull("POWHEG LO vs ggxy LO", pl["sigma_fb"], pl["mc_err_fb"], lo["sigma_fb"], lo["mc_err_fb"], "two codes, same run")
    if a.grids:
        g, r = load(a.grids), REF["grids_nlo"]
        pull("grids NLO integral", g["sigma_fb"], g["mc_err_fb"], r["sigma_fb"], r["mc_err_fb"])
        if a.fo and "nlo" in load(a.fo):
            n = load(a.fo)["nlo"]
            pull("grids NLO vs ggxy NLO", g["sigma_fb"], g["mc_err_fb"], n["sigma_fb"], n["mc_err_fb"], "two codes, same run")
    if a.events:
        ev, r = load(a.events), REF["events_nlo"]
        if "scale_points" in ev:
            points("events NLO", ev["scale_points"], r["scale_points"])
        else:
            pull("events NLO", ev["sigma_fb"], ev["mc_err_fb"], r["sigma_fb"], r["mc_err_fb"])
        if "pdf_hessian_rel" in ev:
            within("events PDF (Hessian)", ev["pdf_hessian_rel"], r["pdf_hessian_rel"], 0.003)
        if "alphas_rel_for_0p002" in ev:
            within("events alpha_s (+-0.002)", ev["alphas_rel_for_0p002"], r["alphas_rel_for_0p002"], 0.005)
        if a.fo and "nlo" in load(a.fo):
            n = load(a.fo)["nlo"]
            pull("events NLO vs ggxy NLO", ev["sigma_fb"], ev["mc_err_fb"], n["sigma_fb"], n["mc_err_fb"], "two codes, same run")
    if a.events_lo:
        # LO (Born-only) events from LO grids: the event weights against the ggxy LO at all 7 points
        points("LO events vs ggxy LO", load(a.events_lo)["scale_points"], REF["fo_ggllh"]["lo"]["points"])
    for opt, key in (("mcfm_stable_nnpdf31", "NNPDF31"), ("mcfm_stable_pdf4lhc21", "PDF4LHC21")):
        p = getattr(a, opt)
        if p:
            m, r = load(p), REF["mcfm_lo_stable"][key]
            pull(f"MCFM stable-Z LO {key}", m["sigma_fb"], m["mc_err_fb"], r["sigma_fb"], r["mc_err_fb"])
    for c in a.cache:
        d = load(c)
        flag(f"cache closure {c.name}", d.get("identical") and d.get("compared", 0) > 0,
             f"{d.get('kinematics_identical')}/{d.get('compared')} events, {d.get('weights_identical')}/{d.get('weights')} weights")

    if not ROWS:
        raise SystemExit("nothing to check: give at least one result file")
    for r in ROWS:
        s = "PASS" if r["pass"] else "FAIL"
        if "pull" in r:
            print(f"{s}  {r['check']:40s} {r['value']:10.4f} +- {r['err']:.4f}  ref {r['ref']:10.4f} +- {r['ref_err']:.4f}  pull {r['pull']:+5.2f}  ({r['what']})")
        elif "tol" in r:
            print(f"{s}  {r['check']:40s} {r['value']:10.4f}  ref {r['ref']:10.4f}  tolerance {r['tol']}  ({r['what']})")
        else:
            print(f"{s}  {r['check']:40s} {r['what']}" + (f"  p = {r['pvalue']:.3g}" if "pvalue" in r else ""))
    nfail = sum(not r["pass"] for r in ROWS)
    print(f"{len(ROWS) - nfail}/{len(ROWS)} checks pass")
    if a.out:
        a.out.write_text(json.dumps(ROWS, indent=1) + "\n")
    return 1 if nfail else 0


if __name__ == "__main__":
    raise SystemExit(main())
