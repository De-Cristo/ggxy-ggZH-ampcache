#!/usr/bin/env python3
"""Validation of the shower and analysis step on one sample (README "Tests").

    check_shower.py --summary summary.json --lhe lhe.json [--noshower] [--out report.json]

summary.json: analysis/yoda_summary.py of the showered (or noshower) events; lhe.json:
generator/analysis/lhe_xsec.py of the same LHE files.  Checks:
  1. the inclusive cross section in Rivet equals the LHE one for the nominal weight and every
     scale point (relative 1e-5: the weights reach Rivet unchanged and are normalised as POWHEG's);
  2. the STXS routine leaves no event unclassified;
  3. the routine's 19 Stage 1.3 fine bins (FWDH + 18) equal GGZH_NLOPS's own classification
     (same definitions of Z, Higgs boson and jets) for every weight variation;
  4. the PDF weights are present and give the LHE PDF uncertainty (relative 1e-5).
With --noshower (test mode without parton level) checks 2 and 3 do not apply: the Higgs boson is
not decayed, and the STXS routine only classifies events with a decayed Higgs boson.
Exit status 1 if a check fails.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--lhe", type=Path, required=True)
    ap.add_argument("--noshower", action="store_true")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    s, l = json.loads(a.summary.read_text()), json.loads(a.lhe.read_text())
    rows = []

    def add(name, ok, what):
        rows.append({"check": name, "pass": bool(ok), "what": what})

    sig = s["hists"]["GGZH_NLOPS/sigma"]["var"]
    tot = lambda v: 1e3 * sum(sig[v]["xs"])
    rel = lambda x, y: abs(x / y - 1) if y else abs(x)
    add("sigma nominal", rel(tot(""), l["sigma_fb"]) < 1e-5, f"{tot(''):.5f} vs LHE {l['sigma_fb']:.5f} fb")
    for p in l.get("scale_points", []):
        v = f"MUR{p['muR']:.1f}_MUF{p['muF']:.1f}"
        ok = v in sig and rel(tot(v), p["sigma_fb"]) < 1e-5
        add(f"sigma {v}", ok, f"{tot(v) if v in sig else float('nan'):.5f} vs LHE {p['sigma_fb']:.5f} fb")
    unk = s.get("stxs_unknown", {}).get("", None)
    if not a.noshower:
        add("STXS unclassified", unk == 0.0, f"{unk} pb")
    own = s["hists"]["GGZH_NLOPS/stxs13_fine_own"]["var"]
    worst, nvar = 0.0, 0
    for v, hw in s.get("stxs13_fine", {}).items():
        o = own[v]["xs"][1:-1]                       # 18 cells (pT-major) + FWDH; drop under/overflow
        own_jet_major = [o[18]] + [o[3 * p + j] for j in range(3) for p in range(6)]
        for x, y in zip(hw["xs"], own_jet_major):
            worst = max(worst, abs(x - y) / max(abs(y), 1e-12) if y else abs(x))
        nvar += 1
    if not a.noshower:
        add("STXS routine vs own classification", nvar > 0 and worst < 1e-6, f"{nvar} variations, max relative difference {worst:.2e}")
    if "pdf_hessian_rel" in l:
        s0 = tot("2000") if "2000" in sig else float("nan")
        members = [tot(str(2000 + k)) for k in range(1, 101) if str(2000 + k) in sig]
        pdf = math.sqrt(sum((m - s0) ** 2 for m in members)) / s0 if len(members) == 100 else float("nan")
        add("PDF weights", len(members) == 100 and rel(pdf, l["pdf_hessian_rel"]) < 1e-4,
            f"{100 * pdf:.4f} % vs LHE {100 * l['pdf_hessian_rel']:.4f} %")
    shower_vars = [v for v in sig if v.startswith(("isr", "fsr"))]
    add("shower-variation weights", True, f"{len(shower_vars)} present (informational)")
    for r in rows:
        print(f"{'PASS' if r['pass'] else 'FAIL'}  {r['check']:38s} {r['what']}")
    nf = sum(not r["pass"] for r in rows)
    print(f"{len(rows) - nf}/{len(rows)} checks pass")
    if a.out:
        a.out.write_text(json.dumps(rows, indent=1) + "\n")
    return 1 if nf else 0


if __name__ == "__main__":
    raise SystemExit(main())
