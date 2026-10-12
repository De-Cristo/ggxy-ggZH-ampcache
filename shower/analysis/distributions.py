#!/usr/bin/env python3
"""Differential distributions of a showered sample with the weight variations needed for bands.

    distributions.py --summary shower-nlo-monash.json [--obs ptz mzh ptzh njets]
                     [--pdf-weight 3000 --pdf-nmem 40] --out dists-nlo.json

Input: a finished summary of analysis/yoda_summary.py.  Per observable (GGZH_NLOPS histograms), in fb
per bin with underflow first and overflow last, as in the summary: the bin edges, the nominal values
and MC errors, and the variations
  scale    the seven (mu_R, mu_F) points (MUR1.0_MUF1.0 first)
  pdf      the central member and the Hessian members (weight ids --pdf-weight ... + --pdf-nmem)
  alphas   the two alpha_s members after the Hessian ones
  isr, fsr the shower renormalisation-scale variations by 1/2 and 2
so that a plot can rebin first and then build the bands (envelope, Hessian sum, quadrature).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

SEVEN = ["MUR1.0_MUF1.0", "MUR0.5_MUF0.5", "MUR2.0_MUF2.0", "MUR2.0_MUF1.0", "MUR0.5_MUF1.0", "MUR1.0_MUF2.0", "MUR1.0_MUF0.5"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--obs", nargs="+", default=["ptz", "mzh", "ptzh", "njets"])
    ap.add_argument("--pdf-weight", type=int, default=3000)
    ap.add_argument("--pdf-nmem", type=int, default=40)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    s = json.loads(a.summary.read_text())
    fb = lambda v: [1e3 * x for x in v]
    out = {"summary": a.summary.name, "unit": "fb per bin, underflow first and overflow last", "obs": {}}
    for o in a.obs:
        h = s["hists"][f"GGZH_NLOPS/{o}"]
        v = h["var"]
        get = lambda name: fb(v[name]["xs"])
        pdf = [str(a.pdf_weight + k) for k in range(a.pdf_nmem + 1)]
        als = [str(a.pdf_weight + a.pdf_nmem + k) for k in (1, 2)]
        missing = [n for n in SEVEN + pdf + als if n not in v]
        if missing:
            raise SystemExit(f"{o}: variations missing: {missing[:5]} ...")
        out["obs"][o] = {"edges": h["edges"], "central": get(""), "err": fb(v[""]["err"]),
                         "scale": [get(n) for n in SEVEN], "pdf": [get(n) for n in pdf], "alphas": [get(n) for n in als],
                         "isr": [get(f"isr.murfac={f}") for f in ("0.5", "2.0") if f"isr.murfac={f}" in v],
                         "fsr": [get(f"fsr.murfac={f}") for f in ("0.5", "2.0") if f"fsr.murfac={f}" in v]}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out) + "\n")
    print(f"{a.out}: {', '.join(a.obs)} from {a.summary.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
