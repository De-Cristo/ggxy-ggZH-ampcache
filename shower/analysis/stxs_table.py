#!/usr/bin/env python3
"""A prediction in the STXS Stage 1.3 fine bins of gg -> Z(ll)H with its uncertainties from the weights.

    stxs_table.py --summary shower-nlo-monash.json [--tunes shower-nlo-cp5.json shower-nlo-a14.json]
                  [--pdf-weight 3000 --pdf-nmem 40] --out stxs.json

Input: finished summaries of analysis/yoda_summary.py.  Per cell (FWDH, then the 18 fine bins in
the routine's order: p_T(V) bins for 0 jets, for 1 jet, for >= 2 jets), in fb:
  xs, mc_err      nominal weight
  scale           seven-point envelope of the MUR/MUF weights, relative [up, down]
  pdf             Hessian: sqrt(sum_k (x_k - x_0)^2) / x_0 over the --pdf-nmem members after the
                  central member (weight id --pdf-weight)
  alphas          half the difference of the two members after the Hessian ones, relative (for
                  PDF4LHC21 delta alpha_s = 0.001)
  pdf_alphas      both in quadrature
  isr, fsr        envelope of the shower renormalisation-scale weights isr.murfac / fsr.murfac in
                  {0.5, 2}, relative [up, down] (informational)
  tunes           each tune given with --tunes relative to this summary (same events, other shower;
                  informational)
The combination of the shower components into one uncertainty is left to the user.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

SEVEN = ["MUR1.0_MUF1.0", "MUR0.5_MUF0.5", "MUR2.0_MUF2.0", "MUR2.0_MUF1.0", "MUR0.5_MUF1.0", "MUR1.0_MUF2.0", "MUR1.0_MUF0.5"]
PTV = ["0-75", "75-150", "150-250", "250-400", "400-600", ">600"]
JETS = ["0J", "1J", ">=2J"]
CELLS = ["FWDH"] + [f"{p} {j}" for j in JETS for p in PTV]


def env(v, names, i, x0):
    vals = [v[n]["xs"][i] for n in names if n in v]
    if not vals or x0 == 0:
        return [float("nan"), float("nan")]
    return [max(vals) / x0 - 1, min(vals) / x0 - 1]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--tunes", type=Path, nargs="*", default=[])
    ap.add_argument("--pdf-weight", type=int, default=3000)
    ap.add_argument("--pdf-nmem", type=int, default=40)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    s = json.loads(a.summary.read_text())
    v = s["stxs13_fine"]
    tunes = {p.stem: json.loads(p.read_text())["stxs13_fine"][""] for p in a.tunes}
    w0, nm = a.pdf_weight, a.pdf_nmem
    pdfs = [str(w0 + k) for k in range(1, nm + 1)]
    if not all(k in v for k in [str(w0)] + pdfs):
        raise SystemExit(f"PDF weights {w0}..{w0 + nm} not in {a.summary}")
    rows = []
    for i, name in enumerate(CELLS):
        x0, e0 = v[""]["xs"][i], v[""]["err"][i]
        c = v[str(w0)]["xs"][i]
        pdf = math.sqrt(sum((v[k]["xs"][i] - c) ** 2 for k in pdfs)) / c if c else float("nan")
        alo, ahi = str(w0 + nm + 1), str(w0 + nm + 2)
        asr = 0.5 * (v[ahi]["xs"][i] - v[alo]["xs"][i]) / c if c and alo in v and ahi in v else float("nan")
        rows.append({"cell": name, "xs_fb": 1e3 * x0, "mc_err_fb": 1e3 * e0,
                     "scale": env(v, SEVEN, i, x0), "pdf": pdf, "alphas": asr, "pdf_alphas": math.hypot(pdf, asr),
                     "isr": env(v, ["isr.murfac=0.5", "isr.murfac=2.0"], i, x0),
                     "fsr": env(v, ["fsr.murfac=0.5", "fsr.murfac=2.0"], i, x0),
                     "tunes": {t: (d["xs"][i] / x0 - 1 if x0 else float("nan")) for t, d in tunes.items()}})
    tot = sum(r["xs_fb"] for r in rows)
    out = {"summary": a.summary.name, "tunes": list(tunes), "pdf_weight": w0, "pdf_nmem": nm, "sigma_fb": tot, "cells": rows}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(f"{'cell':14s} {'xs [fb]':>9s} {'MC %':>6s} {'scale %':>13s} {'PDF+as %':>8s} {'isr %':>12s} {'fsr %':>12s}  tunes %")
    for r in rows:
        print(f"{r['cell']:14s} {r['xs_fb']:9.4f} {100 * r['mc_err_fb'] / r['xs_fb'] if r['xs_fb'] else 0:6.1f} "
              f"{100 * r['scale'][0]:+6.1f}/{100 * r['scale'][1]:+5.1f} {100 * r['pdf_alphas']:8.1f} "
              f"{100 * r['isr'][0]:+5.1f}/{100 * r['isr'][1]:+5.1f} {100 * r['fsr'][0]:+5.1f}/{100 * r['fsr'][1]:+5.1f}  "
              + " ".join(f"{t.replace('shower-nlo-', '')} {100 * x:+.1f}" for t, x in r["tunes"].items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
