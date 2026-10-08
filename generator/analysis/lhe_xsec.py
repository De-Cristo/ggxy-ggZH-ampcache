#!/usr/bin/env python3
"""Inclusive cross section and its scale, PDF and alpha_s variations from POWHEG LHE files.

POWHEG normalises its events so that the cross section is the sum of the event weights
divided by the number of events, for the nominal weight (XWGTUP) and for every alternative
weight of the reweighting list (<rwgt><wgt id=...>, absolute values).  Several files (seeds)
are pooled: sigma_k = sum over all events / total number of events.  The MC error quoted is
sqrt(sum w_k^2) / N.  The meaning of each weight id is read from the <initrwgt> header
(lhapdf=<id> [renscfact=<x> facscfact=<y>]).

  scale:   the 9 (mu_R, mu_F) points at the nominal PDF; 7-point envelope without (1/2, 2), (2, 1/2)
  PDF:     Hessian, sqrt(sum_k (sigma_k - sigma_0)^2) / sigma_0 over the members 1..--pdf-nmem of
           --pdf-id (default NNPDF3.1, LHAPDF 325300: 100 eigenvectors)
  alpha_s: half the difference of the two alpha_s members (325401: 0.116, 325402: 0.120),
           i.e. the change for delta alpha_s = 0.002

    lhe_xsec.py --files out/events-*.lhe.gz --out events.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np

WDEF = re.compile(r"<weight\s+id=['\"]([^'\"]+)['\"]\s*>\s*lhapdf=(\d+)(?:\s+renscfact=(\S+)\s+facscfact=(\S+))?")
WGT = re.compile(r"<wgt\s+id=['\"]([^'\"]+)['\"]\s*>\s*([^<\s]+)\s*</wgt>")


def fnum(s: str) -> float:
    return float(s.replace("d", "e").replace("D", "E"))


def read(path: Path):
    """-> weight definitions {id: (lhapdf, muR, muF)}, list of ids, nominal weights, matrix of alternative weights."""
    opener = gzip.open if path.suffix == ".gz" else open
    wdef, ids, xw, rows = {}, None, [], []
    in_event = header_line = False
    with opener(path, "rt") as f:
        for line in f:
            if not ids and "<weight" in line:
                m = WDEF.search(line)
                if m:
                    wdef[m.group(1)] = (int(m.group(2)), fnum(m.group(3)) if m.group(3) else None,
                                        fnum(m.group(4)) if m.group(4) else None)
            if line.startswith("<event"):
                in_event, header_line, cur = True, True, []
                continue
            if in_event:
                if header_line:
                    xw.append(fnum(line.split()[2])); header_line = False
                elif "<wgt" in line:
                    cur.extend(WGT.findall(line))
                elif line.startswith("</event"):
                    in_event = False
                    if cur:
                        these = [c[0] for c in cur]
                        if ids is None:
                            ids = these
                        if these != ids:
                            raise SystemExit(f"{path}: weight ids change between events")
                        rows.append([fnum(c[1]) for c in cur])
    return wdef, ids or [], np.array(xw), np.array(rows) if rows else np.zeros((len(xw), 0))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--files", nargs="+", type=Path, required=True)
    ap.add_argument("--pdf-id", type=int, default=325300)
    ap.add_argument("--pdf-nmem", type=int, default=100)
    ap.add_argument("--label", default="events")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    wdef, ids, xw, W = {}, None, [], []
    for f in args.files:
        d, i, x, w = read(f)
        wdef.update(d)
        if ids is None:
            ids = i
        if i != ids:
            raise SystemExit(f"{f}: different weight list")
        xw.append(x); W.append(w)
    xw = np.concatenate(xw); W = np.concatenate(W) if ids else np.zeros((len(xw), 0))
    n = len(xw)
    if n == 0:
        raise SystemExit("no events")
    pb_to_fb = 1e3
    out = {"label": args.label, "unit": "fb", "files": len(args.files), "events": n,
           "sigma_fb": pb_to_fb * xw.sum() / n, "mc_err_fb": pb_to_fb * np.sqrt((xw ** 2).sum()) / n,
           "negative_weight_fraction": float((xw < 0).mean())}
    s = W.sum(0) / n * pb_to_fb if ids else np.array([])
    e = np.sqrt((W ** 2).sum(0)) / n * pb_to_fb if ids else np.array([])
    col = {w: k for k, w in enumerate(ids or [])}

    scale = {}
    for w, (pdf, r, f_) in wdef.items():
        if r is not None and pdf == args.pdf_id and w in col:
            scale[(r, f_)] = col[w]
    if (1.0, 1.0) in scale:
        c = s[scale[(1.0, 1.0)]]
        out["scale_points"] = [{"muR": r, "muF": f_, "sigma_fb": float(s[k]), "mc_err_fb": float(e[k]),
                                "rel_to_central": float(s[k] / c - 1)} for (r, f_), k in sorted(scale.items())]
        seven = [k for (r, f_), k in scale.items() if (r, f_) not in ((0.5, 2.0), (2.0, 0.5))]
        out["scale7_band_rel"] = [float(max(s[seven]) / c - 1), float(min(s[seven]) / c - 1)]
        out["scale9_band_rel"] = [float(max(s[list(scale.values())]) / c - 1), float(min(s[list(scale.values())]) / c - 1)]
        out["nominal_over_rwl_central"] = float(out["sigma_fb"] / c)

    member = {pdf: col[w] for w, (pdf, r, _) in wdef.items() if r is None and w in col}
    p0 = member.get(args.pdf_id)
    if p0 is not None and all(args.pdf_id + k in member for k in range(1, args.pdf_nmem + 1)):
        s0 = s[p0]
        out["pdf_hessian_rel"] = float(np.sqrt(sum((s[member[args.pdf_id + k]] - s0) ** 2
                                                    for k in range(1, args.pdf_nmem + 1))) / s0)
        a_lo, a_hi = member.get(args.pdf_id + args.pdf_nmem + 1), member.get(args.pdf_id + args.pdf_nmem + 2)
        if a_lo is not None and a_hi is not None:
            out["alphas_rel_for_0p002"] = float(0.5 * (s[a_hi] - s[a_lo]) / s0)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(f"{args.label}: {len(args.files)} files, {n} events, sigma = {out['sigma_fb']:.4f} +- {out['mc_err_fb']:.4f} fb, "
          f"negative weights {100 * out['negative_weight_fraction']:.2f} %")
    if "scale7_band_rel" in out:
        print(f"  7-point band {100 * out['scale7_band_rel'][0]:+.1f}/{100 * out['scale7_band_rel'][1]:+.1f} %")
    if "pdf_hessian_rel" in out:
        print(f"  PDF (Hessian) {100 * out['pdf_hessian_rel']:.2f} %", end="")
        print(f", alpha_s (+-0.002) {100 * out['alphas_rel_for_0p002']:+.2f} %" if "alphas_rel_for_0p002" in out else "")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
