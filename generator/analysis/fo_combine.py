#!/usr/bin/env python3
"""Combine ggxy fixed-order runs (fo-ggllh / fo-ggzh histogram files, one per seed).

Each file holds, per histogram bin, the cross section for the 7 scale points of the
"mzh2" grid (mu_R, mu_F in {m_ZH/2 x (1/2, 1, 2)}, central first) and their MC errors, in
nb per unit of the observable.  Seeds with the same number of points are averaged; the
error of the mean is sqrt(sum err^2)/n.  fo-ggllh is for one lepton flavour: --nlep 3
gives Z -> l+l- with l = e, mu, tau.

    fo_combine.py --lo out/lo/histo_lo_*.txt --nlo out/nlo/histo_nlo_*.txt --nlep 3 --out fo.json

Output (fb): inclusive cross section and 7-point band at LO and NLO, K = NLO/LO, and the
differential distributions in m_ZH, pT(Z) (bins of the comparison in the reference) and pT(ZH), with
the 7-point band per bin and the differential K factor.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

NB_TO_FB = 1e6
GRID = [(1.0, 1.0), (0.5, 0.5), (2.0, 2.0), (2.0, 1.0), (0.5, 1.0), (1.0, 2.0), (1.0, 0.5)]   # (mu_R, mu_F)/mu_0, ggxy mc_gen order
REBIN = {
    "mzh25": [150, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 700, 800, 1000, 1250],
    "ptz25": [0, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 800, 1000],
    "ptzh": [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360, 420, 480, 540, 600],
}


def read_histo(path: Path, obs: str, scale: str = "mzh2"):
    """-> (edges, values[bin, 7], errors[bin, 7]) in nb per unit of the observable."""
    lines = path.read_text().splitlines()
    scales = lines[1].split("=")[1].split()
    ngrid = int(lines[2].split("=")[1])
    ipos = scales.index(scale)
    for i, line in enumerate(lines):
        tok = line.split()
        if tok and tok[0] == obs:
            edges = np.array([float(x) for x in tok[1:]])
            rows = np.array([[float(x) for x in lines[i + 1 + b].split()] for b in range(len(edges) - 1)])
            val = rows[:, 1 + ipos * ngrid: 1 + (ipos + 1) * ngrid]
            err = rows[:, 1 + (len(scales) + ipos) * ngrid: 1 + (len(scales) + ipos + 1) * ngrid]
            return edges, val, err
    raise KeyError(f"histogram {obs} not in {path}")


def combine(files: list[Path], obs: str, nlep: float):
    """Cross section per bin [fb] for the 7 scale points, averaged over seeds."""
    parts = [read_histo(f, obs) for f in files]
    edges = parts[0][0]
    assert all(np.array_equal(p[0], edges) for p in parts), f"different binning of {obs}"
    w = np.diff(edges)[:, None] * nlep * NB_TO_FB
    val = np.mean([p[1] for p in parts], axis=0) * w
    err = np.sqrt(np.sum([p[2] ** 2 for p in parts], axis=0)) / len(parts) * w
    return edges, val, err


def rebin(edges, val, err, new):
    idx = [int(np.argmin(abs(edges - x))) for x in new]
    assert np.allclose(edges[idx], new), (edges[idx], new)
    v = np.array([val[a:b].sum(0) for a, b in zip(idx[:-1], idx[1:])])
    e = np.array([np.sqrt((err[a:b] ** 2).sum(0)) for a, b in zip(idx[:-1], idx[1:])])
    return np.asarray(new, dtype=float), v, e


def inclusive(files, nlep):
    _, v, e = combine(files, "cross", nlep)
    c = float(v[0, 0])
    return {"sigma_fb": c, "mc_err_fb": float(e[0, 0]), "seeds": len(files),
            "points": [{"muR": r, "muF": f, "sigma_fb": float(v[0, k]), "mc_err_fb": float(e[0, k])} for k, (r, f) in enumerate(GRID)],
            "band_rel": [float(v[0].max() / c - 1), float(v[0].min() / c - 1)]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lo", nargs="*", type=Path, default=[])
    ap.add_argument("--nlo", nargs="*", type=Path, default=[])
    ap.add_argument("--nlep", type=float, default=3.0, help="3 for fo-ggllh (l = e, mu, tau), 1 for fo-ggzh")
    ap.add_argument("--label", default="")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    out = {"label": args.label, "nlep": args.nlep, "unit": "fb"}
    for order, files in (("lo", args.lo), ("nlo", args.nlo)):
        if files:
            out[order] = inclusive(files, args.nlep)
    if args.lo and args.nlo:
        out["K"] = out["nlo"]["sigma_fb"] / out["lo"]["sigma_fb"]
        out["K_err"] = out["K"] * float(np.hypot(out["lo"]["mc_err_fb"] / out["lo"]["sigma_fb"],
                                                 out["nlo"]["mc_err_fb"] / out["nlo"]["sigma_fb"]))
        out["differential"] = {}
        for obs, new in REBIN.items():
            try:
                lo = rebin(*combine(args.lo, obs, args.nlep), new)
                nlo = rebin(*combine(args.nlo, obs, args.nlep), new)
            except KeyError:
                continue
            with np.errstate(divide="ignore", invalid="ignore"):      # p_T(ZH): LO only in the first bin
                k = nlo[1][:, 0] / lo[1][:, 0]
                kerr = k * np.hypot(lo[2][:, 0] / lo[1][:, 0], nlo[2][:, 0] / nlo[1][:, 0])
            fin = lambda a: [float(x) if np.isfinite(x) else None for x in a]
            out["differential"][obs] = {"edges": lo[0].tolist(), "lo_fb": lo[1][:, 0].tolist(), "nlo_fb": nlo[1][:, 0].tolist(),
                                        "lo_err_fb": lo[2][:, 0].tolist(), "nlo_err_fb": nlo[2][:, 0].tolist(),
                                        "K": fin(k), "K_err": fin(kerr),
                                        # seven-point scale envelope per bin, [lowest, highest] in fb
                                        "lo_band_fb": np.stack([lo[1].min(1), lo[1].max(1)], 1).tolist(),
                                        "nlo_band_fb": np.stack([nlo[1].min(1), nlo[1].max(1)], 1).tolist()}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    for order in ("lo", "nlo"):
        if order in out:
            r = out[order]
            print(f"{order.upper():3s} {r['sigma_fb']:.4f} +- {r['mc_err_fb']:.4f} fb  "
                  f"band {100 * r['band_rel'][0]:+.1f}/{100 * r['band_rel'][1]:+.1f} %  ({r['seeds']} seeds)")
    if "K" in out:
        print(f"K = {out['K']:.4f} +- {out['K_err']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
