#!/usr/bin/env python3
"""Inclusive cross sections from POWHEG stage-2 integrations (one directory per scale point).

In parallel stage 2 every seed integrates independently on the common stage-1 grid and
writes its own pwg-st2-<seed>-stat.dat; the result of a point is the error-weighted mean over
its seeds.  The first seed's file (or pwg-stat.dat of a gridpack) alone is ONE seed's estimate,
not the combination.  If stage 3 was run, POWHEG's own combination (pwg-st3-0001-stat.dat) is
reported as well; the two agree.

    powheg_stat.py --points r1f1=out/r1f1 r0p5f0p5=out/r0p5f0p5 ... --out powheg-lo.json

The point named r1f1 is the central one; the band is the envelope of all points.

With Born suppression (bornsuppfact) POWHEG integrates sigma times the suppression factor F and
prints "total (btilde+remnants) cross section times suppression factor"; that integral is not
the cross section (the events carry 1/F, and their sum of weights is).  It is then reported as
"sigma_times_suppression_fb" with "born_suppressed": true, and no "sigma_fb" or band.
"""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

TOTAL = re.compile(r"total \(btilde\+remnants\) cross section\s+(times\s+suppression factor\s+)?in pb\s+(\S+)\s+\+-\s+(\S+)")


def point(d: Path):
    files = sorted(glob.glob(str(d / "pwg-st2-*-stat.dat"))) or sorted(glob.glob(str(d / "pwg-[0-9]*-stat.dat")))
    if not files:
        raise SystemExit(f"no stage-2 stat files in {d}")
    vals, supp = [], set()
    for f in files:
        m = TOTAL.search(Path(f).read_text())
        if not m:
            raise SystemExit(f"no total cross section in {f}")
        supp.add(bool(m.group(1)))
        vals.append((float(m.group(2).replace("D", "E")), float(m.group(3).replace("D", "E"))))
    if len(supp) > 1:
        raise SystemExit(f"{d}: stat files with and without Born suppression")
    w = [1.0 / e ** 2 for _, e in vals]
    mean = sum(v * wi for (v, _), wi in zip(vals, w)) / sum(w)
    res = {"sigma_fb": 1e3 * mean, "mc_err_fb": 1e3 * sum(w) ** -0.5, "seeds": len(vals),
           "per_seed_fb": [1e3 * v for v, _ in vals], "born_suppressed": supp.pop()}
    st3 = sorted(glob.glob(str(d / "pwg-st3-*-stat.dat")))
    m = TOTAL.search(Path(st3[0]).read_text()) if st3 else None
    if m:
        res["powheg_st3_fb"] = [1e3 * float(m.group(2).replace("D", "E")), 1e3 * float(m.group(3).replace("D", "E"))]
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", nargs="+", required=True, help="name=directory; r1f1 is the central point")
    ap.add_argument("--label", default="powheg")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    pts = {}
    for spec in args.points:
        name, d = spec.split("=", 1)
        pts[name] = point(Path(d))
    out = {"label": args.label, "unit": "fb", "points": pts}
    if any([p.pop("born_suppressed") for p in pts.values()]):
        # the integral of sigma x F: never report it under the name of a cross section
        for p in pts.values():
            p["sigma_times_suppression_fb"] = p.pop("sigma_fb")
        out["born_suppressed"] = True
        if "r1f1" in pts:
            out["sigma_times_suppression_fb"] = pts["r1f1"]["sigma_times_suppression_fb"]
            out["mc_err_fb"] = pts["r1f1"]["mc_err_fb"]
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(out, indent=1) + "\n")
        for name, p in pts.items():
            print(f"{name:10s} {p['sigma_times_suppression_fb']:.4f} +- {p['mc_err_fb']:.4f} fb ({p['seeds']} seeds; "
                  "Born-suppressed integral, not the cross section)")
        return 0
    if "r1f1" in pts:
        c = pts["r1f1"]["sigma_fb"]
        out["sigma_fb"] = c
        out["mc_err_fb"] = pts["r1f1"]["mc_err_fb"]
        out["band_rel"] = [max(p["sigma_fb"] for p in pts.values()) / c - 1, min(p["sigma_fb"] for p in pts.values()) / c - 1]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    for name, p in pts.items():
        st3 = f"; POWHEG stage-3 combination {p['powheg_st3_fb'][0]:.4f} +- {p['powheg_st3_fb'][1]:.4f}" if "powheg_st3_fb" in p else ""
        print(f"{name:10s} {p['sigma_fb']:.4f} +- {p['mc_err_fb']:.4f} fb ({p['seeds']} seeds){st3}")
    if "band_rel" in out and len(pts) > 1:
        print(f"central {out['sigma_fb']:.4f} fb, band {100 * out['band_rel'][0]:+.1f}/{100 * out['band_rel'][1]:+.1f} %")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
