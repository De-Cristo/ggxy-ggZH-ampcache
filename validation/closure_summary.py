#!/usr/bin/env python3
"""Summary of the ggxy amplitude-cache closure jobs, and what their 884-weight events tell us.

Each closure job (a directory out-<i>/ under --area) generated the same events twice, same seed,
with the full CMS weight list: R3.lhe with the cache, R4.lhe without (GGXY_AMPCACHE=0), plus
timing.txt (lines 'R3 cache=1 rc=0 wall=<s>s per_event=<s>s', same for R4, and the cache's
hit/miss line) and payload.log ('host <name>', 'seed <n>').  This script

  1. compares R3 and R4 event by event (kinematics with the #rwgt comment removed, and every
     <wgt> string) and collects the timing and cache hit rate per job;
  2. uses the pooled events (identical in R3 and R4) to evaluate the inclusive PDF and alpha_s
     uncertainty of the NLO cross section from the PDF weights (NNPDF3.1 Hessian, LHAPDF
     325300: members 1-100 eigenvectors, 101/102 alpha_s = 0.116/0.120), with a bootstrap
     over events for the statistical error of these ratios, and the nine-point scale envelope.

    closure_summary.py --area <dir with out-<i>/> [--out closure.json]
"""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
AREA = None  # set with --area
EV = re.compile(r"<event>(.*?)</event>", re.S)
WGT = re.compile(r"<wgt id='([^']+)'>\s*([^<]+)</wgt>")
WDEF = re.compile(r"<weight id='([^']+)'\s*>\s*lhapdf=(\d+)(?:\s+renscfact=(\S+)\s+facscfact=(\S+))?")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--area", type=Path, required=True, help="directory with out-<i>/{R3.lhe,R4.lhe,timing.txt,payload.log}")
    ap.add_argument("--out", type=Path, default=Path("closure.json"))
    args = ap.parse_args(argv)
    global AREA
    AREA = args.area
    jobs = []; rows = []; wid_order = None; wdef = {}
    for d in sorted(glob.glob(str(AREA / "out-*")), key=lambda p: int(p.rsplit("-", 1)[1])):
        d = Path(d)
        if not (d / "R3.lhe").is_file() or not (d / "R4.lhe").is_file():
            continue
        t = (d / "timing.txt").read_text()
        m3 = re.search(r"R3 cache=1 rc=0 wall=(\d+)s per_event=([\d.]+)s", t); m4 = re.search(r"R4 cache=0 rc=0 wall=(\d+)s per_event=([\d.]+)s", t)
        hit = re.search(r"enabled=1 hits=(\d+) misses=(\d+)", t)
        pl = (d / "payload.log").read_text() if (d / "payload.log").is_file() else ""
        host = re.search(r"host (\S+)", pl); seed = re.search(r"seed (\d+)", pl)
        s3, s4 = (d / "R3.lhe").read_text(), (d / "R4.lhe").read_text()
        if not wdef:
            for wid, pdf, r, f in WDEF.findall(s3):
                wdef[wid] = {"lhapdf": int(pdf), "mur": r or None, "muf": f or None}
        e3, e4 = EV.findall(s3), EV.findall(s4)
        kin_same = w_same = w_tot = 0
        for a, b in zip(e3, e4):
            ka = re.sub(r"<rwgt>.*?</rwgt>|#rwgt.*", "", a, flags=re.S).strip(); kb = re.sub(r"<rwgt>.*?</rwgt>|#rwgt.*", "", b, flags=re.S).strip()
            kin_same += ka == kb
            wa, wb = WGT.findall(a), WGT.findall(b)
            w_tot += len(wa); w_same += sum(x[1].strip() == y[1].strip() for x, y in zip(wa, wb))
            ids = [x[0] for x in wa]
            if wid_order is None:
                wid_order = ids
            if ids == wid_order:
                rows.append([float(x[1]) for x in wa])
        jobs.append({"job": d.name, "host": host.group(1) if host else None, "seed": int(seed.group(1)) if seed else None,
                     "events": len(e3), "events_r4": len(e4), "kinematics_identical": kin_same, "weights": w_tot, "weights_identical": w_same,
                     "cached_s_per_event": float(m3.group(2)), "uncached_s_per_event": float(m4.group(2)),
                     "speedup": float(m4.group(1)) / float(m3.group(1)),
                     "cache_hit_rate": int(hit.group(1)) / (int(hit.group(1)) + int(hit.group(2))) if hit else None})
    W = np.array(rows)                         # events x weights
    idx = {w: i for i, w in enumerate(wid_order)}
    by_pdf = {wdef[w]["lhapdf"]: idx[w] for w in wid_order if w in wdef and wdef[w]["mur"] is None and 325300 <= wdef[w]["lhapdf"] <= 325402}
    scale_ids = [w for w in wid_order if w in wdef and wdef[w]["mur"] is not None]
    nominal = idx["1001"]

    def estimates(Wb):
        s = Wb.sum(0); s0 = s[nominal]
        pdf = np.sqrt(sum((s[by_pdf[325300 + k]] - s[by_pdf[325300]]) ** 2 for k in range(1, 101))) / s[by_pdf[325300]]
        als = 0.5 * (s[by_pdf[325402]] - s[by_pdf[325401]]) / s[by_pdf[325300]]
        sc = np.array([s[idx[w]] / s0 - 1 for w in scale_ids])
        return pdf, als, sc.max(), sc.min()
    pdf, als, up, dn = estimates(W)
    rng = np.random.default_rng(1)
    boot = np.array([estimates(W[rng.integers(0, len(W), len(W))]) for _ in range(500)])
    err = boot.std(0)
    tot = {k: sum(j[k] for j in jobs) for k in ("events", "kinematics_identical", "weights", "weights_identical")}
    out = {"area": str(AREA), "jobs": jobs, "total": tot,
           "aggregate_speedup": sum(j["uncached_s_per_event"] * j["events"] for j in jobs) / sum(j["cached_s_per_event"] * j["events"] for j in jobs),
           "cached_s_per_event_mean": float(np.mean([j["cached_s_per_event"] for j in jobs])),
           "uncached_s_per_event_mean": float(np.mean([j["uncached_s_per_event"] for j in jobs])),
           "pooled_events": len(W),
           "pdf_hessian_rel": float(pdf), "pdf_hessian_rel_booterr": float(err[0]),
           "alphas_rel_for_delta_0p002": float(als), "alphas_rel_booterr": float(err[1]),
           "scale9_envelope_rel": [float(up), float(dn)], "scale9_envelope_booterr": [float(err[2]), float(err[3])],
           "note": "PDF: NNPDF3.1 NNLO Hessian (LHAPDF 325300), sqrt(sum_k (s_k - s_0)^2)/s_0 over members 1-100; alpha_s: half the difference of members 102 (0.120) and 101 (0.116)"}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(f"{'job':7s} {'seed':>6s} {'host':34s} {'ev':>3s} {'kin':>4s} {'weights':>12s} {'cached':>7s} {'uncached':>8s} {'x':>5s} {'hit':>5s}")
    for j in jobs:
        print(f"{j['job']:7s} {j['seed']!s:>6s} {str(j['host'])[:34]:34s} {j['events']:3d} {j['kinematics_identical']:4d} {j['weights_identical']:6d}/{j['weights']:<6d} {j['cached_s_per_event']:6.2f}s {j['uncached_s_per_event']:7.1f}s {j['speedup']:5.1f} {100*j['cache_hit_rate']:4.1f}%")
    print(f"total: {tot}; aggregate speed-up {out['aggregate_speedup']:.1f}x; mean {out['cached_s_per_event_mean']:.2f} s vs {out['uncached_s_per_event_mean']:.1f} s per event")
    print(f"pooled {len(W)} events: PDF (Hessian) {100*pdf:.2f} ± {100*err[0]:.2f} %, alpha_s(±0.002) {100*als:+.2f} ± {100*err[1]:.2f} %, "
          f"9-point scale {100*up:+.1f}/{100*dn:+.1f} % (± {100*err[2]:.1f}/{100*err[3]:.1f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
