#!/usr/bin/env python3
"""Results of the amplitude-cache timing scan (run/cache_scan.sh), combined over jobs.

    cache_scan.py --jobs <job>/out ... --out cache-scan.json [--extrapolate 884]

Each job directory holds w<N>/timing.txt (run/cache_ab.sh: the same events with the cache on and
off, N weights) and w<N>/compare-<seed>.json.  Per N and job:
  evaluations    one-loop amplitude evaluations per event, from the cache counters of the cached
                 run: with the cache = misses / events; without = (hits + misses) / events (every
                 request is an evaluation).  Independent of the machine.
  time_rel       wall time per event relative to the run without the cache and without weights
                 (N = 0) of the same job, for the uncached and the cached run: ratios on one machine,
                 no absolute times are stored.
  identical      the cached and uncached LHE files are identical (compare_lhe.py).
Combined over jobs: mean and spread (min, max).  A straight line through the evaluations and the
relative times against N (least squares) extrapolates to --extrapolate weights (default 884, the
size of a list used in experiments; stock POWHEG-BOX-V2 allows at most 500).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_timing(path: Path):
    d = {}
    for line in path.read_text().splitlines():
        f = line.split()
        if len(f) < 7 or f[0] != "seed":
            continue
        kv = dict(x.split("=", 1) for x in f if "=" in x)
        mode = kv.get("cache")
        d[mode] = {"events": int(f[f.index("events") + 1]), "wall": float(f[f.index("wall") + 1].rstrip("s")),
                   "hits": int(kv.get("hits", 0)), "misses": int(kv.get("misses", 0))}
    return d


def fit(xs, ys):
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    a = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx if sxx else 0.0
    return a, my - a * mx


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=Path, nargs="+", required=True)
    ap.add_argument("--extrapolate", type=int, default=884)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    per = {}                                       # N -> list of per-job records
    for j in a.jobs:
        base = None
        runs = {}
        for w in sorted(j.glob("w*/timing.txt")):
            n = int(w.parent.name[1:])
            t = parse_timing(w)
            if "0" not in t or "1" not in t:
                continue
            cmp_ = [json.loads(c.read_text()) for c in w.parent.glob("compare-*.json")]
            runs[n] = (t, all(c.get("identical") and c.get("compared", 0) > 0 for c in cmp_) if cmp_ else None)
        if 0 not in runs:
            raise SystemExit(f"{j}: no run without weights (w0)")
        t0 = runs[0][0]["0"]
        base = t0["wall"] / t0["events"]
        for n, (t, ident) in runs.items():
            c, u = t["1"], t["0"]
            req = c["hits"] + c["misses"]
            per.setdefault(n, []).append({
                "job": str(j.parent.name), "events": c["events"],
                "evaluations_uncached": req / c["events"], "evaluations_cached": c["misses"] / c["events"],
                "time_rel_uncached": (u["wall"] / u["events"]) / base, "time_rel_cached": (c["wall"] / c["events"]) / base,
                "identical": ident})
    rows = []
    for n in sorted(per):
        r = per[n]
        m = lambda k: sum(x[k] for x in r) / len(r)
        rows.append({"weights": n, "jobs": len(r), **{k: m(k) for k in ("evaluations_uncached", "evaluations_cached",
                                                                          "time_rel_uncached", "time_rel_cached")},
                     "time_rel_uncached_range": [min(x["time_rel_uncached"] for x in r), max(x["time_rel_uncached"] for x in r)],
                     "time_rel_cached_range": [min(x["time_rel_cached"] for x in r), max(x["time_rel_cached"] for x in r)],
                     "evaluation_ratio": m("evaluations_uncached") / m("evaluations_cached") if m("evaluations_cached") else None,
                     "time_ratio": m("time_rel_uncached") / m("time_rel_cached"),
                     "identical": all(x["identical"] for x in r), "per_job": r})
    xs = [r["weights"] for r in rows]
    ext = {}
    for k in ("evaluations_uncached", "evaluations_cached", "time_rel_uncached", "time_rel_cached"):
        s, i = fit(xs, [r[k] for r in rows])
        ext[k] = {"slope": s, "intercept": i, "at": s * a.extrapolate + i}
    ext["evaluation_ratio_at"] = ext["evaluations_uncached"]["at"] / ext["evaluations_cached"]["at"]
    ext["time_ratio_at"] = ext["time_rel_uncached"]["at"] / ext["time_rel_cached"]["at"]
    out = {"jobs": len(a.jobs), "extrapolate_to": a.extrapolate, "points": rows, "linear_fit": ext}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(f"{'weights':>7s} {'evals/event':>12s} {'(cached)':>9s} {'ratio':>6s} {'time/t0':>8s} {'(cached)':>9s} {'ratio':>6s} identical")
    for r in rows:
        print(f"{r['weights']:7d} {r['evaluations_uncached']:12.1f} {r['evaluations_cached']:9.1f} "
              f"{r['evaluation_ratio'] or float('nan'):6.2f} {r['time_rel_uncached']:8.2f} {r['time_rel_cached']:9.2f} "
              f"{r['time_ratio']:6.2f} {r['identical']}")
    print(f"linear fit at {a.extrapolate} weights: evaluation ratio {ext['evaluation_ratio_at']:.1f}, time ratio {ext['time_ratio_at']:.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
