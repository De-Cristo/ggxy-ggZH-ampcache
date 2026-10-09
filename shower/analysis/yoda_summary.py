#!/usr/bin/env python3
"""Rivet output of ggzh-shower (YODA text) -> compact JSON, per job and combined.

Reads the RAW objects (sums of weights before the analyses' finalize): for every histogram and
every weight variation the per-bin sum of weights and of squared weights, underflow first and
overflow last, and the event counter.  This makes the combination of jobs exact (sums add up)
and keeps the overflow bins, which the finalised Estimate objects of Rivet 4 drop.  Plain
Python: the `yoda` module is not needed.

    yoda_summary.py --raw  <file.yoda> --out job.json          (per job: sums only)
    yoda_summary.py --merge job-*.json --out summary.json      (sums of several jobs, finished)
    yoda_summary.py <file.yoda> ... --out summary.json         (read and finish directly)

Finished output (pb per bin; histograms incl. under/overflow): for each histogram the edges and,
per weight variation (the nominal is ""), "xs" = sum w / N and "err" = sqrt(sum w^2) / N.  From the
STXS routine only the gg -> ZH(l l) bins are kept: "stxs13_fine" (FWDH + 18), "stxs13" (FWDH + 8),
and "stxs_unknown" (events the routine could not classify; must be 0).
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

STXS = {"stxs13_fine": ("STXS_stage1_3_fine_pTjet30", 187, 19),   # GG2ZH offsets in the routine
        "stxs13": ("STXS_stage1_3_pTjet30", 59, 9)}
KEEP_STXS = {h for h, _, _ in STXS.values()}
BEGIN = re.compile(r"^BEGIN (\S+) (\S+)$")
PATH = re.compile(r"^/RAW/(?P<ana>[^/]+)/(?P<name>[^\[]+)(\[(?P<var>.*)\])?$")
CNT = re.compile(r"^/RAW/_EVTCOUNT(\[(?P<var>.*)\])?$")


def parse(path: Path) -> dict:
    """RAW sums of one YODA file."""
    out = {"files": 1, "events": 0, "hists": {}}
    typ = key = None
    body = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            m = BEGIN.match(line)
            if m:
                typ, key, body = m.group(1), m.group(2), []
                continue
            if line.startswith("END "):
                _store(out, typ, key, body)
                typ = key = None
                continue
            if typ is not None:
                body.append(line)
    return out


def _store(out, typ, key, body):
    if typ.startswith("YODA_COUNTER"):
        m = CNT.match(key)
        if m and not m.group("var"):
            row = [l for l in body if l and not l.startswith("#") and ":" not in l and l != "---"]
            out["events"] = int(round(float(row[0].split()[2])))
        return
    if not typ.startswith("YODA_HISTO1D"):
        return
    m = PATH.match(key)
    if not m:
        return
    ana, name, var = m.group("ana"), m.group("name"), m.group("var") or ""
    if ana.startswith("SimplifiedTemplateCrossSections") and name not in KEEP_STXS:
        return
    edges = None
    rows = []
    for l in body:
        if l.startswith("Edges(A1):"):
            edges = [float(x) for x in l.split(":", 1)[1].strip().strip("[]").split(",")]
        elif l and not l.startswith("#") and l != "---" and ":" not in l:
            rows.append(l.split())
    sw = [float(r[0]) for r in rows]
    sw2 = [float(r[1]) for r in rows]
    if ana.startswith("SimplifiedTemplateCrossSections"):
        ana = "STXS"
        for k, (h, off, n) in STXS.items():
            if h == name:                              # keep unknown (bin 0) + the GG2ZH bins
                sw = [sw[1]] + sw[1 + off: 1 + off + n]
                sw2 = [sw2[1]] + sw2[1 + off: 1 + off + n]
                name, edges = k, None
    h = out["hists"].setdefault(f"{ana}/{name}", {"edges": edges, "var": {}})
    h["var"][var] = {"sumw": sw, "sumw2": sw2}


def merge(parts: list[dict]) -> dict:
    out = {"files": 0, "events": 0, "hists": {}}
    for p in parts:
        out["files"] += p["files"]
        out["events"] += p["events"]
        for hk, h in p["hists"].items():
            o = out["hists"].setdefault(hk, {"edges": h["edges"], "var": {}})
            for v, d in h["var"].items():
                ov = o["var"].setdefault(v, {"sumw": [0.0] * len(d["sumw"]), "sumw2": [0.0] * len(d["sumw2"])})
                ov["sumw"] = [a + b for a, b in zip(ov["sumw"], d["sumw"])]
                ov["sumw2"] = [a + b for a, b in zip(ov["sumw2"], d["sumw2"])]
    return out


def finish(raw: dict) -> dict:
    n = raw["events"]
    res = {"files": raw["files"], "events": n, "unit": "pb per bin; histograms incl. underflow (first) and overflow (last)",
           "hists": {}}
    for hk, h in raw["hists"].items():
        e = {"edges": h["edges"], "var": {}}
        for v, d in h["var"].items():
            e["var"][v] = {"xs": [x / n for x in d["sumw"]], "err": [math.sqrt(max(x, 0.0)) / n for x in d["sumw2"]]}
        if hk.startswith("STXS/"):
            k = hk.split("/", 1)[1]
            res[k] = {v: {"xs": d["xs"][1:], "err": d["err"][1:]} for v, d in e["var"].items()}
            res.setdefault("stxs_unknown", {})
            for v, d in e["var"].items():
                res["stxs_unknown"][v] = max(res["stxs_unknown"].get(v, 0.0), d["xs"][0])
        else:
            res["hists"][hk] = e
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", type=Path)
    ap.add_argument("--merge", nargs="*", type=Path, default=[])
    ap.add_argument("--raw", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.merge:
        raw = merge([json.loads(p.read_text()) for p in a.merge])
    elif a.files:
        raw = merge([parse(f) for f in a.files])
    else:
        raise SystemExit("give YODA files or --merge")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    if a.raw:
        a.out.write_text(json.dumps(raw, separators=(",", ":")) + "\n")
        print(f"{raw['files']} files, {raw['events']} events, {len(raw['hists'])} histograms (raw sums)")
        return 0
    res = finish(raw)
    a.out.write_text(json.dumps(res, separators=(",", ":")) + "\n")
    s = res["hists"].get("GGZH_NLOPS/sigma", {}).get("var", {}).get("", {}).get("xs", [0, 0, 0])
    nv = len(res["hists"].get("GGZH_NLOPS/sigma", {}).get("var", {}))
    print(f"{res['files']} files, {res['events']} events, {len(res['hists'])} histograms, {nv} weight variations; "
          f"sigma = {1e3 * sum(s):.4f} fb; STXS unknown {1e3 * res.get('stxs_unknown', {}).get('', 0.0):.4f} fb")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
