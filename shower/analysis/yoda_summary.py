#!/usr/bin/env python3
"""Rivet output (YODA) of ggzh-shower -> JSON with every histogram and every weight variation.

    yoda_summary.py <file.yoda> [more.yoda ...] --out summary.json

Several input files are combined as equivalent runs (same process, independent seeds): each file
is normalised to sum of weights / number of events, so the combination is the average weighted
by the number of events of each file (read from the analyses' event counters).  The RAW objects
are ignored.  Per histogram: bin edges, and per weight variation (the nominal is "") the
contents in pb per bin including underflow (first) and overflow (last) and their errors.
The STXS histograms of SimplifiedTemplateCrossSections are reduced to the gg -> ZH(ll) bins:
"stxs13_fine" (FWDH + 18 fine bins), "stxs13" (FWDH + 8 bins) and "stxs_unknown" (events
the routine could not classify).
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import yoda

# offsets of the production mode in the routine's linearised category histograms (GG2ZH)
STXS = {"stxs13_fine": ("STXS_stage1_3_fine_pTjet30", 187, 19),
        "stxs13": ("STXS_stage1_3_pTjet30", 59, 9)}
VAR = re.compile(r"^(?P<path>[^\[]+)(\[(?P<var>.*)\])?$")


def bins(h):
    """(edges, sumW incl. under/overflow, errW incl. under/overflow) for YODA 1 and 2."""
    try:
        bb = list(h.bins(True))                 # YODA 2: includeOverflows
        edges = list(h.xEdges())
        return edges, [b.sumW() for b in bb], [b.errW() for b in bb]
    except TypeError:
        bb = list(h.bins())
        edges = [b.xMin() for b in bb] + [bb[-1].xMax()]
        uf, of = h.underflow(), h.overflow()
        return (edges, [uf.sumW()] + [b.sumW() for b in bb] + [of.sumW()],
                [uf.errW()] + [b.errW() for b in bb] + [of.errW()])


def read(path: Path):
    """-> {histogram path: {variation: (edges, sumw, err)}}, number of events."""
    out, nev = {}, None
    for key, obj in yoda.read(str(path)).items():
        if key.startswith("/RAW/") or key.startswith("/_"):
            continue
        if "Histo1D" not in type(obj).__name__ and not hasattr(obj, "xEdges"):
            continue
        m = VAR.match(key)
        hp, var = m.group("path"), m.group("var") or ""
        out.setdefault(hp, {})[var] = bins(obj)
    for key, obj in yoda.read(str(path)).items():         # event count: /_EVTCOUNT if present
        if key.startswith("/_EVTCOUNT") and "[" not in key:
            try:
                nev = obj.numEntries()
            except Exception:
                pass
    return out, nev


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    total, ntot, comb = None, 0, {}
    for f in a.files:
        h, nev = read(f)
        nev = nev or 1
        ntot += nev
        for hp, vs in h.items():
            for var, (edges, sw, ew) in vs.items():
                c = comb.setdefault(hp, {}).setdefault(var, {"edges": edges, "sw": [0.0] * len(sw), "e2": [0.0] * len(ew)})
                c["sw"] = [x + nev * y for x, y in zip(c["sw"], sw)]
                c["e2"] = [x + (nev * y) ** 2 for x, y in zip(c["e2"], ew)]
    res = {"files": [str(f) for f in a.files], "events": ntot, "unit": "pb per bin, incl. under/overflow", "hists": {}}
    for hp, vs in comb.items():
        name = hp.split("/")[-1]
        ana = hp.split("/")[1] if hp.count("/") >= 2 else ""
        entry = {"edges": next(iter(vs.values()))["edges"], "variations": {}}
        for var, c in vs.items():
            entry["variations"][var] = {"sumw": [x / ntot for x in c["sw"]], "err": [x ** 0.5 / ntot for x in c["e2"]]}
        res["hists"][f"{ana}/{name}"] = entry
    # gg -> ZH(ll) STXS bins
    for key, (hname, off, n) in STXS.items():
        src = res["hists"].get(f"SimplifiedTemplateCrossSections/{hname}")
        if not src:
            continue
        res[key] = {var: {"sumw": v["sumw"][1 + off: 1 + off + n], "err": v["err"][1 + off: 1 + off + n]}
                    for var, v in src["variations"].items()}          # +1: underflow comes first
        if key == "stxs13_fine":
            res["stxs_unknown"] = {var: v["sumw"][1] for var, v in src["variations"].items()}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res) + "\n")
    nom = res["hists"].get("GGZH_NLOPS/sigma", {}).get("variations", {}).get("", {})
    print(f"{len(a.files)} files, {ntot} events, {len(res['hists'])} histograms, "
          f"{len(next(iter(comb.values())))} weight variations; sigma = {1e3 * sum(nom.get('sumw', [0])):.4f} fb")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
