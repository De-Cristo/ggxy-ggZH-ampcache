#!/usr/bin/env python3
"""Event-by-event comparison of two LHE files generated with the same seed, e.g. with the
amplitude cache on and off (GGXY_AMPCACHE=1 / 0).

Kinematics: every event block with the <rwgt> block and POWHEG's '#rwgt' comment line removed
must be identical as text (the comment carries internal reweighting information and differs
at 1e-15 relative, from dividing out and multiplying back alpha_s).  Weights: every
<wgt id=...> value must be identical as a string.

    compare_lhe.py cached.lhe uncached.lhe [--out compare.json]

Exit status 0 if everything is identical, 1 otherwise.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

WGT = re.compile(r"<wgt\s+id=['\"]([^'\"]+)['\"]\s*>\s*([^<\s]+)\s*</wgt>")


def events(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    cur = None
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("<event"):
                cur = []
            elif cur is not None:
                if line.startswith("</event"):
                    yield cur
                    cur = None
                else:
                    cur.append(line)


def split(ev):
    kin, wgts, in_rwgt = [], [], False
    for line in ev:
        if "<rwgt>" in line:
            in_rwgt = True
        if in_rwgt or line.lstrip().startswith("#rwgt"):
            wgts.extend(WGT.findall(line))
        else:
            kin.append(line.strip())
        if "</rwgt>" in line:
            in_rwgt = False
    return "\n".join(kin), wgts


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", type=Path)
    ap.add_argument("b", type=Path)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args(argv)
    n = kin_same = w_tot = w_same = 0
    na = nb = 0
    ib = events(args.b)
    for ea in events(args.a):
        na += 1
        eb = next(ib, None)
        if eb is None:
            break
        nb += 1
        ka, wa = split(ea)
        kb, wb = split(eb)
        n += 1
        kin_same += ka == kb
        w_tot += max(len(wa), len(wb))
        w_same += sum(x == y for x, y in zip(wa, wb))
    nb += sum(1 for _ in ib)
    res = {"a": str(args.a), "b": str(args.b), "events_a": na, "events_b": nb, "compared": n,
           "kinematics_identical": kin_same, "weights": w_tot, "weights_identical": w_same}
    res["identical"] = na == nb and kin_same == n and w_same == w_tot
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(res, indent=1) + "\n")
    print(f"events {na} / {nb}; kinematics identical {kin_same}/{n}; weights identical {w_same}/{w_tot}"
          f" -> {'IDENTICAL' if res['identical'] else 'DIFFERENT'}")
    return 0 if res["identical"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
