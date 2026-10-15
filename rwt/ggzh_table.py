#!/usr/bin/env python3
"""Merge the per-event tables of the Rivet routine GGZH_EVENTTABLE (shower/rivet/; one table-<seed>.csv.gz per
LHE file, production/paper.sh runs table-*) into one .npz with the column names that ggzh_reweight_common.Sample
and the map, closure and classifier scripts read.

Mapping: ptz, yz, mz, pth, yh, mzh, ptzh, yzh, drzh as written; htxs_njets30 = njets30 (jets pT > 30
GeV without a rapidity cut, as GGZH_NLOPS); njets_isr = njets30eta47; ptj1 = ptj1eta47 (|eta| < 4.7);
cell as written (-1 = |y_H| > 2.5); htxs_cat = 500 for FWDH, 501 + cell-derived jet-major code otherwise
(only "== 500" is used downstream); lhe_mll = mz (last-copy Z mass); lhe_vpt, lhe_pth, lhe_mzh, scalepdf =
nan (not available after the shower); w0 = genw = nominal weight (--nominal 2000: the NNPDF3.1 central
member of a PDF4LHC21 sample); sw = scale weights / nominal in the NanoAOD order of LHEScaleWeight (index 4
= nominal), ones when the sample has no scale weights; run = seed, event = ievt, ifile = index of the
table in the sorted list.

    ggzh_table.py --label lo-cms --tables 'tables/lo-cms/table-*.csv.gz' --out lo-cms.npz [--nominal 2000]
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np

# LHEScaleWeight order of NanoAOD: [0] (0.5,0.5) [1] (0.5,1) [2] (0.5,2) [3] (1,0.5) [4] (1,1) [5] (1,2) [6] (2,0.5) [7] (2,1) [8] (2,2)
NANO_ORDER = ["1009", "1007", "1008", "1003", "1001", "1002", "1006", "1004", "1005"]


def read_table(path: Path):
    with gzip.open(path, "rt") as f:
        header = f.readline().strip().split(",")
    data = np.genfromtxt(path, delimiter=",", skip_header=1, dtype=np.float64)
    data = np.atleast_2d(data)
    return header, data


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", required=True)
    ap.add_argument("--tables", nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--nominal", default="nominal", help="weight column used as nominal: 'nominal' or '2000' (NNPDF3.1 central "
                    "of a PDF4LHC21 sample); scale weights are then taken as ratios to the original nominal")
    args = ap.parse_args(argv)
    files = sorted({f for t in args.tables for f in glob.glob(t)}, key=lambda s: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", s)])
    if not files:
        raise SystemExit("no tables")
    cols, per_file, headers = {}, [], set()
    parts = []
    for i, f in enumerate(files):
        h, d = read_table(Path(f))
        headers.add(tuple(h))
        seed = int(re.search(r"table-(\d+)\.csv", f).group(1))
        parts.append((h, d, i, seed))
        per_file.append(len(d))
    if len(headers) != 1:
        raise SystemExit(f"tables have different headers: {headers}")
    h = list(headers.pop())
    d = np.concatenate([p[1] for p in parts])
    ix = {k: j for j, k in enumerate(h)}
    wnames = [k for k in h if k.startswith("w:")]
    nom0 = d[:, ix["w:nominal"]]
    nom = d[:, ix["w:" + args.nominal]] if args.nominal != "nominal" else nom0
    scale_ids = {}
    for k in wnames:
        m = re.search(r"^w:(100[1-9])$", k)
        if m:
            scale_ids[m.group(1)] = k
    if len(scale_ids) == 9:
        sw = np.stack([d[:, ix[scale_ids[i]]] / np.where(nom0 != 0, nom0, np.nan) for i in NANO_ORDER], axis=1)
    else:
        sw = np.ones((len(d), 9))
    n = len(d)
    cell = d[:, ix["cell"]].astype(np.int8)
    nj = d[:, ix["njets30"]].astype(int)
    htxs_cat = np.where(cell < 0, 500, 501 + 5 * np.minimum(nj, 2) + np.minimum(np.maximum(cell, 0) // 3, 4)).astype(np.int16)
    out = {
        "ptz": d[:, ix["ptz"]], "yz": d[:, ix["yz"]], "mz": d[:, ix["mz"]], "pth": d[:, ix["pth"]], "yh": d[:, ix["yh"]],
        "mzh": d[:, ix["mzh"]], "ptzh": d[:, ix["ptzh"]], "yzh": d[:, ix["yzh"]], "drzh": d[:, ix["drzh"]],
        "etaz": d[:, ix["etaz"]], "phiz": d[:, ix["phiz"]], "etah": d[:, ix["etah"]], "phih": d[:, ix["phih"]],
        "mll_dressed": np.full(n, np.nan), "ptj1": d[:, ix["ptj1eta47"]], "htxs_pth": d[:, ix["pth"]], "htxs_vpt": d[:, ix["ptz"]],
        "lhe_vpt": np.full(n, np.nan), "lhe_pth": np.full(n, np.nan), "lhe_mll": d[:, ix["mz"]], "lhe_mzh": np.full(n, np.nan),
        "scalepdf": np.full(n, np.nan),
        "htxs_njets30": np.minimum(nj, 100).astype(np.int8), "njets_isr": d[:, ix["njets30eta47"]].astype(np.int8), "cell": cell,
        "nlhe": np.zeros(n, np.int8), "htxs_cat": htxs_cat, "run": np.concatenate([np.full(len(p[1]), p[3], np.uint32) for p in parts]),
        "lumi": np.zeros(n, np.uint32), "event": d[:, ix["ievt"]].astype(np.uint64),
        "ifile": np.concatenate([np.full(len(p[1]), p[2], np.int16) for p in parts]),
        "w0": nom, "genw": nom, "sw": sw.astype(np.float32), "nb_from_h": np.full(n, 2, np.int8), "nl_from_z": np.full(n, 2, np.int8),
    }
    for k in ["ptz", "yz", "mz", "pth", "yh", "mzh", "ptzh", "yzh", "drzh", "etaz", "phiz", "etah", "phih", "mll_dressed", "ptj1", "htxs_pth", "htxs_vpt",
              "lhe_vpt", "lhe_pth", "lhe_mll", "lhe_mzh", "scalepdf"]:
        out[k] = out[k].astype(np.float32)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, **out)
    meta = {"label": args.label, "events": n, "files": files, "events_per_file": per_file, "weights_in_tables": wnames,
            "scale_weights_found": sorted(scale_ids), "events_without_z_or_h": int(np.sum(~np.isfinite(out["ptz"]))),
            "fwdh_events": int(np.sum(cell < 0)), "negative_w0_events": int(np.sum(nom < 0)), "sum_w0": float(nom.sum()),
            "unweighted": bool(np.all(nom == nom[0])), "nominal_column": args.nominal}
    args.out.with_suffix(".json").write_text(json.dumps(meta, indent=1) + "\n")
    print(f"{args.label}: {len(files)} tables, {n} events, scale weights {len(scale_ids)}/9, "
          f"FWDH {meta['fwdh_events']}, negative {meta['negative_w0_events']}, mean nominal weight {nom.mean():.6g} -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
