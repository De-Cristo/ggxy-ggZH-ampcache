#!/usr/bin/env python3
"""Validation summary of a gg -> ZH(ll) NanoGEN sample (one pass over the files).

Normalisation (POWHEG): cross section = sum of weights / number of events, with the weight
LHEWeight_originalXWGTUP x LHEScaleWeight[k] for the nine scale points; MC error sqrt(sum w0^2)/N.
NanoAOD order of LHEScaleWeight: [0] (muR, muF) = (0.5, 0.5), [1] (0.5, 1), [2] (0.5, 2), [3] (1, 0.5),
[4] (1, 1), [5] (1, 2), [6] (2, 0.5), [7] (2, 1), [8] (2, 2).

Contents of the output JSON:
  totals       inclusive cross section at the nine scale points (all events), 7-point band,
               fraction of negative weights
  stxs         HTXS Stage 1.2 fine categories for gg -> ZH(ll) (HTXS_stage1_2_fine_cat_pTjet30GeV
               501-515 = 0J/1J/>=2J x pT(V) [0,75,150,250,400,inf), jet-major), the > 400 GeV bins
               split at 600 GeV with the generator Z pT (Stage 1.3 layout, 18 cells; 500 = |y_H| > 2.5
               is counted separately), and the fraction of events in category 0 (must be 0: with
               ProductionMode AUTO every ggxy_ggZH event lands there)
  hists        pT(Z), pT(H), m_ZH, pT(ZH), y_ZH (last-copy generator Z and H), HTXS njets30 and the
               leading generator-jet pT (> 30 GeV, |eta| < 4.7), all nine scale points, with
               under/overflow, inclusive in y_H

    nanogen_summary.py --files 'out/nanogen_*.root' --label nlo --out summary-nlo.json
    nanogen_summary.py --merge chunk-*.json --label nlo --out summary-nlo.json   (combine chunks)
Needs uproot and awkward; about 1 s per 1000-event file.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

SCALES = [(0.5, 0.5), (0.5, 1.0), (0.5, 2.0), (1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (2.0, 0.5), (2.0, 1.0), (2.0, 2.0)]
CENTRAL = 4
SEVEN = [k for k, s in enumerate(SCALES) if s not in ((0.5, 2.0), (2.0, 0.5))]
HISTS = {
    "ptz": np.arange(0.0, 1000.0 + 1e-9, 25.0),
    "pth": np.arange(0.0, 1000.0 + 1e-9, 25.0),
    "mzh": np.arange(150.0, 1250.0 + 1e-9, 25.0),
    "ptzh": np.arange(0.0, 300.0 + 1e-9, 10.0),
    "yzh": np.linspace(-3.0, 3.0, 21),
    "njets": np.array([-0.5, 0.5, 1.5, 2.5, 3.5]),
    "ptj1": np.arange(0.0, 1000.0 + 1e-9, 25.0),
}
CELLS = [f"{p}_{j}" for p in ("0_75", "75_150", "150_250", "250_400", "400_600", "GE600") for j in ("0J", "1J", "GE2J")]
BRANCHES = ["LHEWeight_originalXWGTUP", "LHEScaleWeight", "HTXS_njets30", "HTXS_stage1_2_fine_cat_pTjet30GeV",
            "GenPart_pdgId", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass", "GenPart_statusFlags",
            "GenJet_pt", "GenJet_eta"]


def four_vector(pt, eta, phi, m):
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    return np.stack([np.sqrt(px ** 2 + py ** 2 + pz ** 2 + m ** 2), px, py, pz])


def cells_of(cat, zpt):
    """Stage 1.3 cell 0..17 (pT-major, as CELLS), -1 for |y_H| > 2.5 or other categories."""
    out = np.full(len(cat), -1, dtype=int)
    ok = (cat >= 501) & (cat <= 515)
    k = cat[ok] - 501
    jet, ptb = k // 5, k % 5
    q = np.where(ptb < 4, ptb, np.where(zpt[ok] >= 600.0, 5, 4))
    out[ok] = q * 3 + jet
    return out


def summarise(files: list[str], label: str) -> dict:
    import awkward as ak
    import uproot
    n = 0; nneg = 0; ncat0 = 0; nfwd = 0; nnozh = 0
    tot = np.zeros(9); tot2 = 0.0
    cw = np.zeros((9, 18)); cw2 = np.zeros(18)
    hw = {h: np.zeros((9, len(e) + 1)) for h, e in HISTS.items()}; hw2 = {h: np.zeros(len(e) + 1) for h, e in HISTS.items()}
    for f in files:
        a = uproot.open(f)["Events"].arrays(BRANCHES)
        n += len(a)
        w0 = ak.to_numpy(a["LHEWeight_originalXWGTUP"]); sv = ak.to_numpy(a["LHEScaleWeight"])
        if sv.ndim != 2 or sv.shape[1] != 9:
            raise SystemExit(f"{f}: LHEScaleWeight does not have 9 entries per event")
        wk = w0[:, None] * sv
        tot += wk.sum(0); tot2 += float((w0 ** 2).sum()); nneg += int((w0 < 0).sum())
        cat = ak.to_numpy(a["HTXS_stage1_2_fine_cat_pTjet30GeV"])
        ncat0 += int((cat == 0).sum()); nfwd += int((cat == 500).sum())
        last = (a["GenPart_statusFlags"] & (1 << 13)) > 0
        vec = {}
        for name, pdg in (("z", 23), ("h", 25)):
            sel = (a["GenPart_pdgId"] == pdg) & last
            first = lambda b: ak.to_numpy(ak.fill_none(ak.firsts(a[b][sel]), np.nan))
            vec[name] = four_vector(first("GenPart_pt"), first("GenPart_eta"), first("GenPart_phi"), first("GenPart_mass"))
        nnozh += int(np.sum(np.isnan(vec["z"][0]) | np.isnan(vec["h"][0])))
        pt = lambda v: np.hypot(v[1], v[2])
        zh = vec["z"] + vec["h"]
        zpt = np.nan_to_num(pt(vec["z"]), nan=-1.0)
        cell = cells_of(cat, zpt); keep = cell >= 0
        for k in range(9):
            cw[k] += np.bincount(cell[keep], weights=wk[keep, k], minlength=18)
        cw2 += np.bincount(cell[keep], weights=(w0 ** 2)[keep], minlength=18)
        with np.errstate(invalid="ignore", divide="ignore"):
            mzh = np.sqrt(np.maximum(zh[0] ** 2 - zh[1] ** 2 - zh[2] ** 2 - zh[3] ** 2, 0.0))
            yzh = 0.5 * np.log((zh[0] + zh[3]) / (zh[0] - zh[3]))
        jets = a["GenJet_pt"][(a["GenJet_pt"] > 30.0) & (abs(a["GenJet_eta"]) < 4.7)]
        ptj1 = ak.to_numpy(ak.fill_none(ak.max(jets, axis=1), -1.0))
        obs = {"ptz": pt(vec["z"]), "pth": pt(vec["h"]), "mzh": mzh, "ptzh": pt(zh), "yzh": yzh,
               "njets": np.minimum(ak.to_numpy(a["HTXS_njets30"]), 3).astype(float), "ptj1": ptj1}
        for h, e in HISTS.items():
            x = obs[h]
            ok = np.isfinite(x) & ((x >= 0) | (h == "yzh"))      # ptj1 = -1: no jet, not filled
            idx = np.searchsorted(e, x[ok], side="right")       # 0 = underflow, len(e) = overflow
            for k in range(9):
                hw[h][k] += np.bincount(idx, weights=wk[ok, k], minlength=len(e) + 1)
            hw2[h] += np.bincount(idx, weights=(w0 ** 2)[ok], minlength=len(e) + 1)
    if n == 0:
        raise SystemExit("no events")
    return {"label": label, "files": len(files), "events": n,
            "sums": {"tot": tot.tolist(), "tot2": tot2, "neg": nneg, "cat0": ncat0, "fwdh": nfwd, "nozh": nnozh,
                     "cells": cw.tolist(), "cells2": cw2.tolist(),
                     "hists": {h: hw[h].tolist() for h in HISTS}, "hists2": {h: hw2[h].tolist() for h in HISTS}}}


def finish(raw: dict) -> dict:
    """Sums -> cross sections (fb) and errors."""
    s, n = raw["sums"], raw["events"]
    fb = 1e3 / n
    tot = np.asarray(s["tot"]) * fb; err = np.sqrt(s["tot2"]) * fb
    c = tot[CENTRAL]
    out = {"label": raw["label"], "files": raw["files"], "events": n, "unit": "fb",
           "totals": {"scale_points": [{"muR": r, "muF": f, "sigma_fb": float(tot[k]), "mc_err_fb": float(err)}
                                       for k, (r, f) in enumerate(SCALES)],
                      "sigma_fb": float(c), "mc_err_fb": float(err),
                      "scale7_band_rel": [float(max(tot[SEVEN]) / c - 1), float(min(tot[SEVEN]) / c - 1)],
                      "negative_weight_fraction": s["neg"] / n},
           "stxs": {"cells": CELLS, "category0_fraction": s["cat0"] / n, "fwdh_fraction": s["fwdh"] / n,
                    "events_without_z_or_h": s["nozh"],
                    "sigma_fb": (np.asarray(s["cells"])[CENTRAL] * fb).tolist(),
                    "mc_err_fb": (np.sqrt(np.asarray(s["cells2"])) * fb).tolist(),
                    "scale_points_fb": (np.asarray(s["cells"]) * fb).tolist()},
           "hists": {}}
    for h, e in HISTS.items():
        w = np.asarray(s["hists"][h]) * fb
        out["hists"][h] = {"edges": e.tolist(), "note": "bins incl. underflow (first) and overflow (last), fb per bin",
                           "sigma_fb": w[CENTRAL].tolist(), "mc_err_fb": (np.sqrt(np.asarray(s["hists2"][h])) * fb).tolist(),
                           "scale_points_fb": w.tolist()}
    return out


def merge(parts: list[dict], label: str) -> dict:
    raw = {"label": label, "files": sum(p["files"] for p in parts), "events": sum(p["events"] for p in parts), "sums": {}}
    s0 = parts[0]["sums"]
    for k, v in s0.items():
        if isinstance(v, dict):
            raw["sums"][k] = {h: np.sum([np.asarray(p["sums"][k][h]) for p in parts], axis=0).tolist() for h in v}
        else:
            raw["sums"][k] = np.sum([np.asarray(p["sums"][k]) for p in parts], axis=0).tolist()
    return raw


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--files", nargs="*", default=[], help="NanoGEN files or glob patterns")
    ap.add_argument("--merge", nargs="*", type=Path, default=[], help="raw chunk outputs (--raw) to combine")
    ap.add_argument("--raw", action="store_true", help="write the raw sums (for chunks to be merged)")
    ap.add_argument("--label", default="sample")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.merge:
        raw = merge([json.loads(p.read_text()) for p in a.merge], a.label)
    else:
        files = sorted({f for pat in a.files for f in (glob.glob(pat) or [pat])})
        if not files:
            raise SystemExit("no NanoGEN files")
        raw = summarise(files, a.label)
    res = raw if a.raw else finish(raw)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    if not a.raw:
        t = res["totals"]
        print(f"{a.label}: {res['files']} files, {res['events']} events; sigma = {t['sigma_fb']:.3f} +- {t['mc_err_fb']:.3f} fb, "
              f"7-point band {100 * t['scale7_band_rel'][0]:+.1f}/{100 * t['scale7_band_rel'][1]:+.1f} %, "
              f"negative weights {100 * t['negative_weight_fraction']:.2f} %, HTXS category 0: {100 * res['stxs']['category0_fraction']:.2f} %")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
