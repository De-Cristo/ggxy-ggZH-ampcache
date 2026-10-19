#!/usr/bin/env python3
"""Summary of the reweighting study with the ATLAS-style source (set atlas13600 of production/paper.sh,
rwt/README.md): one JSON with the numbers and figure data, in the layout of ggzh_rwt_summary.py (the CMS-style
study), from the samples, the classifier weights and the closure outputs of rwt/run.sh.

Methods, settings and pass criteria as in the CMS-style study:
  clf   classifier on m(ZH), pT(ZH), pT(Z), N_jets, |y(ZH)| (ggzh_reweight_bdt.py, run rwt-clf-atlas)
  map1  binned pT(Z) x N_jets map (ggzh_reweight_map.py, closure by ggzh_reweight_closure.py)
"atlas_jet": the classifier with the leading-jet pT as a sixth input (run rwt-clf-jet-atlas); its weights are
also a series of the hard-jet data.

    ggzh_rwt_summary_atlas.py --dir <directory with the files below> --out rwt-summary-atlas.json [--grids-xsec grids-xsec.json]

Files in --dir: lo-atlas.npz (ATLAS-style LO+PS: LO events with mu = H_T, PDF4LHC21; A14 as the ATLAS job options
set it, the ATLAS base settings and the ATLAS POWHEG veto), nlo-a14.npz (target: the production NLO events with the
same shower settings, PDF4LHC21 nominal); weights-clf-atlas.npz, weights-clf-jet-atlas.npz, bdt-clf-atlas.json,
bdt-clf-jet-atlas.json; closure-map-atlas.json.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ggzh_reweight_common import Sample, dumps  # noqa: E402
from ggzh_reweight_map import apply_map, build_map  # noqa: E402
from ggzh_reweight_closure import test_pois  # noqa: E402
from ggzh_rwt_summary import HARDJET_EDGES, METHODS, REGIONS, frac_region  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--grids-xsec", type=Path, default=None, help="grids-xsec.json of the LO grids of the source (LO cross section)")
    ap.add_argument("--package", default=None, help="package version for the provenance (default: the first line of VERSION)")
    args = ap.parse_args(argv)
    D = args.dir
    version = args.package
    if version is None:
        vf = Path(__file__).resolve().parents[1] / "VERSION"
        version = vf.read_text().splitlines()[0].strip() if vf.exists() else "unknown"
    src = Sample("lo", D / "lo-atlas.npz"); tgt = Sample("nlo", D / "nlo-a14.npz")
    ss, ts = src.in_window, tgt.in_window
    w = {"lo": np.ones(src.n)}           # every LO event counted once, as in the CMS-style study
    w["map1"] = apply_map(build_map(src, tgt, "1001"), src)[0]

    def load_w(name):
        z = np.load(D / name)
        assert np.array_equal(z["run"], src["run"]) and np.array_equal(z["event"], src["event"]), f"{name}: weights not aligned"
        return z["w_full"]
    w["clf"] = load_w("weights-clf-atlas.npz")
    w6 = load_w("weights-clf-jet-atlas.npz")
    out = {"provenance": {"package": version, "source": "lo-atlas.npz (ATLAS-style LO+PS, mu = H_T, PDF4LHC21, A14, ATLAS veto)",
                          "target": "nlo-a14.npz (production NLO events, Born suppression 250 GeV, PDF4LHC21 weight, A14)",
                          "closure_json": {"map1": "closure-map-atlas.json", "clf": "bdt-clf-atlas.json", "clf6": "bdt-clf-jet-atlas.json"}}}
    grids = json.loads(args.grids_xsec.read_text()) if args.grids_xsec and args.grids_xsec.exists() else None
    out["samples"] = {"lo_events": int(ss.sum()), "nlo_events": int(ts.sum()),
                      "sigma_lo_fb": grids["points"]["r1f1"]["sigma_fb"] if grids else None,
                      "sigma_lo_err_fb": grids["points"]["r1f1"]["mc_err_fb"] if grids else None, "sigma_nlo_fb": tgt.sigma_fb(),
                      "nlo_neff": float(tgt.w0[ts].sum() ** 2 / np.sum(tgt.w0[ts] ** 2)),
                      "lo_events_above_nominal_weight": int(np.sum(src["w0"].astype(float) > 1.0001 * np.median(src["w0"].astype(float))))}
    for t in (400, 600):
        m = ts & (tgt.ptz >= t); out["samples"][f"nlo_neff_ptz{t}"] = float(tgt.w0[m].sum() ** 2 / np.sum(tgt.w0[m] ** 2))
    out["regions"] = {}
    for var, a, b in REGIONS:
        fn, en = frac_region(tgt, tgt.w0, ts, var, a, b)
        row = {"nlo_share": fn}
        for k in ["lo"] + METHODS:
            fl, el = frac_region(src, w[k], ss, var, a, b)
            row[k] = [fl / fn, (fl / fn) * np.sqrt((el / fl) ** 2 + (en / fn) ** 2)]
        out["regions"][f"{var}_{a:g}_{b:g}"] = row
    out["pois"] = {k: test_pois(src, w[k], ss, tgt, tgt.w0, ts)["pulls"] for k in ["lo"] + METHODS}
    out["max_poi_pull"] = {k: float(max(abs(v) for v in out["pois"][k].values())) for k in out["pois"]}
    out["neff"] = {k: float(w[k][ss].sum() ** 2 / np.sum(w[k][ss] ** 2)) for k in METHODS}
    c1 = json.loads((D / "closure-map-atlas.json").read_text())
    c3 = json.loads((D / "bdt-clf-atlas.json").read_text())
    chi = {"lo": {**c3["full"]["unweighted"]["variables"], **c3["full"]["unweighted"]["features"]},
           "map1": {k: v for k, v in c1["full"]["variables"]["reweighted"].items() if not k.startswith("_")},
           "clf": {**c3["full"]["reweighted"]["variables"], **c3["full"]["reweighted"]["features"]}}
    out["chi2"] = {m: {k: [v["chi2"], v["ndof"]] for k, v in chi[m].items()} for m in chi}
    out["criteria"] = {"map1": {"1": c1["summary"]["1_stage13"], "2": c1["summary"]["2_variables"], "3": c1["summary"]["3_stability"]},
                       "clf": {"1": c3["summary"]["1_stage13"], "2": c3["summary"]["2_variables"], "2b": c3["summary"]["2b_features"],
                               "3": c3["summary"]["3_stability"]}}

    def hj(s):
        return s["ptj1"].astype(float) > s["mzh"].astype(float)
    # the six-input classifier (the layout of the CMS-style variants)
    c6 = json.loads((D / "bdt-clf-jet-atlas.json").read_text())
    v = {"source": "ATLAS-style; six inputs (leading-jet pT added)", "closure_json": "bdt-clf-jet-atlas.json", "regions": {}}
    for var, a, b in REGIONS:
        fn, en = frac_region(tgt, tgt.w0, ts, var, a, b); row = {}
        for k, ww in (("lo", w["lo"]), ("clf", w6)):
            fl, el = frac_region(src, ww, ss, var, a, b)
            row[k] = [fl / fn, (fl / fn) * np.sqrt((el / fl) ** 2 + (en / fn) ** 2)]
        v["regions"][f"{var}_{a:g}_{b:g}"] = row
    v["pois"] = {k: test_pois(src, ww, ss, tgt, tgt.w0, ts)["pulls"] for k, ww in (("lo", w["lo"]), ("clf", w6))}
    v["max_poi_pull"] = {k: float(max(abs(x) for x in v["pois"][k].values())) for k in v["pois"]}
    v["neff_clf"] = float(w6[ss].sum() ** 2 / np.sum(w6[ss] ** 2))
    v["chi2"] = {"lo": {k: [x["chi2"], x["ndof"]] for k, x in {**c6["full"]["unweighted"]["variables"], **c6["full"]["unweighted"]["features"]}.items()},
                 "clf": {k: [x["chi2"], x["ndof"]] for k, x in {**c6["full"]["reweighted"]["variables"], **c6["full"]["reweighted"]["features"]}.items()}}
    v["criteria"] = {"1": c6["summary"]["1_stage13"], "2": c6["summary"]["2_variables"], "2b": c6["summary"]["2b_features"], "3": c6["summary"]["3_stability"]}
    v["hardjet"] = {}
    for lab, fs in (("all", ss), ("ptz400", ss & (src.ptz > 400)), ("pth450", ss & (src["pth"].astype(float) > 450))):
        v["hardjet"][lab] = {k: float(ww[fs & hj(src)].sum() / ww[fs].sum()) for k, ww in (("lo", w["lo"]), ("clf", w6))}
    out["atlas_jet"] = v
    out["hardjet"] = {}
    for lab, fs, ft in (("all", ss, ts), ("ptz400", ss & (src.ptz > 400), ts & (tgt.ptz > 400)),
                        ("pth450", ss & (src["pth"].astype(float) > 450), ts & (tgt["pth"].astype(float) > 450))):
        out["hardjet"][lab] = {"lo": float(w["lo"][fs & hj(src)].sum() / w["lo"][fs].sum()), "lo_raw": int((fs & hj(src)).sum()),
                               "clf": float(w["clf"][fs & hj(src)].sum() / w["clf"][fs].sum()),
                               "clf6": float(w6[fs & hj(src)].sum() / w6[fs].sum()),
                               "nlo": float(tgt.w0[ft & hj(tgt)].sum() / tgt.w0[ft].sum())}
    e = np.array(HARDJET_EDGES)
    for var in ("ptz", "pth"):
        for k, s, ww, sel in (("lo", src, w["lo"], ss), ("clf", src, w["clf"], ss), ("clf6", src, w6, ss), ("nlo", tgt, tgt.w0, ts)):
            x = s[var].astype(float); h = hj(s)
            num = np.array([ww[sel & h & (x >= a) & (x < b)].sum() for a, b in zip(e[:-1], e[1:])])
            den = np.array([ww[sel & (x >= a) & (x < b)].sum() for a, b in zip(e[:-1], e[1:])])
            num2 = np.array([np.sum(ww[sel & h & (x >= a) & (x < b)] ** 2) for a, b in zip(e[:-1], e[1:])])
            out["hardjet"].setdefault(f"vs_{var}", {})[k] = {"edges": e.tolist(), "frac": (num / den).tolist(), "err": (np.sqrt(num2) / den).tolist()}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(dumps(out))
    print(f"wrote {args.out}")
    print("max POI pull:", {k: round(x, 2) for k, x in out["max_poi_pull"].items()}, " six inputs:", round(v["max_poi_pull"]["clf"], 2),
          " N_eff:", {k: round(x) for k, x in out["neff"].items()}, round(v["neff_clf"]))
    print("hard jets:", {k: {kk: (round(100 * vv, 2) if isinstance(vv, float) else vv) for kk, vv in d.items()} for k, d in out["hardjet"].items()
                         if not k.startswith("vs_")})
    return 0


if __name__ == "__main__":
    sys.exit(main())
