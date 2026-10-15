#!/usr/bin/env python3
"""Summary of the LO+PS -> NLO+PS reweighting study (rwt/README.md): one JSON with every number the paper quotes
and the data of its figures, from the samples, the classifier weights and the closure outputs of rwt/run.sh.

Methods (settings and pass criteria fixed before running):
  clf   classifier on m(ZH), pT(ZH), pT(Z), N_jets, |y(ZH)| (ggzh_reweight_bdt.py); the only method in the
        paper's figures and tables
  map1  binned pT(Z) x N_jets map (ggzh_reweight_map.py, closure by ggzh_reweight_closure.py); quoted in one
        sentence of the paper, as the reweighting in the STXS variables alone
"default": the same LO events showered with Pythia's default start (no veto), reweighted by a classifier with the
same settings; "default_jet", "cms_jet": the classifier with the leading-jet pT as a sixth input, on the
default-start and the CMS-style source.  The six-input default-start weights are the fourth series of the
hard-jet figure.

    ggzh_rwt_summary.py --dir <directory with the files below> --out rwt-summary.json [--grids-xsec grids-xsec.json]

Files in --dir: samples lo-cms.npz (CMS-style LO+PS, the source), lo-cms-default.npz (default start), nlo-cp5.npz
(NLO+PS target, production NLO events with the NNPDF3.1 weight); classifier weights weights-clf.npz,
weights-clf-default.npz, weights-clf-jet.npz, weights-clf-jet-default.npz and their closure outputs bdt-clf.json,
bdt-clf-default.json, bdt-clf-jet.json, bdt-clf-jet-default.json; closure-map.json (map1).  The paper draws its
figures from this JSON and from bdt-clf.json and bdt-clf-jet-default.json (tools/plot_rwt.py of the paper).
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

METHODS = ["map1", "clf"]
REGIONS = [("ptz", 250, 400), ("ptz", 400, 600), ("ptz", 600, np.inf), ("pth", 200, 450), ("pth", 450, np.inf),
           ("ptj1", 200, 400), ("ptj1", 450, np.inf), ("drzh", 1.0, 2.4), ("ptzh", 150, 300), ("lhe_mll", 15, 40)]
HARDJET_EDGES = [0, 75, 150, 250, 400, 600, 1000.0]


def frac_region(s, w, sel, var, a, b):
    x = s[var].astype(float); m = sel & (x >= a) & (x < b)
    return w[m].sum() / w[sel].sum(), np.sqrt(np.sum(w[m] ** 2)) / w[sel].sum()


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
    src = Sample("lo", D / "lo-cms.npz"); tgt = Sample("nlo", D / "nlo-cp5.npz")
    ss, ts = src.in_window, tgt.in_window
    # weights of the methods (nominal), all on the same source events
    w = {"lo": np.ones(src.n)}
    w["map1"] = apply_map(build_map(src, tgt, "1001"), src)[0]

    def load_w(name, s_):
        z = np.load(D / name)
        assert np.array_equal(z["run"], s_["run"]) and np.array_equal(z["event"], s_["event"]), f"{name}: weights not aligned"
        return z["w_full"]
    w["clf"] = load_w("weights-clf.npz", src)
    dsrc = Sample("lo", D / "lo-cms-default.npz"); dss = dsrc.in_window
    wd = {"lo": np.ones(dsrc.n), "clf": load_w("weights-clf-default.npz", dsrc)}
    out = {"provenance": {"package": version, "source": "lo-cms.npz (CMS-style LO+PS)",
                          "target": "nlo-cp5.npz (production NLO events, Born suppression 250 GeV, NNPDF3.1 weight)",
                          "closure_json": {"map1": "closure-map.json", "clf": "bdt-clf.json"}}}
    grids = json.loads(args.grids_xsec.read_text()) if args.grids_xsec and args.grids_xsec.exists() else None
    out["samples"] = {"lo_events": int(ss.sum()), "nlo_events": int(ts.sum()),
                      "sigma_lo_fb": grids["points"]["r1f1"]["sigma_fb"] if grids else None,
                      "sigma_lo_err_fb": grids["points"]["r1f1"]["mc_err_fb"] if grids else None, "sigma_nlo_fb": tgt.sigma_fb(),
                      "nlo_neff": float(tgt.w0[ts].sum() ** 2 / np.sum(tgt.w0[ts] ** 2))}
    for t in (400, 600):
        m = ts & (tgt.ptz >= t); out["samples"][f"nlo_neff_ptz{t}"] = float(tgt.w0[m].sum() ** 2 / np.sum(tgt.w0[m] ** 2))
    # region ratios, reweighted (or unweighted) LO / NLO, each normalised to itself
    out["regions"] = {}
    for var, a, b in REGIONS:
        fn, en = frac_region(tgt, tgt.w0, ts, var, a, b)
        row = {"nlo_share": fn}
        for k in ["lo"] + METHODS:
            fl, el = frac_region(src, w[k], ss, var, a, b)
            row[k] = [fl / fn, (fl / fn) * np.sqrt((el / fl) ** 2 + (en / fn) ** 2)]
        out["regions"][f"{var}_{a:g}_{b:g}"] = row
    # POI pulls
    out["pois"] = {k: test_pois(src, w[k], ss, tgt, tgt.w0, ts)["pulls"] for k in ["lo"] + METHODS}
    out["max_poi_pull"] = {k: float(max(abs(v) for v in out["pois"][k].values())) for k in out["pois"]}
    out["neff"] = {k: float(w[k][ss].sum() ** 2 / np.sum(w[k][ss] ** 2)) for k in METHODS}
    # pass/fail and chi2 from the closure runs (criteria fixed before running)
    c1 = json.loads((D / "closure-map.json").read_text())
    c3 = json.loads((D / "bdt-clf.json").read_text())
    chi = {"lo": {**c3["full"]["unweighted"]["variables"], **c3["full"]["unweighted"]["features"]},   # same binning as clf
           "map1": {k: v for k, v in c1["full"]["variables"]["reweighted"].items() if not k.startswith("_")},
           "clf": {**c3["full"]["reweighted"]["variables"], **c3["full"]["reweighted"]["features"]}}
    out["chi2"] = {m: {k: [v["chi2"], v["ndof"]] for k, v in chi[m].items()} for m in chi}
    out["criteria"] = {"map1": {"1": c1["summary"]["1_stage13"], "2": c1["summary"]["2_variables"], "3": c1["summary"]["3_stability"]},
                       "clf": {"1": c3["summary"]["1_stage13"], "2": c3["summary"]["2_variables"], "2b": c3["summary"]["2b_features"],
                               "3": c3["summary"]["3_stability"]}}

    # hard-jet configurations: leading jet pT > m(ZH)
    def hj(s):
        return s["ptj1"].astype(float) > s["mzh"].astype(float)

    def variant(s_, ssel, wts, cjson, desc):
        c = json.loads((D / cjson).read_text())
        v = {"source": desc, "closure_json": cjson, "regions": {}}
        for var, a, b in REGIONS:
            fn, en = frac_region(tgt, tgt.w0, ts, var, a, b); row = {}
            for k in ("lo", "clf"):
                fl, el = frac_region(s_, wts[k], ssel, var, a, b)
                row[k] = [fl / fn, (fl / fn) * np.sqrt((el / fl) ** 2 + (en / fn) ** 2)]
            v["regions"][f"{var}_{a:g}_{b:g}"] = row
        v["pois"] = {k: test_pois(s_, wts[k], ssel, tgt, tgt.w0, ts)["pulls"] for k in ("lo", "clf")}
        v["max_poi_pull"] = {k: float(max(abs(x) for x in v["pois"][k].values())) for k in v["pois"]}
        v["neff_clf"] = float(wts["clf"][ssel].sum() ** 2 / np.sum(wts["clf"][ssel] ** 2))
        v["chi2"] = {"lo": {k: [x["chi2"], x["ndof"]] for k, x in {**c["full"]["unweighted"]["variables"], **c["full"]["unweighted"]["features"]}.items()},
                     "clf": {k: [x["chi2"], x["ndof"]] for k, x in {**c["full"]["reweighted"]["variables"], **c["full"]["reweighted"]["features"]}.items()}}
        v["criteria"] = {"1": c["summary"]["1_stage13"], "2": c["summary"]["2_variables"], "2b": c["summary"]["2b_features"], "3": c["summary"]["3_stability"]}
        v["hardjet"] = {}
        for lab, fs in (("all", ssel), ("ptz400", ssel & (s_.ptz > 400)), ("pth450", ssel & (s_["pth"].astype(float) > 450))):
            v["hardjet"][lab] = {k: float(wts[k][fs & hj(s_)].sum() / wts[k][fs].sum()) for k in ("lo", "clf")}
        return v
    out["default"] = variant(dsrc, dss, wd, "bdt-clf-default.json",
                             "lo-cms-default.npz (same LO events, shower card lo-default: no veto, pTmaxMatch 0); five inputs")
    wd6 = load_w("weights-clf-jet-default.npz", dsrc)   # six inputs: also the default-start series of the hard-jet figure
    out["default_jet"] = variant(dsrc, dss, {"lo": wd["lo"], "clf": wd6}, "bdt-clf-jet-default.json",
                                 "default start; six inputs (leading-jet pT added)")
    out["cms_jet"] = variant(src, ss, {"lo": w["lo"], "clf": load_w("weights-clf-jet.npz", src)}, "bdt-clf-jet.json",
                             "CMS-style; six inputs (leading-jet pT added)")
    out["hardjet"] = {}
    for lab, fs, ft in (("all", ss, ts), ("ptz400", ss & (src.ptz > 400), ts & (tgt.ptz > 400)),
                        ("pth450", ss & (src["pth"].astype(float) > 450), ts & (tgt["pth"].astype(float) > 450))):
        out["hardjet"][lab] = {"lo": float(src.w0[fs & hj(src)].sum() / src.w0[fs].sum()), "lo_raw": int((fs & hj(src)).sum()),
                               "clf": float(w["clf"][fs & hj(src)].sum() / w["clf"][fs].sum()),
                               "nlo": float(tgt.w0[ft & hj(tgt)].sum() / tgt.w0[ft].sum())}
    # data of the hard-jet figure: share of events with p_T(j1) > m_ZH per bin of p_T(Z) and p_T(H)
    e = np.array(HARDJET_EDGES)
    for var in ("ptz", "pth"):
        for k, s, ww, sel in (("lo", src, w["lo"], ss), ("clf", src, w["clf"], ss), ("dclf", dsrc, wd6, dss), ("nlo", tgt, tgt.w0, ts)):
            x = s[var].astype(float); h = hj(s)
            num = np.array([ww[sel & h & (x >= a) & (x < b)].sum() for a, b in zip(e[:-1], e[1:])])
            den = np.array([ww[sel & (x >= a) & (x < b)].sum() for a, b in zip(e[:-1], e[1:])])
            num2 = np.array([np.sum(ww[sel & h & (x >= a) & (x < b)] ** 2) for a, b in zip(e[:-1], e[1:])])
            out["hardjet"].setdefault(f"vs_{var}", {})[k] = {"edges": e.tolist(), "frac": (num / den).tolist(), "err": (np.sqrt(num2) / den).tolist()}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(dumps(out))
    print(f"wrote {args.out}")
    print("max POI pull:", {k: round(v, 2) for k, v in out["max_poi_pull"].items()}, " N_eff:", {k: round(v) for k, v in out["neff"].items()})
    print("hard jets:", {k: {kk: (round(100 * vv, 2) if isinstance(vv, float) else vv) for kk, vv in v.items()} for k, v in out["hardjet"].items()
                         if not k.startswith("vs_")})
    return 0


if __name__ == "__main__":
    sys.exit(main())
