#!/usr/bin/env python3
"""Reusable LO+PS -> NLO+PS reweighting models of gg -> ZH (rwt/README.md): the classifier, trained once on all
source and target events with the settings of ggzh_reweight_bdt.py, exported to ONNX, and the binned p_T(Z) x N_jets
map (ggzh_reweight_map.py) as a correctionlib table; each for the nominal target and its seven scale variations
(SCALE7 of ggzh_reweight_common.py).  The per-event weights of the paper are out of fold (five folds) and are not a
single model; these are.  Validity: the source definition below; the paper's section on reweighting gives the tests
these weights pass and fail (a vetoed LO+PS sample lacks events with a jet harder than m_ZH, which no weight recovers).

    ggzh_rwt_export.py --source lo-cms.npz --target nlo-cp5.npz --out-dir model [--compare-weights weights-clf.npz]

Outputs in --out-dir:
  clf5-<var>.onnx, clf6-<var>.onnx   classifier with five inputs (m_ZH, p_T(ZH), p_T(Z), N_jets, |y_ZH|) and with six
                                     (and the leading-jet p_T); <var> = nominal or MUR<r>_MUF<f>
  map-ptz-njets.json                 correctionlib: weight(systematic, ptz, njets)
  model-card.json                    inputs (order, definitions, units), weight formula and normalisations, training
                                     samples, settings, software versions, checks (ONNX against scikit-learn,
                                     correctionlib against the map code)
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ggzh_reweight_common import FIT_PT_MIN, MU_OF_ID, N_JET, PT_EDGES, SCALE7, Sample, dumps  # noqa: E402
from ggzh_reweight_map import FIRST_FIT_BIN, apply_map, build_map  # noqa: E402
from ggzh_reweight_bdt import CLF, fit, prepare  # noqa: E402

INPUTS = {
    "mzh": ("invariant mass of the Z + H system", "GeV"),
    "ptzh": ("transverse momentum of the Z + H system", "GeV"),
    "ptz": ("transverse momentum of the Z boson", "GeV"),
    "njets": ("number of jets with p_T > 30 GeV, no rapidity requirement, capped at 3", ""),
    "ayzh": ("|rapidity| of the Z + H system", ""),
    "ptj1": ("transverse momentum of the leading jet with |eta| < 4.7; -1 without such a jet", "GeV"),
}
OBJECTS = ("Z and H: the last copies of the bosons in the event record (Z from its decay leptons' parent); jets: anti-kT "
           "R = 0.4 from all stable particles except the descendants of the Z and the Higgs boson")
MODELS = {"clf5": ["mzh", "ptzh", "ptz", "njets", "ayzh"], "clf6": ["mzh", "ptzh", "ptz", "njets", "ayzh", "ptj1"]}


def varname(wid):
    if wid == "1001":
        return "nominal"
    r, f = MU_OF_ID[wid]
    return f"MUR{float(r):.1f}_MUF{float(f):.1f}"


def to_onnx(clf, nfeat):
    from skl2onnx import convert_sklearn
    from skl2onnx.common.data_types import FloatTensorType
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return convert_sklearn(clf, initial_types=[("X", FloatTensorType([None, nfeat]))],
                               options={id(clf): {"zipmap": False}}, target_opset=17)


def snap_thresholds(clf):
    """ONNX stores the tree thresholds in float32.  Each threshold t becomes the largest float32 value not above t: for
    float32 inputs (the card's input type) x <= t is unchanged, and scikit-learn and ONNX take the same branches."""
    for trees in clf._predictors:
        for tree in trees:
            t = tree.nodes["num_threshold"]
            t32 = t.astype(np.float32)
            up = t32.astype(np.float64) > t
            t32[up] = np.nextafter(t32[up], np.float32(-np.inf))
            tree.nodes["num_threshold"] = t32.astype(np.float64)


def onnx_p(model_bytes, X):
    import onnxruntime as ort
    so = ort.SessionOptions()   # one thread: onnxruntime ignores OMP_NUM_THREADS (run.sh) and takes every core otherwise
    so.intra_op_num_threads = so.inter_op_num_threads = 1
    s = ort.InferenceSession(model_bytes, sess_options=so, providers=["CPUExecutionProvider"])
    names = [o.name for o in s.get_outputs()]
    return s.run(None, {"X": X.astype(np.float32)})[names.index("probabilities")][:, 1].astype(np.float64)


# the samples of each model set: "cms" (set 13600, the paper's reweighting study) and "atlas" (set atlas13600)
STYLES = {
    "cms": {"label": "CMS-style",
            "source": "CMS-style LO+PS at 13.6 TeV: LO gg -> l+l-H events of ggxy_ggZH (POWHEG-BOX-V2, bornonly), NNPDF3.1 NNLO, "
                      "fixed scale mu_R = mu_F = m_Z + m_H, 15 < m_ll < 150 GeV, showered with Pythia 8 (CP5, the CMS common "
                      "settings and the CMS POWHEG veto at SCALUP = m_ZH); production runs table-lo-cms",
            "target": "NLO+PS at 13.6 TeV: the production NLO events (Born suppression 250 GeV) with the NNPDF3.1 central member "
                      "as nominal weight, showered with CP5 and the CMS common settings; production run table-nlo-cp5; the "
                      "variations use the target's scale weights"},
    "atlas": {"label": "ATLAS-style",
              "source": "ATLAS-style LO+PS at 13.6 TeV: LO gg -> l+l-H events of ggxy_ggZH (POWHEG-BOX-V2, bornonly), PDF4LHC21, "
                        "dynamic scale mu_R = mu_F = H_T = m_T(H) + p_T(l1) + p_T(l2) (dynscale 3), 15 < m_ll < 150 GeV, showered "
                        "with Pythia 8 (A14 as the ATLAS job options set it, the ATLAS base settings and the ATLAS POWHEG veto at "
                        "SCALUP = m_ZH); set atlas13600, runs events-lo and table-lo-atlas",
              "target": "NLO+PS at 13.6 TeV: the production NLO events (Born suppression 250 GeV, PDF4LHC21 nominal weight), "
                        "showered with the same A14 settings, the ATLAS base settings and the NLO POWHEG veto; set atlas13600, "
                        "run table-nlo-a14; the variations use the target's scale weights"},
}


def map_correction(maps, label="CMS-style"):
    """correctionlib (schema v2) of the binned maps: per systematic, a binning in njets (0, 1, >=2) of a binning in
    p_T(Z) (25 GeV bins to 250 GeV: the binned ratio; above: exp(polynomial in log p_T(Z)) of the jet bin's fit),
    divided by the map's normalisation (mean weight 1 over the source events in the m_ll window)."""
    edges = [float(x) for x in PT_EDGES[:FIRST_FIT_BIN + 1]] + [1e6]   # last bin: p_T >= 250 GeV ("clamp" above 1e6)
    systs = []
    for wid, m in maps.items():
        jets = []
        for j in range(N_JET):
            content = []
            for i in range(FIRST_FIT_BIN):
                r = m.r[i * N_JET + j]
                content.append(float((1.0 if not np.isfinite(r) else r) / m.norm))
            coef = list(m.fits[j].coef)                       # highest power first, of log(r) vs log(pT)
            deg = len(coef) - 1
            terms = [f"[{k}]*pow(log(x),{deg - k})" if deg - k > 1 else (f"[{k}]*log(x)" if deg - k == 1 else f"[{k}]")
                     for k in range(len(coef))]
            pars = [float(c) for c in coef]
            pars[-1] -= float(np.log(m.norm))                 # the constant term carries the normalisation
            content.append({"nodetype": "formula", "expression": "exp(" + "+".join(terms) + ")", "parser": "TFormula",
                            "variables": ["ptz"], "parameters": pars})
            jets.append({"nodetype": "binning", "input": "ptz", "edges": edges, "content": content, "flow": "clamp"})
        systs.append({"key": varname(wid), "value": {"nodetype": "binning", "input": "njets", "edges": [0.0, 1.0, 2.0, 100.0],
                                                      "content": jets, "flow": "clamp"}})
    corr = {"name": "ggzh_lo_to_nlo_ptz_njets", "version": 1,
            "description": f"gg -> Z(ll)H at 13.6 TeV: weight of {label} LO+PS events towards NLO+PS in the generator-level "
                           "p_T(Z) x N_jets binning (rwt/README.md); shape only, mean 1 over the source sample inside "
                           "15 < m_ll < 150 GeV; events outside that window: weight 1",
            "inputs": [{"name": "systematic", "type": "string", "description": "nominal or MUR<r>_MUF<f>: the NLO+PS target "
                        "at that renormalisation and factorisation scale"},
                       {"name": "ptz", "type": "real", "description": "p_T of the Z boson (last copy in the record), GeV"},
                       {"name": "njets", "type": "real", "description": "number of jets, p_T > 30 GeV, no rapidity requirement"}],
            "output": {"name": "weight", "type": "real"},
            "data": {"nodetype": "category", "input": "systematic", "content": systs}}
    return {"schema_version": 2, "description": "ggxy-ggZH-ampcache rwt/: LO+PS -> NLO+PS reweighting map", "corrections": [corr]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--target", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--compare-weights", type=Path, default=None, help="per-event out-of-fold weights of rwt-clf, for a comparison")
    ap.add_argument("--package", default=None, help="package version for the card (default: the first line of VERSION)")
    ap.add_argument("--style", choices=sorted(STYLES), default="cms", help="the source and target descriptions of the card")
    args = ap.parse_args(argv)
    sty = STYLES[args.style]
    t0 = time.time()
    out = args.out_dir; out.mkdir(parents=True, exist_ok=True)
    version = args.package
    if version is None:
        vf = Path(__file__).resolve().parents[1] / "VERSION"
        version = vf.read_text().splitlines()[0].strip() if vf.exists() else "unknown"
    import onnx
    import onnxruntime
    import sklearn
    import skl2onnx
    src, tgt = Sample("source", args.source), Sample("target", args.target)
    prepare(src); prepare(tgt)                                # adds ayzh and njets (capped at 3)
    ss, ts = src.in_window, tgt.in_window
    wids = ["1001"] + list(SCALE7)
    card = {"package": version,
            "software": {"scikit-learn": sklearn.__version__, "skl2onnx": skl2onnx.__version__, "onnx": onnx.__version__,
                         "onnxruntime": onnxruntime.__version__},
            "source": sty["source"], "target": sty["target"],
            "objects": OBJECTS, "source_events_in_window": int(ss.sum()), "target_events_in_window": int(ts.sum()),
            "validity": "shape weights for samples defined as the source; outside that definition untested.  The paper's "
                        "reweighting section gives the tests (STXS fractions, distributions, split-sample stability) and their "
                        "outcome; a vetoed LO+PS sample lacks events with a jet harder than m_ZH, which no weight recovers.",
            "classifier": {"settings": CLF, "training": "all source events vs all target events in the window, balanced "
                           "classes, target events with their weights (negative weights set to zero)",
                           "input_dtype": "float32",
                           "thresholds": "each tree threshold stored as the largest float32 value not above the trained "
                                         "one, so that float32 inputs take the same branches in ONNX as in scikit-learn",
                           "output": "'probabilities'[:, 1] = p (probability of NLO+PS)",
                           "weight": "w = p / (1 - p) / c, with c the mean of p / (1 - p) over the source events in the "
                                     "window (below, per model and variation); events outside 15 < m_ll < 150 GeV: weight 1; "
                                     "for another sample, renormalise to mean 1 over its events in the window",
                           "models": {}},
            "map": {"file": "map-ptz-njets.json", "correction": "ggzh_lo_to_nlo_ptz_njets", "checks": {}},
            "target_sigma_fb": {varname(w): 1e3 * float(np.sum(tgt.weight(w)[ts])) / tgt.n for w in wids}}
    for name, feats in MODELS.items():
        Xs = np.stack([src[k].astype(float) for k in feats], axis=1)
        Xt = np.stack([tgt[k].astype(float) for k in feats], axis=1)
        Xs32 = Xs.astype(np.float32).astype(float)            # the source inputs as the models take them (float32)
        mc = {"inputs": [{"name": k, "definition": INPUTS[k][0], "unit": INPUTS[k][1]} for k in feats], "variations": {}}
        for wid in wids:
            t1 = time.time()
            clf = fit(Xs[ss], Xt[ts], tgt.weight(wid)[ts])
            snap_thresholds(clf)
            p = clf.predict_proba(Xs32[ss])[:, 1]
            c = float(np.mean(p / (1.0 - p)))
            model = to_onnx(clf, len(feats)).SerializeToString()
            fname = f"{name}-{varname(wid)}.onnx"
            (out / fname).write_bytes(model)
            po = onnx_p(model, Xs32[ss])
            wo, wk = po / (1 - po) / c, p / (1 - p) / c
            mc["variations"][varname(wid)] = {"onnx": fname, "c": c, "n_iter": int(clf.n_iter_),
                                              "max_abs_dp_onnx_vs_sklearn": float(np.max(np.abs(po - p))),
                                              "max_rel_dw_onnx_vs_sklearn": float(np.max(np.abs(wo / wk - 1))),
                                              "n_eff_source": float(wk.sum() ** 2 / np.sum(wk ** 2)), "seconds": round(time.time() - t1, 1)}
            print(f"{fname}: {clf.n_iter_} iterations, c = {c:.5f}, ONNX max |dp| {mc['variations'][varname(wid)]['max_abs_dp_onnx_vs_sklearn']:.2e}, "
                  f"{time.time() - t1:.0f} s", flush=True)
            if wid == "1001" and name == "clf5" and args.compare_weights and args.compare_weights.exists():
                z = np.load(args.compare_weights)
                assert np.array_equal(z["run"], src["run"]) and np.array_equal(z["event"], src["event"]), "weights not aligned"
                wf = z["w_full"][ss]
                card["classifier"]["comparison_with_out_of_fold_weights"] = {
                    "file": args.compare_weights.name, "mean_ratio": float(np.mean(wk / wf)),
                    "correlation": float(np.corrcoef(wk, wf)[0, 1]), "rms_rel_diff": float(np.sqrt(np.mean((wk / wf - 1) ** 2)))}
        card["classifier"]["models"][name] = mc
    # the binned map
    maps = {wid: build_map(src, tgt, wid) for wid in wids}
    cset = map_correction(maps, sty["label"])
    (out / "map-ptz-njets.json").write_text(json.dumps(cset, indent=1) + "\n")
    import correctionlib
    ev = correctionlib.CorrectionSet.from_file(str(out / "map-ptz-njets.json"))["ggzh_lo_to_nlo_ptz_njets"]
    nj = np.minimum(src["htxs_njets30"].astype(float), 2)
    for wid, m in maps.items():
        w_code = apply_map(m, src)[0][ss]
        w_cl = ev.evaluate(varname(wid), src.ptz[ss], nj[ss])
        card["map"]["checks"][varname(wid)] = {"max_abs_diff_correctionlib_vs_code": float(np.max(np.abs(w_cl - w_code))),
                                               "norm": m.norm, "fit_degrees": [f.degree for f in m.fits],
                                               "fits_pass_p_gt_0p05": [bool(f.passed) for f in m.fits]}
    card["map"]["fit_from_ptz_gev"] = FIT_PT_MIN
    card["map"]["correctionlib"] = correctionlib.__version__
    card["seconds"] = round(time.time() - t0, 1)
    (out / "model-card.json").write_text(dumps(card))
    worst = max(v["max_abs_dp_onnx_vs_sklearn"] for m in card["classifier"]["models"].values() for v in m["variations"].values())
    worstmap = max(v["max_abs_diff_correctionlib_vs_code"] for v in card["map"]["checks"].values())
    print(f"model card: {out / 'model-card.json'}; ONNX vs scikit-learn max |dp| {worst:.2e}; correctionlib vs code max |dw| {worstmap:.2e}; "
          f"{card['seconds']} s")
    return 0 if worst < 1e-5 and worstmap < 1e-9 else 1


if __name__ == "__main__":
    sys.exit(main())
