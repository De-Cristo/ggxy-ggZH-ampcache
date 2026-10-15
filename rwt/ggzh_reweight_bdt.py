#!/usr/bin/env python3
"""Multivariate (classifier) LO+PS -> NLO+PS reweighting of gg -> ZH at generator level, with the settings
and pass criteria fixed before this script was run on the full samples (rwt/README.md).

Series compared everywhere: "CMS-style LO+PS", "CMS-style LO x rwt + PS", "NLO+PS".
Features: m(ZH), pT(ZH), pT(Z), N_jets (0,1,2,>=3), |y(ZH)|.  HistGradientBoostingClassifier, balanced
classes, LO weights p/(1-p) out-of-fold (5 folds), mean 1.  Criteria: 1 POI pulls < 2; 2 pT(H), leading-jet
pT, dR(Z,H), y(H), m_ll; 2b the features; 3 even/odd NLO halves both ways.  Run by rwt/run.sh (scikit-learn
of the LCG_106 view).

    ggzh_reweight_bdt.py --source lo-cms.npz --target nlo-cp5.npz --out bdt.json --fig-dir figs --weights-out w.npz

--jet-feature (added to the plan 2026-10-11, before its run) adds the leading-jet pT as a sixth input; it then moves from
the test variables (criterion 2) to the features (2b).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ggzh_reweight_common as C  # noqa: E402
from ggzh_reweight_common import SCALE7, CLOSURE_VARS, OTHER_VARS, Sample, dumps  # noqa: E402
from ggzh_reweight_2d import shapes, test_vars  # noqa: E402
from ggzh_reweight_closure import test_pois  # noqa: E402

TEST_VARS = {k: CLOSURE_VARS[k] for k in ("pth", "ptj1", "drzh", "yh", "lhe_mll")}
FEATURE_VARS = {"ptz": OTHER_VARS["ptz"], "njets": OTHER_VARS["njets"], "mzh": CLOSURE_VARS["mzh"], "ptzh": CLOSURE_VARS["ptzh"],
                "ayzh": (np.arange(0.0, 3.0 + 1e-9, 0.25), "|y(ZH)|", "rapidity of the Z+H system")}
FEATURES = ["mzh", "ptzh", "ptz", "njets", "ayzh"]
CLF = dict(learning_rate=0.05, max_iter=400, max_leaf_nodes=31, min_samples_leaf=100, l2_regularization=1.0,
           early_stopping=True, validation_fraction=0.2, n_iter_no_change=20, random_state=1)
NFOLD = 5
LABELS = {"official": "CMS-style LO+PS", "reweighted": "CMS-style LO x rwt + PS", "nlo": "NLO+PS"}


def prepare(s: Sample):
    s.cols["ayzh"] = np.abs(s["yzh"].astype(float))
    s.cols["njets"] = np.minimum(s["htxs_njets30"].astype(float), 3)
    return np.stack([s[k].astype(float) for k in FEATURES], axis=1)


def fit(Xs, Xt, wt):
    from sklearn.ensemble import HistGradientBoostingClassifier
    wt = np.clip(wt, 0.0, None)
    X = np.vstack([Xs, Xt]); y = np.concatenate([np.zeros(len(Xs)), np.ones(len(Xt))])
    n = len(X) / 2.0
    sw = np.concatenate([np.full(len(Xs), n / len(Xs)), wt * n / wt.sum()])
    clf = HistGradientBoostingClassifier(**CLF)
    clf.fit(X, y, sample_weight=sw)
    return clf


def lo_weights(src: Sample, Xs, ssel, tgt: Sample, Xt, tsel, log, tag):
    idx = np.flatnonzero(ssel)
    fold = np.arange(len(idx)) % NFOLD
    w = np.ones(src.n); info = []
    for k in range(NFOLD):
        t0 = time.time()
        clf = fit(Xs[idx[fold != k]], Xt[tsel], tgt.w0[tsel])
        p = clf.predict_proba(Xs[idx[fold == k]])[:, 1]
        w[idx[fold == k]] = p / (1.0 - p)
        info.append({"fold": k, "n_iter": int(clf.n_iter_), "seconds": round(time.time() - t0, 1)})
        log(f"  [{tag}] fold {k}: {clf.n_iter_} iterations, {time.time() - t0:.0f} s")
    w[ssel] /= w[ssel].mean()
    return w, info


def evaluate(src, w, ssel, tgt, ttest, tag, log):
    allv = {**TEST_VARS, **FEATURE_VARS}
    band = [tgt.weight(wid) for wid in SCALE7]
    target = shapes(tgt, tgt.w0, ttest, allv, band)
    res = {"tag": tag, "n_eff": float(w[ssel].sum() ** 2 / np.sum(w[ssel] ** 2)),
           "weight_quantiles": dict(zip(["0", "0.1", "1", "50", "99", "99.9", "100"], np.percentile(w[ssel], [0, 0.1, 1, 50, 99, 99.9, 100]).tolist()))}
    for name, ws in (("reweighted", w), ("unweighted", np.ones(src.n))):
        sh = shapes(src, ws, ssel, allv)
        t2, ok2 = test_vars({k: sh[k] for k in TEST_VARS}, {k: target[k] for k in TEST_VARS}, TEST_VARS)
        t2b, ok2b = test_vars({k: sh[k] for k in FEATURE_VARS}, {k: target[k] for k in FEATURE_VARS}, FEATURE_VARS)
        st = test_pois(src, ws, ssel, tgt, tgt.w0, ttest)
        res[name] = {"stage13": st, "variables": t2, "features": t2b, "pass_1": st["pass"], "pass_2": ok2, "pass_2b": ok2b}
        log(f"  [{tag}] {name:10s} crit 1 {'PASS' if st['pass'] else 'FAIL'} (max |pull| {st['max_abs_pull']:.2f}); "
            f"crit 2 {'PASS' if ok2 else 'FAIL'}: " + ", ".join(f"{k} {t2[k]['chi2']:.0f}/{t2[k]['ndof']} out {t2[k]['bins_outside_band_2sigma']}" for k in TEST_VARS)
            + f"; crit 2b {'PASS' if ok2b else 'FAIL'}: " + ", ".join(f"{k} {t2b[k]['chi2']:.0f}/{t2b[k]['ndof']} out {t2b[k]['bins_outside_band_2sigma']}" for k in FEATURE_VARS))
    res["pass"] = bool(res["reweighted"]["pass_1"] and res["reweighted"]["pass_2"] and res["reweighted"]["pass_2b"])
    return res


def figure(path, res, title):
    from ggzh_reweight_plots import figure as fig_grid, plt, shape_panel
    saved = dict(C.LABELS); C.LABELS.update(LABELS)
    names = list(FEATURE_VARS) + list(TEST_VARS)
    fig, pairs = fig_grid(4, 3)
    for (ax, rax), name in zip(pairs, names):
        sec = "features" if name in FEATURE_VARS else "variables"
        v = res["reweighted"][sec][name]; u = res["unweighted"][sec][name]
        ser = {"nlo": (v["f_target"], v["e_target"]), "official": (u["f_test"], u["e_test"]), "reweighted": (v["f_test"], v["e_test"])}
        label = (FEATURE_VARS.get(name) or TEST_VARS[name])[1]
        shape_panel(ax, rax, v["edges"], ser, "nlo", band=(v["band_lo"], v["band_hi"]), xlabel=label, logy=name not in ("yh", "drzh", "njets", "ayzh"),
                    title=f"{name}{' (feature)' if name in FEATURE_VARS else ''}: p = {v['pvalue']:.2f}", flags=v["flags"], legend=(name == names[0]))
        rax.set_ylabel("/ NLO+PS")
    for ax, rax in pairs[len(names):]:
        ax.axis("off"); rax.axis("off")
    fig.suptitle(title, x=0.01, y=0.92, ha="left")
    fig.savefig(path, bbox_inches="tight"); plt.close(fig)
    C.LABELS.clear(); C.LABELS.update(saved)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--target", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--fig-dir", type=Path, default=None)
    ap.add_argument("--weights-out", type=Path, default=None)
    ap.add_argument("--jet-feature", action="store_true", help="leading-jet pT as a sixth input (step 4)")
    args = ap.parse_args(argv)
    if args.jet_feature:
        FEATURES.append("ptj1"); FEATURE_VARS["ptj1"] = TEST_VARS.pop("ptj1")
    t0 = time.time(); lines = []
    def log(s):
        print(s, flush=True); lines.append(s)
    import sklearn
    src, tgt = Sample("source", args.source), Sample("target", args.target)
    Xs, Xt = prepare(src), prepare(tgt)
    ss, ts = src.in_window, tgt.in_window
    idx = np.arange(tgt.n); A, B = ts & (idx % 2 == 0), ts & (idx % 2 == 1)
    log(f"sklearn {sklearn.__version__}; source {ss.sum()} events, target {ts.sum()} (half A {A.sum()}, half B {B.sum()})")
    out = {"plan": "rwt/README.md", "features": FEATURES, "classifier": CLF, "nfold": NFOLD,
           "inputs": {"source": str(args.source), "target": str(args.target)}}
    weights = {}
    for tag, train, test in (("full", ts, ts), ("A_to_B", A, B), ("B_to_A", B, A)):
        w, info = lo_weights(src, Xs, ss, tgt, Xt, train, log, tag)
        res = evaluate(src, w, ss, tgt, test, tag, log); res["folds"] = info
        out[tag] = res; weights[f"w_{tag}"] = w
    out["summary"] = {"1_stage13": out["full"]["reweighted"]["pass_1"], "2_variables": out["full"]["reweighted"]["pass_2"],
                      "2b_features": out["full"]["reweighted"]["pass_2b"], "3_stability": bool(out["A_to_B"]["pass"] and out["B_to_A"]["pass"]),
                      "failing_full": [k for k in TEST_VARS if not out["full"]["reweighted"]["variables"][k]["pass"]]
                                      + [k for k in FEATURE_VARS if not out["full"]["reweighted"]["features"][k]["pass"]]}
    out["summary"]["all_pass"] = bool(out["summary"]["1_stage13"] and out["summary"]["2_variables"] and out["summary"]["2b_features"] and out["summary"]["3_stability"])
    log(f"summary: {out['summary']}  ({time.time() - t0:.0f} s)")
    out["log"] = lines
    args.out.parent.mkdir(parents=True, exist_ok=True); args.out.write_text(dumps(out))
    if args.weights_out:
        np.savez_compressed(args.weights_out, run=src["run"], event=src["event"], in_window=ss, **weights)
    if args.fig_dir:
        args.fig_dir.mkdir(parents=True, exist_ok=True)
        figure(args.fig_dir / "closure-bdt.png", out["full"], "classifier reweighting: CMS-style LO+PS reweighted to NLO+PS (NLO training = test events)")
        figure(args.fig_dir / "closure-bdt-halfA-to-B.png", out["A_to_B"], "classifier trained with NLO half A, tested against half B")
        figure(args.fig_dir / "closure-bdt-halfB-to-A.png", out["B_to_A"], "classifier trained with NLO half B, tested against half A")
        print(f"figures in {args.fig_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
