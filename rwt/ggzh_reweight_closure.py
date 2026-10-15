#!/usr/bin/env python3
"""Closure tests of the binned pT(V) x N_jets reweighting map (ggzh_reweight_map.py), with the pass criteria
fixed before any result was looked at (rwt/README.md):

  1. Stage 1.3: the 8 POI fractions of the reweighted source sample (FWDH excluded, 15 < m_ll < 150)
     agree with NLO+PS, every pull |f_rw - f_NLO| / sqrt(df_rw^2 + df_NLO^2) < 2 (MC statistics only).
  2. Variables the map does not use - pT(H), pT(ZH), m_ZH, leading-jet pT, y_H, dR(Z,H), m_ll -
     bin by bin vs NLO+PS: chi2 probability > 5 % in every distribution (MC statistics; bins merged
     to >= 20 raw events in both samples), and no bin outside the NLO scale band (min/max over the
     nominal and the seven variations of SCALE7) by more than 2 sigma of the MC statistics.
  3. Stability: map built on one half of NLO+PS (even events), tests 1-2 against the other half, and
     the same with the halves exchanged; both directions must pass.
The same tests are also evaluated for the unweighted source sample and, if given, a reference LO+PS
sample.  Also used as a library: test_pois (criterion 1) by the classifier reweighting.

    ggzh_reweight_closure.py --official <source>.npz --nlo <target>.npz --out closure.json [--fig-dir figs]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ggzh_reweight_common import (SCALE7, CELL_LABELS, CLOSURE_VARS, COLORS, LABELS, N_JET, N_PT, NTUPLE_DIR, POIS, PT_EDGES,  # noqa: E402
                                  Sample, cell_sums, chi2_shape, dumps, fractions, hist,
                                  merge_sparse_bins, poi_sums, regroup)
from ggzh_reweight_map import FIRST_FIT_BIN, apply_map, build_map  # noqa: E402

MIN_RAW = 20


# ----------------------------------------------------------------------------------------------- tests
def shapes(sample: Sample, w, sel, band_weights=None):
    """Per closure variable: fractions, errors, raw counts (in range), and the band if band_weights
    (list of weight arrays for the variations) is given."""
    out = {}
    for name, (edges, label, note) in CLOSURE_VARS.items():
        x = sample[name].astype(float)
        s, s2, c = hist(x[sel], w[sel], edges)
        f, e = fractions(s, s2)
        d = {"f": f, "e": e, "raw": c}
        if band_weights is not None:
            fs = [f] + [fractions(*hist(x[sel], bw[sel], edges)[:2])[0] for bw in band_weights]
            fs = np.array(fs); d["variations"] = fs; d["lo"], d["hi"] = fs.min(axis=0), fs.max(axis=0)
        out[name] = d
    return out


def test_variables(test: dict, target: dict):
    """Criterion 2 for one test sample against the target (which carries the band)."""
    res, all_pass = {}, True
    for name in CLOSURE_VARS:
        t, g = test[name], target[name]
        groups = merge_sparse_bins(CLOSURE_VARS[name][0], [t["raw"], g["raw"]], MIN_RAW)
        ft, fg = regroup(t["f"], groups), regroup(g["f"], groups)
        et, eg = np.sqrt(regroup(t["e"] ** 2, groups)), np.sqrt(regroup(g["e"] ** 2, groups))
        fv = np.array([regroup(v, groups) for v in g["variations"]])      # band of the merged bins
        lo, hi = fv.min(axis=0), fv.max(axis=0)
        c2, ndof, p = chi2_shape(ft, et, fg, eg)
        sig = np.sqrt(et ** 2 + eg ** 2)
        with np.errstate(invalid="ignore"):
            flags = (ft > hi + 2 * sig) | (ft < lo - 2 * sig)
        pulls = np.where(sig > 0, (ft - fg) / np.where(sig > 0, sig, 1), 0.0)
        ok = (p > 0.05) and not flags.any()
        all_pass &= ok
        res[name] = {"chi2": c2, "ndof": ndof, "pvalue": p, "n_bins_merged": len(groups), "bins_outside_band_2sigma": int(flags.sum()),
                     "max_abs_pull": float(np.max(np.abs(pulls))), "pass": bool(ok), "groups": groups,
                     "f_test": ft, "e_test": et, "f_target": fg, "e_target": eg, "band_lo": lo, "band_hi": hi, "flags": flags}
    return res, bool(all_pass)


def test_pois(sample: Sample, w, sel, target: Sample, wt, selt):
    s, s2, _ = cell_sums(sample, w, sel); f, e = fractions(s, s2)
    st, st2, _ = cell_sums(target, wt, selt); ft, et = fractions(st, st2)
    pf, pft = poi_sums(f), poi_sums(ft)
    pe = {p: float(np.sqrt(sum(e[q * 3 + j] ** 2 for j in jets))) for p, (q, jets) in POIS.items()}
    pet = {p: float(np.sqrt(sum(et[q * 3 + j] ** 2 for j in jets))) for p, (q, jets) in POIS.items()}
    pulls = {p: (pf[p] - pft[p]) / np.sqrt(pe[p] ** 2 + pet[p] ** 2) if pe[p] ** 2 + pet[p] ** 2 > 0 else 0.0 for p in POIS}
    cells_pull = (f - ft) / np.sqrt(e ** 2 + et ** 2 + 1e-30)
    ok = all(abs(v) < 2 for v in pulls.values())
    return {"poi_frac": pf, "poi_err": pe, "target_poi_frac": pft, "target_poi_err": pet, "pulls": pulls, "max_abs_pull": float(max(abs(v) for v in pulls.values())),
            "cell_frac": f, "cell_err": e, "target_cell_frac": ft, "target_cell_err": et, "cell_pulls": cells_pull, "pass": bool(ok)}


# ----------------------------------------------------------------------------------------------- main
def run_tests(official: Sample, nlo: Sample, lo: Sample | None, nlo_map_sel, nlo_test_sel, tag: str, log):
    """Build the 8 maps on nlo_map_sel, reweight the official sample, run criteria 1 and 2 against nlo_test_sel."""
    maps = {wid: build_map(official, nlo, wid, nlo_sel=nlo_map_sel) for wid in ["1001"] + SCALE7}
    w = {wid: apply_map(m, official)[0] for wid, m in maps.items()}
    dw = apply_map(maps["1001"], official)[1]
    sel_o = official.in_window
    target = shapes(nlo, nlo.w0, nlo_test_sel, band_weights=[nlo.weight(wid) for wid in SCALE7])
    res = {"tag": tag, "maps": {wid: m.to_json() for wid, m in maps.items()},
           "nlo_events_map": int(nlo_map_sel.sum()), "nlo_events_test": int(nlo_test_sel.sum()),
           "n_eff_in_window": float(w["1001"][sel_o].sum() ** 2 / np.sum(w["1001"][sel_o] ** 2)),
           "weight_quantiles": dict(zip(["0", "1", "5", "50", "95", "99", "100"], np.percentile(w["1001"][sel_o], [0, 1, 5, 50, 95, 99, 100]).tolist())),
           "mean_rel_map_stat_err": float(np.mean(dw[sel_o] / w["1001"][sel_o]))}
    samples = {"reweighted": (official, w["1001"], sel_o), "official": (official, np.ones(official.n), sel_o)}
    if lo is not None:
        samples["lo"] = (lo, lo.w0, lo.in_window)
    res["stage13"], res["variables"] = {}, {}
    for name, (s, ws, sel) in samples.items():
        res["stage13"][name] = test_pois(s, ws, sel, nlo, nlo.w0, nlo_test_sel)
        sh = shapes(s, ws, sel)
        res["variables"][name], ok = test_variables(sh, target)
        res["variables"][name]["_all_pass"] = ok
    res["target_shapes"] = target
    res["pass_1_stage13"] = res["stage13"]["reweighted"]["pass"]
    res["pass_2_variables"] = res["variables"]["reweighted"]["_all_pass"]
    log(f"[{tag}] maps on {nlo_map_sel.sum()} NLO events, test on {nlo_test_sel.sum()}; N_eff {res['n_eff_in_window']:.0f}; "
        f"criterion 1 (POI pulls < 2): {'PASS' if res['pass_1_stage13'] else 'FAIL'} (max |pull| {res['stage13']['reweighted']['max_abs_pull']:.2f}); "
        f"criterion 2: {'PASS' if res['pass_2_variables'] else 'FAIL'}")
    for name in CLOSURE_VARS:
        v = res["variables"]["reweighted"][name]; o = res["variables"]["official"][name]
        log(f"    {name:8s} reweighted: p = {v['pvalue']:.3f} (chi2 {v['chi2']:.1f}/{v['ndof']}), bins outside band {v['bins_outside_band_2sigma']}, max|pull| {v['max_abs_pull']:.1f} {'ok' if v['pass'] else 'FAIL'}"
            f"   | unweighted official: p = {o['pvalue']:.1e}, outside band {o['bins_outside_band_2sigma']}")
    return res, maps, w


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--official", type=Path, required=True, help="the source LO+PS sample (.npz)")
    ap.add_argument("--nlo", type=Path, required=True, help="the NLO+PS target (.npz)")
    ap.add_argument("--lo", type=Path, default=None, help="optional reference LO+PS sample (.npz), tests only")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--fig-dir", type=Path, default=None)
    ap.add_argument("--apply-to", type=Path, default=None,
                    help="also apply the nominal map built on the source to this sample and run criteria 1-2 on it against the same target")
    ap.add_argument("--source-label", default=None, help="legend label of the source sample")
    ap.add_argument("--target-label", default=None, help="legend label of the target sample")
    ap.add_argument("--reweighted-label", default=None, help="legend label of the reweighted source")
    args = ap.parse_args(argv)
    import ggzh_reweight_common as _c
    if args.source_label:
        _c.LABELS["official"] = args.source_label
        _c.LABELS["reweighted"] = args.source_label + ", reweighted"
    if args.target_label:
        _c.LABELS["nlo"] = args.target_label
    if getattr(args, "reweighted_label", None):
        _c.LABELS["reweighted"] = args.reweighted_label
    t0 = time.time()
    lines = []

    def log(s):
        print(s, flush=True); lines.append(s)

    official, nlo = Sample("official", args.official), Sample("nlo", args.nlo)
    lo = Sample("lo", args.lo) if args.lo is not None and args.lo.exists() else None
    out = {"inputs": {"official": str(args.official), "nlo": str(args.nlo), "lo": str(args.lo) if lo else None},
           "criteria": {"1": "all 8 POI pulls < 2 (MC stat)", "2": "chi2 p > 0.05 in every closure variable and no merged bin outside the NLO 8-point band by > 2 sigma",
                        "3": "criteria 1-2 with the map from one NLO half and the other half as target, both directions"},
           "min_raw_events_per_bin": MIN_RAW}
    all_sel = nlo.in_window
    # full map
    full, maps, w = run_tests(official, nlo, lo, all_sel, all_sel, "full", log)
    out["full"] = full
    # deliverable test (defined 2026-10-09 before any result of the public samples was seen): the map of
    # the source applied to another LO sample, criteria 1-2 against the same target
    if args.apply_to is not None:
        other = Sample("applied", args.apply_to)
        wo, _ = apply_map(maps["1001"], other)
        sel_x = other.in_window
        wo = np.where(sel_x, wo / np.mean(wo[sel_x]), 1.0)          # shape only: mean 1 in the window
        target = full["target_shapes"]
        app = {"path": str(args.apply_to), "events": other.n, "events_in_window": int(sel_x.sum()),
               "n_eff_in_window": float(wo[sel_x].sum() ** 2 / np.sum(wo[sel_x] ** 2)), "stage13": {}, "variables": {}}
        for name, ws in (("reweighted", wo), ("unweighted", np.ones(other.n))):
            app["stage13"][name] = test_pois(other, ws, sel_x, nlo, nlo.w0, all_sel)
            app["variables"][name], ok = test_variables(shapes(other, ws, sel_x), target)
            app["variables"][name]["_all_pass"] = ok
        app["pass_1_stage13"] = app["stage13"]["reweighted"]["pass"]
        app["pass_2_variables"] = app["variables"]["reweighted"]["_all_pass"]
        out["applied"] = app
        log(f"[applied to {args.apply_to.name}] {other.n} events; N_eff {app['n_eff_in_window']:.0f}; criterion 1: "
            f"{'PASS' if app['pass_1_stage13'] else 'FAIL'} (max |pull| {app['stage13']['reweighted']['max_abs_pull']:.2f}, "
            f"unweighted {app['stage13']['unweighted']['max_abs_pull']:.2f}); criterion 2: {'PASS' if app['pass_2_variables'] else 'FAIL'}")
        for name in CLOSURE_VARS:
            v = app["variables"]["reweighted"][name]; o = app["variables"]["unweighted"][name]
            log(f"    {name:8s} reweighted: p = {v['pvalue']:.3f} (chi2 {v['chi2']:.1f}/{v['ndof']}), bins outside band {v['bins_outside_band_2sigma']} "
                f"{'ok' if v['pass'] else 'FAIL'}   | unweighted: p = {o['pvalue']:.1e}, outside band {o['bins_outside_band_2sigma']}")

    # stability
    idx = np.arange(nlo.n)
    half_a, half_b = all_sel & (idx % 2 == 0), all_sel & (idx % 2 == 1)
    ab, _, _ = run_tests(official, nlo, None, half_a, half_b, "map A -> test B", log)
    ba, _, _ = run_tests(official, nlo, None, half_b, half_a, "map B -> test A", log)
    out["stability"] = {"A_to_B": ab, "B_to_A": ba}
    out["pass_3_stability"] = bool(ab["pass_1_stage13"] and ab["pass_2_variables"] and ba["pass_1_stage13"] and ba["pass_2_variables"])
    out["summary"] = {"1_stage13": full["pass_1_stage13"], "2_variables": full["pass_2_variables"], "3_stability": out["pass_3_stability"],
                      "variables_failing": [n for n in CLOSURE_VARS if not full["variables"]["reweighted"][n]["pass"]],
                      "all_pass": bool(full["pass_1_stage13"] and full["pass_2_variables"] and out["pass_3_stability"])}
    out["log"] = lines; out["seconds"] = round(time.time() - t0, 1)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(dumps(out))
    log(f"summary: {out['summary']}  ({out['seconds']} s)")
    if args.fig_dir:
        make_figures(args.fig_dir, official, nlo, lo, full, ab, maps, w, out.get("applied"))
    return 0


def variables_figure(path, res_rw, res_ref, title, labels=None):
    """Closure variables: target, unweighted reference and reweighted series with the band and the flags."""
    from ggzh_reweight_plots import figure, plt, shape_panel
    import ggzh_reweight_common as _c
    saved = dict(_c.LABELS)
    if labels:
        _c.LABELS.update(labels)
    fig, pairs = figure(2, 4)
    names = list(CLOSURE_VARS)
    for (ax, rax), name in zip(pairs, names):
        v = res_rw[name]; o = res_ref[name]
        groups = v["groups"]; edges = CLOSURE_VARS[name][0]
        e = np.array([edges[lo_] for lo_, _ in groups] + [edges[groups[-1][1]]])
        ser = {"nlo": (v["f_target"], v["e_target"]), "official": (o["f_test"], o["e_test"]), "reweighted": (v["f_test"], v["e_test"])}
        shape_panel(ax, rax, e, ser, "nlo", band=(v["band_lo"], v["band_hi"]), xlabel=CLOSURE_VARS[name][1], logy=name not in ("yh", "drzh"),
                    title=f"{name}: p = {v['pvalue']:.2f}", flags=v["flags"], legend=(name == names[0]))
    pairs[-1][0].axis("off"); pairs[-1][1].axis("off")
    fig.suptitle(title, x=0.01, y=0.95, ha="left")
    fig.savefig(path, bbox_inches="tight"); plt.close(fig)
    _c.LABELS.clear(); _c.LABELS.update(saved)


def make_figures(fig_dir: Path, official, nlo, lo, full, ab, maps, w, applied=None):
    from ggzh_reweight_plots import figure, plt, shape_panel
    fig_dir.mkdir(parents=True, exist_ok=True)
    src, tgt = LABELS["official"], LABELS["nlo"]
    for tag, res, fname in (("full map", full, "closure-variables.png"), ("map on target half A, test on half B", ab, "closure-variables-halfA-to-B.png")):
        variables_figure(fig_dir / fname, res["variables"]["reweighted"], res["variables"]["official"],
                         f"closure: {src} reweighted to {tgt} ({tag}); bins merged to >= {MIN_RAW} raw events")
    if applied is not None:
        variables_figure(fig_dir / "closure-variables-applied.png", applied["variables"]["reweighted"], applied["variables"]["unweighted"],
                         f"deliverable test: the map of {src} applied to {Path(applied['path']).stem}, vs {tgt}",
                         labels={"official": f"{Path(applied['path']).stem} (unweighted)", "reweighted": f"{Path(applied['path']).stem}, reweighted"})
    # the map
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    m = maps["1001"]
    centers = np.array([0.5 * (PT_EDGES[i] + PT_EDGES[i + 1]) if np.isfinite(PT_EDGES[i + 1]) else 700.0 for i in range(N_PT)])
    for j, ax in enumerate(axes):
        b = np.arange(N_PT) * N_JET + j
        ax.errorbar(centers, m.r[b], yerr=m.dr[b], fmt="o", ms=3, color=COLORS["nlo"], label="binned ratio f_NLO / f_official")
        for wid in SCALE7:
            mk = maps[wid]
            x = np.linspace(np.log(250.0), np.log(1000.0), 100)
            ax.plot(np.exp(x), np.exp(np.polyval(mk.fits[j].coef, x)), color=COLORS["band"], alpha=0.25, lw=0.8)
        x = np.linspace(np.log(250.0), np.log(1000.0), 100)
        f = m.fits[j]
        ax.plot(np.exp(x), np.exp(np.polyval(f.coef, x)), color=COLORS["reweighted"], lw=1.6, label=f"fit above 250 GeV: degree {f.degree}, p = {f.pvalue:.2f}")
        ax.axvline(250.0, color="#52514e", lw=0.6, ls="--")
        ax.set_yscale("log"); ax.set_xlim(0, 1000); ax.set_xlabel("pT(Z) [GeV]"); ax.set_title(["0 jets", "1 jet", ">= 2 jets"][j] + " (HTXS N_jets30)", loc="left")
        ax.legend(frameon=False, fontsize=7, loc="lower left")
    axes[0].set_ylabel("weight (before mean-1 scaling)")
    fig.suptitle("reweighting map: NLO+PS / official per pT(Z) x N_jets bin; thin lines = the seven scale-variation maps", x=0.01, y=1.02, ha="left")
    fig.savefig(fig_dir / "map.png", bbox_inches="tight"); plt.close(fig)
    # POI pulls and cells
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 3.6), gridspec_kw={"width_ratios": [1, 2]})
    st = full["stage13"]
    for i, k in enumerate(("reweighted", "official")):
        ax1.bar(np.arange(8) + (i - 0.5) * 0.38, [st[k]["pulls"][p] for p in POIS], width=0.36, color=COLORS[k], label=LABELS[k])
    ax1.axhline(2, color="#e34948", lw=0.8, ls="--"); ax1.axhline(-2, color="#e34948", lw=0.8, ls="--")
    ax1.set_xticks(np.arange(8)); ax1.set_xticklabels(list(POIS), rotation=60, ha="right", fontsize=7); ax1.set_ylabel("pull vs NLO+PS (MC stat)")
    ax1.set_ylim(-12, 12); ax1.legend(frameon=False, fontsize=7); ax1.set_title("Stage 1.3 POI fractions", loc="left")
    xx = np.arange(18)
    for i, k in enumerate(("nlo", "official", "reweighted")):
        d = st["reweighted"] if k == "nlo" else st[k]
        f = d["target_cell_frac"] if k == "nlo" else d["cell_frac"]; e = d["target_cell_err"] if k == "nlo" else d["cell_err"]
        ax2.errorbar(xx + (i - 1) * 0.22, f, yerr=e, fmt="o", ms=3, color=COLORS[k], label=LABELS[k], elinewidth=0.8)
    ax2.set_yscale("log"); ax2.set_xticks(xx); ax2.set_xticklabels(CELL_LABELS, rotation=60, ha="right", fontsize=7); ax2.set_ylabel("fraction (non-FWDH, in window)")
    ax2.legend(frameon=False, fontsize=7, loc="lower left"); ax2.set_title("Stage 1.3 cells", loc="left")
    fig.savefig(fig_dir / "closure-stage13.png", bbox_inches="tight"); plt.close(fig)
    # weights
    fig, ax = plt.subplots(figsize=(5, 3.4))
    ws = w["1001"][official.in_window]
    ax.hist(np.clip(ws, 0, 3), bins=60, range=(0, 3), color=COLORS["reweighted"], histtype="step", lw=1.2)
    ax.set_xlabel("event weight (nominal map)"); ax.set_ylabel("official events"); ax.set_yscale("log"); ax.set_title(f"weights: N_eff = {full['n_eff_in_window']:.0f} of {official.in_window.sum()}", loc="left")
    fig.savefig(fig_dir / "weights.png", bbox_inches="tight"); plt.close(fig)
    print(f"figures in {fig_dir}")


if __name__ == "__main__":
    sys.exit(main())
