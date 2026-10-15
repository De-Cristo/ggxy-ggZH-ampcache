#!/usr/bin/env python3
"""Plot helpers for the ggZH reweighting study: normalised shape comparisons with a ratio panel.
Series colours are fixed by identity (ggzh_reweight_common.COLORS), marks are thin step lines,
the NLO seven-point band is a light fill, errors are MC statistics."""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ggzh_reweight_common import COLORS, LABELS  # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8,
                     "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 120, "savefig.dpi": 150})


def _centers(edges):
    e = np.array(edges, dtype=float)
    c = 0.5 * (e[:-1] + e[1:])
    if np.isinf(e[-1]):
        c[-1] = e[-2] * 1.25
    return c


def shape_panel(ax, rax, edges, series, ref, band=None, xlabel="", logy=True, title="", flags=None, legend=True):
    """series: dict name -> (f, err); ref: name of the reference (ratio denominator);
    band: (lo, hi) arrays for the reference's scale band; flags: bool array of bins failing a test."""
    e = np.array(edges, dtype=float)
    if np.isinf(e[-1]):
        e[-1] = e[-2] + (e[-2] - e[-3])
    x = np.append(e[:-1], e[-1])
    fr, er = series[ref]
    with np.errstate(divide="ignore", invalid="ignore"):
        if band is not None:
            lo, hi = band
            ax.fill_between(x, np.append(lo, lo[-1]), np.append(hi, hi[-1]), step="post", color=COLORS["band"], alpha=0.15, lw=0,
                            label="NLO+PS 7-point scale band")
            rax.fill_between(x, np.append(lo / fr, (lo / fr)[-1]), np.append(hi / fr, (hi / fr)[-1]), step="post", color=COLORS["band"], alpha=0.15, lw=0)
        for name, (f, err) in series.items():
            c = COLORS[name]
            ax.step(x, np.append(f, f[-1]), where="post", color=c, lw=1.4 if name == ref else 1.1, label=LABELS[name])
            ax.errorbar(_centers(e), f, yerr=err, fmt="none", ecolor=c, elinewidth=0.6, capsize=0)
            if name != ref:
                r = f / fr
                rerr = np.sqrt(err ** 2 + (f * er / fr) ** 2) / fr
                rax.errorbar(_centers(e), r, yerr=rerr, fmt="o", ms=2.2, color=c, elinewidth=0.6, capsize=0, label=f"{LABELS[name]} / {LABELS[ref]}")
        if flags is not None and flags.any():
            cx = _centers(e)[flags]
            rax.plot(cx, np.full(len(cx), 1.0), marker="x", ls="none", color="#e34948", ms=5, mew=1.2, label="outside band (> 2 sigma)")
    rax.axhline(1.0, color="#52514e", lw=0.8)
    ax.set_xlim(e[0], e[-1]); rax.set_xlim(e[0], e[-1])
    if logy:
        ax.set_yscale("log")
        pos = np.concatenate([f[f > 0] for f, _ in series.values()])
        ax.set_ylim(max(pos.min() * 0.5, 1e-7), pos.max() * 3)
    rax.set_ylim(0.0, 2.5)
    rax.set_xlabel(xlabel); ax.set_ylabel("fraction of events / bin"); rax.set_ylabel("ratio")
    ax.text(0.02, 0.97, title, transform=ax.transAxes, ha="left", va="top", fontsize=10, fontweight="bold")
    if legend:
        ax.legend(loc="lower left", frameon=False)


def figure(nrows, ncols, width=4.2, height=3.8):
    """nrows x ncols panels, each a (main, ratio) pair sharing x; a clear gap between rows."""
    fig = plt.figure(figsize=(width * ncols, height * nrows))
    outer = fig.add_gridspec(nrows, ncols, hspace=0.45, wspace=0.3)
    pairs = []
    for r in range(nrows):
        for c in range(ncols):
            inner = outer[r, c].subgridspec(2, 1, height_ratios=[3, 1.2], hspace=0.07)
            ax = fig.add_subplot(inner[0]); rax = fig.add_subplot(inner[1], sharex=ax)
            ax.tick_params(labelbottom=False)
            pairs.append((ax, rax))
    return fig, pairs
