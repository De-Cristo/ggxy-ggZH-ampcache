"""Figure style for Physical Review D, shared by every figure script of the paper.

Rules (APS Journals Style Guide for Authors, February 2026, figures): at the printed size the lettering
height and symbol width are at least 2 mm, data points at least 1 mm in diameter, lines at least 0.5 pt;
axis units in parentheses, e.g. "p_T(Z) (GeV)"; integer axis ticks where possible; avoid small open
symbols.  Column widths 8.6 cm (3.375 in) and 17.8 cm (7.0 in).

How they are met here:
  * figures are made at their printed width and included at natural size (single column) or at
    \\textwidth (double column, figure*), so the sizes below are the printed sizes;
  * STIX capitals are 0.662 em high, so 2 mm needs at least 8.6 pt: every text element is 9 pt
    (capitals 2.1 mm); sub- and superscripts are smaller by construction;
  * markers at least MS = 3.5 pt (1.2 mm), filled; lines and ticks at least 0.5 pt;
  * fonts embedded as TrueType (pdf.fonttype 42), no Type 3.

    import prd_style; prd_style.apply(plt)
"""
SINGLE, DOUBLE = 3.375, 7.0   # inches
FONT = 9.0                    # pt, every text element
MS = 3.5                      # pt, smallest marker diameter
LW = 1.0                      # pt, data lines


def apply(plt):
    plt.rcParams.update({
        "font.family": "STIXGeneral", "mathtext.fontset": "stix",
        "font.size": FONT, "axes.labelsize": FONT, "axes.titlesize": FONT, "legend.fontsize": FONT,
        "legend.title_fontsize": FONT, "xtick.labelsize": FONT, "ytick.labelsize": FONT, "figure.titlesize": FONT,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": 0.6, "lines.linewidth": LW, "lines.markersize": MS, "errorbar.capsize": 0,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.size": 3.5, "ytick.major.size": 3.5, "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
        "legend.frameon": False, "legend.handlelength": 2.4, "savefig.dpi": 600})


def unit(label, u):
    """Axis label with the unit in parentheses (APS style): unit("$p_T(Z)$", "GeV") -> "$p_T(Z)$ (GeV)"."""
    return f"{label} ({u})" if u else label
