#!/usr/bin/env python3
"""PDF and alpha_s uncertainties at fixed order from one ggxy run per PDF member (run/fo_members.sh).

    fo_pdf.py --members <job>/out/m* ... [--nlep 3] [--pdf-nmem 40] [--label ...] --out fo-pdf.json

Each m<k>/ holds the histogram file of member k (one seed, the same for all members).  Inclusive
cross section at the central scale per member (analysis/fo_combine.py); then, as for the event
weights (analysis/lhe_xsec.py):
  PDF      Hessian: sqrt(sum_{k=1..nmem} (sigma_k - sigma_0)^2) / sigma_0
  alpha_s  half the difference of the two members after the Hessian ones, relative to sigma_0
           (PDF4LHC21_40_pdfas: alpha_s(m_Z) = 0.117 / 0.119, i.e. delta = 0.001)
The members are separate Monte Carlo integrations, so their statistical errors add to the Hessian
sum: "mc_noise" = sqrt(sum_k (err_k^2 + err_0^2)) / sigma_0, an upper bound of that contribution
(it enters in quadrature).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fo_combine import inclusive  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--members", type=Path, nargs="+", required=True)
    ap.add_argument("--nlep", type=float, default=3.0)
    ap.add_argument("--pdf-nmem", type=int, default=40)
    ap.add_argument("--mode", default="lo")
    ap.add_argument("--label", default="")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    sig = {}
    for d in a.members:
        k = int(d.name[1:])
        files = sorted(d.glob(f"histo_{a.mode}_*.txt"))
        if not files:
            raise SystemExit(f"{d}: no histo_{a.mode}_*.txt")
        r = inclusive(files, a.nlep)
        sig[k] = (r["sigma_fb"], r["mc_err_fb"])
    nm = a.pdf_nmem
    need = [0] + list(range(1, nm + 3))
    missing = [k for k in need if k not in sig]
    if missing:
        raise SystemExit(f"members missing: {missing}")
    s0, e0 = sig[0]
    pdf = math.sqrt(sum((sig[k][0] - s0) ** 2 for k in range(1, nm + 1))) / s0
    noise = math.sqrt(sum(sig[k][1] ** 2 + e0 ** 2 for k in range(1, nm + 1))) / s0
    asr = 0.5 * (sig[nm + 2][0] - sig[nm + 1][0]) / s0
    out = {"label": a.label, "mode": a.mode, "nlep": a.nlep, "unit": "fb", "sigma_fb": s0, "mc_err_fb": e0,
           "pdf_hessian_rel": pdf, "pdf_mc_noise_rel": noise, "alphas_rel_half_difference": asr,
           "members": {str(k): {"sigma_fb": v[0], "mc_err_fb": v[1]} for k, v in sorted(sig.items())}}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(f"{a.label}: sigma_0 {s0:.4f} +- {e0:.4f} fb; PDF (Hessian) {100 * pdf:.2f} % (MC noise <= {100 * noise:.2f} %), "
          f"alpha_s (+-0.001) {100 * asr:+.2f} %")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
