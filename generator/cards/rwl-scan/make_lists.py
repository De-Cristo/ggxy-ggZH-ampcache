#!/usr/bin/env python3
"""Weight lists of increasing size for the amplitude-cache timing scan (production/paper.sh run
cache-scan).  Every list is the first N weights of one master list, so a longer list only adds
weights:

   9   scale points (mu_R, mu_F) with PDF4LHC21                      ids 1001-1009
  52   + PDF4LHC21_40_pdfas, members 0-42 (LHAPDF 93300-93342)         ids 3000-3042
 155   + NNPDF3.1 NNLO Hessian + alpha_s (325300-325402)              ids 2000-2102  (= ../pwg-rwl-pdf4lhc21.dat)
 256   + NNPDF4.0 NNLO, members 0-100 (331100-331200)                  ids 4000-4100
 315   + CT18NNLO, members 0-58 (14000-14058)                           ids 5000-5058
 380   + MSHT20nnlo_as118, members 0-64 (27400-27464)                   ids 6000-6064
 439   + CT18ZNNLO, members 0-58 (14100-14158)                          ids 7000-7058
 500   + NNPDF3.1 NNLO replicas (303600-...)                            ids 8000-...
Stock POWHEG-BOX-V2 accepts at most 500 weights.  The lists beyond 155 are for timing only: a list
cut inside a PDF set (e.g. 300) is not meant for PDF uncertainties; those groups say combine='none'.

    python3 make_lists.py 9 52 155 300 450      -> pwg-rwl-scan-<N>.dat in this directory
"""
import sys
from pathlib import Path

SCALE = [(1.0, 1.0), (1.0, 2.0), (1.0, 0.5), (2.0, 1.0), (2.0, 2.0), (2.0, 0.5), (0.5, 1.0), (0.5, 2.0), (0.5, 0.5)]
SETS = [("PDF4LHC21", 3000, 93300, 43, "hessian"), ("NNPDF31_hessian", 2000, 325300, 103, "hessian"),
        ("NNPDF40", 4000, 331100, 101, "none"), ("CT18NNLO", 5000, 14000, 59, "none"),
        ("MSHT20", 6000, 27400, 65, "none"), ("CT18ZNNLO", 7000, 14100, 59, "none"),
        ("NNPDF31_replicas", 8000, 303600, 101, "none")]
f = lambda x: {1.0: "1d0", 2.0: "2d0", 0.5: "0.5d0"}[x]


def master():
    out = [("scale_variation", "envelope", f"<weight id='{1001 + i}'> lhapdf=93300 renscfact={f(r)} facscfact={f(m)} </weight>")
           for i, (r, m) in enumerate(SCALE)]
    for name, wid, lha, n, comb in SETS:
        out += [(name, comb, f"<weight id='{wid + k}'> lhapdf={lha + k} </weight>") for k in range(n)]
    return out


def write(n):
    lines, group = ["<initrwgt>"], None
    for g, comb, w in master()[:n]:
        if g != group:
            if group:
                lines.append("</weightgroup>")
            lines.append(f"<weightgroup name='{g}' combine='{comb}' >"); group = g
        lines.append(w)
    lines += ["</weightgroup>", "</initrwgt>"]
    p = Path(__file__).resolve().parent / f"pwg-rwl-scan-{n}.dat"
    p.write_text("\n".join(lines) + "\n")
    return p


if __name__ == "__main__":
    for a in sys.argv[1:]:
        n = int(a)
        if not 0 < n <= 500:
            raise SystemExit("1 <= N <= 500 (stock POWHEG-BOX-V2)")
        print(write(n).name, n)
