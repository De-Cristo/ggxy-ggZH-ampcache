#!/usr/bin/env python3
"""Cross section and seven-point band of one MCFM run (mcfm/run_mcfm.sh output directory).

Reads the last 'Value of integral is X +- Y fb' of stdout.log and the scale envelope from the
*total_cross_scalevar_maximum/minimum.txt files (absolute shifts in fb).

    mcfm_xsec.py <workdir>/out/ggzh-lo-stable-nnpdf31 [--out mcfm-stable-nnpdf31.json]
"""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path


def total(path: str) -> float:
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            return float(line.split()[2])
    raise SystemExit(f"no value in {path}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir", type=Path)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    vals = re.findall(r"Value of integral is\s+(\S+)\s+\S+\s+(\S+)\s+fb", (a.dir / "stdout.log").read_text(errors="replace"))
    if not vals:
        raise SystemExit(f"no 'Value of integral' in {a.dir}/stdout.log")
    sig, err = float(vals[-1][0]), float(vals[-1][1])
    up = glob.glob(str(a.dir / "*total_cross_scalevar_maximum.txt")); dn = glob.glob(str(a.dir / "*total_cross_scalevar_minimum.txt"))
    out = {"dir": str(a.dir), "sigma_fb": sig, "mc_err_fb": err}
    if up and dn:
        out["band_rel"] = [total(up[0]) / sig, total(dn[0]) / sig]
    br = re.findall(r"brzee\s*=\s*(\S+)", (a.dir / "stdout.log").read_text(errors="replace"))
    if br:
        out["brzee"] = float(br[-1])
    if a.out:
        a.out.write_text(json.dumps(out, indent=1) + "\n")
    band = f", band {100 * out['band_rel'][0]:+.1f}/{100 * out['band_rel'][1]:+.1f} %" if "band_rel" in out else ""
    print(f"{a.dir.name}: {sig:.3f} +- {err:.3f} fb{band}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
