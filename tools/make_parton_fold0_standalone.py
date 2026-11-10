"""Embed an existing, checksum-verified single-fold report; read no event data."""

import argparse
import base64
import hashlib
import html
import json
from pathlib import Path


def checked(directory, name):
    path = directory / name
    if path.stat().st_size > 8 * 1024**2:
        raise ValueError(f"Not a small report artifact: {name}")
    data = path.read_bytes()
    expected = Path(str(path) + ".sha256").read_text().split()[0]
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"Checksum mismatch: {name}")
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists() or Path(str(args.out) + ".sha256").exists():
        raise FileExistsError("Use a new report name; existing products are immutable")

    document = checked(args.input_dir, "index.html").decode("utf-8")
    raw_result = checked(args.input_dir, "result.json")
    result = json.loads(raw_result)
    if result["status"] != "SINGLE_FOLD_EXPLORATORY":
        raise ValueError("Expected the single-fold exploratory report")
    if result["checkpoint_reload"] != "PASS":
        raise ValueError("The plotted checkpoint must have passed its reload check")

    for name in ("pth", "ptj1", "ptz", "lhe_mll"):
        reference = f'src="{name}.png"'
        if document.count(reference) != 1:
            raise ValueError(f"Unexpected image reference: {name}")
        encoded = base64.b64encode(checked(args.input_dir, name + ".png")).decode("ascii")
        document = document.replace(reference, f'alt="{name} closure" src="data:image/png;base64,{encoded}"')
    download = base64.b64encode(raw_result).decode("ascii")
    document = document.replace('href="result.json"',
        f'download="result.json" href="data:application/json;base64,{download}"')

    stage = result["parton"]["reweighted"]["stage13"]
    rows = []
    for name, pull in stage["pulls"].items():
        rows.append(
            f"<tr><td>{html.escape(name)}</td>"
            f'<td>{stage["poi_frac"][name]:.6g}</td>'
            f'<td>{stage["target_poi_frac"][name]:.6g}</td>'
            f'<td>{pull:.4f}</td><td>{"PASS" if abs(pull) < 2 else "FAIL"}</td></tr>'
        )
    mass = result["mll_15_40"]["parton"]
    explanation = (
        '<section class="context"><h2>What failed?</h2>'
        '<p><strong>These are physics-closure tests of an accepted checkpoint.</strong> '
        'The checkpoint reload check passed. Separate training attempts that failed '
        'their reload checks are excluded from this report.</p>'
        '<ul><li>Criterion 1 requires every one of the eight STXS POI fractions '
        f'to have |pull| &lt; 2. The largest |pull| here is {stage["max_abs_pull"]:.4f}.</li>'
        '<li>Criterion 2 requires each tested shape to have chi-squared p &gt; 0.05 '
        '<em>and</em> zero merged bins beyond the NLO scale envelope by more than '
        'two MC standard deviations. Either failure is enough; zero bins outside '
        'the band alone does not imply a pass.</li>'
        '<li>Criterion 3 tests the two independent target-half directions. '
        '<strong>Pending</strong> here; it has not been evaluated.</li></ul>'
        f'<p>Source: {result["source_test_events"]:,} held-out LO events, '
        f'effective sample size {result["parton"]["n_eff"]:,.1f}. '
        f'Target: {result["target_events"]:,} NLO events. '
        f'{html.escape(result["target_scope"])}.</p>'
        f'<p>{html.escape(result["interpretation"])}</p>'
        '<h3>Eight STXS POI fractions</h3><table><tr><th>Region</th>'
        '<th>LO × ParT fraction</th><th>NLO fraction</th><th>Pull</th><th>Pass</th></tr>'
        + ''.join(rows) + '</table>'
        f'<p>Integrated m_ll = 15–40 GeV fraction ratio: {mass["ratio"]:.4f} '
        f'± {mass["error"]:.4f} (MC). This integrated ratio does not replace '
        'the differential mass-shape test below.</p>'
        '<p>The red crosses on the ratio plots mark merged bins that exceed '
        'the target scale band by more than two MC standard deviations. '
        'All plots below are embedded in this file.</p></section>'
    )
    document = document.replace('<h2>ParT</h2>', explanation + '<h2>ParT</h2>', 1)
    style = ('<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<style>body{font:16px/1.5 system-ui,sans-serif;max-width:1050px;'
        'margin:2em auto;padding:0 1em;color:#18212b;background:white}'
        'img{max-width:100%;height:auto}table{border-collapse:collapse;'
        'font-size:.9em}td,th{border:1px solid #bbb;padding:.35em .6em;'
        'text-align:left}pre{overflow:auto;background:#f3f5f7;padding:1em}'
        '.context{border-left:4px solid #b86814;padding-left:1em}'
        'h1,h2,h3{line-height:1.2}a{color:#075fa8}</style>')
    document = document.replace('<body>', style + '<body>', 1)
    data = document.encode("utf-8")
    with args.out.open("xb") as output:
        output.write(data)
    checksum = hashlib.sha256(data).hexdigest()
    with Path(str(args.out) + ".sha256").open("x") as output:
        output.write(checksum + "\n")
    print(json.dumps(dict(output=str(args.out), bytes=len(data), sha256=checksum,
                         source_selection=result["selection_sha256"],
                         status=result["status"])))


if __name__ == "__main__":
    main()
