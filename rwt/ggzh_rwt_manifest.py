#!/usr/bin/env python3
"""Freeze and verify the available-table selection authorized on 2026-11-01.

This module reads small inventories/bookkeeping JSONs only. Event reads and training belong
to rwt/run.sh on workers. A manifest is immutable once the analysis starts; all inputs and
criteria are inherited from the original rwt10m definitions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

PACKAGE = Path(__file__).resolve().parents[1]
SCHEMA = "ggzh-rwt-available-v1"
PRODUCTION = "3161a6224c1c354da81671c6b54fe66a27f303e8"
SAMPLES = ("lo-cms", "nlo-cp5")


def digest(manifest):
    value = {k:v for k,v in manifest.items() if k != "selection_sha256"}
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load(path):
    return validate(json.loads(Path(path).read_text()))


def validate(d):
    if d["schema"] != SCHEMA or d["production_commit"] != PRODUCTION:
        raise ValueError("wrong selection schema or production pin")
    if d["kind"] not in ("production", "pilot"):
        raise ValueError("wrong selection kind")
    if d["nominal"] != nominal(d["kind"] == "pilot"):
        raise ValueError("nominal sample differs from the public production definitions")
    missing = []
    for sample in SAMPLES:
        n = d["nominal"][sample]
        full = set(range(n["first_seed"], n["first_seed"] + n["files"]))
        selected = d["included"][sample]
        if not selected or selected != sorted(set(selected)) or not set(selected) <= full:
            raise ValueError("invalid, duplicated or unordered seeds: " + sample)
        missing.extend((sample, s) for s in sorted(full - set(selected)))
        if d["expected_events"][sample] != len(selected)*n["events_per_file"]:
            raise ValueError("event-count declaration does not match selection")
    if len(missing) > (0 if d["kind"] == "pilot" else 2):
        raise ValueError("selection exceeds the authorized missing-file limit")
    if sorted(missing) != sorted((x["sample"],x["seed"]) for x in d["omitted"]):
        raise ValueError("omission record does not match selection")
    if d["selection_sha256"] != digest(d):
        raise ValueError("selection manifest hash mismatch")
    return d


def products(manifest):
    for sample in SAMPLES:
        for seed in manifest["included"][sample]:
            name = f"table-{seed}.csv.gz"
            yield f"tables/{sample}/{name}", f"tables/{sample}-rwt10m-13600/{name}"


def nominal(pilot=False):
    setname = "rwt10m-pilot" if pilot else "rwt10m"
    paper = PACKAGE/"production/paper.sh"
    lines = subprocess.check_output(["bash", str(paper), "inputs", setname, "rwt-clf-a-10m", "0"], text=True)
    seeds = {s:[] for s in SAMPLES}
    for line in lines.splitlines():
        local, _ = line.split()
        m = re.fullmatch(r"tables/(lo-cms|nlo-cp5)/table-(\d+)\.csv\.gz", local)
        if m:
            seeds[m[1]].append(int(m[2]))
    out = {}
    for sample, run in (("lo-cms","events-lo-cms-10m"), ("nlo-cp5","events-nlo-10m")):
        s = sorted(seeds[sample])
        if not s or s != list(range(s[0],s[0]+len(s))):
            raise ValueError("original run does not define a consecutive nominal sample")
        shown = subprocess.check_output(["bash",str(paper),"show",setname,run,"0"], text=True)
        n = int(re.search(r"run_job.sh events (\d+) ",shown)[1])
        out[sample] = dict(first_seed=s[0],files=len(s),events_per_file=n)
    return out


def freeze(inventory, output, pilot=False):
    if output.exists():
        raise FileExistsError("refusing to replace a frozen manifest: " + str(output))
    inv = json.loads(inventory.read_text())
    present = set(inv["products"])
    d = dict(schema=SCHEMA, kind="pilot" if pilot else "production",
             production_commit=PRODUCTION, snapshot_utc=inv["checked_utc"],
             authorization="User amendment 2026-11-01: use available tables if at most two files are missing; record actual counts. A/B/C and criteria 1-3 unchanged.",
             nominal=nominal(pilot),included={},omitted=[],expected_events={})
    reasons = inv.get("omission_reasons",{})
    for sample in SAMPLES:
        n = d["nominal"][sample]
        selected = []
        for seed in range(n["first_seed"], n["first_seed"]+n["files"]):
            product = f"tables/{sample}-rwt10m-13600/table-{seed}.csv.gz"
            if product in present and product+".sha256" in present:
                selected.append(seed)
            else:
                d["omitted"].append(dict(sample=sample,seed=seed,
                    reason=reasons.get(sample,{}).get(str(seed),"not available as a complete product/checksum pair at snapshot")))
        d["included"][sample] = selected
        d["expected_events"][sample] = len(selected)*n["events_per_file"]
    if len(d["omitted"]) > (0 if pilot else 2):
        raise ValueError(f"not ready: {len(d['omitted'])} table files unavailable")
    d["selection_sha256"] = digest(d)
    validate(d)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("x") as f:
        f.write(json.dumps(d,indent=2)+"\n")
    return d


def verify(manifest, directory):
    counts = {}
    for sample in SAMPLES:
        meta = json.loads((directory/(sample+".json")).read_text())
        names = [Path(x).name for x in meta["files"]]
        expected = [f"table-{s}.csv.gz" for s in manifest["included"][sample]]
        if names != expected:
            raise ValueError("loaded files differ from frozen selection: " + sample)
        n = manifest["nominal"][sample]["events_per_file"]
        if meta["events_per_file"] != [n]*len(expected):
            raise ValueError("wrong per-file event count: " + sample)
        if meta["events"] != manifest["expected_events"][sample]:
            raise ValueError("loaded event count differs from frozen selection: " + sample)
        if meta["nominal_column"] != ("2000" if sample == "nlo-cp5" else "nominal"):
            raise ValueError("wrong nominal weight: " + sample)
        counts[sample] = meta["events"]
    return counts


def annotate(manifest, directory, result, model=None):
    counts = verify(manifest,directory)
    data = json.loads(result.read_text())
    data["input_selection"] = manifest
    provenance = data.setdefault("provenance",{})
    provenance.update(package=(PACKAGE/"VERSION").read_text().strip(),
                      study="rwt10m-available"+("-pilot" if manifest["kind"]=="pilot" else ""),
                      source_events=counts["lo-cms"],target_events=counts["nlo-cp5"],
                      input_selection_sha256=manifest["selection_sha256"])
    result.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    if model is not None:
        cardpath = model/"model-card.json"
        card = json.loads(cardpath.read_text())
        card["input_selection"] = manifest
        card["training_sample_counts"] = counts
        card["analysis_package"] = (PACKAGE/"VERSION").read_text().strip()
        cardpath.write_text(json.dumps(card,indent=2,allow_nan=False)+"\n")
        (model/"selection-manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
        (directory/"model-card.json").write_text(cardpath.read_text())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command",required=True)
    f = sub.add_parser("freeze")
    f.add_argument("--inventory",type=Path,required=True)
    f.add_argument("--out",type=Path,required=True)
    f.add_argument("--pilot",action="store_true")
    i = sub.add_parser("inputs"); i.add_argument("manifest",type=Path)
    v = sub.add_parser("verify"); v.add_argument("manifest",type=Path);v.add_argument("--dir",type=Path,required=True)
    a = sub.add_parser("annotate");a.add_argument("manifest",type=Path);a.add_argument("--dir",type=Path,required=True)
    a.add_argument("--result",type=Path,required=True);a.add_argument("--model",type=Path)
    args = ap.parse_args()
    if args.command == "freeze":
        d = freeze(args.inventory,args.out,args.pilot)
        print(json.dumps(dict(expected_events=d["expected_events"],omitted=d["omitted"],selection_sha256=d["selection_sha256"])))
    else:
        d = load(args.manifest)
        if args.command == "inputs":
            for local, product in products(d): print(local,product)
        elif args.command == "verify": print(json.dumps(verify(d,args.dir)))
        else: annotate(d,args.dir,args.result,args.model)


if __name__ == "__main__":
    main()
