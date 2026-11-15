"""Freeze calibration, evaluation and ensemble-member selections of scalar-table shards (ATLAS-style inputs).

Reads the delivered-shard catalog of the scalar production (atlas-<ch>-material.json)
and the small per-shard manifests and validation reports only.  In a salted SHA256 order of the file seeds, the
first 200 LO / 400 NLO eligible shards form the calibration set, the next 200 / 400 the evaluation set, and every
member draws 200 / 400 from the remaining pool in its own salted order.  The selections use the sample keys of the
CMS bundles ('lo-cms' = LO+PS source, 'nlo' = NLO+PS target) so the same member and assessment code reads them.
"""
import argparse
import hashlib
import json
from pathlib import Path

COUNTS = (('lo-atlas', 'lo-cms', 200, 5000), ('nlo-a14', 'nlo', 400, 2500))
SEED_BASE = {'atlas-zll': 770100, 'atlas-zvv': 780100}

def digest(data):
    return hashlib.sha256(data).hexdigest()

def ch_of(storage):
    return storage.split('-', 1)[1]

def eligible(key, x, storage, order, rows):
    manifest = Path(x['manifest']); assert manifest.stat().st_size < 100000
    b = manifest.read_bytes(); assert digest(b) == x['manifest_sha256']
    m = json.loads(b); assert m['primary_included'] and m['events'] == rows
    assert m['channel'] == storage and m['order'] == order and m['seed'] == x['seed']
    files = {}
    for f in m['files']:
        path = Path(f['path']); assert f['bytes'] > 0 and path.stat().st_size == f['bytes']
        files[f['role']] = dict(role=f['role'], product=f['product'], sha256=f['sha256'], bytes=f['bytes'])
    assert set(files) == {'legacy_table', 'schema', 'validation'}
    v = Path(m['files'][[f['role'] for f in m['files']].index('validation')]['path']).read_bytes()
    assert digest(v) == files['validation']['sha256']; v = json.loads(v)
    if v.get('schema') == 'ggzh-particle-v1':
        # Shard reused from the particle production: its own validation format (table checked there against the
        # legacy writer; the manifest carries the table checksum).
        meta = v['metadata']
        assert v['events'] == v['expected'] == rows and not v['missing_ordinals'] and v['legacy_comparison']['pass_']
        assert meta['channel'] in (storage, ch_of(storage)) and meta['order'] == order and meta['seed'] == x['seed']
    else:
        assert v['accepted'] and v['primary_included'] and v['events'] == v['expected'] == rows
        assert v['channel'] == storage and v['order'] == order and v['seed'] == x['seed']
        assert v['table_sha256'] == files['legacy_table']['sha256']
    return dict(key=key, seed=x['seed'], events=rows, producer=m['package'], primary_included=True,
                manifest_sha256=x['manifest_sha256'], lhe_sha256=m.get('input_lhe_sha256'), sumw=v['sumw'],
                files=list(files.values()))

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', required=True, type=Path)
    p.add_argument('--storage', required=True, choices=sorted(SEED_BASE))
    p.add_argument('--members', type=int, default=20)
    p.add_argument('--out', required=True, type=Path, help='directory, e.g. production/ensemble')
    a = p.parse_args()
    st = a.storage; ch = st.split('-')[1]; salt = f'ggzh-scalar-{st}-ensemble-20261112-v2'
    raw = a.catalog.read_bytes(); catalog = json.loads(raw)
    names = ['calib', 'eval']+[f'member-{i:02d}' for i in range(a.members)]
    for n in names: assert not (a.out/f'{st}-v2-{n}.json').exists(), n
    pools = {}; rejected = []
    for order, key, count, rows in COUNTS:
        good = []
        for k, x in catalog['delivered'].items():
            if x['channel'] != st or x['order'] != order or not x['primary_included']: continue
            try: good.append(eligible(k, x, st, order, rows))
            except (AssertionError, OSError, ValueError, KeyError) as exc: rejected.append(dict(key=k, reason=type(exc).__name__))
        good.sort(key=lambda e: digest(f'{salt}|roles|{order}|{e["seed"]}'.encode()))
        assert len(good) >= 3*count, (order, len(good), 'need calibration + evaluation + one member pool')
        pools[key] = good
    selections = {}
    for order, key, count, rows in COUNTS:
        good = pools[key]
        selections.setdefault('calib', {})[key] = good[:count]
        selections.setdefault('eval', {})[key] = good[count:2*count]
        rest = good[2*count:]
        for i in range(a.members):
            n = f'member-{i:02d}'
            selections.setdefault(n, {})[key] = sorted(rest, key=lambda e: digest(f'{salt}|{n}|{order}|{e["seed"]}'.encode()))[:count]
        pools[key] = rest
    a.out.mkdir(parents=True, exist_ok=True)
    summary = dict(schema='ggzh-scalar-ensemble-index-v1', storage=st, channel=ch, style='atlas', salt=salt,
                   catalog=a.catalog.name, catalog_checked=catalog.get('checked'), catalog_sha256=digest(raw),
                   source_order='lo-atlas', target_order='nlo-a14', nominal_target_weight='w:nominal (PDF4LHC21 central = id 1001)',
                   eligibility_rejections=rejected, pool_shards={k: len(v) for k, v in pools.items()}, selections={})
    for i, n in enumerate(names):
        samples = {k: sorted(xs, key=lambda e: e['seed']) for k, xs in selections[n].items()}
        result = dict(schema='ggzh-scalar-ensemble-selection-v1', study=f'ensemble-{st}-v2-{n}', channel=ch, style='atlas',
                      storage=st, role=n, selection_rule=f'SHA256 {salt}|...|order|seed; scalar tables only',
                      catalog_checked=catalog.get('checked'), catalog_sha256=digest(raw), events_per_sample=1000000,
                      sample_orders={'lo-cms': 'lo-atlas', 'nlo': 'nlo-a14'},
                      initialization_seed=SEED_BASE[st]+i-2 if n.startswith('member') else None, samples=samples)
        path = a.out/f'{st}-v2-{n}.json'; path.write_text(json.dumps(result, indent=2)+'\n')
        summary['selections'][n] = dict(path=str(path), sha256=digest(path.read_bytes()),
                                         seeds={k: [e['seed'] for e in xs] for k, xs in samples.items()})
    fixed = {k: set(summary['selections']['calib']['seeds'][k]) | set(summary['selections']['eval']['seeds'][k]) for k in ('lo-cms', 'nlo')}
    assert not set(summary['selections']['calib']['seeds']['nlo']) & set(summary['selections']['eval']['seeds']['nlo'])
    for n in names[2:]:
        for k in ('lo-cms', 'nlo'): assert not set(summary['selections'][n]['seeds'][k]) & fixed[k]
    (a.out/f'{st}-v2-ensemble.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(dict(storage=st, pool_shards=summary['pool_shards'], rejected=len(rejected),
                          selections={n: s['sha256'][:12] for n, s in summary['selections'].items()})))

if __name__ == '__main__':
    main()
