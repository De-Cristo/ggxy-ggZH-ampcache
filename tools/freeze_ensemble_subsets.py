"""Freeze the evaluation and ensemble-member selections from small manifests/validation only.

The calibration set is the existing frozen production/parton-<channel>-1m-v1.json.  The
evaluation set takes the first 200 LO / 400 NLO eligible shards outside it in a salted SHA256
order; every member then draws 200 LO / 400 NLO shards from the remaining pool in its own
salted order.  Members overlap each other but never the calibration or evaluation shards.
"""
import argparse
import hashlib
import json
from pathlib import Path

COUNTS = (('lo-cms', 200, 5000), ('nlo', 400, 2500))
SEED_BASE = {'zll': 750000, 'zvv': 760000}

def digest(data):
    return hashlib.sha256(data).hexdigest()

def eligible(key, x, channel, order, rows, producer, parent):
    manifest = Path(x['manifest']); assert manifest.stat().st_size < 100000
    b = manifest.read_bytes(); assert digest(b) == x['manifest_sha256']
    m = json.loads(b); assert m['primary_included'] and m['events'] == rows
    assert m['channel'] == channel and m['order'] == order and m['seed'] == x['seed']
    files = []; validation = None
    for f in m['files']:
        path = Path(f['path']); assert f['bytes'] > 0 and path.stat().st_size == f['bytes']
        product = 'particle-v1/'+str(path).split('/particle-v1/', 1)[1]
        files.append(dict(role=f['role'], product=product, sha256=f['sha256'], bytes=f['bytes']))
        if f['role'] == 'validation':
            assert path.stat().st_size < 100000
            b = path.read_bytes(); assert digest(b) == f['sha256']; validation = json.loads(b)
    assert {f['role'] for f in files} == {'parton_inputs', 'closure', 'legacy_table', 'schema', 'validation'}
    v = validation; meta = v['metadata']
    assert meta['package'].split()[0] == producer and meta['frozen_manifest'] == parent
    assert meta['primary_included'] and meta['channel'] == channel and meta['order'] == order and meta['seed'] == x['seed']
    assert v['events'] == v['expected'] == rows and not v['missing_ordinals']
    assert v['legacy_comparison']['pass_'] and v['views']['partons-v1']['overflow_events'] == 0
    assert v['forward_validation']['partons-v1']['finite']
    entry = dict(key=key, seed=x['seed'], events=rows, producer=producer, primary_included=True,
                 manifest_sha256=x['manifest_sha256'], lhe_sha256=meta['lhe_sha256'], files=files)
    return entry, v['views']['partons-v1']['max_count']

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', required=True, type=Path)
    p.add_argument('--channel', required=True, choices=('zll', 'zvv'))
    p.add_argument('--members', type=int, default=20)
    p.add_argument('--out', required=True, type=Path, help='directory, e.g. production/ensemble')
    a = p.parse_args()
    ch = a.channel; salt = f'ggzh-parton-{ch}-ensemble-20261111-v1'
    calibration_path = Path(f'production/parton-{ch}-1m-v1.json')
    calibration_raw = calibration_path.read_bytes(); calibration = json.loads(calibration_raw)
    parent = calibration['parent_selection_sha256']
    producers = {x['producer'] for xs in calibration['samples'].values() for x in xs}; assert len(producers) == 1
    producer = producers.pop()
    raw = a.catalog.read_bytes(); catalog = json.loads(raw)
    a.out.mkdir(parents=True, exist_ok=True)
    names = ['eval']+[f'member-{i:02d}' for i in range(a.members)]
    for n in names: assert not (a.out/f'{ch}-{n}.json').exists(), n
    taken = {o: {x['seed'] for x in calibration['samples'][o]} for o, _, _ in COUNTS}
    pools = {}; rejected = []; maximum = {}
    for order, count, rows in COUNTS:
        candidates = [(k, x) for k, x in catalog['delivered'].items()
                      if x['channel'] == ch and x['order'] == order and x['primary_included']
                      and x['package'] == producer and x['seed'] not in taken[order]]
        good = []
        for key, x in candidates:
            try:
                good.append(eligible(key, x, ch, order, rows, producer, parent))
            except (AssertionError, OSError, ValueError, KeyError) as exc:
                rejected.append(dict(key=key, reason=type(exc).__name__))
        pools[order] = good
    selections = {}
    for order, count, rows in COUNTS:
        good = sorted(pools[order], key=lambda e: digest(f'{salt}|eval|{order}|{e[0]["seed"]}'.encode()))
        assert len(good) >= 2*count, (order, len(good))
        selections.setdefault('eval', {})[order] = good[:count]
        pools[order] = good[count:]
        for i in range(a.members):
            n = f'member-{i:02d}'
            ranked = sorted(pools[order], key=lambda e: digest(f'{salt}|{n}|{order}|{e[0]["seed"]}'.encode()))
            selections.setdefault(n, {})[order] = ranked[:count]
    tokens = max(m for s in selections.values() for xs in s.values() for _, m in xs)
    length = max(4+tokens, calibration['maximum_tensor_length'])
    summary = dict(schema='ggzh-parton-ensemble-index-v1', channel=ch, salt=salt, catalog_checked=catalog['checked'],
                   catalog_sha256=digest(raw), calibration_selection=str(calibration_path),
                   calibration_selection_sha256=digest(calibration_raw), producer=producer,
                   parent_selection_sha256=parent, maximum_tensor_length=length, eligibility_rejections=rejected,
                   pool_shards={o: len(pools[o]) for o, _, _ in COUNTS}, selections={})
    for i, n in enumerate(names):
        samples = {o: sorted((e for e, _ in xs), key=lambda e: e['seed']) for o, xs in selections[n].items()}
        result = dict(schema='ggzh-parton-ensemble-selection-v1', study=f'ensemble-{ch}-v1-{n}', channel=ch, role=n,
                      selection_rule=f'SHA256 {salt}|{n}|order|seed; first 200 LO/400 NLO eligible producer shards '
                                     'outside the calibration selection'+('' if n == 'eval' else ' and the evaluation selection'),
                      catalog_checked=catalog['checked'], catalog_sha256=digest(raw), parent_selection_sha256=parent,
                      calibration_selection_sha256=digest(calibration_raw), events_per_sample=1000000,
                      maximum_radiation_tokens={o: max(m for _, m in xs) for o, xs in selections[n].items()},
                      maximum_tensor_length=length,
                      initialization_seed=None if n == 'eval' else SEED_BASE[ch]+i-1, samples=samples)
        path = a.out/f'{ch}-{n}.json'; path.write_text(json.dumps(result, indent=2)+'\n')
        summary['selections'][n] = dict(path=str(path), sha256=digest(path.read_bytes()),
                                         seeds={o: [e['seed'] for e in xs] for o, xs in samples.items()})
    ev = {o: set(summary['selections']['eval']['seeds'][o]) for o, _, _ in COUNTS}
    for n in names[1:]:
        for o, _, _ in COUNTS:
            s = set(summary['selections'][n]['seeds'][o])
            assert len(s) == len(summary['selections'][n]['seeds'][o]) and not s & ev[o] and not s & taken[o]
    (a.out/f'{ch}-ensemble.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(dict(channel=ch, maximum_tensor_length=length, pool_shards=summary['pool_shards'],
                          rejected=len(rejected), selections={n: s['sha256'][:12] for n, s in summary['selections'].items()})))

if __name__ == '__main__':
    main()
