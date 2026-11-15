"""Worker-only ATLAS Born scale validation and frozen calibration companions."""
import concurrent.futures
import json
import math
import re
import sys
import tarfile
import time
from pathlib import Path
import scale_inputs as base
from scale_atlas_control import ht

np = base.np
REFERENCE_SHA = '0d1f704d6ac72a7c534832d284e2707e026c6bb5b2fa74f8ef7fb4300c2b4ad1'
SELECTION_SHA = {
    'zll': 'c84a4ed1df7d3f822b9230261e0db9b6282a8c1eb90ffc8775949a131f033f4f',
    'zvv': '0730348768f677516fab0b8919324e07b419bd0f1b9f083f1844f765fa6d73ef',
}
# Reuse the passed control at its immutable package and product fingerprint.
CONTROL_PACKAGE = 'bf4a9f53f8fec05141e6784c2c643f052dd730d5'
CONTROL_SHA = '7a03ee81f93478f53b96d7207183a27007b5faa8a874496bea41483aae88ddf4'
RTOL, ATOL = 1e-3, 1e-8


def card_check(card, channel):
    required = dict(lhans1=93300, lhans2=93300, dynscale=3, bornonly=1,
                    loevents=1, ebeam1=6800, ebeam2=6800, mtscheme=0,
                    alphas_from_pdf=1, renscfact=1, facscfact=1,
                    zmass=91.1876, hmass=125, topmass=172.5,
                    vdecaymode=11 if channel == 'zll' else 12)
    assert card.get('fixedscale', 0) != 1
    for key, value in required.items():
        assert math.isclose(card[key], value, rel_tol=1e-12), (key, card.get(key), value)
    return {key: card[key] for key in required}


def ratios(event, pdf):
    return np.asarray(base.born_ratio(event, pdf, pdf, base.POINTS, False, mu0=ht(event)))


def validate():
    path = Path('atlas-scale-reference.lhe.gz')
    assert base.sha(path) == REFERENCE_SHA
    old = json.loads(Path('previous-control.json').read_text())
    assert old['accepted'] and old['lhe_sha256'] == REFERENCE_SHA and old['events'] == 5000
    card, defs = base.header(path)
    checked = card_check(card, 'zll')
    for key, (r, f) in zip(base.IDS, base.POINTS):
        assert 'lhapdf=93300' in defs[key]
        assert base.num(re.search(r'renscfact=(\S+)', defs[key])[1]) == r
        assert base.num(re.search(r'facscfact=(\S+)', defs[key])[1]) == f
    pdf = base.lhapdf.mkPDF(base.PDFS[93300], 0)
    audit = dict(schema='ggzh-atlas-born-control-v2', accepted=False,
                 package=Path('VERSION').read_text().strip(), reference_cluster=1445174,
                 lhe_sha256=REFERENCE_SHA, previous_control_sha256=base.sha('previous-control.json'),
                 pdf=base.PDFS[93300], lhapdf_id=93300, card=checked,
                 scale='HT=mT(H)+pT(daughter1)+pT(daughter2)',
                 scale_ids=base.IDS, scale_points=base.POINTS,
                 rtol=RTOL, atol=ATOL, compared='stored weight / nominal LHE weight',
                 events=0, max_relative_error_by_id={}, failed_by_id={})
    worst = np.zeros(9); failed = np.zeros(9, dtype=np.int64)
    try:
        for event in base.events(path):
            exact = np.asarray([event['weights'][key] / event['weight'] for key in base.IDS])
            predicted = ratios(event, pdf)
            assert np.isfinite(exact).all() and np.all(exact != 0)
            worst = np.maximum(worst, np.abs(predicted - exact) / np.abs(exact))
            failed += ~np.isclose(predicted, exact, rtol=RTOL, atol=ATOL)
            audit['events'] += 1
        audit.update(max_relative_error_by_id=dict(zip(base.IDS, worst.tolist())),
                     failed_by_id=dict(zip(base.IDS, failed.tolist())),
                     max_relative_error=float(worst.max()), failed_comparisons=int(failed.sum()),
                     comparisons=audit['events'] * 9)
        assert audit['events'] == 5000
        assert not failed.any(), audit['failed_by_id']
        audit['accepted'] = True
    except Exception as exc:
        audit['error'] = repr(exc)
        raise
    finally:
        Path('scale-audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print('ATLAS_CONTROL_PASS', json.dumps(audit), flush=True)


def calibration(channel):
    start = time.time()
    selection = Path(f'production/scale-selections/atlas-{channel}-v2-calib.json')
    assert base.sha(selection) == SELECTION_SHA[channel]
    spec = json.loads(selection.read_text())
    assert spec['role'] == 'calib' and spec['style'] == 'atlas' and spec['channel'] == channel
    entries = sorted(spec['samples']['lo-cms'], key=lambda x: x['seed'])
    assert len(entries) == len({x['seed'] for x in entries}) == 200
    assert sum(x['events'] for x in entries) == 1000000
    control = json.loads(Path('control.json').read_text())
    assert control['accepted'] and control['failed_comparisons'] == 0 and control['events'] == 5000
    assert control['rtol'] == RTOL and control['atol'] == ATOL and control['lhe_sha256'] == REFERENCE_SHA
    assert control['package'].split()[0] == CONTROL_PACKAGE
    assert base.sha('control.json') == CONTROL_SHA
    assert base.sha('calib-prepared.tgz') == Path('calib-prepared.tgz.sha256.ref').read_text().split()[0]
    # Read only the LO calibration array and its small audit from the prepared bundle.
    with tarfile.open('calib-prepared.tgz') as archive:
        members = archive.getmembers()
        assert all((m.name == 'parton-data' or m.name.startswith('parton-data/'))
                   and '..' not in Path(m.name).parts and not m.issym() and not m.islnk()
                   for m in members)
        wanted = ['parton-data/lo-cms.npz', 'parton-data/audit.json', 'parton-data/SHA256SUMS']
        for name in wanted: archive.extract(name, '.')
    checks = dict((name.lstrip('./'), digest) for digest, name in
                  (line.split(maxsplit=1) for line in Path('parton-data/SHA256SUMS').read_text().splitlines()))
    for name in ('lo-cms.npz', 'audit.json'):
        assert base.sha(Path('parton-data') / name) == checks[name]
    prepared_audit = json.loads(Path('parton-data/audit.json').read_text())
    assert prepared_audit['selection_sha256'] == SELECTION_SHA[channel]
    with np.load('parton-data/lo-cms.npz', allow_pickle=False) as original:
        seed = original['run'].copy(); row = original['event'].copy()
        nominal_table = original['w0'].copy(); mz_table = original['mz'].copy()
    n = len(seed); assert n == 1000000
    data = dict(seed=seed, table_row=row, nominal_lhe=np.empty(n), scale_born=np.empty((n, 9)))
    pdf = base.lhapdf.mkPDF(base.PDFS[93300], 0)
    audit = dict(schema='ggzh-atlas-scale-calibration-v2', accepted=False, channel=channel,
                 package=Path('VERSION').read_text().strip(), selection_role='calib',
                 selection_sha256=SELECTION_SHA[channel], prepared_sha256=base.sha('calib-prepared.tgz'),
                 control_sha256=base.sha('control.json'), control=control,
                 pdf=base.PDFS[93300], lhapdf_id=93300, scale=control['scale'],
                 scale_ids=base.IDS, scale_points=base.POINTS, shards=[],
                 unit='pb', fixed_shower=True, row_order='exact prepared LO (run,event)',
                 alignment='complete 5000-row LHE/table shards; nominal weight and Z mass checked in order')
    def fetch(entry):
        sample = 'lo-atlas-parton10m-13600' if channel == 'zll' else 'lo-atlas-zvv-bdt-v1'
        product = f'lhe/{sample}/events-{entry["seed"]}.lhe.gz'
        return base.download(product, Path('lhe') / str(entry['seed']) / 'events.lhe.gz', entry['lhe_sha256'])
    try:
        offset = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            futures = {i: pool.submit(fetch, entry) for i, entry in enumerate(entries[:4])}
            for j, entry in enumerate(entries):
                path = futures.pop(j).result()
                if j + 4 < len(entries): futures[j + 4] = pool.submit(fetch, entries[j + 4])
                rows = entry['events']; assert rows == 5000
                sl = slice(offset, offset + rows)
                assert np.array_equal(seed[sl], np.full(rows, entry['seed']))
                assert np.array_equal(row[sl], np.arange(rows))
                card, defs = base.header(path); checked = card_check(card, channel)
                nominal = []; grid = []; masses = []
                for event in base.events(path):
                    assert not event['weights'], 'Unexpected stored weight block; review input'
                    nominal.append(event['weight']); grid.append(event['weight'] * ratios(event, pdf))
                    daughters = [p['v'] for p in event['particles'] if p['status'] == 1 and abs(p['id']) in (11,12,13,14,15,16)]
                    assert len(daughters) == 2
                    v = np.sum(daughters, axis=0)
                    masses.append(math.sqrt(max(0., v[3]**2 - np.dot(v[:3], v[:3]))))
                nominal = np.asarray(nominal); grid = np.asarray(grid); masses = np.asarray(masses)
                assert len(nominal) == rows and np.isfinite(grid).all()
                np.testing.assert_allclose(nominal, nominal_table[sl], rtol=1e-6, atol=1e-10)
                np.testing.assert_allclose(masses, mz_table[sl], rtol=1e-5, atol=1e-6)
                np.testing.assert_allclose(grid[:, 0], nominal, rtol=1e-12, atol=1e-14)
                data['nominal_lhe'][sl] = nominal; data['scale_born'][sl] = grid
                audit['shards'].append(dict(seed=entry['seed'], rows=rows, lhe_sha256=entry['lhe_sha256'], card=checked,
                    max_relative_mz_difference=float(np.max(abs(masses-mz_table[sl])/masses))))
                offset += rows; path.unlink(); Path(str(path)+'.sha256').unlink()
                print('ATLAS_SCALE_SHARD', channel, j+1, len(entries), entry['seed'], flush=True)
        assert offset == n
        out = Path('scale-inputs'); out.mkdir()
        np.savez_compressed(out/'lo-cms.npz', **data, scale_ids=np.asarray(base.IDS), scale_points=np.asarray(base.POINTS))
        with np.load(out/'lo-cms.npz', allow_pickle=False) as saved:
            assert np.array_equal(saved['seed'], seed) and np.array_equal(saved['table_row'], row)
            assert saved['scale_born'].shape == (n, 9) and list(saved['scale_ids']) == base.IDS
        audit.update(accepted=True, rows=n, joined=True, output_sha256=base.sha(out/'lo-cms.npz'),
                     elapsed_seconds=time.time()-start, sumw=data['scale_born'].sum(axis=0).tolist())
        (out/'README.md').write_text('''# ATLAS LO Born scale companions\n\nPDF4LHC21_40_pdfas, LHAPDF 93300. mu0 = H_T = mT(H) + pT(Z daughter 1) + pT(Z daughter 2), POWHEG dynscale 3, no factor 1/2.\nscale-inputs/lo-cms.npz is the compatibility filename for the ATLAS LO+PS source.\nscale_born contains absolute signed LHE weights in pb, columns 1001--1009; scale_points gives (muR/mu0, muF/mu0).\nseed and table_row exactly match the frozen prepared calibration bundle. No evaluation events.\nBorn hard-scale variations at fixed shower; no shower-kernel variation or extra normalization.\n''')
        (out/'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
        paths = sorted(out.iterdir()); (out/'SHA256SUMS').write_text(''.join(f'{base.sha(p)}  {p.name}\n' for p in paths))
        with tarfile.open('scale-inputs.tgz', 'w:gz', compresslevel=1) as archive: archive.add(out, arcname='scale-inputs')
    except Exception as exc:
        audit['error'] = repr(exc)
        raise
    finally:
        Path('scale-audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    print('ATLAS_SCALE_CALIBRATION_COMPLETE', channel, n, flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'validate': validate()
    else:
        assert sys.argv[1] == 'calibration'
        calibration(sys.argv[2])
