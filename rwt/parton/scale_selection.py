"""Worker-only Born companions for frozen selections, keyed by original LHE ordinal."""
import csv
import gzip
import json
from pathlib import Path
import sys
import tarfile
import time

import scale_inputs as base


def main():
    start = time.time()
    selection = Path(sys.argv[1])
    spec = json.loads(selection.read_text())
    ch = spec['channel']
    assert ch in ('zll', 'zvv')
    entries = sorted(spec['samples']['lo-cms'], key=lambda x: x['seed'])
    assert len(entries) == len({e['seed'] for e in entries}) == 200
    assert sum(e['events'] for e in entries) == 1000000
    pdfs = {i: base.lhapdf.mkPDF(name, 0) for i, name in base.PDFS.items()}
    audit = dict(schema='ggzh-scale-selection-v1', channel=ch,
                 package=Path('VERSION').read_text().strip(),
                 selection_sha256=base.sha(selection), selection_role=spec['role'],
                 control=base.control(pdfs), scale_points=base.POINTS,
                 join_keys=['seed', 'lhe_ordinal'], ordinal_base=0, unit='pb',
                 fixed_shower=True, shards=[], accepted=False)
    out = Path('scale-inputs'); out.mkdir()
    try:
        suffix = 'rwt10m-13600' if ch == 'zll' else 'vv-rwt10m-13600'
        total = 0
        with gzip.open(out/'lo-cms.csv.gz', 'wt', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['seed', 'lhe_ordinal', 'w:nominal'] +
                            ['nnpdf:'+k for k in base.IDS] + ['pdf4lhc:'+k for k in base.IDS])
            for entry in entries:
                seed = entry['seed']
                product = f'lhe/lo-cms-{suffix}/events-{seed}.lhe.gz'
                path = base.download(product, Path('lhe')/str(seed)/'events.lhe.gz', entry['lhe_sha256'])
                card, defs = base.header(path)
                checked = base.check_card(card, 'lo-cms', 325300, True, ch)
                count = 0
                for ordinal, event in enumerate(base.events(path)):
                    assert not event['weights']
                    nominal = event['weight']
                    weights = [nominal*x for x in base.born_ratio(event, pdfs[325300], pdfs[325300], base.POINTS, True)]
                    other = [nominal*x for x in base.born_ratio(event, pdfs[325300], pdfs[93300], base.POINTS, True)]
                    assert base.math.isclose(weights[0], nominal, rel_tol=1e-12, abs_tol=1e-14)
                    assert all(base.math.isfinite(x) for x in [nominal]+weights+other)
                    writer.writerow([seed, ordinal, nominal]+weights+other)
                    count += 1
                assert count == entry['events'] == 5000
                total += count
                audit['shards'].append(dict(seed=seed, rows=count, lhe_sha256=entry['lhe_sha256'], card=checked))
                path.unlink(); Path(str(path)+'.sha256').unlink()
                print('SCALE_SHARD', ch, seed, count, flush=True)
        assert total == 1000000
        audit.update(accepted=True, rows=total, output_sha256=base.sha(out/'lo-cms.csv.gz'),
                     elapsed_seconds=time.time()-start,
                     limitation='Born hard-scale variation at fixed shower; no NLO or shower-kernel approximation. Join via the existing original LHE ordinal map.')
    except Exception as exc:
        audit['error'] = repr(exc)
        raise
    finally:
        Path('scale-audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    (out/'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    (out/'SHA256SUMS').write_text(''.join(f'{base.sha(p)}  {p.name}\n' for p in sorted(out.iterdir())))
    with tarfile.open('scale-inputs.tgz', 'w:gz', compresslevel=1) as t:
        t.add(out, arcname='scale-inputs')


if __name__ == '__main__':
    main()
