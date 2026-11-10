#!/usr/bin/env python3
"""Worker-only native weight replay: preserve event identities and check old-weight controls."""
import collections
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tarfile

SEED = 101654
EXPECTED = 'e9a007aafe0f1cfab8e434e7bd3f6335c5a96400b45c7bd7283189b9f95c9d47'
SCALES = [(1,1),(1,2),(1,.5),(2,1),(2,2),(2,.5),(.5,1),(.5,2),(.5,.5)]
RTOL, ATOL = 1e-4, 1e-8


def parse(text):
    result = []
    for block in re.findall(r'<event>\s*(.*?)</event>', text, re.S):
        lines = block.splitlines()
        n = int(lines[0].split()[0])
        core = [[float(x.replace('D','E')) for x in ln.split()] for ln in lines[:n+1]]
        replay = [ln.split()[1:] for ln in lines if ln.strip().startswith('#rwgt ')]
        assert len(replay) == 1 and len(replay[0]) == 6, 'missing replay record'
        r = replay[0]
        assert int(r[0]) in (1,2,3) and int(r[1]) > 0
        assert math.isfinite(float(r[2])) and float(r[2]) != 0
        assert int(r[3]) == SEED and min(int(r[4]),int(r[5])) >= 0
        weights = {k: float(v.replace('D','E')) for k,v in
                   re.findall(r"<wgt\s+id=['\"]([^'\"]+)['\"]\s*>\s*([^<]+)</wgt>", block)}
        assert all(math.isfinite(v) for row in core for v in row)
        assert all(math.isfinite(v) for v in weights.values())
        result.append((core, r, weights))
    assert len(result) == 2500, f'event count {len(result)}'
    return result


def prepare():
    blob = Path('original.lhe.gz').read_bytes()
    assert hashlib.sha256(blob).hexdigest() == EXPECTED
    text = gzip.decompress(blob).decode()
    events = parse(text)
    required = {str(i) for i in range(1001,1010)} | {'2000'}
    assert all(required <= set(e[2]) for e in events)
    run = Path('topup-run'); run.mkdir()
    with tarfile.open('grids.tgz') as tf:
        for m in tf.getmembers():
            if m.isfile() and '/' not in m.name and m.name.startswith('pwg') and m.name.endswith('.dat'):
                (run/m.name).write_bytes(tf.extractfile(m).read())
    assert list(run.glob('pwggrid-*.dat')) and (run/'pwgubound.dat').is_file()
    card = text.split('Input file powheg.input contained:',1)[1].split('End of powheg.input content',1)[0]
    settings = {'rwl_add':'1', 'compute_rwgt':'0', 'rwl_file':"'topup-weights.xml'",
                'rwl_format_rwgt':'1', 'rwl_group_events':'2000', 'lhapdf6maxsets':'50',
                'manyseeds':'1', 'parallelstage':'4', 'numevts':'2500', 'iseed':str(SEED),
                'compress_lhe':'0', 'clobberlhe':'0'}
    lines = [ln for ln in card.splitlines() if not ln.split() or ln.split()[0] not in settings]
    lines += [f'{k} {v}' for k,v in settings.items()]
    (run/'powheg.input').write_text('\n'.join(lines)+'\n')
    (run/'pwgseeds.dat').write_text(str(SEED)+'\n')
    (run/'seed-index.txt').write_text('1\n')
    (run/'pwgevents.lhe').write_text(text)
    xml = ['<initrwgt>']
    for name,base,pdf in [('NNPDF31_scale_topup',5100,325300),('PDF4LHC_replay_control',6100,93300)]:
        xml.append(f"<weightgroup name='{name}' combine='envelope'>")
        for i,(r,f) in enumerate(SCALES,1):
            xml.append(f"<weight id='{base+i}'> lhapdf={pdf} renscfact={r} facscfact={f} </weight>")
        xml.append('</weightgroup>')
    xml.append('</initrwgt>')
    (run/'topup-weights.xml').write_text('\n'.join(xml)+'\n')
    audit = dict(seed=SEED,events=len(events),input_sha256=EXPECTED,
                 replay_types=dict(collections.Counter(e[1][0] for e in events)),
                 package=Path('VERSION').read_text().strip(),accepted=False,stage='metadata_verified',
                 rtol=RTOL,atol=ATOL)
    Path('topup-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit),flush=True)


def validate():
    audit = json.loads(Path('topup-audit.json').read_text())
    try:
        old = parse(gzip.decompress(Path('original.lhe.gz').read_bytes()).decode())
        out = Path('topup-run/pwgevents-rwgt.lhe')
        new = parse(out.read_text())
        maxima = collections.defaultdict(float); bad = 0
        sums = collections.defaultdict(float)
        for a,b in zip(old,new):
            assert a[0] == b[0], 'event record changed'
            assert a[1][:2]+a[1][3:] == b[1][:2]+b[1][3:], 'replay identity changed'
            assert all(k in b[2] and v == b[2][k] for k,v in a[2].items()), 'old weight changed'
            assert all(str(k) in b[2] for k in list(range(5101,5110))+list(range(6101,6110)))
            for control,ref in [(str(6100+i),str(1000+i)) for i in range(1,10)]+[('5101','2000')]:
                diff = abs(b[2][control]-a[2][ref]); rel = diff/max(abs(a[2][ref]),1e-300)
                maxima[f'{control}:{ref}'] = max(maxima[f'{control}:{ref}'],rel)
                maxima[f'type_{a[1][0]}'] = max(maxima[f'type_{a[1][0]}'],rel)
                bad += diff > ATOL+RTOL*abs(a[2][ref])
            sums['original_nominal'] += a[0][0][2]
            for k,v in b[2].items(): sums[k] += v
        audit.update(control_max_relative=dict(maxima),failed_controls=bad,
                     sigma_pb={k:v/2500 for k,v in sums.items()},records_preserved=True,
                     old_weights_preserved=True)
        assert bad == 0, f'{bad} control comparisons failed'
        times = Path('topup-time.txt').read_text().splitlines()[-1].split()
        cpu = (float(times[0])+float(times[1]))/3600
        audit.update(replay_cpu_hours=cpu,replay_wall_seconds=float(times[2]))
        assert cpu <= .90, 'replay cost exceeds twice generation baseline'
        with out.open('rb') as src, gzip.open('topup.lhe.gz','wb') as dst:
            import shutil
            shutil.copyfileobj(src,dst)
        audit.update(accepted=True,stage='validated',output_sha256=hashlib.sha256(Path('topup.lhe.gz').read_bytes()).hexdigest())
    except Exception as e:
        audit.update(accepted=False,stage='STOP',error=str(e))
        raise
    finally:
        Path('topup-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
        print(json.dumps(audit),flush=True)


if __name__ == '__main__':
    assert os.environ.get('_CONDOR_SCRATCH_DIR') and os.uname().nodename.split('.')[0] != 'cmscon'
    {'prepare':prepare,'validate':validate}[sys.argv[1]]()
