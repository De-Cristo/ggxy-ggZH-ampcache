#!/usr/bin/env python3
"""Apply the accepted native replay to all original Zll shards; export keyed weights."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sys

import nlo_topup as replay


def configure():
    assert os.environ.get('_CONDOR_SCRATCH_DIR') and os.uname().nodename.split('.')[0] != 'cmscon'
    replay.SEED = int(os.environ['NLO_TOPUP_SEED'])
    assert 100001 <= replay.SEED <= 104000
    replay.EXPECTED = Path('original.lhe.gz.sha256.ref').read_text().split()[0]
    assert re.fullmatch('[0-9a-f]{64}', replay.EXPECTED)


def prepare():
    replay.prepare()
    if os.environ.get('NLO_TOPUP_REUSE') == '1':
        assert replay.SEED == 101654
        a = json.loads(Path('reused-audit.json').read_text())
        assert a['accepted'] and a['failed_controls'] == 0 and a['input_sha256'] == replay.EXPECTED
        assert a['package'].split()[0] == '41d66c5ece3d433c7778b985017515e0669d6098'
        blob = Path('reused-topup.lhe.gz').read_bytes()
        assert hashlib.sha256(blob).hexdigest() == a['output_sha256']
        Path('topup-run/pwgevents-rwgt.lhe').write_bytes(gzip.decompress(blob))
        Path('topup-time.txt').write_text('0 0 0\n')


def validate():
    replay.validate()
    a = json.loads(Path('topup-audit.json').read_text())
    a.update(schema='ggzh-nlo-scale-topup-v1', weight_unit='pb', ordinal_base=0,
             scale_points=replay.SCALES, source_pdf=93300, added_pdf=325300,
             input_production_package='3161a6224c1c354da81671c6b54fe66a27f303e8',
             join_keys=['seed','lhe_ordinal'])
    if os.environ.get('NLO_TOPUP_REUSE') == '1':
        old = json.loads(Path('reused-audit.json').read_text())
        a['reused_calculation'] = dict(package=old['package'], seed=replay.SEED,
                                      output_sha256=old['output_sha256'],
                                      previously_accounted_cpu_hours=old['replay_cpu_hours'])
    events = replay.parse(gzip.decompress(Path('topup.lhe.gz').read_bytes()).decode())
    ids = ['2000']+[str(i) for i in range(1001,1010)]+[str(i) for i in range(5101,5110)]
    with gzip.open('topup-weights.csv.gz','wt',newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['seed','lhe_ordinal','w:nominal']+['w:'+k for k in ids])
        for ordinal,(core,_,weights) in enumerate(events):
            writer.writerow([replay.SEED,ordinal,core[0][2]]+[weights[k] for k in ids])
    a['weights_sha256'] = hashlib.sha256(Path('topup-weights.csv.gz').read_bytes()).hexdigest()
    a['companion_rows'] = len(events)
    a['companion_columns'] = ['seed','lhe_ordinal','w:nominal']+['w:'+k for k in ids]
    Path('topup-audit.json').write_text(json.dumps(a,indent=2)+'\n')
    print('TOPUP_COMPANION_ACCEPTED', replay.SEED, len(events), flush=True)


if __name__ == '__main__':
    configure()
    {'prepare':prepare,'validate':validate}[sys.argv[1]]()
