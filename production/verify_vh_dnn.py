"""Worker-only integrity/checkpoint verification; physics failure is retained."""
import json
import os
from pathlib import Path
import platform
import sys

if platform.node().split('.')[0] == 'cmscon' or not os.environ.get('_CONDOR_SCRATCH_DIR'):
    raise SystemExit('Compute worker only')

import numpy as np
import torch
sys.path[:0] = ['vh', 'rwt']
from train import network, predict_logits
from ggzh_reweight_common import Sample, FEATURE_SETS
from ggzh_reweight_bdt import prepare

torch.set_num_threads(int(os.environ.get('RWT_THREADS', '4')))
out = Path('dnn-out')
report = json.loads((out/'dnn-clf-c.json').read_text())
assert report['features'] == FEATURE_SETS['C']
assert (out/'COMPLETED').is_file()
source = Sample('source', Path('rwt-out/lo-cms.npz'))
x = prepare(source, report['features']).astype('float32')
indices = np.flatnonzero(source.in_window)
folds = np.arange(len(indices)) % 5
weights = np.load(out/'weights-clf-c.npz', allow_pickle=False)
assert np.array_equal(weights['run'], source['run'])
assert np.array_equal(weights['event'], source['event'])
assert np.array_equal(weights['in_window'], source.in_window)
checked = []
for tag in ('full', 'A_to_B', 'B_to_A'):
    w = weights['w_'+tag]
    assert w.shape == (source.n,) and np.isfinite(w).all() and (w > 0).all()
    assert abs(float(w[source.in_window].mean()) - 1) < 1e-10
    assert (w[~source.in_window] == 1).all()
    for k in range(5):
        card = torch.load(out/f'{tag}-fold-{k}.pt', map_location='cpu', weights_only=True)
        assert card['features'] == report['features'] and card['tag'] == tag and card['fold'] == k
        net = network(card['n_features'], card['hidden']).cuda()
        net.load_state_dict(card['state'])
        ii = indices[folds == k][:128]
        assert len(ii)
        z = predict_logits(net, x[ii], card, 'cuda')
        expected = np.exp(z - report[tag]['log_normalization'])
        np.testing.assert_allclose(w[ii], expected, rtol=5e-5, atol=1e-9)
        checked.append(f'{tag}-fold-{k}')
        del net
    assert all(key in report[tag] for key in ('reweighted', 'unweighted', 'n_eff', 'pass'))
assert all(key in report['diagnostics'] for key in ('regions', 'hardjet', 'mll_x_ptz', 'neff'))
assert set(report['summary']) == {'1_stage13', '2_variables', '2b_features', '3_stability', 'all_pass'}
verification = dict(status='PASS', device='cuda', reloaded_checkpoints=checked,
                    source_events=source.n, selection_sha256=report['input_selection']['selection_sha256'],
                    physics_summary=report['summary'], trainer_commit=os.environ['VH_TRAINER_COMMIT'],
                    package=Path('VERSION').read_text().strip())
(out/'verification.json').write_text(json.dumps(verification, indent=2)+'\n')
print(json.dumps(verification), flush=True)
