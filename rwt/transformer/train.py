#!/usr/bin/env python3
"""Fresh dense classifier ratios for ggZH. Run on a compute worker only.

Public ggxy-ggZH-ampcache supplies feature definitions and unchanged closure tests.
No nlorw, CMS pretrained model, binary extractor or original repository data is used.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

os.environ.setdefault('PYTHONDONTWRITEBYTECODE', '1')
from protocol import validate_config, source_fold, target_half


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def network(n, hidden):
    import torch
    from torch import nn
    from vendor.ParticleTransformer import ParticleTransformer
    class FeatureParT(nn.Module):
        def __init__(self):
            super().__init__()
            self.scale=nn.Parameter(torch.randn(n,64)/8)
            self.identity=nn.Parameter(torch.randn(n,64)*.02)
            self.backbone=ParticleTransformer(input_dim=64,num_classes=1,embed_dims=(),
                pair_input_dim=0,pair_embed_dims=None,num_heads=4,num_layers=2,num_cls_layers=2,
                trim=False,use_amp=False)
        def forward(self,x):
            tokens=x.unsqueeze(-1)*self.scale+self.identity
            mask=torch.ones((x.shape[0],1,n),dtype=torch.bool,device=x.device)
            return self.backbone(tokens,mask=mask)
    return FeatureParT()


def fit(xs, xt, wt, config, seed, device):
    import numpy as np
    import torch
    from torch.nn.functional import binary_cross_entropy_with_logits
    torch.manual_seed(seed)
    if device.startswith('cuda'):
        torch.cuda.manual_seed_all(seed)
    rng = np.random.default_rng(seed)
    splits = []
    for x in (xs, xt):
        order = rng.permutation(len(x))
        nval = max(1, int(len(x) * config['validation_fraction']))
        if len(x) - nval < 2:
            raise ValueError('Insufficient events for train/validation split')
        splits.append((order[nval:], order[:nval]))
    # The held-out source fold and opposite target half never reach this function.
    train = np.concatenate([xs[splits[0][0]], xt[splits[1][0]]])
    mu = train.mean(0, dtype=np.float64).astype('float32')
    sd = train.std(0, dtype=np.float64).astype('float32')
    sd = np.where(sd > 1e-8, sd, 1).astype('float32')
    del train

    def dataset(which):
        si, ti = splits[0][which], splits[1][which]
        x = np.concatenate([xs[si], xt[ti]])
        x = ((x - mu) / sd).astype('float32')
        y = np.concatenate([np.zeros(len(si)), np.ones(len(ti))]).astype('float32')
        tw = np.maximum(wt[ti], 0).astype('float64')
        if tw.sum() <= 0:
            raise ValueError('Target class has no positive weight')
        # Each class has total loss coefficient 1/2, as in the BDT control.
        w = np.concatenate([np.full(len(si), .5 / len(si)), .5 * tw / tw.sum()])
        return x, y, w.astype('float32')

    x, y, w = dataset(0)
    xv, yv, wv = dataset(1)
    net = network(xs.shape[1], config['hidden']).to(device)
    opt = torch.optim.AdamW(net.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay'])
    best, best_state, bad, history = float('inf'), None, 0, []
    bs = config['batch_size']
    for epoch in range(config['epochs']):
        net.train()
        order = rng.permutation(len(x))
        for start in range(0, len(x), bs):
            ii = order[start:start + bs]
            xb = torch.as_tensor(x[ii], device=device)
            yb = torch.as_tensor(y[ii], device=device)
            wb = torch.as_tensor(w[ii], device=device)
            loss = (binary_cross_entropy_with_logits(net(xb).squeeze(1), yb, reduction='none') * wb).sum() * len(x) / len(ii)
            if not torch.isfinite(loss):
                raise ValueError('Non-finite training loss')
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
        net.eval()
        score = 0.
        with torch.no_grad():
            for start in range(0, len(xv), bs):
                sl = slice(start, start + bs)
                logits = net(torch.as_tensor(xv[sl], device=device)).squeeze(1)
                score += float((binary_cross_entropy_with_logits(logits, torch.as_tensor(yv[sl], device=device), reduction='none') * torch.as_tensor(wv[sl], device=device)).sum())
        if not np.isfinite(score):
            raise ValueError('Non-finite validation loss')
        history.append(score)
        print(f'epoch {epoch + 1}: validation BCE {score:.8g}', flush=True)
        if score < best - 1e-6:
            best, bad = score, 0
            best_state = copy.deepcopy({k: v.detach().cpu() for k, v in net.state_dict().items()})
        else:
            bad += 1
            if bad >= config['patience']:
                break
    net.load_state_dict(best_state)
    return net, dict(state=best_state, mean=mu.tolist(), std=sd.tolist(), hidden=config['hidden'],
                     n_features=xs.shape[1], seed=seed, validation_loss=history, best_validation_loss=best)


def predict_logits(net, x, card, device, batch=8192):
    import numpy as np
    import torch
    out = np.empty(len(x), dtype='float64')
    mu, sd = np.asarray(card['mean'], dtype='float32'), np.asarray(card['std'], dtype='float32')
    net.eval()
    with torch.no_grad():
        for start in range(0, len(x), batch):
            sl = slice(start, start + batch)
            z = ((x[sl] - mu) / sd).astype('float32')
            out[sl] = net(torch.as_tensor(z, device=device)).squeeze(1).cpu().numpy()
    if not np.isfinite(out).all():
        raise ValueError('Non-finite predicted log ratio')
    return out


def normalized_ratios(logits):
    import numpy as np
    from scipy.special import logsumexp
    if len(logits) == 0:
        raise ValueError('No source events in the mass window')
    lognorm = float(logsumexp(logits) - np.log(len(logits)))
    weights = np.exp(logits - lognorm)
    if not np.isfinite(weights).all() or np.any(weights == 0):
        raise ValueError('Invalid/underflowed weights; no silent clipping is allowed')
    return weights, lognorm


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--target', type=Path, required=True)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--paper-package', type=Path, required=True)
    ap.add_argument('--config', type=Path, default=Path(__file__).with_name('config.json'))
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--device', default='cuda')
    ap.add_argument('--worker', action='store_true', help='Explicit acknowledgement: this is a compute worker')
    a = ap.parse_args()
    if not a.worker or platform.node().split('.')[0] == 'cmscon':
        ap.error('Training is allowed only on a compute worker; pass --worker there')
    if a.out.exists():
        ap.error('Output already exists; use a new directory')
    import numpy as np
    import torch
    sys.path.insert(0, str(a.paper_package.resolve() / 'rwt'))
    import ggzh_reweight_bdt as BDT
    import ggzh_reweight_common as C
    import ggzh_rwt_manifest as M
    from ggzh_rwt10m import diagnostics
    config = json.loads(a.config.read_text()); validate_config(config)
    manifest = M.load(a.manifest)
    if a.source.parent != a.target.parent or a.source.name != 'lo-cms.npz' or a.target.name != 'nlo-cp5.npz':
        ap.error('Use paired lo-cms.npz and nlo-cp5.npz in the same directory with their JSON bookkeeping')
    M.verify(manifest, a.source.parent)
    torch.set_num_threads(int(os.environ.get('RWT_THREADS', '4')))
    BDT.configure(C.FEATURE_SETS[config['features']])
    src, tgt = C.Sample('source', a.source), C.Sample('target', a.target)
    xs, xt = BDT.prepare(src).astype('float32'), BDT.prepare(tgt).astype('float32')
    if not np.isfinite(xs[src.in_window]).all() or not np.isfinite(xt[tgt.in_window]).all():
        raise ValueError('Non-finite in-window feature; fix input provenance, do not drop rows')
    a.out.mkdir(parents=True)
    source_idx = np.flatnonzero(src.in_window)
    fold = source_fold(np.arange(len(source_idx)), config['folds'])
    target_idx = np.arange(tgt.n)
    half_a = tgt.in_window & (target_half(target_idx) == 0)
    half_b = tgt.in_window & (target_half(target_idx) == 1)
    result = dict(method='fresh-dense-classifier', config=config, features=BDT.FEATURES,
                  input_selection=manifest, provenance=dict(paper_package=(a.paper_package / 'VERSION').read_text().strip(),
                  training_script_sha256=digest(__file__), config_sha256=digest(a.config),
                  source_sha256=digest(a.source), target_sha256=digest(a.target),
                  source_events=src.n, target_events=tgt.n, numpy=np.__version__, torch=torch.__version__,
                  negative_target_training='clip to zero; signed weights retained in closure'))
    weights = dict(run=src['run'], event=src['event'], in_window=src.in_window)
    for tag_i, (tag, train, test) in enumerate([('full', tgt.in_window, tgt.in_window), ('A_to_B', half_a, half_b), ('B_to_A', half_b, half_a)]):
        logits = np.empty(len(source_idx), dtype='float64')
        cards = []
        for k in range(5):
            si = source_idx[fold != k]
            net, card = fit(xs[si], xt[train], tgt.w0[train], config, config['seed'] + tag_i * 100 + k, a.device)
            logits[fold == k] = predict_logits(net, xs[source_idx[fold == k]], card, a.device, config['batch_size'])
            card.update(features=BDT.FEATURES, tag=tag, fold=k, config=config, input_selection=manifest)
            torch.save(card, a.out / f'{tag}-fold-{k}.pt')
            cards.append({key: value for key, value in card.items() if key != 'state'})
            del net
        w = np.ones(src.n)
        w[source_idx], lognorm = normalized_ratios(logits)
        weights['w_' + tag] = w
        result[tag] = BDT.evaluate(src, w, src.in_window, tgt, test, tag, print)
        result[tag]['models'] = cards
        result[tag]['log_normalization'] = lognorm
    r = result['full']['reweighted']
    result['summary'] = {'1_stage13': r['pass_1'], '2_variables': r['pass_2'], '2b_features': r['pass_2b'],
                         '3_stability': result['A_to_B']['pass'] and result['B_to_A']['pass']}
    result['summary']['all_pass'] = all(result['summary'].values())
    result['diagnostics'] = diagnostics(src, weights['w_full'], tgt)
    result['claim_reweighting_works'] = bool(result['summary']['all_pass'] and result['diagnostics']['hardjet']['pth450']['within_band'])
    name = config['features'].lower()
    np.savez_compressed(a.out / f'weights-clf-{name}.npz', **weights)
    (a.out / f'dnn-clf-{name}.json').write_text(C.dumps(result))
    (a.out / 'COMPLETED').write_text('All three directions evaluated; see JSON for physics pass/fail.\n')


if __name__ == '__main__':
    main()
