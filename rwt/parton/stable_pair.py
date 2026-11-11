"""Versioned rapidity from retained mass/pT; no changes to stored four-vectors."""
import torch
from ParticleTransformer import delta_r2, to_m2

CONVENTION = 'transverse-mass-rapidity-v1'

def rapidity(features, vectors, mask):
    # Feature 5 is log1p(m/GeV), evaluated before float32 E/pz rounding.
    mass = torch.expm1(features[:, 5:6].double())
    pt2 = vectors[:, :2].double().square().sum(dim=1, keepdim=True)
    mt2 = pt2 + mass.square()
    # Padding has no physical rapidity. Active undefined tokens fail explicitly.
    if not bool(torch.all(torch.isfinite(mt2) & ((mt2 > 0) | ~mask))):
        raise ValueError('Nonfinite or zero active transverse mass')
    y = torch.asinh(vectors[:, 2:3].double() / torch.sqrt(mt2.masked_fill(~mask, 1)))
    y = y.masked_fill(~mask, 0).to(vectors.dtype)
    if not bool(torch.isfinite(y).all()):
        raise ValueError('Nonfinite transverse-mass rapidity')
    return y

def pairwise(xi, xj):
    # Fifth coordinate is separately computed rapidity, not a fifth momentum.
    vi, vj = xi[:, :4], xj[:, :4]
    pti = vi[:, :2].square().sum(dim=1, keepdim=True).sqrt()
    ptj = vj[:, :2].square().sum(dim=1, keepdim=True).sqrt()
    phii, phij = torch.atan2(vi[:, 1:2], vi[:, :1]), torch.atan2(vj[:, 1:2], vj[:, :1])
    delta = delta_r2(xi[:, 4:5], phii, xj[:, 4:5], phij).sqrt()
    eps = 1e-8
    ptmin = torch.minimum(pti, ptj)
    return torch.cat((torch.log((ptmin*delta).clamp(min=eps)),
        torch.log((ptmin/(pti+ptj).clamp(min=eps)).clamp(min=eps)),
        torch.log(delta.clamp(min=eps)), torch.log(to_m2(vi+vj, eps=eps))), dim=1)
