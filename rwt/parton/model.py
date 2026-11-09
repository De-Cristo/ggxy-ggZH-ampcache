"""Fresh particle-attention ZH classifier, with event/shower context after class attention."""
import sys
from pathlib import Path
import torch
from torch import nn
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'transformer/vendor'))
from ParticleTransformer import ParticleTransformer

class ZHParT(nn.Module):
    def __init__(self,config,normalization):
        super().__init__();cfg=config['model']
        self.register_buffer('particle_mean',torch.tensor(normalization['particle_mean'],dtype=torch.float32).reshape(1,7,1))
        self.register_buffer('particle_std',torch.tensor(normalization['particle_std'],dtype=torch.float32).reshape(1,7,1))
        self.register_buffer('global_mean',torch.tensor(normalization['global_mean'],dtype=torch.float32).reshape(1,24))
        self.register_buffer('global_std',torch.tensor(normalization['global_std'],dtype=torch.float32).reshape(1,24))
        self.backbone=ParticleTransformer(input_dim=7,num_classes=None,fc_params=None,pair_input_dim=4,
            embed_dims=cfg['embed_dims'],pair_embed_dims=cfg['pair_embed_dims'],num_heads=cfg['num_heads'],
            num_layers=cfg['num_layers'],num_cls_layers=cfg['num_cls_layers'],trim=False,use_amp=False)
        self.head=nn.Sequential(nn.Linear(cfg['embed_dims'][-1]+24,128),nn.GELU(),nn.Dropout(.1),nn.Linear(128,1))
    def forward(self,x,v,mask,global_features):
        x=((x-self.particle_mean)/self.particle_std).masked_fill(~mask,0)
        embedding=self.backbone(x,v=v,mask=mask)
        context=(global_features-self.global_mean)/self.global_std
        return self.head(torch.cat([embedding,context],dim=1)).squeeze(1)
