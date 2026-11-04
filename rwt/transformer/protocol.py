"""Small deterministic split/config helpers; no data or ML imports."""
import math


def source_fold(window_index, folds=5):
    return window_index % folds


def target_half(original_index):
    return original_index % 2


def validate_config(c):
    if c['features'] not in ('A', 'B', 'C') or c['folds'] != 5:
        raise ValueError('Use the registered A/B/C inputs and five source folds')
    for k in ('epochs', 'patience', 'batch_size'):
        if type(c[k]) is not int or c[k] < 1:
            raise ValueError(k)
    if not c['hidden'] or any(type(n) is not int or n < 1 for n in c['hidden']):
        raise ValueError('hidden')
    for k in ('learning_rate', 'weight_decay', 'validation_fraction'):
        if not math.isfinite(c[k]) or c[k] < 0:
            raise ValueError(k)
    if not 0 < c['validation_fraction'] < .5 or c['learning_rate'] == 0:
        raise ValueError('learning_rate/validation_fraction')
