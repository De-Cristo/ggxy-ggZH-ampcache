#!/usr/bin/env python3
"""Worker-only pilot: compare appended event-table columns to independent GGZH_NLOPS objects."""
import json
import sys
import numpy as np

a = np.genfromtxt(sys.argv[1], delimiter=',', names=True)
b = np.genfromtxt(sys.argv[2], delimiter=',', names=True)
assert len(a) == len(b) and len(a) > 0
out = {'events': len(a), 'checks': {}}
for key in b.dtype.names:
    name = 'ptj1eta47' if key == 'ptj1' else key
    assert np.allclose(a[name], b[key], rtol=6e-6, atol=6e-6, equal_nan=True), key
    out['checks'][key] = float(np.nanmax(np.abs(a[name]-b[key])))
assert np.all((a['ptj2'] == -1) == (a['njets30eta47'] < 2))
assert np.all(a['ptl1'] >= a['ptl2'])
assert np.all(a['ptl2'] > 0), 'missing charged Z child'
out['pass'] = True
print(json.dumps(out, indent=1))
