#!/bin/bash
# Source on a worker after the production LCG view. The public wheel bundles ZSTD.
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || return 2
mkdir -p particle-logs
PARTICLE_WHEEL=$(python3 -c 'import json;print(json.load(open("shower/particle/pyarrow-wheel.json"))["filename"])')
if [[ ! -f particle-python/.ready ]]; then
  python3 -m pip install --no-index --no-deps --disable-pip-version-check --target particle-python "$PARTICLE_WHEEL" >particle-logs/pip-arrow.log 2>&1 || return
  touch particle-python/.ready
fi
export PYTHONPATH=$PWD/particle-python${PYTHONPATH:+:$PYTHONPATH}
export LD_LIBRARY_PATH=$PWD/particle-python/pyarrow${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
python3 - <<'PY' >particle-logs/python-preflight.json
import json,pyarrow as pa,pyarrow.parquet as pq
assert pa.__version__=='20.0.0' and pa.Codec.is_available('zstd'), 'Required ZSTD codec unavailable'
print(json.dumps(dict(pyarrow=pa.__version__,module=pa.__file__,zstd=pa.Codec.is_available('zstd')),indent=2))
PY
