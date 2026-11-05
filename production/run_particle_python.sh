#!/bin/bash
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
source shower/versions.sh
source "$SHOWER_LCG_VIEW/setup.sh" >/dev/null 2>&1
python3 - <<'PY'
import hashlib,json,urllib.request
from pathlib import Path
m=json.loads(Path('shower/particle/pyarrow-wheel.json').read_text())
path=Path(m['filename']);urllib.request.urlretrieve(m['url'],path)
assert path.stat().st_size==m['size']
assert hashlib.sha256(path.read_bytes()).hexdigest()==m['sha256']
PY
source shower/particle/python_env.sh
cp particle-logs/python-preflight.json python-codec.json
