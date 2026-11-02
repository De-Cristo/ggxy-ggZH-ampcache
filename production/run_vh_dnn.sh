#!/bin/bash
# Worker-only GPU training; run/resources/inputs are defined in paper.sh.
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=${RWT_THREADS:-4}
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt
source "$VIEW/setup.sh" >/dev/null 2>&1
mkdir -p dnn-logs
# Keep the LCG libraries for its Python executable but isolate Python packages.
unset PYTHONPATH PYTHONHOME
python3 -m venv dnn-env
source dnn-env/bin/activate
python -m pip install --disable-pip-version-check --no-cache-dir \
    torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121 >dnn-logs/pip-torch.log 2>&1
python -m pip install --disable-pip-version-check --no-cache-dir \
    numpy==1.26.4 scipy==1.13.1 matplotlib==3.8.4 >dnn-logs/pip-analysis.log 2>&1
python -m pip freeze >dnn-logs/pip-freeze.txt
tar xzf vh-code.tgz
test "$(cat VH-COMMIT)" = "$VH_TRAINER_COMMIT"
sha256sum -c VH-SHA256SUMS
python - <<'PY'
import json, os, torch
assert torch.cuda.is_available(), 'CUDA unavailable: no CPU fallback'
assert torch.cuda.device_count() == 1, 'Expected exactly one assigned GPU'
p = torch.cuda.get_device_properties(0)
torch.ones(1, device='cuda').sum().item()
record = dict(torch=str(torch.__version__), cuda=torch.version.cuda,
              device=p.name, total_memory=p.total_memory,
              visible_devices=os.getenv('CUDA_VISIBLE_DEVICES'))
print(json.dumps(record), flush=True)
open('dnn-logs/gpu.json', 'w').write(json.dumps(record, indent=2)+'\n')
PY
VH_WORKER_TEST=1 python -m unittest discover -s vh/tests -v 2>&1 | tee dnn-logs/smoke.log
/usr/bin/time -p -o dnn-logs/train-time.txt python vh/train.py --worker --device cuda \
    --paper-package "$PWD" --manifest "$RWT_INPUT_MANIFEST" \
    --source rwt-out/lo-cms.npz --target rwt-out/nlo-cp5.npz \
    --config vh/config.json --out dnn-out 2>&1 | tee dnn-logs/train.log
python production/verify_vh_dnn.py
cp dnn-logs/gpu.json dnn-logs/pip-freeze.txt VH-COMMIT VH-SHA256SUMS dnn-out/
tar czf dnn-models.tgz dnn-out/*.pt dnn-out/gpu.json dnn-out/pip-freeze.txt \
    dnn-out/VH-COMMIT dnn-out/VH-SHA256SUMS vh/config.json dnn-out/COMPLETED
tar czf out-rwt-dnn-c-10m.tgz dnn-out/*.json dnn-logs VH-COMMIT VH-SHA256SUMS
