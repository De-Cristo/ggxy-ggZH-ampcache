#!/bin/bash
# Real 1M/class study only; all reads/training/checks execute on Condor workers.
set -eo pipefail
[[ $(hostname -s) != cmscon && -n ${_CONDOR_SCRATCH_DIR:-} ]] || exit 2
TASK=${1:?}; INDEX=${2:-0}
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=${RWT_THREADS:-4} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh >/dev/null 2>&1
unset PYTHONPATH PYTHONHOME
mkdir -p parton-logs
python3 -m venv parton-env
source parton-env/bin/activate
python -m pip install --no-cache-dir --disable-pip-version-check numpy==1.26.4 scipy==1.13.1 matplotlib==3.8.4 scikit-learn==1.2.2 >parton-logs/pip-analysis.log 2>&1
if [[ $TASK = fit || $TASK = ensemble-fit ]]; then
  python -m pip install --no-cache-dir --disable-pip-version-check torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121 >parton-logs/pip-torch.log 2>&1
fi
if [[ $TASK = bdt-deliver ]]; then
  python -m pip install --no-cache-dir --disable-pip-version-check skl2onnx==1.16.0 onnxconverter-common==1.14.0 onnx==1.16.2 onnxruntime==1.19.2 protobuf==3.20.2 >parton-logs/pip-onnx.log 2>&1
fi
python -m pip freeze >parton-logs/environment.txt
extract() {
  mkdir -p "$2"
  python - "$1" "$2" <<'PY'
from pathlib import Path
import sys,tarfile
with tarfile.open(sys.argv[1]) as t:
    assert all((m.name=='parton-data' or m.name.startswith('parton-data/')) and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers())
    t.extractall(sys.argv[2])
PY
  (cd "$2/parton-data" && sha256sum -c --quiet SHA256SUMS)
}
if [[ $TASK != prepare && $TASK != ensemble-prepare && $TASK != ensemble-assess && $TASK != bdt-assess && $TASK != pdf-transport && $TASK != bdt-deliver ]]; then extract parton-prepared.tgz .; fi
if [[ $TASK = bdt-deliver ]]; then extract eval-prepared.tgz eval; fi
if [[ $TASK = ensemble-fit || $TASK = ensemble-assess || $TASK = bdt-member || $TASK = bdt-assess ]]; then
  extract calib-prepared.tgz calib
  extract eval-prepared.tgz eval
fi
case $TASK in
  prepare) SCRIPT=rwt/parton/prepare.py; ARGS=() ;;
  fit) SCRIPT=rwt/parton/fit.py; ARGS=("$INDEX") ;;
  ensemble-prepare) SCRIPT=rwt/parton/prepare.py; ARGS=() ;;
  ensemble-fit) SCRIPT=rwt/parton/ensemble_fit.py; ARGS=() ;;
  ensemble-assess) SCRIPT=rwt/parton/ensemble_assess.py; ARGS=() ;;
  bdt-member) SCRIPT=rwt/parton/bdt_member.py; ARGS=() ;;
  bdt-assess) SCRIPT=rwt/parton/ensemble_assess.py; ARGS=() ;;
  pdf-transport) SCRIPT=rwt/parton/pdf_transport.py; ARGS=() ;;
  bdt-deliver) SCRIPT=rwt/parton/ensemble_deliver.py; ARGS=() ;;
  bdt) SCRIPT=rwt/parton/control.py; ARGS=() ;;
  jet-compare) SCRIPT=rwt/parton/jet_compare.py; ARGS=() ;;
  report) SCRIPT=rwt/parton/report.py; ARGS=() ;;
  *) exit 2 ;;
esac
/usr/bin/time -v -o parton-logs/resources.txt python "$SCRIPT" "${ARGS[@]}" 2>&1 | tee parton-logs/task.log
