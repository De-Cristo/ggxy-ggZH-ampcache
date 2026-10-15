#!/bin/bash
# Job steps of the LO+PS -> NLO+PS reweighting study (rwt/README.md), run in a job directory of
# production/paper.sh (runs rwt-*, set 13600) with the per-event tables of the table-* runs as inputs:
#
#   rwt/run.sh clf <name> <source> <target> [ggzh_reweight_bdt.py options, e.g. --jet-feature]
#        classifier weights of <source> towards <target>
#        -> weights-<name>.npz, rwt-out/bdt-<name>.json (products), out-rwt-<name>.tgz (result)
#   rwt/run.sh summary
#        closure of the binned map (ggzh_reweight_closure.py) and the summary (ggzh_rwt_summary.py), with the
#        weights and classifier outputs of the four rwt-clf* runs in rwt-in/
#        -> out-rwt-summary.tgz (rwt-summary.json, closure-map.json, the classifier outputs the paper draws)
#   rwt/run.sh model
#        the reusable models (ggzh_rwt_export.py): classifier with five and six inputs as ONNX and the binned map as a
#        correctionlib table, nominal and seven scale variations; the ONNX converter from the wheels in wheels/
#        (rwt/requirements-export.txt, checked hashes) -> model.tgz (product), out-rwt-model.tgz (card and log)
#
# Samples: tables/<sample>/table-<seed>.csv.gz -> rwt-out/<sample>.npz (ggzh_table.py; the NLO+PS target
# nlo-cp5 with the NNPDF3.1 central member as nominal weight, --nominal 2000).  Environment: the LCG_106 view
# (scikit-learn 1.2.2, the version the published weights were made with), single-threaded.
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
VIEW=${RWT_LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt}
[ -f "$VIEW/setup.sh" ] || { echo "FATAL: no $VIEW" >&2; exit 64; }
source "$VIEW/setup.sh" > /dev/null 2>&1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
mkdir -p rwt-out
npz() {  # npz <sample>
  [ -f "rwt-out/$1.npz" ] && return 0
  local nom=nominal; [ "$1" = nlo-cp5 ] && nom=2000
  ls tables/"$1"/table-*.csv.gz > /dev/null 2>&1 || { echo "FATAL: no tables/$1/table-*.csv.gz" >&2; return 2; }
  python3 "$HERE/ggzh_table.py" --label "$1" --tables "tables/$1/table-*.csv.gz" --out "rwt-out/$1.npz" --nominal "$nom"
}
WHAT=${1:?clf|summary|model}; shift
case $WHAT in
  clf)
    NAME=${1:?name}; SRC=${2:?source sample}; TGT=${3:?target sample}; shift 3
    npz "$SRC" && npz "$TGT" || exit 2
    python3 "$HERE/ggzh_reweight_bdt.py" --source "rwt-out/$SRC.npz" --target "rwt-out/$TGT.npz" --out "rwt-out/bdt-$NAME.json" \
        --weights-out "weights-$NAME.npz" "$@" 2>&1 | tee "rwt-out/bdt-$NAME.log"
    rc=${PIPESTATUS[0]}
    [ -s "rwt-out/bdt-$NAME.json" ] && [ -s "weights-$NAME.npz" ] || rc=10
    tar czf "out-rwt-$NAME.tgz" rwt-out/*.json rwt-out/*.log
    exit $rc ;;
  summary)
    for s in lo-cms lo-cms-default nlo-cp5; do npz "$s" || exit 2; done
    for f in rwt-in/*; do ln -sf "$PWD/$f" "rwt-out/$(basename "$f")"; done
    python3 "$HERE/ggzh_reweight_closure.py" --official rwt-out/lo-cms.npz --nlo rwt-out/nlo-cp5.npz --out rwt-out/closure-map.json \
        2>&1 | tee rwt-out/closure-map.log
    rc=${PIPESTATUS[0]}
    [ $rc = 0 ] && python3 "$HERE/ggzh_rwt_summary.py" --dir rwt-out --out rwt-out/rwt-summary.json --grids-xsec rwt-out/grids-xsec.json \
        2>&1 | tee rwt-out/summary.log
    [ $rc = 0 ] && rc=${PIPESTATUS[0]}
    [ -s rwt-out/rwt-summary.json ] || rc=10
    tar czhf out-rwt-summary.tgz rwt-out/rwt-summary.json rwt-out/closure-map.json rwt-out/bdt-clf*.json rwt-out/grids-xsec.json \
        rwt-out/lo-cms.json rwt-out/lo-cms-default.json rwt-out/nlo-cp5.json rwt-out/*.log
    exit $rc ;;
  model)
    for s in lo-cms nlo-cp5; do npz "$s" || exit 2; done
    mkdir -p pylib
    python3 -m pip install --no-index --no-deps --require-hashes --find-links wheels --target pylib -r "$HERE/requirements-export.txt" \
        > rwt-out/pip.log 2>&1 || { cat rwt-out/pip.log; exit 2; }
    export PYTHONPATH=$PWD/pylib${PYTHONPATH:+:$PYTHONPATH}
    python3 "$HERE/ggzh_rwt_export.py" --source rwt-out/lo-cms.npz --target rwt-out/nlo-cp5.npz --out-dir model \
        --compare-weights rwt-in/weights-clf.npz 2>&1 | tee rwt-out/export.log
    rc=${PIPESTATUS[0]}
    [ -s model/model-card.json ] || rc=10
    tar czf model.tgz model; cp -p model/model-card.json rwt-out/
    tar czf out-rwt-model.tgz rwt-out/model-card.json rwt-out/export.log rwt-out/pip.log rwt-out/lo-cms.json rwt-out/nlo-cp5.json
    exit $rc ;;
  *) echo "unknown step $WHAT" >&2; exit 2 ;;
esac
