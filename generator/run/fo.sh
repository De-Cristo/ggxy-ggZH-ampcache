#!/bin/bash
# Fixed-order LO or NLO with the ggxy Monte Carlo (bin/fo-ggllh or bin/fo-ggzh): <nseeds>
# independent seeds in parallel on this node, seeds <first>, <first>+1, ...
#
#   GGZH_PREFIX=<prefix> run/fo.sh <lo|nlo> <outdir> <first seed> <nseeds> [np] [np_opt] [pdf] [program] [sqrt(s) GeV]
#                                  [mtscheme] [mu_t/m_ZH]
#
# Defaults: np = 2e6 (LO) / 1e6 (NLO), np_opt = 2e5 / 1.5e5, NNPDF3.1, fo-ggllh, 13600 GeV.
# Reference statistics: LO 8 seeds (2001-2008), NLO 40 seeds (3001-3040); about 950 s per LO
# seed and 5000 s per NLO seed on one core.  Each seed writes <outdir>/histo_<mode>_<seed>.txt
# with the 7-point scale grid; combine with analysis/fo_combine.py.
# The stable-Z cross check: program fo-ggzh, LO, 2 seeds (5001, 5002), with NNPDF3.1 and with
# PDF4LHC21_40_pdfas (combine with --nlep 1).
# Top-mass scheme: mtscheme 0 (on-shell, default) or 1 (MS-bar, m_t(mu_t) with mu_t = <mu_t/m_ZH> x m_ZH;
# 0 or empty: mu_t = m_t(m_t)); passed to the program only when given.
set -o pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
MODE=${1:?mode lo|nlo}; OUT=${2:?outdir}; FIRST=${3:?first seed}; NS=${4:?nseeds}
case $MODE in lo) NP=${5:-2000000}; NPOPT=${6:-200000};; nlo) NP=${5:-1000000}; NPOPT=${6:-150000};; *) die "mode lo|nlo";; esac
source "$GEN/versions.sh"
PDF=${7:-$PDF_MAIN}; PROG=${8:-fo-ggllh}; ECM=${9:-13600}
MT=(); [ -n "${10:-}" ] && MT=("${10}" "${11:-0}")
load_env
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
manifest "$OUT/manifest.txt" "fo.sh $PROG $MODE seeds $FIRST..+$NS np $NP np_opt $NPOPT pdf $PDF sqrt(s) $ECM${MT:+ mtscheme ${MT[0]} mu_t/m_ZH ${MT[1]}}"
pids=()
for k in $(seq 0 $((NS - 1))); do
  s=$((FIRST + k))
  ( t=$SECONDS; "$PROG" "$MODE" "$s" "$OUT/" "$NP" "$NPOPT" "$PDF" "$ECM" "${MT[@]}" > "$OUT/run-$MODE-$s.log" 2>&1; rc=$?
    echo "seed $s rc $rc wall $((SECONDS - t)) s" >> "$OUT/timing.txt"; exit $rc ) &
  pids+=($!)
done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
n=$(ls "$OUT"/histo_"$MODE"_*.txt 2>/dev/null | wc -l)
echo "$(stamp) $PROG $MODE: $n histogram files in $OUT"
[ "$rc" = 0 ] && [ "$n" -ge "$NS" ] || die "some seeds failed, see $OUT/run-*.log"
