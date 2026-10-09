#!/bin/bash
# Combine the job results of one set of production/paper.sh into the result files the paper
# quotes, and run the checks.  Needs bash, GNU tar and Python 3 with numpy; minutes, one core.
#
#   production/analyze.sh <set> <results> <outdir>
#
# <results>/<run>/<p>/ holds what job p of <run> left in its job directory: out-*.tgz,
# <run>-<p>.raw.json (showers), VERSION.  run_local.sh's <workdir>/jobs/<set> has this layout;
# only successful jobs belong there.  A run is used only when all its jobs are present;
# <outdir>/STATUS.txt lists every run with its jobs and the package commits they ran with.
#
# Output in <outdir> (fb unless stated otherwise):
#   fo-llh.json                 ggxy l-l+H (l = e, mu, tau) LO, NLO, K and differential K
#   fo-zh.json                  ggxy stable Z (13600, 13000: LO and NLO)
#   fo-zh-nnpdf31.json, fo-zh-pdf4lhc21.json      (verify: stable-Z LO with each PDF)
#   powheg-lo.json              POWHEG LO at the seven scale points
#   grids-<name>.json           the integral of each set of grids
#   mcfm-<card>.json            MCFM LO, every card
#   events-nlo.json, events-lo.json   cross sections, scale points, PDF and alpha_s from the
#                               event weights (PDF4LHC21; verify: NNPDF3.1)
#   events-nlo-nnpdf31.json, events-lo-nnpdf31.json   (13600, 13000: PDF uncertainty with the
#                               NNPDF3.1 weights of the same events)
#   events-nlo-bs250.json       (verify: weighted events with Born suppression)
#   cache/compare-<seed>.json   cache closure
#   shower-<variant>.json       showers, exact sums over all jobs (pb, analysis/yoda_summary.py)
#   check-shower-<variant>.json shower/check_shower.py of each shower against its events
#   lops-<tune>.json            rescaled LO+PS against NLO+PS (analysis/compare_orders.py)
#   bornsupp-shower.json        (verify: NLO+PS from weighted against unweighted events)
#   check-generator.json        generator/check.py (13600, 13000: --paper-energy)
#   SHA256SUMS                  of every file above
# Exit status 1 if a step or a check failed (listed at the end of STATUS.txt).
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); PKG=$(cd "$HERE/.." && pwd)
SET=${1:?set}; RES=${2:?results directory}; OUT=${3:?output directory}
RES=$(cd "$RES" && pwd) || exit 2
mkdir -p "$OUT/work" "$OUT/cache" || exit 2; OUT=$(cd "$OUT" && pwd); W=$OUT/work
GA=$PKG/generator/analysis; SA=$PKG/shower/analysis
if [ "$SET" = verify ]; then PDFID=325300; NMEM=100; PDFW=2000      # pwg-rwl.dat: NNPDF3.1
else PDFID=93300; NMEM=40; PDFW=3000; fi                             # pwg-rwl-pdf4lhc21.dat: PDF4LHC21

declare -A NJOBS; RUNS=()
while read -r r j _; do NJOBS[$r]=$j; RUNS+=("$r"); done < <("$HERE/paper.sh" list "$SET")
[ ${#RUNS[@]} -gt 0 ] || { echo "FATAL: no runs in set $SET" >&2; exit 2; }
present() { local n=0 p; for p in $(seq 0 $((NJOBS[$1] - 1))); do [ -d "$RES/$1/$p" ] && n=$((n + 1)); done; echo $n; }

# ---------------------------------------------------------------- status of the results
{
  echo "production/analyze.sh set $SET, package $(head -1 "$PKG/VERSION")"
  echo "results $RES"
  for r in "${RUNS[@]}"; do
    n=$(present "$r")
    c=$(cat "$RES/$r"/*/VERSION 2>/dev/null | awk '{print substr($1, 1, 12)}' | sort -u | tr '\n' ' ')
    printf '%-34s %4d/%-4d %-10s package %s\n' "$r" "$n" "${NJOBS[$r]}" \
      "$([ "$n" = "${NJOBS[$r]}" ] && echo complete || echo INCOMPLETE)" "${c:-?}"
  done
} > "$OUT/STATUS.txt"
cat "$OUT/STATUS.txt"

FAILED=()
step() {  # step <name> <command ...>: run it, remember a failure
  local name=$1; shift
  echo "--- $name"
  "$@"; local rc=$?
  [ $rc = 0 ] || { echo "*** $name failed (exit $rc)"; FAILED+=("$name"); }
  return $rc
}
# have <run>: the run is in the set and all its jobs are present
have() {
  [ -n "${NJOBS[$1]:-}" ] || return 1
  [ "$(present "$1")" = "${NJOBS[$1]}" ] || { echo "--- $1: skipped, incomplete"; return 1; }
}
# unpack <run> <member patterns ...>: the needed files of every job's out-*.tgz into work/<run>/<p>
unpack() {
  local r=$1 p t; shift
  [ -f "$W/$r/.unpacked" ] && return 0
  for p in $(seq 0 $((NJOBS[$r] - 1))); do
    mkdir -p "$W/$r/$p"
    for t in "$RES/$r/$p"/out-*.tgz; do [ -f "$t" ] && tar xzf "$t" -C "$W/$r/$p" --wildcards "$@" 2>/dev/null; done
  done
  touch "$W/$r/.unpacked"
}

# ---------------------------------------------------------------- fixed order (ggxy)
fo() {  # fo <output> <nlep> <LO run> [<NLO run>]
  local out=$1 nlep=$2 args=()
  if have "$3"; then unpack "$3" 'out/histo_*'; args+=(--lo "$W/$3"/*/out/histo_lo_*.txt); fi
  if [ -n "${4:-}" ] && have "$4"; then unpack "$4" 'out/histo_*'; args+=(--nlo "$W/$4"/*/out/histo_nlo_*.txt); fi
  [ ${#args[@]} = 0 ] || step "$out" python3 "$GA/fo_combine.py" "${args[@]}" --nlep "$nlep" --label "$SET $out" --out "$OUT/$out.json"
}
fo fo-llh 3 fo-llh-lo fo-llh-nlo                      # fo-ggllh is for one lepton flavour
if [ "$SET" = verify ]; then
  fo fo-zh-nnpdf31 1 fo-zh-lo
  fo fo-zh-pdf4lhc21 1 fo-zh-lo-pdf4lhc21
else
  fo fo-zh 1 fo-zh-lo fo-zh-nlo
fi

# ---------------------------------------------------------------- POWHEG LO, grids, MCFM
if have powheg-lo; then
  unpack powheg-lo 'out/*'
  pts=(); for d in "$W"/powheg-lo/0/out/*/; do d=${d%/}; pts+=("$(basename "$d")=$d"); done
  step powheg-lo python3 "$GA/powheg_stat.py" --points "${pts[@]}" --label "$SET powheg-lo" --out "$OUT/powheg-lo.json"
fi
for r in grids-nlo grids-lo grids-nlo-bs250; do
  have "$r" || continue
  unpack "$r" 'grids/grids-xsec.json' 'grids/manifest.txt'
  step "$r" cp "$W/$r/0/grids/grids-xsec.json" "$OUT/$r.json"
done
if have mcfm; then unpack mcfm 'out/*.json'; step mcfm cp "$W"/mcfm/0/out/mcfm-*.json "$OUT/"; fi

# ---------------------------------------------------------------- events, cache closure
events() {  # events <run> <output> <LHAPDF id> <Hessian members>
  have "$1" || return 0
  unpack "$1" 'out/xsec-raw-*.json'
  step "$2" python3 "$GA/lhe_xsec.py" --merge "$W/$1"/*/out/xsec-raw-*.json --pdf-id "$3" --pdf-nmem "$4" \
    --label "$SET $1" --out "$OUT/$2.json"
}
events events-nlo events-nlo $PDFID $NMEM
events events-lo events-lo $PDFID $NMEM
if [ "$SET" = verify ]; then
  events events-nlo-bs250 events-nlo-bs250 $PDFID $NMEM
else
  events events-nlo events-nlo-nnpdf31 325300 100
  events events-lo events-lo-nnpdf31 325300 100
fi
if have cache; then unpack cache 'out/compare-*.json'; step cache cp "$W"/cache/*/out/compare-*.json "$OUT/cache/"; fi

# ---------------------------------------------------------------- showers
for r in "${RUNS[@]}"; do
  case $r in shower-*) ;; *) continue ;; esac
  have "$r" || continue
  x=${r#shower-}
  step "$r" python3 "$SA/yoda_summary.py" --merge "$RES/$r"/*/"$r"-*.raw.json --out "$OUT/shower-$x.json" || continue
  case $x in nlo-bs250-*) l=events-nlo-bs250 ;; nlo-*) l=events-nlo ;; *) l=events-lo ;; esac
  [ -f "$OUT/$l.json" ] && step "check-shower-$x" python3 "$PKG/shower/check_shower.py" --summary "$OUT/shower-$x.json" \
    --lhe "$OUT/$l.json" --pdf-weight $PDFW --pdf-nmem $NMEM --out "$OUT/check-shower-$x.json"
done
# rescaled LO+PS against NLO+PS, per tune, with every LO+PS shower start present
for f in "$OUT"/shower-nlo-*.json; do
  [ -f "$f" ] || continue
  t=${f##*/shower-nlo-}; t=${t%.json}
  case $t in bs250-*) continue ;; esac
  los=(); for g in "$OUT"/shower-lo-*-"$t".json; do [ -f "$g" ] && los+=("$g"); done
  [ ${#los[@]} = 0 ] || step "lops-$t" python3 "$SA/compare_orders.py" --nlo "$f" --lo "${los[@]}" --out "$OUT/lops-$t.json"
done
if [ -f "$OUT/shower-nlo-bs250-monash.json" ] && [ -f "$OUT/shower-nlo-monash.json" ]; then
  step bornsupp-shower python3 "$SA/compare_samples.py" --a "$OUT/shower-nlo-bs250-monash.json" \
    --b "$OUT/shower-nlo-monash.json" --out "$OUT/bornsupp-shower.json"
fi

# ---------------------------------------------------------------- generator checks
C=()
opt() { [ -f "$OUT/$2" ] && C+=("$1" "$OUT/$2"); }
opt --fo fo-llh.json; opt --powheg-lo powheg-lo.json; opt --events events-nlo.json; opt --events-lo events-lo.json
if [ "$SET" = verify ]; then
  opt --stable-nnpdf31 fo-zh-nnpdf31.json; opt --stable-pdf4lhc21 fo-zh-pdf4lhc21.json; opt --grids grids-nlo.json
  opt --events-bornsupp events-nlo-bs250.json
  opt --mcfm-stable-nnpdf31 mcfm-ggzh-lo-stable-nnpdf31.json; opt --mcfm-stable-pdf4lhc21 mcfm-ggzh-lo-stable-pdf4lhc21.json
else
  # the production NLO grids are integrated with Born suppression; the events give the cross section
  opt --stable-pdf4lhc21 fo-zh.json
fi
cc=("$OUT"/cache/compare-*.json); [ -f "${cc[0]}" ] && C+=(--cache "${cc[@]}")
if [ ${#C[@]} -gt 0 ]; then
  [ "$SET" = verify ] || C+=(--paper-energy "$SET")
  step check-generator python3 "$PKG/generator/check.py" "${C[@]}" --out "$OUT/check-generator.json"
fi

rm -rf "$W"
( cd "$OUT" && find . -name '*.json' | sort | xargs -r sha256sum > SHA256SUMS )
echo "failed: ${FAILED[*]:-none}" >> "$OUT/STATUS.txt"
echo "=== analyze.sh $SET: ${#FAILED[@]} failed step(s)${FAILED:+: ${FAILED[*]}}"
[ ${#FAILED[@]} = 0 ]
