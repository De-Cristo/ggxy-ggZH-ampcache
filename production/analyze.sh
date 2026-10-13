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
#   cache/compare-<seed>.json   cache closure (file names relative to the job)
#   fo-llh-mt<f>.json, fo-zh-mt<f>.json   fixed order in the MS-bar top-mass scheme, mu_t = m_ZH / f, f = 1, 2, 4
#                               (fo-zh-mt<f>: 13600 only); ratios to fo-llh.json and fo-zh.json give the scheme band
#   fo-llh-pdf.json, fo-zh-pdf.json   (13600) PDF and alpha_s at fixed order (LO), one run per PDF
#                               member: check of the event-weight uncertainties (analysis/fo_pdf.py)
#   cache-scan.json             (13600) the cache timing scan: amplitude evaluations per event and
#                               relative times against the number of weights (analysis/cache_scan.py)
#   cache/speedup.json          per seed, same events with and without the amplitude cache:
#                               "evaluations" = one-loop amplitude evaluations without / with the cache
#                               (from the cache counters: (hits + misses) / misses; independent of the
#                               machine), "wall" = wall time without / with (both runs in the same job;
#                               depends on the load of the machine)
#   shower-<variant>.json       showers, exact sums over all jobs (pb, analysis/yoda_summary.py)
#   check-shower-<variant>.json shower/check_shower.py of each shower against its events
#   lops-<tune>.json            rescaled LO+PS against NLO+PS (analysis/compare_orders.py)
#   dists-nlo.json              (13600, 13000) NLO+PS (Monash) p_T(Z), m_ZH, p_T(ZH), N_jets with the
#                               variations for bands (analysis/distributions.py)
#   stxs-nlo.json               (13600, 13000) NLO+PS (Monash) in the STXS fine bins with scale, PDF,
#                               alpha_s, shower-weight and tune uncertainties (analysis/stxs_table.py)
#   parton-nlo.json, parton-lo.json   the NLO and LO events at parton level (no shower; summaries)
#   check-parton-<order>.json   shower/check_shower.py --noshower of each against its events
#   effect-<tune>.json          fixed order vs parton level vs shower at each order, and the high-pT
#                               test of the NLO events (analysis/shower_effect.py)
#   bornsupp-shower.json        (verify: NLO+PS from weighted against unweighted events)
#   bornsupp-parton.json        (verify: the same at parton level)
#   stxs-nlo-<v>.json, veto-<v>-shower.json   (13600) POWHEG-veto variations v = pthard1, ptdef2
#   stxs-nlo-hdampinf.json, hdamp-shower.json, effect-hdampinf-monash.json   (13600) the matching
#                               variation h_damp = infinity: STXS table, against h_damp = 250 GeV bin by
#                               bin, parton level against fixed order
#   effect-bs250-monash.json    (verify: the weighted NLO events against fixed order; effect-<tune>.json
#                               holds the unweighted ones)
#   check-generator.json        generator/check.py (13600, 13000: --paper-energy)
#   SHA256SUMS                  of every file above
# The outputs record no absolute paths (file names only), so they can be published as they are.
# Exit status 1 if a step or a check failed (listed at the end of STATUS.txt).
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); PKG=$(cd "$HERE/.." && pwd)
SET=${1:?set}; RES=${2:?results directory}; OUT=${3:?output directory}
RES=$(cd "$RES" && pwd) || exit 2
mkdir -p "$OUT/work" "$OUT/cache" || exit 2; OUT=$(cd "$OUT" && pwd); W=$OUT/work
cd "$OUT" || exit 2                     # inputs of the comparison scripts by relative name
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
# top-mass scheme: MS-bar with mu_t = m_ZH / f (l l H at both energies, stable Z at 13.6 TeV)
for f in 1 2 4; do
  fo "fo-llh-mt$f" 3 "fo-llh-lo-mt$f" "fo-llh-nlo-mt$f"
  fo "fo-zh-mt$f" 1 "fo-zh-lo-mt$f" "fo-zh-nlo-mt$f"
done

# ---------------------------------------------------------------- POWHEG LO, grids, MCFM
if have powheg-lo; then
  unpack powheg-lo 'out/*'
  pts=(); for d in "$W"/powheg-lo/0/out/*/; do d=${d%/}; pts+=("$(basename "$d")=$d"); done
  step powheg-lo python3 "$GA/powheg_stat.py" --points "${pts[@]}" --label "$SET powheg-lo" --out "$OUT/powheg-lo.json"
fi
for r in grids-nlo grids-lo grids-nlo-bs250; do
  have "$r" || continue
  unpack "$r" 'grids/grids-xsec.json' 'grids/manifest.txt'
  if [ ! -f "$W/$r/0/grids/grids-xsec.json" ] && grep -q bornsuppfact "$W/$r/0/grids/manifest.txt" 2>/dev/null; then
    # before 2026-10-10 analysis/powheg_stat.py did not read the Born-suppressed total, so those grid
    # jobs wrote no grids-xsec.json; that integral (sigma x F) is not a cross section and not needed
    echo "--- $r: Born-suppressed grids without grids-xsec.json (integral of sigma x F, not used)"
    continue
  fi
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
  events events-nlo-hdampinf events-nlo-hdampinf $PDFID $NMEM
  events events-nlo events-nlo-nnpdf31 325300 100
  events events-lo events-lo-nnpdf31 325300 100
fi
if have cache; then
  unpack cache 'out/compare-*.json' 'out/timing.txt'
  for f in "$W"/cache/*/out/compare-*.json; do
    sed -E 's#"(a|b)": "[^"]*/(cached|uncached)/#"\1": "\2/#' "$f" > "$OUT/cache/${f##*/}" || FAILED+=(cache)
  done
  # timing.txt: "seed <s> cache=<0|1> events <n> wall <t>s per_event <x>s ggxy ampcache: ... hits=<h> misses=<m> ..."
  cat "$W"/cache/*/out/timing.txt 2>/dev/null | awk '
    $3 == "cache=1" { c[$2] = $7 + 0; for (i = 1; i <= NF; i++) { split($i, kv, "="); if (kv[1] == "hits") h[$2] = kv[2]; if (kv[1] == "misses") m[$2] = kv[2] } }
    $3 == "cache=0" { u[$2] = $7 + 0 }
    END { printf "{"; n = 0
          for (s in c) if (s in u && c[s] > 0 && m[s] > 0)
            printf "%s\n \"%s\": {\"evaluations\": %.3f, \"wall\": %.3f}", (n++ ? "," : ""), s, (h[s] + m[s]) / m[s], u[s] / c[s]
          print "\n}" }' > "$OUT/cache/speedup.json"
fi

for r in fo-llh-pdf fo-zh-pdf; do
  have "$r" || continue
  unpack "$r" 'out/m*/histo_*'
  step "$r" python3 "$GA/fo_pdf.py" --members "$W/$r"/*/out/m* --nlep "$([ "$r" = fo-llh-pdf ] && echo 3 || echo 1)" \
    --pdf-nmem $NMEM --label "$SET $r" --out "$r.json"
done
if have cache-scan; then
  unpack cache-scan 'out/w*/timing.txt' 'out/w*/compare-*.json'
  step cache-scan python3 "$GA/cache_scan.py" --jobs "$W"/cache-scan/*/out --out cache-scan.json
fi

# ---------------------------------------------------------------- showers
for r in "${RUNS[@]}"; do
  case $r in shower-nlo-*|shower-lo-*|parton-*) ;; *) continue ;; esac
  have "$r" || continue
  case $r in
    parton-*) sum=$r; chk=check-$r; ns=(--noshower); l=events-${r#parton-} ;;
    *) x=${r#shower-}; sum=shower-$x; chk=check-shower-$x; ns=()
       case $x in nlo-bs250-*) l=events-nlo-bs250 ;; nlo-hdampinf-*) l=events-nlo-hdampinf ;; nlo-*) l=events-nlo ;; *) l=events-lo ;; esac ;;
  esac
  step "$r" python3 "$SA/yoda_summary.py" --merge "$RES/$r"/*/"$r"-*.raw.json --out "$OUT/$sum.json" || continue
  [ -f "$OUT/$l.json" ] && step "$chk" python3 "$PKG/shower/check_shower.py" --summary "$OUT/$sum.json" \
    --lhe "$OUT/$l.json" "${ns[@]}" --pdf-weight $PDFW --pdf-nmem $NMEM --out "$OUT/$chk.json"
done
# rescaled LO+PS against NLO+PS, per tune, with every LO+PS shower start present
for f in "$OUT"/shower-nlo-*.json; do
  [ -f "$f" ] || continue
  t=${f##*/shower-nlo-}; t=${t%.json}
  case $t in bs250-*) continue ;; esac
  los=(); for g in "$OUT"/shower-lo-*-"$t".json; do [ -f "$g" ] && los+=("$g"); done
  [ ${#los[@]} = 0 ] || step "lops-$t" python3 "$SA/compare_orders.py" --nlo "${f##*/}" --lo "${los[@]##*/}" --out "lops-$t.json"
done
# differential distributions of the production NLO+PS (Monash) with the variations for bands
if [ "$SET" != verify ] && [ -f shower-nlo-monash.json ]; then
  step dists-nlo python3 "$SA/distributions.py" --summary shower-nlo-monash.json --pdf-weight $PDFW --pdf-nmem $NMEM \
    --out dists-nlo.json
fi
# the STXS table of the production: NLO+PS with Monash, the other tunes for comparison
if [ "$SET" != verify ] && [ -f shower-nlo-monash.json ]; then
  tn=(); for t in cp5 a14; do [ -f "shower-nlo-$t.json" ] && tn+=("shower-nlo-$t.json"); done
  step stxs-nlo python3 "$SA/stxs_table.py" --summary shower-nlo-monash.json --tunes "${tn[@]}" \
    --pdf-weight $PDFW --pdf-nmem $NMEM --out stxs-nlo.json
fi
# fixed order, parton level and shower at each order, per tune (needs the parton-level runs)
if [ -f parton-nlo.json ] && [ -f parton-lo.json ] && [ -f fo-llh.json ]; then
  for f in shower-nlo-*.json; do
    t=${f#shower-nlo-}; t=${t%.json}
    case $t in bs250-*|hdampinf-*|pthard1-*|ptdef2-*) continue ;; esac
    los=(); for g in shower-lo-*-"$t".json; do [ -f "$g" ] && los+=("$g"); done
    step "effect-$t" python3 "$SA/shower_effect.py" --fo fo-llh.json --parton-nlo parton-nlo.json --parton-lo parton-lo.json \
      --nlo-ps "$f" --lo-ps "${los[@]}" --out "effect-$t.json"
  done
fi
# matching variation h_damp = infinity (13600): STXS table, parton level against fixed order, and the
# showered sample bin by bin against the nominal h_damp = 250 GeV (expected to differ: report only)
# matching variations of the POWHEG veto (13600): STXS tables and comparison with the nominal shower
for v in pthard1 ptdef2; do
  [ -f "shower-nlo-$v-monash.json" ] || continue
  step "stxs-nlo-$v" python3 "$SA/stxs_table.py" --summary "shower-nlo-$v-monash.json" \
    --pdf-weight $PDFW --pdf-nmem $NMEM --out "stxs-nlo-$v.json"
  [ -f shower-nlo-monash.json ] && step "veto-$v-shower" python3 "$SA/compare_samples.py" --a "shower-nlo-$v-monash.json" \
    --b shower-nlo-monash.json --report-only --out "veto-$v-shower.json"
done
if [ -f shower-nlo-hdampinf-monash.json ]; then
  step stxs-nlo-hdampinf python3 "$SA/stxs_table.py" --summary shower-nlo-hdampinf-monash.json \
    --pdf-weight $PDFW --pdf-nmem $NMEM --out stxs-nlo-hdampinf.json
  [ -f shower-nlo-monash.json ] && step hdamp-shower python3 "$SA/compare_samples.py" --a shower-nlo-hdampinf-monash.json \
    --b shower-nlo-monash.json --report-only --out hdamp-shower.json
  if [ -f parton-nlo-hdampinf.json ] && [ -f parton-lo.json ] && [ -f fo-llh.json ]; then
    step effect-hdampinf-monash python3 "$SA/shower_effect.py" --fo fo-llh.json --parton-nlo parton-nlo-hdampinf.json \
      --parton-lo parton-lo.json --nlo-ps shower-nlo-hdampinf-monash.json --out effect-hdampinf-monash.json
  fi
fi
if [ -f parton-nlo-bs250.json ] && [ -f parton-nlo.json ]; then
  step bornsupp-parton python3 "$SA/compare_samples.py" --a parton-nlo-bs250.json --b parton-nlo.json --out bornsupp-parton.json
  if [ -f parton-lo.json ] && [ -f fo-llh.json ] && [ -f shower-nlo-bs250-monash.json ]; then
    step effect-bs250-monash python3 "$SA/shower_effect.py" --fo fo-llh.json --parton-nlo parton-nlo-bs250.json \
      --parton-lo parton-lo.json --nlo-ps shower-nlo-bs250-monash.json --out effect-bs250-monash.json
  fi
fi
if [ -f "$OUT/shower-nlo-bs250-monash.json" ] && [ -f "$OUT/shower-nlo-monash.json" ]; then
  step bornsupp-shower python3 "$SA/compare_samples.py" --a shower-nlo-bs250-monash.json \
    --b shower-nlo-monash.json --out bornsupp-shower.json
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
