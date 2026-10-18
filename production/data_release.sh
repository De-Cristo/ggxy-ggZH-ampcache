#!/bin/bash
# The data record of the paper (Zenodo): every product of production/paper.sh except the compiled installations,
# and the analysis outputs of production/analyze.sh, packed into one archive per sample or category (Zenodo takes at
# most 100 files per record), with a checksum manifest and production/DATA-README.md.
#
#   production/data_release.sh list
#        the record: "<archive> <product>" for every product, and "<archive> analysis:<set>" for the analysis outputs
#   production/data_release.sh stage <outdir> --stage <xrootd url of the products> [--stage-for <glob>=<url> ...]
#                                  --analysis <set>=<dir> [...]
#        fetch every product (checked against its <product>.sha256 when that exists), pack the archives, add the
#        analysis outputs (each <dir> = the analysis/ directory of analyze.sh <set>), README.md and MANIFEST.sha256
#        (sha256 of every archive, and of every file inside them).  Needs xrdcp and about 4 GB in <outdir>; run on a
#        machine with access to the products, not on a shared login node.  --stage-for: products whose path
#        matches the shell glob are fetched from <url> instead (first match wins), e.g. the Z -> nu nubar products
#        of the vv sets from their own staging area: --stage-for '*-vv-*'=<url>.
#
# Archives: lhe-<sample>.tar (LHE events), grids.tar (POWHEG grids and the LO grid cross section of the CMS-style
# sample), tables-<sample>.tar (per-event tables of the reweighting study), rwt-weights.tar (per-event classifier
# weights and the classifier outputs), rwt-model.tgz (the reusable models and their card), analysis-<set>.tar.
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SETS="common 13600 13000 verify vv13600 vv13000"
archive() {  # archive <product>: the archive a product goes to (empty: not released)
  case $1 in
    install/*)              echo "" ;;
    lhe/*)                  p=${1#lhe/}; echo "lhe-${p%%/*}.tar" ;;
    grids/*)                echo "grids.tar" ;;
    tables/*)               p=${1#tables/}; echo "tables-${p%%/*}.tar" ;;
    rwt/model-*)            echo "rwt-model.tgz" ;;
    rwt/*)                  echo "rwt-weights.tar" ;;
    *)                      echo "other.tar" ;;
  esac
}
list() {
  local s r n p f prod a
  for s in $SETS; do
    while read -r r n _; do
      for p in $(seq 0 $((n - 1))); do
        while read -r f prod; do a=$(archive "$prod"); [ -n "$a" ] && echo "$a $prod"; done < <("$HERE/paper.sh" outputs "$s" "$r" "$p")
      done
    done < <("$HERE/paper.sh" list "$s")
    [ "$s" = common ] || echo "analysis-$s.tar analysis:$s"
  done
}
CMD=${1:?list|stage}; shift
case $CMD in
  list) list ;;
  stage)
    OUTD=${1:?outdir}; shift; STAGE=""; declare -A AN; FORG=(); FORU=()
    while [ $# -gt 0 ]; do
      case $1 in
        --stage) STAGE=$2; shift 2 ;;
        --stage-for) FORG+=("${2%%=*}"); FORU+=("${2#*=}"); shift 2 ;;
        --analysis) AN[${2%%=*}]=${2#*=}; shift 2 ;;
        *) echo "unknown option $1" >&2; exit 2 ;;
      esac
    done
    src() {  # src <product>: the staging url the product is fetched from
      local i; for i in "${!FORG[@]}"; do [[ $1 == ${FORG[$i]} ]] && { echo "${FORU[$i]}"; return; }; done; echo "$STAGE"
    }
    [ -n "$STAGE" ] || { echo "--stage <xrootd url of the products> needed" >&2; exit 2; }
    command -v xrdcp > /dev/null || { echo "xrdcp not found" >&2; exit 2; }
    mkdir -p "$OUTD/files" "$OUTD/record"; OUTD=$(cd "$OUTD" && pwd)
    list > "$OUTD/record.list"
    bad=0
    while read -r a prod; do
      case $prod in analysis:*) continue ;; esac
      d=$OUTD/files/$a; mkdir -p "$d/$(dirname "$prod")"
      [ -s "$d/$prod" ] && continue
      u=$(src "$prod")
      xrdcp -s -f "$u/$prod" "$d/$prod" || { echo "MISSING $prod" >&2; bad=1; continue; }
      if xrdcp -s -f "$u/$prod.sha256" "$d/$prod.sha256.tmp" 2>/dev/null; then
        [ "$(cut -d' ' -f1 "$d/$prod.sha256.tmp")" = "$(sha256sum "$d/$prod" | cut -d' ' -f1)" ] || { echo "CHECKSUM $prod" >&2; bad=1; }
        rm -f "$d/$prod.sha256.tmp"
      fi
    done < "$OUTD/record.list"
    for s in "${!AN[@]}"; do
      [ -d "${AN[$s]}" ] || { echo "no analysis directory ${AN[$s]}" >&2; bad=1; continue; }
      mkdir -p "$OUTD/files/analysis-$s.tar/analysis-$s"; cp -rp "${AN[$s]}"/. "$OUTD/files/analysis-$s.tar/analysis-$s/"
    done
    for s in $SETS; do [ "$s" = common ] || [ -n "${AN[$s]:-}" ] || { echo "no --analysis $s=<dir>" >&2; bad=1; }; done
    [ $bad = 0 ] || { echo "stopped: missing or wrong files" >&2; exit 3; }
    ( cd "$OUTD/files" && find . -type f | sort | sed 's|^\./||' | while read -r f; do sha256sum "$f"; done ) > "$OUTD/record/FILES.sha256"
    for a in $(awk '{print $1}' "$OUTD/record.list" | sort -u); do
      case $a in
        *.tgz) ( cd "$OUTD/files/$a" && tar czf "$OUTD/record/$a" . ) ;;
        *)     ( cd "$OUTD/files/$a" && tar cf "$OUTD/record/$a" . ) ;;
      esac
    done
    cp -p "$HERE/DATA-README.md" "$OUTD/record/README.md"
    ( cd "$OUTD/record" && sha256sum *.tar *.tgz README.md FILES.sha256 > MANIFEST.sha256 )
    echo "record in $OUTD/record: $(ls "$OUTD/record" | wc -l) files, $(du -sh "$OUTD/record" | cut -f1)" ;;
  *) echo "unknown command $CMD" >&2; exit 2 ;;
esac
