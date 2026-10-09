#!/bin/bash
# Run every job of a set of production/paper.sh one after another on this machine, with a local
# product store.  The reference implementation of the production (slow: the sets need hundreds of
# CPU hours; on a batch system run the same jobs in parallel, see README.md).
#
#   production/run_local.sh <set> <workdir> [runs ...]
#
# <workdir>/products/<product> is the store; job p of a run works in <workdir>/jobs/<set>/<run>/<p>
# with a fresh copy of the package and its inputs, and leaves its results there.  Runs of the set
# "common" (the builds) must have been run first.  A job that already finished (.done) is skipped.
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); PKG=$(cd "$HERE/.." && pwd)
SET=${1:?set}; WD=${2:?workdir}; shift 2
mkdir -p "$WD/products"; WD=$(cd "$WD" && pwd)
RUNS=("$@"); [ ${#RUNS[@]} -gt 0 ] || RUNS=($("$HERE/paper.sh" list "$SET" | awk '{print $1}'))
for run in "${RUNS[@]}"; do
  jobs=$("$HERE/paper.sh" list "$SET" | awk -v r="$run" '$1 == r {print $2}')
  [ -n "$jobs" ] || { echo "FATAL: no run $run in set $SET"; exit 2; }
  for p in $(seq 0 $((jobs - 1))); do
    J=$WD/jobs/$SET/$run/$p; [ -f "$J/.done" ] && continue
    rm -rf "$J"; mkdir -p "$J"
    ( cd "$PKG" && git archive HEAD VERSION production generator shower 2>/dev/null || tar cf - VERSION production generator shower ) | tar xf - -C "$J"
    while read -r f prod; do
      mkdir -p "$J/$(dirname "$f")"; cp -p "$WD/products/$prod" "$J/$f" || { echo "FATAL: product $prod missing"; exit 2; }
    done < <("$HERE/paper.sh" inputs "$SET" "$run" "$p")
    echo "$(date -u +%FT%TZ) $SET $run $p"
    ( cd "$J" && production/paper.sh run "$SET" "$run" "$p" > job.log 2>&1 ) || { echo "FATAL: $SET $run $p failed, see $J/job.log"; exit 3; }
    while read -r f prod; do
      mkdir -p "$WD/products/$(dirname "$prod")"; cp -p "$J/$f" "$WD/products/$prod"
    done < <("$HERE/paper.sh" outputs "$SET" "$run" "$p")
    touch "$J/.done"
  done
done
