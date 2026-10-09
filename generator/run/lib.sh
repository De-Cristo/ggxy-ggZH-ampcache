# Shared helpers of the run scripts (sourced).

GEN=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
stamp() { date -u +%FT%TZ; }
die() { echo "FATAL: $*" >&2; exit 1; }

# The installation made by build.sh: GGZH_PREFIX=<prefix>, or env.sh already sourced.
load_env() {
  [ -n "${GGZH_PREFIX:-}" ] || die "set GGZH_PREFIX to the build.sh prefix"
  [ -f "$GGZH_PREFIX/env.sh" ] || die "no $GGZH_PREFIX/env.sh; run build.sh first"
  source "$GGZH_PREFIX/env.sh"
  command -v pwhg_main >/dev/null || die "pwhg_main not in PATH"
  [ -n "${NCPU:-}" ] || NCPU=$(nproc)
}

# make_card <base card> <output> key=value ...   (removes every line starting with a key, appends the new values)
make_card() {
  local base=$1 out=$2; shift 2
  local expr=() kv
  for kv in "$@"; do expr+=(-e "/^[[:space:]]*${kv%%=*}[[:space:]]/d"); done
  sed "${expr[@]}" "$base" > "$out"
  for kv in "$@"; do echo "${kv%%=*} ${kv#*=}" >> "$out"; done
}

# run_parallel <n> <log prefix>: POWHEG with manyseeds, seed index 1..n on this node (n <= NCPU is sensible)
run_parallel() {
  local n=$1 log=$2 i rc=0 pids=()
  for i in $(seq 1 "$n"); do
    ( echo "$i" | pwhg_main > "$log-$i.log" 2>&1 ) &
    pids+=($!)
  done
  for i in "${pids[@]}"; do wait "$i" || rc=1; done
  return $rc
}

# commit of this package: VERSION as written by git archive (export-subst), else the git checkout
package_version() {
  local v; v=$(head -1 "$GEN/../VERSION" 2>/dev/null)
  case $v in ''|*Format*) v="$(git -C "$GEN" rev-parse HEAD 2>/dev/null || echo unknown)$(git -C "$GEN" status --porcelain 2>/dev/null | grep -q . && echo ' modified')";; esac
  echo "$v"
}

# record what ran where
manifest() {  # manifest <file> <text...>
  { echo "$(stamp) $*"; echo "host $(hostname)"; echo "package $(package_version)"
    head -6 "$GGZH_PREFIX/BUILD-MANIFEST.txt" 2>/dev/null; } >> "$1"
}
