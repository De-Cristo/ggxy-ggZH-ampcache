#!/bin/bash
# Clone the three pinned repositories into <dir> and pack them as <dir>.tgz, for building on
# machines without access to gitlab.com (build.sh with SRC_DIR=<dir>, or condor/run_job.sh,
# which uses sources.tgz when it is transferred).
#
#   ./fetch_sources.sh <dir>
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/versions.sh"
D=${1:?usage: fetch_sources.sh <dir>}
mkdir -p "$D"
clone() { [ -d "$1/.git" ] || git clone -q "$2" "$1"; git -C "$1" checkout -q "$3"; echo "$1 $(git -C "$1" rev-parse HEAD)"; }
clone "$D/ggxy-release" "$GGXY_URL" "$GGXY_COMMIT"
clone "$D/POWHEG-BOX-V2" "$POWHEG_URL" "$POWHEG_COMMIT"
clone "$D/POWHEG-BOX-V2/ggxy_ggZH" "$GGZH_URL" "$GGZH_COMMIT"
tar czf "$D.tgz" -C "$(dirname "$D")" "$(basename "$D")"
ls -la "$D.tgz"
