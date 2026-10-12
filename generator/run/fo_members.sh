#!/bin/bash
# The ggxy fixed-order cross section for every member of a PDF set, one run per member with the same
# seed: an independent check of the PDF and alpha_s uncertainties obtained from event weights.
#
#   GGZH_PREFIX=<prefix> run/fo_members.sh <lo|nlo> <outdir> <seed> <np> <np_opt> <PDF set> <first member> <n members> [program] [sqrt(s)]
#
# ggxy loads member 0 of the set it is given (LHAPDF::mkPDF(name, 0)).  For member k this script
# makes a one-member copy of the set, <set>_m<k> (the set's .info with NumMembers 1 and the member's
# grid as member 0; a member's own header, e.g. its alpha_s(m_Z), is kept), in a directory put on
# LHAPDF_DATA_PATH, and runs run/fo.sh with it.  NCPU members at a time.
# Output: <outdir>/m<k>/ (fo.sh output: histo_<mode>_<seed>.txt, log, manifest) for each member.
set -o pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
MODE=${1:?lo|nlo}; OUT=${2:?outdir}; SEED=${3:?seed}; NP=${4:?np}; NPOPT=${5:?np_opt}; SET=${6:?PDF set}
FIRST=${7:?first member}; NM=${8:?members}; PROG=${9:-fo-ggllh}; ECM=${10:-13600}
load_env
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
src=""; IFS=: read -ra dirs <<< "${LHAPDF_DATA_PATH:-}"
for d in "${dirs[@]}" "$(lhapdf-config --datadir 2>/dev/null)"; do [ -n "$d" ] && [ -f "$d/$SET/$SET.info" ] && { src=$d/$SET; break; }; done
[ -n "$src" ] || die "PDF set $SET not found on LHAPDF_DATA_PATH"
PDFDIR=$OUT/lhapdf; mkdir -p "$PDFDIR"; export LHAPDF_DATA_PATH=$PDFDIR:${LHAPDF_DATA_PATH:-}
manifest "$OUT/manifest.txt" "fo_members.sh $PROG $MODE seed $SEED np $NP np_opt $NPOPT set $SET members $FIRST..+$NM sqrt(s) $ECM ($src)"
NPAR=${NCPU:-1}; running=0; rc=0
for k in $(seq "$FIRST" $((FIRST + NM - 1))); do
  m=$(printf '%s_m%02d' "$SET" "$k"); f=$(printf '%s/%s_%04d.dat' "$src" "$SET" "$k")
  [ -f "$f" ] || die "no member $k of $SET ($f)"
  mkdir -p "$PDFDIR/$m"
  sed 's/^NumMembers:.*/NumMembers: 1/' "$src/$SET.info" > "$PDFDIR/$m/$m.info"
  cp -p "$f" "$PDFDIR/$m/${m}_0000.dat"
  "$GEN/run/fo.sh" "$MODE" "$OUT/m$k" "$SEED" 1 "$NP" "$NPOPT" "$m" "$PROG" "$ECM" > "$OUT/fo-m$k.log" 2>&1 &
  running=$((running + 1))
  if [ "$running" -ge "$NPAR" ]; then wait -n || rc=1; running=$((running - 1)); fi
done
while [ "$running" -gt 0 ]; do wait -n || rc=1; running=$((running - 1)); done
rm -rf "$PDFDIR"
echo "$(stamp) fo_members.sh: $(ls "$OUT"/m*/histo_"$MODE"_*.txt 2>/dev/null | wc -l) of $NM members"
exit $rc
