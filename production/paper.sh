#!/bin/bash
# The production of the gg -> ZH NLO+PS paper: every run, its jobs, seeds, statistics and
# settings.  This file is the only definition of them.  The results are made by executing these
# jobs at the package commit recorded in VERSION, on any HTCondor pool or one after another
# (production/run_local.sh); the outputs of the jobs are independent of where they run, except
# that NLO integrations agree only statistically across compilers and CPUs (generator/README.md).
#
#   production/paper.sh list    <set>               runs of the set in execution order:
#                                                   <run> <jobs> <cpus> <memory GB> <hours>
#   production/paper.sh inputs  <set> <run> <p>     files job p needs:    <file in the job directory> <product>
#   production/paper.sh outputs <set> <run> <p>     products job p makes: <file in the job directory> <product>
#   production/paper.sh show    <set> <run> <p>     the environment and command of job p
#   production/paper.sh run     <set> <run> <p>     execute job p in the current directory
#
# Sets
#   common   the generator and shower builds (used by every other set)
#   13600    PDF4LHC21, sqrt(s) = 13.6 TeV: fixed order, POWHEG LO, grids, events, cache closure,
#            MCFM, showers (NLO+PS and five LO+PS shower starts, tunes Monash, CP5, A14), and the
#            same NLO and LO events at parton level (no shower: the shower effect at each order and
#            the high-pT test of the Born suppression against the fixed-order NLO)
#            and, at 13.6 TeV only, the variations (PDF members, cache scan, h_damp, POWHEG veto, damped
#            LO+PS start, MS-bar top mass) and the reweighting study (rwt/README.md: CMS-style LO+PS with
#            public CMS settings, per-event tables, classifier weights, summary)
#   13000    the same runs, seeds and statistics at sqrt(s) = 13 TeV
#   verify   NNPDF3.1, 13.6 TeV: the verification runs of generator/reference (48 checks), the
#            unweighted validation samples and their showers, the Born-suppression test (also at
#            parton level: unweighted and weighted NLO events against the fixed-order NLO)
#   vv13600  Z -> nu nubar (vdecaymode 12, three flavours, 15 < m_nunu < 150 GeV), otherwise the 13600
#            settings: fixed order (fo-ggvvh), LO/NLO events and showers, matching and reweighting;
#            seeds independent of the l+l- runs (fixed order 12001/13001, NLO events 30001)
#   vv13000  Z -> nu nubar at 13 TeV: fixed order, LO/NLO events and showers
#   atlas13600  the reweighting study with an ATLAS-style LO+PS source (rwt/README.md): LO events with the
#            ATLAS scale H_T = m_T(H) + p_T(l1) + p_T(l2) (generator/patches/ggxy_ggZH-scale-ht-atlas.patch,
#            card dynscale 3) and PDF4LHC21, showered with A14 as athena sets it, the ATLAS base settings and the
#            ATLAS POWHEG veto; target: the 13600 production NLO events showered with the same A14 settings;
#            POWHEG LO at the seven scale points with the same card (check of the scale option); seeds 70001
# Jobs p = 0 .. <jobs>-1; seeds are fixed functions of p (below), identical in 13600 and 13000.
#
# A job directory contains the unpacked package (VERSION, production/, generator/, shower/) and
# the files listed by "inputs".  Products have storage-independent names: install/<name>.tgz,
# grids/<name>.tgz, lhe/<sample>/events-<seed>.lhe.gz.  Everything else a job writes (out-*.tgz,
# *.json, logs) is its result and stays in the job directory.
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PKG=$(cd "$HERE/.." && pwd)
CMD_=${1:?list|inputs|outputs|show|run}; SET=${2:?set: common|13600|13000|verify|vv13600|vv13000|atlas13600}

TUNES=(monash cp5 a14)
LOMODES=(default scalup scalup-half scalup-quarter power)
MTFAC=(1 2 4)                       # MS-bar top-mass scale mu_t = m_ZH / f

S0NLO=10001                        # first seed of the production NLO events
S0LO=20001; S0HD=10001; S0CMS=40001  # first seeds of the LO, h_damp = inf and CMS-style LO events
VDM=""                             # Z decay setting added to the grids of the CMS-style LO (vv sets: vdecaymode=12)
RWT_AVAILABLE_MANIFEST=""           # only the user-authorized available-table analysis sets
RWT_CACHED_INPUTS=""                # parsed tables are staged once through EOS
RWT_RETUNED=""                      # only the registered conditional R6 capacity test
case $SET in
  common) ;;
  13600|13000)
    E=$SET; EB=$([ "$E" = 13600 ] && echo 6800d0 || echo 6500d0)
    PDF=PDF4LHC21_40_pdfas
    CARD="lhans1=93300 lhans2=93300 ebeam1=$EB ebeam2=$EB"
    RWL=pwg-rwl-pdf4lhc21.dat ;;
  vv13600|vv13000)                  # Z -> nu nubar: the card's vdecaymode 11 replaced by 12
    E=${SET#vv}; EB=$([ "$E" = 13600 ] && echo 6800d0 || echo 6500d0)
    PDF=PDF4LHC21_40_pdfas
    CARD="lhans1=93300 lhans2=93300 ebeam1=$EB ebeam2=$EB vdecaymode=12"
    RWL=pwg-rwl-pdf4lhc21.dat
    S0NLO=30001; S0LO=25001; S0HD=35001; S0CMS=45001; VDM=vdecaymode=12 ;;   # seeds independent of the l+l- samples
  atlas13600)                       # ATLAS-style LO: dynamic scale H_T (dynscale 3), otherwise the 13600 settings
    E=13600; EB=6800d0
    PDF=PDF4LHC21_40_pdfas
    CARD="lhans1=93300 lhans2=93300 ebeam1=$EB ebeam2=$EB dynscale=3"
    RWL=pwg-rwl-pdf4lhc21.dat
    S0LO=70001 ;;                   # seeds independent of every other LO sample
  rwt10m|rwt10m-pilot)
    E=13600; EB=6800d0; PDF=PDF4LHC21_40_pdfas
    CARD="lhans1=93300 lhans2=93300 ebeam1=$EB ebeam2=$EB"
    RWL=pwg-rwl-rwt.dat; S0NLO=100001; S0CMS=110001
    NNLO=4000; NCMS=2000; NREP=1; PILOT=""
    if [ "$SET" = rwt10m-pilot ]; then NNLO=2; NCMS=2; NREP=2; PILOT=1; fi ;;
  rwt10m-available|rwt10m-available-pilot)
    E=13600; EB=6800d0; PDF=PDF4LHC21_40_pdfas
    CARD="lhans1=93300 lhans2=93300 ebeam1=$EB ebeam2=$EB"
    RWL=pwg-rwl-rwt.dat; S0NLO=100001; S0CMS=110001
    NNLO=4000; NCMS=2000; NREP=1; PILOT=""
    RWT_AVAILABLE_MANIFEST=production/rwt10m-available.json
    if [ "$SET" = rwt10m-available-pilot ]; then
      NNLO=2; NCMS=2; NREP=2; PILOT=1
      RWT_AVAILABLE_MANIFEST=production/rwt10m-available-pilot.json
    fi ;;
  rwt10m-available-cached|rwt10m-available-cached-pilot|rwt10m-available-cached-r6|rwt10m-available-cached-r6-pilot)
    E=13600; EB=6800d0; PDF=PDF4LHC21_40_pdfas
    CARD="lhans1=93300 lhans2=93300 ebeam1=$EB ebeam2=$EB"
    RWL=pwg-rwl-rwt.dat; S0NLO=100001; S0CMS=110001
    NNLO=4000; NCMS=2000; NREP=1; PILOT=""; RWT_CACHED_INPUTS=1
    RWT_AVAILABLE_MANIFEST=production/rwt10m-available.json
    case $SET in *-r6|*-r6-pilot) RWT_RETUNED=1 ;; esac
    if [[ "$SET" = *-pilot ]]; then
      NNLO=2; NCMS=2; NREP=2; PILOT=1
      RWT_AVAILABLE_MANIFEST=production/rwt10m-available-pilot.json
    fi ;;
  rwt10m-dnn-c|rwt10m-dnn-c-pilot)
    E=13600; NREP=1; PILOT=""; RWT_CACHED_INPUTS=1
    RWT_AVAILABLE_MANIFEST=production/rwt10m-available.json
    if [[ "$SET" = *-pilot ]]; then
      NREP=2; PILOT=1; RWT_AVAILABLE_MANIFEST=production/rwt10m-available-pilot.json
    fi ;;
  vv-rwt10m)
    E=13600; EB=6800d0; PDF=PDF4LHC21_40_pdfas
    CARD="lhans1=93300 lhans2=93300 ebeam1=$EB ebeam2=$EB vdecaymode=12"
    RWL=pwg-rwl-rwt.dat; S0NLO=300001; S0CMS=310001; VDM=vdecaymode=12
    NNLO=4000; NCMS=2000; NREP=1; PILOT="" ;;
  parton-zll-jets-1m-v1|parton-zvv-jets-1m-v1|parton-zll-1m-v1|parton-zvv-1m-v1|parton-zll-1m-recovery-v2) E=13600 ;;
  numerical-zll-1m-v1|numerical-zvv-1m-v1|stable-zll-1m-v1|stable-zvv-1m-v1|parton-zvv-stable-1m-v1|parton-zvv-jets-stable-1m-v1) E=13600 ;;
  scale-cms-zll-eval-v1|scale-cms-zvv-eval-v1|nlo-topup-zvv-1m-v1) E=13600 ;;
  nlo-topup-zll-v1) E=13600 ;;
  nlo-topup-zll-10m-v1|nlo-topup-zll-10m-v2) E=13600 ;;
  uncertainty-zll-1m-v1|uncertainty-zvv-1m-v1) E=13600 ;;
  ensemble-zll-v1|ensemble-zvv-v1|ensemble-zll-v2|ensemble-zvv-v2|ensemble-atlas-zll-v2|ensemble-atlas-zvv-v2)
    # v2: new evaluation/member selections (<channel>-v2-*.json), training domain pT(j1) < m_ZH (RWT_DOMAIN=rho1)
    # atlas-*: scalar tables only (scalar-prepare), own calibration selection; ECHAN = decay channel zll/zvv
    E=13600; ECH=${SET#ensemble-}; EVER=${ECH##*-}; ECH=${ECH%-v*}; ESEL=$ECH; [ "$EVER" = v1 ] || ESEL=$ECH-$EVER
    ECHAN=${ECH#atlas-}; ECALIB=rwt/parton-$ECH-1m-v1/prepared.tgz
    [[ $ECH = atlas-* ]] && ECALIB=rwt/$SET/prepared/calib.tgz ;;
  tables-atlas-zll-v1|tables-atlas-zvv-v1) E=13600 ;;
  bdt17-tables-zll|bdt17-tables-zvv)
    # Table-only copies of the CMS-style BDT-17 inputs for the data record (the prepared bundles also hold particle tensors).
    E=13600; ECH=${SET#bdt17-tables-} ;;
  scale-atlas-control-v2|scale-atlas-zll-calib-v2|scale-atlas-zvv-calib-v2) E=13600 ;;
  particle-atlas-zll)
    E=13600; EB=6800d0; PDF=PDF4LHC21_40_pdfas
    CARD="lhans1=93300 lhans2=93300 ebeam1=$EB ebeam2=$EB dynscale=3"
    NNLO=4000; NATLAS=2000; S0NLO=100001; S0ATLAS=500001 ;;
  particle-zll|particle-zvv)
    E=13600; NNLO=4000; NCMS=2000
    if [ "$SET" = particle-zll ]; then S0NLO=100001; S0CMS=110001; else S0NLO=300001; S0CMS=310001; fi ;;
  rwt10m-transformer-c)
    E=13600; NREP=1; PILOT=""; RWT_CACHED_INPUTS=1
    RWT_AVAILABLE_MANIFEST=production/rwt10m-available.json ;;
  verify)
    E=13600; PDF=NNPDF31_nnlo_as_0118_mc_hessian_pdfas
    CARD=""                       # the card's own settings: NNPDF3.1, 13.6 TeV
    RWL=pwg-rwl.dat ;;
  *) echo "unknown set $SET" >&2; exit 2 ;;
esac

# ---------------------------------------------------------------- runs of each set
runs() {
  case $SET in
    common) echo build shower-build ;;
    parton-zll-jets-1m-v1|parton-zvv-jets-1m-v1) echo parton-fit-1m parton-jet-compare-1m ;;
    parton-zvv-1m-v1) echo parton-prepare-1m parton-fit-1m parton-bdt-1m parton-early-report-1m parton-report-1m ;;
    parton-zll-1m-recovery-v2) echo parton-fit-1m parton-early-report-1m parton-report-1m ;;
    parton-zll-1m-v1) echo parton-prepare-1m parton-fit-1m parton-bdt-1m parton-early-report-1m parton-report-1m parton-reload-diagnose-1m parton-fold0-report-1m ;;
    numerical-zll-1m-v1|numerical-zvv-1m-v1) echo numerical-input-audit ;;
    stable-zll-1m-v1|stable-zvv-1m-v1) echo stable-input-audit ;;
    parton-zvv-stable-1m-v1) echo parton-fit-1m parton-early-report-1m parton-report-1m ;;
    parton-zvv-jets-stable-1m-v1) echo parton-fit-1m parton-jet-compare-1m ;;
    scale-cms-zll-eval-v1|scale-cms-zvv-eval-v1) echo scale-selection ;;
    scale-atlas-control-v2) echo atlas-scale-validate ;;
    scale-atlas-zll-calib-v2|scale-atlas-zvv-calib-v2) echo atlas-scale-calibration ;;
    nlo-topup-zvv-1m-v1) echo nlo-weight-topup-all ;;
    nlo-topup-zll-v1) echo nlo-weight-topup ;;
    nlo-topup-zll-10m-v1|nlo-topup-zll-10m-v2) echo nlo-weight-topup-all ;;
    uncertainty-zll-1m-v1|uncertainty-zvv-1m-v1) echo scale-inputs-1m ;;
    ensemble-zll-v1|ensemble-zvv-v1) echo ensemble-prepare ensemble-smoke ensemble-fit ensemble-assess-smoke ensemble-assess bdt-fit bdt-assess pdf-transport bdt-deliver bdt-assess-variants bdt-assess-variants2 ;;
    ensemble-zll-v2|ensemble-zvv-v2) echo ensemble-prepare bdt-fit bdt-assess bdt-deliver bdt-assess-variants2 bdt-assess-tail bdt-cache bdt-scaleband bdt-tables2 bdt-fit3 bdt-assess3 bdt-confirm bdt-deliver4 bdt-cache4 bdt-deliver5 ;;
    ensemble-atlas-zll-v2|ensemble-atlas-zvv-v2) echo scalar-prepare bdt-fit bdt-assess bdt-deliver bdt-cache bdt-scaleband bdt-tables2 ;;
    bdt17-tables-zll|bdt17-tables-zvv) echo table-only ;;
    tables-atlas-zll-v1) echo scalar-lo-atlas scalar-nlo-a14 parity-lo-atlas parity-nlo-a14 atlas-born-scale-control ;;
    tables-atlas-zvv-v1) echo grids-lo-atlas-vv events-lo-atlas-vv scalar-lo-atlas scalar-nlo-a14 ;;
    particle-atlas-zll) echo events-lo-atlas-10m particle-nlo-a14-10m particle-lo-atlas-10m ;;
    particle-zll|particle-zvv) echo particle-build particle-python particle-splits particle-nlo particle-lo-cms ;;
    vv-rwt10m) echo events-nlo-10m events-lo-cms-10m ;;
    rwt10m-transformer-c) echo rwt-transformer-fold-10m rwt-transformer-merge-10m ;;
    rwt10m-dnn-c|rwt10m-dnn-c-pilot) echo rwt-dnn-c-10m ;;
    rwt10m|rwt10m-pilot)
      echo shower-build-rwt10m events-nlo-10m events-lo-cms-10m table-nlo-cp5-10m table-lo-cms-10m \
           rwt-clf-a-10m rwt-clf-b-10m rwt-clf-c-10m rwt-summary-10m rwt-transfer-10m rwt-model-10m ;;
    rwt10m-available|rwt10m-available-pilot)
      echo rwt-clf-a-10m rwt-clf-b-10m rwt-clf-c-10m rwt-summary-10m rwt-transfer-10m rwt-model-10m ;;
    rwt10m-available-cached-r6|rwt10m-available-cached-r6-pilot)
      echo rwt-clf-a-10m rwt-clf-b-10m rwt-clf-c-10m rwt-summary-10m rwt-transfer-10m rwt-model-10m ;;
    rwt10m-available-cached|rwt10m-available-cached-pilot)
      echo rwt-prepare-10m rwt-clf-a-10m
      [ -z "$PILOT" ] || echo rwt-cache-check-10m
      echo rwt-clf-b-10m rwt-clf-c-10m rwt-summary-10m rwt-transfer-10m rwt-model-10m ;;
    13600|13000)
      echo fo-llh-lo fo-llh-nlo fo-zh-lo fo-zh-nlo powheg-lo grids-lo grids-nlo mcfm \
           events-nlo events-lo cache
      for t in "${TUNES[@]}"; do echo "shower-nlo-$t"; done
      for t in "${TUNES[@]}"; do for m in "${LOMODES[@]}"; do echo "shower-lo-$m-$t"; done; done
      echo parton-nlo parton-lo
      for f in "${MTFAC[@]}"; do echo "fo-llh-lo-mt$f fo-llh-nlo-mt$f"; done
      # (if, not "[ ] &&": the list must end with status 0 for known() under pipefail)
      if [ "$SET" = 13600 ]; then
        for f in "${MTFAC[@]}"; do echo "fo-zh-lo-mt$f fo-zh-nlo-mt$f"; done
        echo cache-scan fo-llh-pdf fo-zh-pdf \
           grids-nlo-hdampinf events-nlo-hdampinf parton-nlo-hdampinf shower-nlo-hdampinf-monash \
           shower-nlo-pthard1-monash shower-nlo-ptdef2-monash shower-lo-damped-monash
        echo grids-lo-cms events-lo-cms table-lo-cms table-lo-cms-default table-nlo-cp5 \
             rwt-clf rwt-clf-default rwt-clf-jet rwt-clf-jet-default rwt-summary rwt-model
      fi ;;
    verify)
      echo fo-llh-lo fo-llh-nlo fo-zh-lo fo-zh-lo-pdf4lhc21 powheg-lo grids-nlo grids-lo \
           grids-nlo-bs250 mcfm events-nlo events-lo events-nlo-bs250 cache
      for t in "${TUNES[@]}"; do echo "shower-nlo-$t"; done
      for m in "${LOMODES[@]}"; do echo "shower-lo-$m-monash"; done
      echo shower-lo-default-cp5 shower-lo-scalup-half-cp5 shower-lo-default-a14 shower-lo-scalup-half-a14
      echo shower-nlo-bs250-monash
      echo parton-nlo parton-lo parton-nlo-bs250 ;;
    vv13600|vv13000)                # Z -> nu nubar: the l+l- runs of the paper except the m_t scheme (one check at 13.6 TeV),
                                    # the stable-Z, MCFM, POWHEG-LO, cache and PDF-member runs (decay independent or code checks)
      [ "$SET" = vv13600 ] && echo build      # ggzh-install-v3 (v2 plus bin/fo-ggvvh), used by both vv sets
      echo fo-vvh-lo fo-vvh-nlo grids-nlo events-nlo grids-lo events-lo parton-nlo parton-lo
      for t in "${TUNES[@]}"; do echo "shower-nlo-$t"; done
      for t in "${TUNES[@]}"; do for m in "${LOMODES[@]}"; do echo "shower-lo-$m-$t"; done; done
      if [ "$SET" = vv13600 ]; then
        echo grids-nlo-hdampinf events-nlo-hdampinf parton-nlo-hdampinf shower-nlo-hdampinf-monash \
             shower-nlo-pthard1-monash shower-nlo-ptdef2-monash shower-lo-damped-monash \
             grids-lo-cms events-lo-cms table-lo-cms table-lo-cms-default table-nlo-cp5 \
             rwt-clf rwt-clf-default rwt-clf-jet rwt-clf-jet-default rwt-summary rwt-model fo-vvh-lo-mt1 fo-vvh-nlo-mt1
      fi ;;
    atlas13600)                     # ggzh-install-v4 (v3 plus the H_T scale option); LO with the card's dynscale 3
      echo build powheg-lo grids-lo events-lo table-lo-atlas table-nlo-a14 \
           rwt-clf-atlas rwt-clf-jet-atlas rwt-summary-atlas rwt-model-atlas ;;
  esac
}

# ---------------------------------------------------------------- one job
# def <run> <p>: sets JOBS CPUS MEM HOURS, ENVS (VAR=value ...), STEP (command), IN and OUTP
# ("<file> <product>" lines).  G = generator/, S = shower/ of the job directory.
G=generator; S=shower
# ggzh-install-v2: the fixed-order programs take the top-mass scheme and mu_t (generator/fo, 2026-10-13); the
# production runs of 2026-10-09..12 used install/ggzh-install, identical apart from these two options
# shower-install-v2: the Rivet plugin also has GGZH_EVENTTABLE (the per-event tables of the reweighting study,
# 2026-10-15); the shower runs of 2026-10-09..13 used install/shower-install, otherwise identical
# ggzh-install-v3: v2 plus bin/fo-ggvvh (Z -> nu nubar), for the vv sets; built by run build of vv13600
# ggzh-install-v4: v3 plus the scale option dynscale 3 (H_T, generator/patches), for set atlas13600; built by its run build
INSTALL=install/ggzh-install-v2; SHOWER_INSTALL=install/shower-install-v2
case $SET in rwt10m*) SHOWER_INSTALL=install/shower-install-rwt10m-0 ;; esac
case $SET in vv*) INSTALL=install/ggzh-install-v3 ;; atlas*) INSTALL=install/ggzh-install-v4 ;; esac
need_install() { IN+=("ggzh-install.tgz $INSTALL.tgz"); }
need_grids()   { IN+=("grids.tgz grids/$1.tgz"); }
# shower job p of a sample: <nper> LHE files of seeds first + nper*p ...
need_lhe() {  # need_lhe <sample> <first seed> <files per job> <p>
  local k s
  for k in $(seq 0 $(($3 - 1))); do s=$(($2 + $3 * $4 + k)); IN+=("lhe/events-$s.lhe.gz lhe/$1/events-$s.lhe.gz"); done
  IN+=("shower-install.tgz $SHOWER_INSTALL.tgz")
}
# reweighting study: the per-event tables of <n> consecutive seeds of a sample (production run table-<sample>)
need_tables() {  # need_tables <sample> <first seed> <n>
  local k sample; sample=$(sname "$1")
  for k in $(seq 0 $(($3 - 1))); do IN+=("tables/$1/table-$(($2 + k)).csv.gz tables/$sample/table-$(($2 + k)).csv.gz"); done
}
# the same with an explicit product sample name (set atlas13600: tables of another set's samples, tables/<product>/)
need_tables_as() {  # need_tables_as <sample in the job> <product sample> <first seed> <n>
  local k; for k in $(seq 0 $(($4 - 1))); do IN+=("tables/$1/table-$(($3 + k)).csv.gz tables/$2/table-$(($3 + k)).csv.gz"); done
}
res() { JOBS=$1; CPUS=$2; MEM=$3; HOURS=$4; }

# sample names: production lhe/<nlo|lo>-<sqrt(s)>, verification lhe/<...>-nnpdf31-13600, Z -> nu nubar <...>-vv-<sqrt(s)>,
# ATLAS-style <...>-atlas-<sqrt(s)>
sname() { case $SET in vv-rwt10m) echo "$1-vv-rwt10m-$E" ;; rwt10m*) echo "$1-rwt10m-$E" ;; verify) echo "$1-nnpdf31-13600" ;; vv*) echo "$1-vv-$E" ;; atlas*) echo "$1-atlas-$E" ;; *) echo "$1-$E" ;; esac; }

campaign_product() {
  local recovery=""
  # B computed but its first upload failed; preserve every original/partial product.
  if [ "$SET" = rwt10m-available-cached-r6 ]; then
    case $1 in weights-clf-b|bdt-clf-b|ensemble-clf-b) recovery=-recovery1 ;; esac
  fi
  echo "rwt/$1-rwt10m${RWT_AVAILABLE_MANIFEST:+-available}${RWT_CACHED_INPUTS:+-cached}${RWT_RETUNED:+-r6}${recovery}${PILOT:+-pilot-$2}.$3"
}
campaign_data_product() { echo "rwt/$1-rwt10m-available-cached${PILOT:+-pilot-$2}.$3"; }
campaign_raw_tables() {
  if [ -n "$RWT_AVAILABLE_MANIFEST" ]; then
    # Resource listing is independent of the selection file and also works from git-show stdin.
    [ "$CMD_" != list ] || return 0
    local selection f prod
    selection=$(python3 "$PKG/rwt/ggzh_rwt_manifest.py" inputs "$PKG/$RWT_AVAILABLE_MANIFEST") || exit 2
    while read -r f prod; do IN+=("$f $prod"); done <<< "$selection"
  else
    need_tables lo-cms "$S0CMS" "$NCMS"; need_tables nlo-cp5 "$S0NLO" "$NNLO"
  fi
}
campaign_tables() {
  if [ -n "$RWT_CACHED_INPUTS" ]; then
    local sample
    for sample in lo-cms nlo-cp5; do
      IN+=("rwt-out/$sample.npz $(campaign_data_product sample-$sample "$p" npz)"
           "rwt-out/$sample.json $(campaign_data_product sample-$sample "$p" json)")
    done
  else
    campaign_raw_tables
  fi
}
campaign_checks() {
  local a
  for a in a b c; do
    IN+=("rwt-in/weights-clf-$a.npz $(campaign_product weights-clf-$a "$1" npz)"
         "rwt-in/bdt-clf-$a.json $(campaign_product bdt-clf-$a "$1" json)")
  done
}

def() {
  local run=$1 p=$2 m t
  ENVS=(); STEP=(); IN=(); OUTP=()
  case $run in
    scale-selection)
                  res 1 2 8 12
                  channel=${SET#scale-cms-}; channel=${channel%-eval-v1}
                  IN=("reference.lhe.gz lhe/lo-13600/events-20001.lhe.gz")
                  STEP=(bash production/run_scale_selection.sh "production/scale-selections/$channel-eval.json")
                  OUTP=("scale-inputs.tgz rwt/$SET/scale-inputs.tgz" "scale-audit.json rwt/$SET/audit.json") ;;
    nlo-weight-topup-all)
                  res 4000 2 8 2
                  if [[ $SET = nlo-topup-zvv-1m-v1 ]]; then
                    s=$((300001+p))
                    IN=("ggzh-install.tgz @vv/install/ggzh-install-v3.tgz"
                        "grids.tgz @vv/grids/grids-nlo-vv-13600.tgz"
                        "original.lhe.gz @vv/lhe/nlo-vv-rwt10m-13600/events-$s.lhe.gz")
                    ENVS=("NLO_TOPUP_SEED=$s" "NLO_TOPUP_CHANNEL=zvv" "NLO_TOPUP_CONTROL_GATE=user-20261111-rtol1e-3")
                  else
                    need_install; need_grids grids-nlo-13600; s=$((100001+p))
                    IN+=("original.lhe.gz lhe/nlo-rwt10m-13600/events-$s.lhe.gz")
                    ENVS=("NLO_TOPUP_SEED=$s")
                  fi
                  if [[ $SET = nlo-topup-zll-10m-v2 ]]; then
                    ENVS+=("NLO_TOPUP_CONTROL_GATE=user-20261111-rtol1e-3")
                  fi
                  if [ "$s" = 101654 ]; then
                    IN+=("reused-topup.lhe.gz rwt/nlo-topup-zll-v1/events-101654-topup.lhe.gz"
                         "reused-audit.json rwt/nlo-topup-zll-v1/audit-101654.json")
                    ENVS+=("NLO_TOPUP_REUSE=1")
                  fi
                  STEP=(bash production/run_nlo_topup_all.sh)
                  OUTP=("topup.lhe.gz rwt/$SET/events-$s-topup.lhe.gz"
                         "topup-audit.json rwt/$SET/audit-$s.json"
                         "topup-weights.csv.gz rwt/$SET/weights-$s.csv.gz") ;;
    nlo-weight-topup)
                  res 1 2 8 2; need_install
                  need_grids grids-nlo-13600
                  IN+=("original.lhe.gz lhe/nlo-rwt10m-13600/events-101654.lhe.gz")
                  STEP=(bash production/run_nlo_topup.sh)
                  OUTP=("topup.lhe.gz rwt/$SET/events-101654-topup.lhe.gz"
                         "topup-audit.json rwt/$SET/audit-101654.json") ;;
    stable-input-audit)
                  res 1 1 8 8
                  channel=${SET#stable-}; channel=${channel%-1m-v1}
                  IN=("parton-prepared.tgz rwt/parton-$channel-1m-v1/prepared.tgz")
                  if [[ $channel = zvv ]]; then
                    IN+=("native.parquet rwt/numerical-recovery-inputs-v1/311108.parquet"
                         "native-inputs.npz particle-v1/zvv/lo-cms/inputs/partons-v1/311108.npz"
                         "native-objects.parquet particle-v1/zvv/lo-cms/objects/311108.parquet")
                  fi
                  STEP=(bash production/run_stable_audit.sh "$channel")
                  OUTP=("stable-audit.json rwt/$SET/audit.json") ;;
    numerical-input-audit)
                  res 1 1 8 8
                  channel=${SET#numerical-}; channel=${channel%-1m-v1}
                  IN=("parton-prepared.tgz rwt/parton-$channel-1m-v1/prepared.tgz")
                  STEP=(bash production/run_numerical_audit.sh "$channel")
                  OUTP=("numerical-audit.json rwt/$SET/audit.json"
                         "numerical-examples.json rwt/$SET/examples.json") ;;
    scale-inputs-1m)
                  res 1 1 8 12
                  channel=${SET#uncertainty-}; channel=${channel%-1m-v1}
                  IN=("parton-prepared.tgz rwt/parton-$channel-1m-v1/prepared.tgz"
                      "reference.lhe.gz lhe/lo-13600/events-20001.lhe.gz")
                  STEP=(bash production/run_scale_inputs.sh "$channel")
                  OUTP=("scale-inputs.tgz rwt/$SET/scale-inputs.tgz"
                         "scale-audit.json rwt/$SET/audit.json") ;;
    rwt-transformer-fold-10m)
                  res 15 8 18 12; campaign_tables
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1")
                  STEP=(bash production/run_transformer.sh fold "$p")
                  OUTP=("transformer-part.tgz rwt/transformer-c/fold-$p.tgz"
                        "transformer-part.json rwt/transformer-c/fold-$p.json") ;;
    rwt-transformer-merge-10m)
                  res 1 8 18 6; campaign_tables
                  for x in $(seq 0 14); do IN+=("parts/fold-$x.tgz rwt/transformer-c/fold-$x.tgz"); done
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1")
                  STEP=(bash production/run_transformer.sh merge 0)
                  OUTP=("transformer-out/result.json rwt/transformer-c/result.json"
                        "transformer-out/weights.npz rwt/transformer-c/weights.npz"
                        "transformer-models.tgz rwt/transformer-c/models.tgz") ;;
    rwt-dnn-c-10m)
                  res "$NREP" 8 18 12; [ -z "$PILOT" ] || res "$NREP" 4 8 2
                  campaign_tables
                  IN+=("vh-code.tgz sources/vh-dnn-14e64cafe1af51a4f6261703cdf0790fd24aa698.tgz")
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1" "RWT_GPUS=1"
                        "VH_TRAINER_COMMIT=14e64cafe1af51a4f6261703cdf0790fd24aa698")
                  STEP=(bash production/run_vh_dnn.sh)
                  OUTP=("dnn-out/dnn-clf-c.json rwt/dnn-clf-c-$SET-$p.json"
                        "dnn-out/weights-clf-c.npz rwt/dnn-weights-clf-c-$SET-$p.npz"
                        "dnn-models.tgz rwt/dnn-models-clf-c-$SET-$p.tgz"
                        "dnn-out/verification.json rwt/dnn-verification-clf-c-$SET-$p.json") ;;
    rwt-prepare-10m)
                  res "$NREP" 1 16 8; campaign_raw_tables
                  ENVS=("RWT_THREADS=1" "PYTHONDONTWRITEBYTECODE=1")
                  STEP=(rwt/run.sh prepare-10m)
                  OUTP=()
                  for t in lo-cms nlo-cp5; do
                    OUTP+=("rwt-out/$t.npz $(campaign_product sample-$t "$p" npz)"
                           "rwt-out/$t.json $(campaign_product sample-$t "$p" json)")
                  done
                  OUTP+=("rwt-out/rwt10m.json $(campaign_product prepared "$p" json)") ;;
    rwt-cache-check-10m)
                  res "$NREP" 1 2 2
                  IN+=("rwt-in/raw-weights.npz rwt/weights-clf-a-rwt10m-available-pilot-$p.npz"
                       "rwt-in/raw-check.json rwt/bdt-clf-a-rwt10m-available-pilot-$p.json"
                       "rwt-in/cached-weights.npz $(campaign_product weights-clf-a "$p" npz)"
                       "rwt-in/cached-check.json $(campaign_product bdt-clf-a "$p" json)")
                  ENVS=("RWT_THREADS=1" "PYTHONDONTWRITEBYTECODE=1")
                  STEP=(rwt/run.sh cache-check-10m)
                  OUTP=("rwt-out/rwt10m.json $(campaign_product cache-check "$p" json)") ;;
    particle-python)
                  res 1 1 2 1
                  STEP=(bash production/run_particle_python.sh)
                  OUTP=("pyarrow-20.0.0-cp312-cp312-manylinux_2_28_x86_64.whl particle-v1/install/pyarrow-20.0.0-cp312-cp312-manylinux_2_28_x86_64.whl"
                         "python-codec.json particle-v1/install/python-codec.json") ;;
    particle-splits)
                  res 1 1 4 2
                  IN=("split-input/lo-cms.npz rwt/sample-lo-cms-rwt10m-available-cached.npz"
                      "split-input/nlo-cp5.npz rwt/sample-nlo-cp5-rwt10m-available-cached.npz")
                  STEP=(bash -c 'source shower/versions.sh && source "$SHOWER_LCG_VIEW/setup.sh" >/dev/null 2>&1 && PYTHONDONTWRITEBYTECODE=1 python3 shower/particle/make_splits.py')
                  OUTP=("split-map.json particle-v1/splits/zll.json") ;;
    particle-build)
                  res 1 1 4 2
                  IN=("shower-install.tgz install/shower-install-rwt10m-0.tgz")
                  STEP=(bash production/run_particle_build.sh)
                  OUTP=("particle-install.tgz particle-v1/install/particle-install.tgz") ;;
    parton-prepare-1m)
                  res 1 4 24 8
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1")
                  # Bulk inputs are the frozen storage-independent artifact/SHA catalog in the package.
                  STEP=(bash production/run_parton_study.sh prepare)
                  OUTP=("parton-prepared.tgz rwt/parton-zll-1m-v1/prepared.tgz"
                         "parton-prepare.json rwt/parton-zll-1m-v1/prepare-audit.json") ;;
    scalar-prepare)
                  # Scalar tables only (ATLAS-style). Job 0: calibration, 1: evaluation, k>1: member k-2.
                  res 22 2 8 6
                  case $p in 0) n=calib ;; 1) n=eval ;; *) n=member-$(printf %02d $((p - 2))) ;; esac
                  ENVS=("RWT_THREADS=2" "PYTHONDONTWRITEBYTECODE=1" "RWT_PARTON_CHANNEL=$ECHAN" "RWT_PARTON_SELECTION=production/ensemble/$ESEL-$n.json")
                  STEP=(bash production/run_parton_study.sh scalar-prepare)
                  OUTP=("parton-prepared.tgz rwt/$SET/prepared/$n.tgz"
                         "parton-prepare.json rwt/$SET/prepared/$n.json") ;;
    ensemble-prepare)
                  # Job 0: evaluation set; job k>0: member k-1 (production/ensemble/<channel>-*.json).
                  res 21 4 8 8
                  if [ "$p" = 0 ]; then n=eval; else n=member-$(printf %02d $((p - 1))); fi
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1" "RWT_PARTON_CHANNEL=$ECHAN"
                        "RWT_PARTON_RAPIDITY=transverse-mass-rapidity-v1" "RWT_PARTON_SELECTION=production/ensemble/$ESEL-$n.json")
                  STEP=(bash production/run_parton_study.sh ensemble-prepare)
                  OUTP=("parton-prepared.tgz rwt/$SET/prepared/$n.tgz"
                         "parton-prepare.json rwt/$SET/prepared/$n.json") ;;
    ensemble-fit|ensemble-smoke)
                  res 20 4 24 36
                  n=member-$(printf %02d "$p")
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1" "RWT_GPUS=1" "RWT_PARTON_CHANNEL=$ECHAN"
                        "RWT_PARTON_RAPIDITY=transverse-mass-rapidity-v1" "RWT_PARTON_SELECTION=production/ensemble/$ESEL-$n.json")
                  # Calibration = the original frozen 1M study bundle; evaluation = prepared job 0.
                  IN=("parton-prepared.tgz rwt/$SET/prepared/$n.tgz"
                      "calib-prepared.tgz $ECALIB"
                      "eval-prepared.tgz rwt/$SET/prepared/eval.tgz")
                  STEP=(bash production/run_parton_study.sh ensemble-fit)
                  OUTP=("ensemble-fit.tgz rwt/$SET/fits/$n.tgz"
                         "ensemble-fit.json rwt/$SET/fits/$n.json")
                  if [ "$run" = ensemble-smoke ]; then
                    # One job on member 00 with 20000 rows per class and one epoch: chain test only.
                    res 1 4 24 4; ENVS+=("RWT_ENSEMBLE_SMOKE_ROWS=20000")
                    OUTP=("${OUTP[@]/\/fits\//\/smoke\/}")
                  fi ;;
    bdt-fit3)
                  # v3 member k: bundles k, k+7, k+14 (mod 20) = 3M/class, minimum leaf size 100.
                  res 20 8 24 12
                  n=member-$(printf %02d "$p"); n2=member-$(printf %02d $(( (p + 7) % 20 ))); n3=member-$(printf %02d $(( (p + 14) % 20 )))
                  ENVS=("RWT_THREADS=8" "PYTHONDONTWRITEBYTECODE=1" "RWT_PARTON_CHANNEL=$ECHAN" "RWT_PARTON_SELECTION=production/ensemble/$ESEL-$n.json"
                        "RWT_DOMAIN=rho1" "RWT_BDT_MIN_LEAF=100" "RWT_SEED_OFFSET=1000")
                  IN=("parton-prepared.tgz rwt/$SET/prepared/$n.tgz" "extra1-prepared.tgz rwt/$SET/prepared/$n2.tgz" "extra2-prepared.tgz rwt/$SET/prepared/$n3.tgz"
                      "calib-prepared.tgz $ECALIB" "eval-prepared.tgz rwt/$SET/prepared/eval.tgz")
                  STEP=(bash production/run_parton_study.sh bdt-member)
                  OUTP=("bdt-member.tgz rwt/$SET/bdt3/$n.tgz" "bdt-member.json rwt/$SET/bdt3/$n.json") ;;
    table-only)
                  # Job 0: calibration bundle; 1: v1 evaluation set (confirmation); 2: v2 evaluation set; 3-22: v2 members.
                  res 23 1 4 2
                  case $p in
                    0) n=calib; src=rwt/parton-$ECH-1m-v1/prepared.tgz ;;
                    1) n=eval-v1; src=rwt/ensemble-$ECH-v1/prepared/eval.tgz ;;
                    2) n=eval-v2; src=rwt/ensemble-$ECH-v2/prepared/eval.tgz ;;
                    *) n=member-$(printf %02d $((p - 3))); src=rwt/ensemble-$ECH-v2/prepared/$n.tgz ;;
                  esac
                  ENVS=("RWT_THREADS=1" "PYTHONDONTWRITEBYTECODE=1" "RWT_PARTON_CHANNEL=$ECH")
                  IN=("parton-prepared.tgz $src")
                  STEP=(bash production/run_parton_study.sh table-only)
                  OUTP=("parton-tables.tgz rwt/$SET/$n.tgz" "table-only.json rwt/$SET/$n.json") ;;
    bdt-fit)
                  # BDT-17 member k on the frozen member-k subset; CPU only, threads = requested CPUs.
                  res 20 8 16 12
                  n=member-$(printf %02d "$p")
                  ENVS=("RWT_THREADS=8" "PYTHONDONTWRITEBYTECODE=1" "RWT_PARTON_CHANNEL=$ECHAN"
                        "RWT_PARTON_RAPIDITY=transverse-mass-rapidity-v1" "RWT_PARTON_SELECTION=production/ensemble/$ESEL-$n.json")
                  IN=("parton-prepared.tgz rwt/$SET/prepared/$n.tgz"
                      "calib-prepared.tgz $ECALIB"
                      "eval-prepared.tgz rwt/$SET/prepared/eval.tgz")
                  [ "$EVER" = v1 ] || ENVS+=("RWT_DOMAIN=rho1")
                  STEP=(bash production/run_parton_study.sh bdt-member)
                  OUTP=("bdt-member.tgz rwt/$SET/bdt/$n.tgz"
                         "bdt-member.json rwt/$SET/bdt/$n.json") ;;
    pdf-transport)
                  # PDF4LHC21 ratios from the paper's 5e5-event 13.6 TeV Z->ll NLO+PS production (bulk inputs
                  # fetched by the script itself with checksum verification).
                  res 1 1 4 6
                  ENVS=("RWT_THREADS=1" "PYTHONDONTWRITEBYTECODE=1" "RWT_PARTON_CHANNEL=$ECHAN")
                  STEP=(bash production/run_parton_study.sh pdf-transport)
                  OUTP=("pdf-transport.npz rwt/$SET/pdf-transport.npz"
                         "pdf-transport.json rwt/$SET/pdf-transport.json") ;;
    bdt-deliver|bdt-deliver4|bdt-deliver5)
                  # ONNX members + calibration/scale/PDF tables + apply/PCA scripts + test vectors + model card.
                  # PDF ratios: the Z->ll 13.6 TeV production for both channels (model card).
                  res 1 4 16 4
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1" "RWT_PARTON_CHANNEL=$ECHAN"
                        "RWT_PARTON_RAPIDITY=transverse-mass-rapidity-v1" "RWT_PARTON_SELECTION=production/ensemble/$ESEL-eval.json")
                  T2=bdt-tables.json; [ "$EVER" = v1 ] || T2=bdt-tables-r3.json
                  IN=("eval-prepared.tgz rwt/$SET/prepared/eval.tgz" "tables.json rwt/$SET/$T2"
                      "pdf-transport.npz rwt/ensemble-zll-v1/pdf-transport.npz")
                  for k in $(seq 0 19); do n=member-$(printf %02d "$k"); IN+=("members/$n.tgz rwt/$SET/bdt/$n.tgz"); done
                  [ "$EVER" = v1 ] || ENVS+=("RWT_DOMAIN=rho1")
                  STEP=(bash production/run_parton_study.sh bdt-deliver)
                  DR=; [ "$EVER" = v1 ] || DR=-r3
                  # r4: r3 tables and members, with the 1 % per-bin residual of the recipe (model card, calibration.json).
                  { [ "$run" = bdt-deliver4 ] || [ "$run" = bdt-deliver5 ] || [[ $ECH = atlas-* ]]; } && { DR=-r4; ENVS+=("RWT_RESIDUAL=0.01"); }
                  # r5: r4 with the scale and PDF factors also applied to the unit weight outside the domain (apply.py).
                  { [ "$run" = bdt-deliver5 ] || [[ $ECH = atlas-* ]]; } && DR=-r5
                  OUTP=("bdt-deliver.tgz rwt/$SET/deliver/ggzh-bdt17-$ECH$DR.tgz"
                         "bdt-deliver.json rwt/$SET/deliver/ggzh-bdt17-$ECH$DR.json") ;;
    bdt-assess|bdt-assess-variants|bdt-assess-variants2|bdt-assess-tail|bdt-cache|bdt-scaleband|bdt-tables2|bdt-assess3|bdt-confirm|bdt-cache4)
                  res 1 4 16 8
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1" "RWT_PARTON_CHANNEL=$ECHAN"
                        "RWT_PARTON_RAPIDITY=transverse-mass-rapidity-v1" "RWT_PARTON_SELECTION=production/ensemble/$ESEL-eval.json")
                  IN=("calib-prepared.tgz $ECALIB"
                      "eval-prepared.tgz rwt/$SET/prepared/eval.tgz")
                  for k in $(seq 0 19); do n=member-$(printf %02d "$k"); IN+=("members/$n.tgz rwt/$SET/bdt/$n.tgz"); done
                  [ "$EVER" = v1 ] || ENVS+=("RWT_DOMAIN=rho1")
                  # ATLAS-style sets: the frozen recipe of the CMS-style confirmation (1 % per-bin residual) from the start.
                  [[ $ECH = atlas-* ]] && ENVS+=("RWT_RESIDUAL=0.01")
                  STEP=(bash production/run_parton_study.sh bdt-assess)
                  OUTP=("ensemble-assess.json rwt/$SET/bdt-assess.json"
                         "ensemble-tables.json rwt/$SET/bdt-tables.json")
                  [ "$run" = bdt-assess ] && [ "$EVER" != v1 ] && OUTP+=("histogram-cache.json rwt/$SET/bdt-histogram-cache.json")
                  if [ "$run" = bdt-scaleband ]; then
                    # Calibration set: LO+PS scale band (Born-formula companion) against the reweighted bands.
                    # Revision r2: outside-domain variation binned in max(pT(Z), pT(H)).
                    ENVS+=("RWT_SCALEBAND=1"); IN+=("scale-inputs.tgz rwt/uncertainty-$ECH-1m-v1/scale-inputs.tgz")
                    # ATLAS-style: Born-formula companion at mu0 = H_T with PDF4LHC21 (Codex, calibration selection only).
                    [[ $ECH = atlas-* ]] && IN[-1]="scale-inputs.tgz rwt/scale-$ECH-calib-v2/scale-inputs.tgz"
                    OUTP=("scaleband-cache.json rwt/$SET/bdt-scaleband-cache-r3.json")
                  fi
                  if [ "$run" = bdt-assess3 ]; then
                    # v3 members; confirmatory test on the (once inspected) v2 evaluation set; cache included.
                    for k in $(seq 0 19); do n=member-$(printf %02d "$k"); IN=("${IN[@]/rwt\/$SET\/bdt\/$n.tgz/rwt\/$SET\/bdt3\/$n.tgz}"); done
                    OUTP=("ensemble-assess.json rwt/$SET/bdt3-assess.json" "ensemble-tables.json rwt/$SET/bdt3-tables.json"
                          "histogram-cache.json rwt/$SET/bdt3-histogram-cache.json")
                  fi
                  if [ "$run" = bdt-confirm ]; then
                    # Frozen recipe with a 1 % per-bin residual, tested on the v1 evaluation set: disjoint from the v2
                    # calibration, evaluation and member selections; member logits recomputed from model.pkl.
                    ENVS=("${ENVS[@]/production\/ensemble\/$ESEL-eval.json/production/ensemble/$ECH-eval.json}")
                    ENVS+=("RWT_CONFIRM=1" "RWT_RESIDUAL=0.01"); IN[1]="eval-prepared.tgz rwt/ensemble-$ECH-v1/prepared/eval.tgz"
                    OUTP=("ensemble-assess.json rwt/$SET/bdt-confirm-assess.json" "ensemble-tables.json rwt/$SET/bdt-confirm-tables.json"
                          "histogram-cache.json rwt/$SET/bdt-confirm-histogram-cache.json")
                  fi
                  if [ "$run" = bdt-tables2 ]; then
                    # Calibration tables only, revision r2 (outside-domain variation in max(pT(Z), pT(H))).
                    ENVS+=("RWT_TABLES_ONLY=1"); OUTP=("ensemble-tables.json rwt/$SET/bdt-tables-r3.json")
                  fi
                  if [ "$run" = bdt-cache4 ]; then
                    # Histogram cache with the raw LO+PS weights ('raw', the LO+PS x K shape) for the closure figures.
                    OUTP=("histogram-cache.json rwt/$SET/bdt-histogram-cache-r4.json")
                  fi
                  if [ "$run" = bdt-cache ]; then
                    # Same computation as bdt-assess; publishes only the fine-binned histogram cache.
                    OUTP=("histogram-cache.json rwt/$SET/bdt-histogram-cache.json")
                  fi
                  if [ "$run" = bdt-assess-tail ]; then
                    # Post-hoc: sparse tail bins merged from the high end to >= 500/1000/2000 raw LO+PS events.
                    ENVS+=("RWT_ASSESS_VARIANTS=1" "RWT_ASSESS_VARIANT_SET=tail")
                    OUTP=("ensemble-assess.json rwt/$SET/bdt-assess-tail.json"
                           "ensemble-tables.json rwt/$SET/bdt-tables-tail.json")
                  fi
                  if [ "$run" = bdt-assess-variants2 ]; then
                    # Post-hoc: bins of >= 2000 raw LO+PS events, VH(bb)-like domain (Delta phi(Z,H), pT(H)/pT(Z), m_ll window).
                    ENVS+=("RWT_ASSESS_VARIANTS=1" "RWT_ASSESS_VARIANT_SET=v2")
                    OUTP=("ensemble-assess.json rwt/$SET/bdt-assess-variants2.json"
                           "ensemble-tables.json rwt/$SET/bdt-tables-variants2.json")
                  fi
                  if [ "$run" = bdt-assess-variants ]; then
                    # Post-hoc diagnostics (wider bins, m_ll > 40 GeV) on the inspected evaluation set.
                    ENVS+=("RWT_ASSESS_VARIANTS=1")
                    OUTP=("ensemble-assess.json rwt/$SET/bdt-assess-variants.json"
                           "ensemble-tables.json rwt/$SET/bdt-tables-variants.json")
                  fi ;;
    ensemble-assess|ensemble-assess-smoke)
                  res 1 4 16 8
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1" "RWT_PARTON_CHANNEL=$ECHAN"
                        "RWT_PARTON_RAPIDITY=transverse-mass-rapidity-v1" "RWT_PARTON_SELECTION=production/ensemble/$ESEL-eval.json")
                  IN=("calib-prepared.tgz $ECALIB"
                      "eval-prepared.tgz rwt/$SET/prepared/eval.tgz")
                  d=fits; [ "$run" = ensemble-assess ] || d=smoke
                  for k in $(seq 0 19); do
                    [ "$d" = fits ] || [ "$k" = 0 ] || continue
                    n=member-$(printf %02d "$k"); IN+=("members/$n.tgz rwt/$SET/$d/$n.tgz")
                  done
                  STEP=(bash production/run_parton_study.sh ensemble-assess)
                  OUTP=("ensemble-assess.json rwt/$SET/assess-$d.json"
                         "ensemble-tables.json rwt/$SET/tables-$d.json") ;;
    parton-fit-1m)
                  res 15 4 24 36
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1" "RWT_GPUS=1")
                  IN=("parton-prepared.tgz rwt/parton-zll-1m-v1/prepared.tgz")
                  STEP=(bash production/run_parton_study.sh fit "$p")
                  OUTP=("parton-fit.tgz rwt/parton-zll-1m-v1/fits/fit-$p.tgz"
                         "parton-fit.json rwt/parton-zll-1m-v1/fits/fit-$p.json") ;;
    parton-jet-compare-1m)
                  res 1 4 24 8
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1")
                  IN=("parton-prepared.tgz rwt/parton-zll-1m-v1/prepared.tgz")
                  for k in $(seq 0 4); do
                    IN+=("baseline/fit-$k.tgz rwt/parton-zll-1m-v1/fits/fit-$k.tgz"
                         "jets/fit-$k.tgz rwt/$SET/fits/fit-$k.tgz")
                  done
                  STEP=(bash production/run_parton_study.sh jet-compare)
                  OUTP=("parton-jet-report.tgz rwt/$SET/comparison.tgz"
                         "parton-jet-summary.json rwt/$SET/summary.json"
                         "parton-jet-report/index.html rwt/$SET/index.html") ;;
    parton-bdt-1m)
                  res 1 8 24 24
                  ENVS=("RWT_THREADS=8" "PYTHONDONTWRITEBYTECODE=1")
                  IN=("parton-prepared.tgz rwt/parton-zll-1m-v1/prepared.tgz")
                  STEP=(bash production/run_parton_study.sh bdt)
                  OUTP=("parton-bdt.tgz rwt/parton-zll-1m-v1/bdt/control.tgz"
                         "parton-bdt.json rwt/parton-zll-1m-v1/bdt/result.json") ;;
    parton-fold0-report-1m)
                  res 1 4 16 4
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1")
                  IN=("parton-prepared.tgz rwt/parton-zll-1m-v1/prepared.tgz"
                      "reference-fit.tgz rwt/parton-zll-1m-v1/fits/fit-0.tgz"
                      "reference-bdt.tgz rwt/parton-zll-1m-v1/bdt/control.tgz")
                  STEP=(bash production/run_parton_fold0_report.sh)
                  OUTP=("parton-fold0-report.tgz rwt/parton-zll-1m-v1/fold0-report.tgz"
                         "parton-fold0-summary.json rwt/parton-zll-1m-v1/fold0-summary.json"
                         "parton-fold0-report/index.html rwt/parton-zll-1m-v1/fold0/index.html"
                         "parton-fold0-report/result.json rwt/parton-zll-1m-v1/fold0/result.json")
                  for v in pth ptj1 ptz lhe_mll; do OUTP+=("parton-fold0-report/$v.png rwt/parton-zll-1m-v1/fold0/$v.png" "parton-fold0-report/$v.pdf rwt/parton-zll-1m-v1/fold0/$v.pdf"); done ;;
    parton-reload-diagnose-1m)
                  res 1 2 8 2
                  ENVS=("RWT_THREADS=2" "PYTHONDONTWRITEBYTECODE=1" "RWT_GPUS=1")
                  IN=("reference-fit.tgz rwt/parton-zll-1m-v1/fits/fit-0.tgz"
                      "reference-input.npz particle-v1/zll/lo-cms/inputs/partons-v1/110011.npz")
                  STEP=(bash production/run_parton_reload_diagnose.sh)
                  OUTP=("parton-reload-diagnosis.json rwt/parton-zll-1m-v1/reload-diagnosis.json") ;;
    parton-early-report-1m)
                  res 1 4 16 4
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1")
                  IN=("parton-prepared.tgz rwt/parton-zll-1m-v1/prepared.tgz")
                  for k in $(seq 0 4); do IN+=("parts/fit-$k.tgz rwt/parton-zll-1m-v1/fits/fit-$k.tgz"); done
                  STEP=(bash production/run_parton_early_report.sh)
                  OUTP=("parton-early-report.tgz rwt/parton-zll-1m-v1/early-report.tgz"
                         "parton-early-summary.json rwt/parton-zll-1m-v1/early-summary.json"
                         "parton-early-report/index.html rwt/parton-zll-1m-v1/early/index.html"
                         "parton-early-report/result.json rwt/parton-zll-1m-v1/early/result.json")
                  for v in pth ptj1 ptz lhe_mll; do
                    OUTP+=("parton-early-report/$v.png rwt/parton-zll-1m-v1/early/$v.png"
                           "parton-early-report/$v.pdf rwt/parton-zll-1m-v1/early/$v.pdf")
                  done ;;
    parton-report-1m)
                  res 1 4 24 8
                  ENVS=("RWT_THREADS=4" "PYTHONDONTWRITEBYTECODE=1")
                  IN=("parton-prepared.tgz rwt/parton-zll-1m-v1/prepared.tgz"
                      "bdt/control.tgz rwt/parton-zll-1m-v1/bdt/control.tgz")
                  for k in $(seq 0 14); do IN+=("parts/fit-$k.tgz rwt/parton-zll-1m-v1/fits/fit-$k.tgz"); done
                  STEP=(bash production/run_parton_study.sh report)
                  OUTP=("parton-report.tgz rwt/parton-zll-1m-v1/report.tgz"
                         "parton-models.tgz rwt/parton-zll-1m-v1/models.tgz"
                         "parton-summary.json rwt/parton-zll-1m-v1/summary.json") ;;
    atlas-scale-validate)
                  res 1 2 4 2
                  IN=("atlas-scale-reference.lhe.gz @atlas/lhe/lo-atlas-13600/events-70001.lhe.gz"
                      "previous-control.json bdt-tables-v1/atlas-zll/scale-control/audit.json")
                  STEP=(bash production/run_atlas_scale_calibration.sh validate)
                  OUTP=("scale-audit.json rwt/$SET/control.json") ;;
    atlas-scale-calibration)
                  res 1 2 8 12
                  channel=${SET#scale-atlas-}; channel=${channel%-calib-v2}
                  IN=("control.json rwt/scale-atlas-control-v2/control.json"
                      "calib-prepared.tgz rwt/ensemble-atlas-$channel-v2/prepared/calib.tgz")
                  STEP=(bash production/run_atlas_scale_calibration.sh calibration "$channel")
                  OUTP=("scale-inputs.tgz rwt/$SET/scale-inputs.tgz" "scale-audit.json rwt/$SET/audit.json") ;;
    atlas-born-scale-control)
                  res 1 2 4 4
                  IN=("atlas-scale-reference.lhe.gz @atlas/lhe/lo-atlas-13600/events-70001.lhe.gz")
                  STEP=(bash production/run_atlas_scale_control.sh)
                  OUTP=("scale-audit.json bdt-tables-v1/atlas-zll/scale-control/audit.json") ;;
    grids-lo-atlas-vv)
                  res 1 8 16 12
                  IN=("ggzh-install.tgz @atlas/install/ggzh-install-v4.tgz")
                  ENVS=("GRID_SETTINGS=lhans1=93300 lhans2=93300 ebeam1=6800d0 ebeam2=6800d0 dynscale=3 bornonly=1 LOevents=1 vdecaymode=12")
                  STEP=(bash production/run_atlas_vv_grid.sh)
                  OUTP=("grids.tgz grids/lo-atlas-zvv-bdt-v1.tgz" "grid-check.json grids/lo-atlas-zvv-bdt-v1.json") ;;
    events-lo-atlas-vv)
                  res 2000 2 4 12; s=$((800001+p))
                  IN=("ggzh-install.tgz @atlas/install/ggzh-install-v4.tgz" "grids.tgz grids/lo-atlas-zvv-bdt-v1.tgz")
                  STEP=(bash production/run_atlas_vv_events.sh "$s")
                  OUTP=("events-$s.lhe.gz lhe/lo-atlas-zvv-bdt-v1/events-$s.lhe.gz") ;;
    scalar-lo-atlas|scalar-nlo-a14|parity-lo-atlas|parity-nlo-a14)
                  channel=zll; [[ $SET != tables-atlas-zvv-v1 ]] || channel=zvv
                  order=${run#*-}; scope=production
                  if [[ $order = lo-atlas ]]; then
                    res 2000 2 4 6; s=$((500001+p)); ss=$((600001+p)); expected=5000
                    lheprod="lhe/lo-atlas-parton10m-13600/events-$s.lhe.gz"
                    if [[ $channel = zvv ]]; then s=$((800001+p)); ss=$((900001+p)); lheprod="lhe/lo-atlas-zvv-bdt-v1/events-$s.lhe.gz"; fi
                  else
                    res 4000 2 4 6; s=$((100001+p)); ss=$((700001+p)); expected=2500
                    lheprod="lhe/nlo-rwt10m-13600/events-$s.lhe.gz"
                    if [[ $channel = zvv ]]; then s=$((300001+p)); ss=$((1000001+p)); lheprod="@vv/lhe/nlo-vv-rwt10m-13600/events-$s.lhe.gz"; fi
                  fi
                  IN=("input.lhe.gz $lheprod" "reference-install.tgz install/shower-install-rwt10m-0.tgz")
                  if [[ $run = parity-* ]]; then
                    res 5 2 4 6; scope=parity
                    IN+=("parity-reference.csv.gz bdt-tables-v1/parity-reference/atlas-zll/$order/$s.csv.gz")
                  fi
                  STEP=(bash production/run_atlas_scalar.sh "$channel" "$order" "$s" "$expected" "$ss")
                  OUTP=("table.csv.gz bdt-tables-v1/atlas-$channel/$scope/$order/closure/$s.csv.gz"
                         "table-schema.json bdt-tables-v1/atlas-$channel/$scope/$order/validation/$s.schema.json"
                         "table-validation.json bdt-tables-v1/atlas-$channel/$scope/$order/validation/$s.json") ;;
    events-lo-atlas-10m)
                  res "$NATLAS" 2 4 12; s=$((S0ATLAS+p))
                  IN=("ggzh-install.tgz @atlas/install/ggzh-install-v4.tgz"
                      "grids.tgz @atlas/grids/grids-lo-atlas-13600.tgz")
                  STEP=(bash production/run_atlas_events.sh "$s")
                  OUTP=("events-$s.lhe.gz lhe/lo-atlas-parton10m-13600/events-$s.lhe.gz") ;;
    particle-nlo-a14-10m|particle-lo-atlas-10m)
                  if [ "$run" = particle-nlo-a14-10m ]; then
                    res "$NNLO" 2 8 12; s=$((S0NLO+p)); ss=$((700001+p)); order=nlo-a14; expected=2500
                    lheprod="lhe/nlo-rwt10m-13600/events-$s.lhe.gz"
                  else
                    res "$NATLAS" 2 8 12; s=$((S0ATLAS+p)); ss=$((600001+p)); order=lo-atlas; expected=5000
                    lheprod="lhe/lo-atlas-parton10m-13600/events-$s.lhe.gz"
                  fi
                  IN=("pyarrow-20.0.0-cp312-cp312-manylinux_2_28_x86_64.whl particle-v1/install/pyarrow-20.0.0-cp312-cp312-manylinux_2_28_x86_64.whl"
                      "input.lhe.gz $lheprod"
                      "particle-install.tgz particle-v1/install/particle-install.tgz"
                      "reference-install.tgz install/shower-install-rwt10m-0.tgz")
                  STEP=(bash production/run_particle_atlas.sh "$order" "$s" "$expected" "$ss")
                  OUTP=("particle-out/records.parquet particle-v1/atlas-zll/$order/records/$s.parquet"
                         "particle-out/objects.parquet particle-v1/atlas-zll/$order/objects/$s.parquet"
                         "particle-out/jets-v1.npz particle-v1/atlas-zll/$order/inputs/jets-v1/$s.npz"
                         "particle-out/partons-v1.npz particle-v1/atlas-zll/$order/inputs/partons-v1/$s.npz"
                         "particle-out/closure.npz particle-v1/atlas-zll/$order/closure/$s.npz"
                         "particle-out/legacy.csv.gz particle-v1/atlas-zll/$order/closure/$s.csv.gz"
                         "particle-out/schema.json particle-v1/atlas-zll/$order/validation/$s.schema.json"
                         "particle-out/validation.json particle-v1/atlas-zll/$order/validation/$s.json") ;;
    particle-nlo|particle-lo-cms)
                  channel=zll; prefix=""; suffix=rwt10m
                  if [ "$SET" = particle-zvv ]; then channel=zvv; prefix="@vv/"; suffix=vv-rwt10m; fi
                  if [ "$run" = particle-nlo ]; then res "$NNLO" 1 8 12; s=$((S0NLO+p)); order=nlo; expected=2500; legacy=nlo-cp5
                  else res "$NCMS" 1 8 12; s=$((S0CMS+p)); order=lo-cms; expected=5000; legacy=lo-cms; fi
                  IN=("pyarrow-20.0.0-cp312-cp312-manylinux_2_28_x86_64.whl particle-v1/install/pyarrow-20.0.0-cp312-cp312-manylinux_2_28_x86_64.whl"
                      "input.lhe.gz ${prefix}lhe/$order-$suffix-13600/events-$s.lhe.gz"
                      "particle-install.tgz particle-v1/install/particle-install.tgz")
                  if [ "$channel" = zll ]; then
                    IN+=("split-map.json particle-v1/splits/zll.json")
                    IN+=("legacy.csv.gz tables/$legacy-rwt10m-13600/table-$s.csv.gz")
                  else
                    IN+=("reference-install.tgz install/shower-install-rwt10m-0.tgz")
                  fi
                  STEP=(bash production/run_particle.sh "$channel" "$order" "$s" "$expected")
                  OUTP=("particle-out/records.parquet particle-v1/$channel/$order/records/$s.parquet"
                         "particle-out/objects.parquet particle-v1/$channel/$order/objects/$s.parquet"
                         "particle-out/jets-v1.npz particle-v1/$channel/$order/inputs/jets-v1/$s.npz"
                         "particle-out/partons-v1.npz particle-v1/$channel/$order/inputs/partons-v1/$s.npz"
                         "particle-out/closure.npz particle-v1/$channel/$order/closure/$s.npz"
                         "particle-out/legacy.csv.gz particle-v1/$channel/$order/closure/$s.csv.gz"
                         "particle-out/schema.json particle-v1/$channel/$order/validation/$s.schema.json"
                         "particle-out/validation.json particle-v1/$channel/$order/validation/$s.json") ;;
    # ---- Pre-registered 10M test; pilot uses the first two seeds/products, which production must reuse.
    shower-build-rwt10m)
                  res "$NREP" 2 4 2
                  STEP=(bash -c "$S/build.sh shower-install && tar czf shower-install.tgz -C shower-install bin lib env.sh BUILD-MANIFEST.txt")
                  OUTP=("shower-install.tgz install/shower-install-rwt10m-$p.tgz") ;;
    events-nlo-10m|events-lo-cms-10m)
                  need_install
                  if [ "$run" = events-nlo-10m ]; then
                    res "$NNLO" 1 4 12; m=2500; s=$((S0NLO+p)); t=nlo
                    if [ "$SET" = vv-rwt10m ]; then need_grids grids-nlo-vv-13600; else need_grids grids-nlo-13600; fi; n="$G/cards/$RWL"
                  else
                    res "$NCMS" 1 4 12; m=5000; s=$((S0CMS+p)); t=lo-cms
                    if [ "$SET" = vv-rwt10m ]; then need_grids grids-lo-cms-vv-13600; else need_grids grids-lo-cms-13600; fi; n=none
                  fi
                  STEP=("$G/condor/run_job.sh" events "$m" "$s" "$n")
                  OUTP=("events-$s.lhe.gz lhe/$(sname "$t")/events-$s.lhe.gz") ;;
    table-nlo-cp5-10m|table-lo-cms-10m)
                  if [ "$run" = table-nlo-cp5-10m ]; then
                    res "$NNLO" 1 4 6; s=$((S0NLO+p)); need_lhe "$(sname nlo)" "$S0NLO" 1 "$p"; m=nlo; t=nlo-cp5
                  else
                    res "$NCMS" 1 4 6; s=$((S0CMS+p)); need_lhe "$(sname lo-cms)" "$S0CMS" 1 "$p"; m=cms-lo-veto; t=lo-cms
                  fi
                  ENVS=("GGZH_SHOWER_TABLE=1" "GGZH_SHOWER_TABLE_DIR=tables" "PYTHONDONTWRITEBYTECODE=1")
                  [ -z "$PILOT" ] || ENVS+=("GGZH_TABLE_CHECK=1")
                  STEP=("$S/run/shower_batch.sh" shower-out "$m" cp5 cms-common -- lhe/*.lhe.gz)
                  OUTP=("tables/table-$s.csv.gz tables/$(sname "$t")/table-$s.csv.gz") ;;
    rwt-clf-a-10m|rwt-clf-b-10m|rwt-clf-c-10m)
                  res "$NREP" 8 16 12; campaign_tables
                  n=${run#rwt-clf-}; n=${n%-10m}; t=${n^^}
                  ENVS=("RWT_THREADS=8" "PYTHONDONTWRITEBYTECODE=1")
                  STEP=(rwt/run.sh clf "clf-$n" lo-cms nlo-cp5 --features "$t" --model-out ensemble)
                  [ -z "$RWT_RETUNED" ] || STEP+=(--retuned)
                  OUTP=("weights-clf-$n.npz $(campaign_product weights-clf-$n "$p" npz)"
                        "rwt-out/bdt-clf-$n.json $(campaign_product bdt-clf-$n "$p" json)"
                        "ensemble.tgz $(campaign_product ensemble-clf-$n "$p" tgz)") ;;
    rwt-summary-10m)
                  res "$NREP" 8 16 6; campaign_tables; campaign_checks "$p"
                  ENVS=("RWT_THREADS=8" "PYTHONDONTWRITEBYTECODE=1")
                  STEP=(rwt/run.sh summary-10m)
                  OUTP=("rwt-out/rwt10m.json $(campaign_product summary "$p" json)") ;;
    rwt-transfer-10m)
                  res "$NREP" 8 16 8; campaign_tables
                  need_tables_as old-lo-cms lo-cms-13600 40001 100
                  need_tables_as old-nlo-cp5 nlo-cp5-13600 10001 200
                  IN+=("rwt-in/ensemble-clf-a.tgz $(campaign_product ensemble-clf-a "$p" tgz)"
                       "rwt-in/model-baseline.tgz rwt/model-13600.tgz")
                  ENVS=("RWT_THREADS=8" "PYTHONDONTWRITEBYTECODE=1")
                  STEP=(rwt/run.sh transfer-10m)
                  OUTP=("rwt-out/rwt10m.json $(campaign_product transfer "$p" json)") ;;
    rwt-model-10m)
                  res "$NREP" 8 16 12; campaign_tables
                  IN+=("rwt-in/summary.json $(campaign_product summary "$p" json)")
                  for x in skl2onnx-1.16.0-py2.py3-none-any.whl onnxconverter_common-1.14.0-py2.py3-none-any.whl; do
                    IN+=("wheels/$x sources/$x")
                  done
                  ENVS=("RWT_THREADS=8" "PYTHONDONTWRITEBYTECODE=1")
                  STEP=(rwt/run.sh model-10m)
                  OUTP=("model.tgz $(campaign_product model "$p" tgz)"
                        "rwt-out/rwt10m.json $(campaign_product model-status "$p" json)") ;;
    # ---- builds
    build)        res 1 8 24 6;  STEP=("$G/condor/run_job.sh" build); OUTP=("ggzh-install.tgz $INSTALL.tgz") ;;
    shower-build) res 1 2 4 2;   STEP=(bash -c "$S/build.sh shower-install && tar czf shower-install.tgz -C shower-install bin lib env.sh BUILD-MANIFEST.txt")
                  OUTP=("shower-install.tgz $SHOWER_INSTALL.tgz") ;;
    # ---- fixed order, ggxy Monte Carlo: LO 8 seeds x 2e6 points (2e5 for the optimisation), NLO 40 x 1e6 (1.5e5)
    fo-llh-lo)    res 2 4 8 6;   need_install; STEP=("$G/condor/run_job.sh" fo lo  $((2001 + 4 * p)) 4 2000000 200000 "$PDF" fo-ggllh "$E") ;;
    fo-llh-nlo)   res 10 4 8 16; need_install; STEP=("$G/condor/run_job.sh" fo nlo $((3001 + 4 * p)) 4 1000000 150000 "$PDF" fo-ggllh "$E") ;;
    fo-zh-lo)     res 1 2 4 6;   need_install; STEP=("$G/condor/run_job.sh" fo lo  $((5001 + 2 * p)) 2 2000000 200000 "$PDF" fo-ggzh "$E") ;;
    fo-zh-lo-pdf4lhc21)          # verification: the stable-Z LO of arXiv:2603.15762, Table 1
                  res 1 2 4 6;   need_install; STEP=("$G/condor/run_job.sh" fo lo  $((5001 + 2 * p)) 2 2000000 200000 PDF4LHC21_40_pdfas fo-ggzh 13600) ;;
    fo-zh-nlo)    res 10 4 8 16; need_install; STEP=("$G/condor/run_job.sh" fo nlo $((6001 + 4 * p)) 4 1000000 150000 "$PDF" fo-ggzh "$E") ;;
    # Z -> nu nubar (vv sets): the statistics of fo-llh, other seeds
    fo-vvh-lo)    res 2 4 8 6;   need_install; STEP=("$G/condor/run_job.sh" fo lo  $((12001 + 4 * p)) 4 2000000 200000 "$PDF" fo-ggvvh "$E") ;;
    fo-vvh-nlo)   res 10 4 8 16; need_install; STEP=("$G/condor/run_job.sh" fo nlo $((13001 + 4 * p)) 4 1000000 150000 "$PDF" fo-ggvvh "$E") ;;
    # ---- POWHEG fixed-order LO at the 7 scale points, 8 seeds
    powheg-lo)    res 1 8 8 12;  need_install; ENVS=("GRID_SETTINGS=$CARD"); STEP=("$G/condor/run_job.sh" lo 8) ;;
    # ---- grids: NLO 16 seeds (production with Born suppression L = 250 GeV), LO 8 seeds; 3 x-grid iterations
    grids-nlo)    res 1 16 32 20; need_install; ENVS=("GRID_SETTINGS=$CARD${CARD:+ }$([ "$SET" = verify ] || echo bornsuppfact=250)")
                  # vv13000: stage checkpoints on the staging area (two evictions after 4.7 and 5.7 h, 2026-10-16)
                  [ "$SET" = vv13000 ] && ENVS+=("GGZH_GRID_CKPT=$(sname grids-nlo)")
                  STEP=("$G/condor/run_job.sh" grids 16 3); OUTP=("grids.tgz grids/$(sname grids-nlo).tgz") ;;
    grids-nlo-bs250)             # verification: Born suppression with NNPDF3.1, against the unweighted grids-nlo
                  res 1 16 32 20; need_install; ENVS=("GRID_SETTINGS=bornsuppfact=250")
                  STEP=("$G/condor/run_job.sh" grids 16 3); OUTP=("grids.tgz grids/$(sname grids-nlo-bs250).tgz") ;;
    grids-lo)     res 1 8 16 12; need_install; ENVS=("GRID_SETTINGS=$CARD${CARD:+ }bornonly=1 LOevents=1")
                  STEP=("$G/condor/run_job.sh" grids 8 3); OUTP=("grids.tgz grids/$(sname grids-lo).tgz")
                  # atlas13600: also the LO cross section of the grids (the reweighting summary quotes it)
                  if [ "$SET" = atlas13600 ]; then OUTP+=("grids/grids-xsec.json grids/$(sname grids-lo).xsec.json"); fi ;;
    # ---- MCFM LO, stable Z, all cards of generator/mcfm/cards
    mcfm)         res 1 8 8 4;   ENVS=("MCFM_SQRTS=$E"); STEP=("$G/condor/run_job.sh" mcfm) ;;
    # ---- events, one seed per job; production NLO 200 x 2500 (weighted), LO 100 x 5000;
    #      verification NLO 50 x 1000 (unweighted and with Born suppression, same seeds), LO 25 x 2000
    events-nlo)   if [ "$SET" = verify ]; then res 50 1 4 8; m=1000; s0=1001; else res 200 1 4 12; m=2500; s0=$S0NLO; fi
                  need_install; need_grids "$(sname grids-nlo)"
                  STEP=("$G/condor/run_job.sh" events $m $((s0 + p)) "$G/cards/$RWL")
                  OUTP=("events-$((s0 + p)).lhe.gz lhe/$(sname nlo)/events-$((s0 + p)).lhe.gz") ;;
    events-nlo-bs250)
                  res 50 1 4 8; need_install; need_grids "$(sname grids-nlo-bs250)"
                  STEP=("$G/condor/run_job.sh" events 1000 $((1001 + p)) "$G/cards/$RWL")
                  OUTP=("events-$((1001 + p)).lhe.gz lhe/$(sname nlo-bs250)/events-$((1001 + p)).lhe.gz") ;;
    events-lo)    if [ "$SET" = verify ]; then res 25 1 4 8; m=2000; s0=2001; else res 100 1 4 12; m=5000; s0=$S0LO; fi
                  need_install; need_grids "$(sname grids-lo)"
                  STEP=("$G/condor/run_job.sh" events $m $((s0 + p)) "$G/cards/$RWL")
                  OUTP=("events-$((s0 + p)).lhe.gz lhe/$(sname lo)/events-$((s0 + p)).lhe.gz") ;;
    # ---- cache closure: the same 50 events with the cache on and off, 3 seeds
    cache)        res 3 1 4 12; need_install; need_grids "$(sname grids-nlo)"
                  s0=$([ "$SET" = verify ] && echo 51000 || echo 52000)
                  STEP=("$G/condor/run_job.sh" cache 50 $((s0 + p)) "$G/cards/$RWL") ;;
    # ---- PDF and alpha_s at fixed order, member by member (13600; independent check of the event
    #      weights): LO, one seed (7001), 2e6 points, PDF4LHC21 members 0-42, 4 per job
    # top-mass scheme (sec. m_t scheme): MS-bar with mu_t = m_ZH / f, the seeds and statistics of the on-shell run
    fo-llh-lo-mt*|fo-llh-nlo-mt*|fo-zh-lo-mt*|fo-zh-nlo-mt*|fo-vvh-lo-mt*|fo-vvh-nlo-mt*)
                  f=${run##*-mt}; base=${run%-mt*}
                  case $base in
                    fo-llh-lo)  res 2 4 8 6;   s0=2001; ns=4; np=2000000; npo=200000; prog=fo-ggllh; mode=lo ;;
                    fo-llh-nlo) res 10 4 8 16; s0=3001; ns=4; np=1000000; npo=150000; prog=fo-ggllh; mode=nlo ;;
                    fo-zh-lo)   res 1 2 4 6;   s0=5001; ns=2; np=2000000; npo=200000; prog=fo-ggzh;  mode=lo ;;
                    fo-zh-nlo)  res 10 4 8 16; s0=6001; ns=4; np=1000000; npo=150000; prog=fo-ggzh;  mode=nlo ;;
                    fo-vvh-lo)  res 2 4 8 6;   s0=12001; ns=4; np=2000000; npo=200000; prog=fo-ggvvh; mode=lo ;;   # vv sets
                    fo-vvh-nlo) res 10 4 8 16; s0=13001; ns=4; np=1000000; npo=150000; prog=fo-ggvvh; mode=nlo ;;
                  esac
                  need_install
                  STEP=("$G/condor/run_job.sh" fo $mode $((s0 + ns * p)) $ns $np $npo "$PDF" $prog "$E" 1 "$(awk -v f="$f" 'BEGIN { print 1 / f }')") ;;
    fo-llh-pdf|fo-zh-pdf)
                  res 11 4 8 4; need_install; m=$((4 * p)); nm=$(( p == 10 ? 3 : 4 ))
                  STEP=("$G/condor/run_job.sh" fo-members lo 7001 2000000 200000 "$PDF" "$m" "$nm"
                        "$([ "$run" = fo-llh-pdf ] && echo fo-ggllh || echo fo-ggzh)" "$E") ;;
    # ---- matching variation h_damp = infinity (13600): POWHEG damps nothing for hdamp <= 0 (no remnant
    #      events), so hdamp=-1 replaces the card's hdamp 250; otherwise as the production (Born
    #      suppression, 200 x 2500 NLO events, NLO+PS with Monash, parton level)
    grids-nlo-hdampinf)
                  res 1 16 32 20; need_install; ENVS=("GRID_SETTINGS=$CARD bornsuppfact=250 hdamp=-1")
                  STEP=("$G/condor/run_job.sh" grids 16 3); OUTP=("grids.tgz grids/$(sname grids-nlo-hdampinf).tgz") ;;
    events-nlo-hdampinf)
                  res 200 1 4 12; need_install; need_grids "$(sname grids-nlo-hdampinf)"
                  STEP=("$G/condor/run_job.sh" events 2500 $((S0HD + p)) "$G/cards/$RWL")
                  OUTP=("events-$((S0HD + p)).lhe.gz lhe/$(sname nlo-hdampinf)/events-$((S0HD + p)).lhe.gz") ;;
    parton-nlo-hdampinf)
                  res 5 4 8 4; need_lhe "$(sname nlo-hdampinf)" $S0HD 40 "$p"
                  STEP=("$S/run/shower_batch.sh" shower-out noshower monash -- lhe/*.lhe.gz) ;;
    shower-nlo-hdampinf-monash)
                  res 10 8 16 6; need_lhe "$(sname nlo-hdampinf)" $S0HD 20 "$p"
                  STEP=("$S/run/shower_batch.sh" shower-out nlo monash ps-variations -- lhe/*.lhe.gz) ;;
    # ---- matching variations of the POWHEG veto in Pythia (13600): the production NLO events
    #      re-showered with POWHEG:pThard = 1 or POWHEG:pTdef = 2 (cards veto-*.cmnd after nlo.cmnd)
    shower-nlo-pthard1-monash|shower-nlo-ptdef2-monash)
                  m=${run#shower-nlo-}; m=${m%-monash}
                  res 10 8 16 6; need_lhe "$(sname nlo)" $S0NLO 20 "$p"
                  STEP=("$S/run/shower_batch.sh" shower-out nlo monash ps-variations "veto-$m" -- lhe/*.lhe.gz) ;;
    # ---- cache timing scan (13600): the same 50 events with 0 to 450 weights, cache on and off,
    #      3 jobs (machines); lists in generator/cards (rwl-scan/: prefixes of one list)
    cache-scan)   res 3 1 4 6; need_install; need_grids "$(sname grids-nlo)"
                  STEP=("$G/condor/run_job.sh" cache-scan 50 $((53000 + p)) none "$G/cards/rwl-scan/pwg-rwl-scan-9.dat"
                        "$G/cards/rwl-scan/pwg-rwl-scan-52.dat" "$G/cards/pwg-rwl.dat" "$G/cards/pwg-rwl-pdf4lhc21.dat"
                        "$G/cards/rwl-scan/pwg-rwl-scan-300.dat" "$G/cards/rwl-scan/pwg-rwl-scan-450.dat") ;;
    # ---- showers: NLO+PS (POWHEG veto) with shower-scale variation weights; LO+PS shower starts
    shower-nlo-bs250-monash)
                  res 1 8 16 4; need_lhe "$(sname nlo-bs250)" 1001 50 "$p"
                  STEP=("$S/run/shower_batch.sh" shower-out nlo monash ps-variations -- lhe/*.lhe.gz) ;;
    shower-nlo-*) t=${run#shower-nlo-}
                  if [ "$SET" = verify ]; then res 1 8 16 4; need_lhe "$(sname nlo)" 1001 50 "$p"
                  else res 10 8 16 6; need_lhe "$(sname nlo)" $S0NLO 20 "$p"; fi
                  STEP=("$S/run/shower_batch.sh" shower-out nlo "$t" ps-variations -- lhe/*.lhe.gz) ;;
    shower-lo-*)  t=${run##*-}; m=${run#shower-lo-}; m=${m%-"$t"}
                  if [ "$SET" = verify ]; then res 1 8 16 4; need_lhe "$(sname lo)" 2001 25 "$p"
                  else res 5 8 16 6; need_lhe "$(sname lo)" $S0LO 20 "$p"; fi
                  STEP=("$S/run/shower_batch.sh" shower-out "lo-$m" "$t" -- lhe/*.lhe.gz) ;;
    # ---- the production events at parton level: Rivet on the LHE kinematics (shower/run/shower.sh
    #      mode noshower: no shower, MPI or hadronisation), all files
    parton-nlo-bs250)
                  res 1 4 8 4; need_lhe "$(sname nlo-bs250)" 1001 50 "$p"
                  STEP=("$S/run/shower_batch.sh" shower-out noshower monash -- lhe/*.lhe.gz) ;;
    parton-nlo)   if [ "$SET" = verify ]; then res 1 4 8 4; need_lhe "$(sname nlo)" 1001 50 "$p"
                  else res 5 4 8 4; need_lhe "$(sname nlo)" $S0NLO 40 "$p"; fi
                  STEP=("$S/run/shower_batch.sh" shower-out noshower monash -- lhe/*.lhe.gz) ;;
    parton-lo)    if [ "$SET" = verify ]; then res 1 4 8 4; need_lhe "$(sname lo)" 2001 25 "$p"
                  else res 5 4 8 4; need_lhe "$(sname lo)" $S0LO 20 "$p"; fi
                  STEP=("$S/run/shower_batch.sh" shower-out noshower monash -- lhe/*.lhe.gz) ;;
    # ---- the reweighting study (rwt/README.md; 13600 only).  CMS-style LO+PS with public CMS settings: LO events
    #      at NNPDF3.1 (the card's PDF) with the fixed scale m_Z + m_H, showered with CP5, the CMS common settings and
    #      the CMS POWHEG veto; the same LO events with Pythia's default start; the production NLO events showered with
    #      CP5 and the CMS common settings (the target).  One row per event (Rivet GGZH_EVENTTABLE), staged as products.
    grids-lo-cms) res 1 8 16 12; need_install; ENVS=("GRID_SETTINGS=${VDM:+$VDM }bornonly=1 LOevents=1 fixedscale=1")
                  STEP=("$G/condor/run_job.sh" grids 8 3)
                  OUTP=("grids.tgz grids/$(sname grids-lo-cms).tgz" "grids/grids-xsec.json grids/$(sname grids-lo-cms).xsec.json") ;;
    events-lo-cms) res 100 1 4 12; need_install; need_grids "$(sname grids-lo-cms)"
                  STEP=("$G/condor/run_job.sh" events 5000 $((S0CMS + p)) none)
                  OUTP=("events-$((S0CMS + p)).lhe.gz lhe/$(sname lo-cms)/events-$((S0CMS + p)).lhe.gz") ;;
    table-lo-cms|table-lo-cms-default|table-nlo-cp5)
                  case $run in
                    table-nlo-cp5) res 200 1 4 6; s=$((S0NLO + p)); need_lhe "$(sname nlo)" $S0NLO 1 "$p"; m=nlo ;;
                    table-lo-cms)  res 100 1 4 6; s=$((S0CMS + p)); need_lhe "$(sname lo-cms)" $S0CMS 1 "$p"; m=cms-lo-veto ;;
                    *)             res 100 1 4 6; s=$((S0CMS + p)); need_lhe "$(sname lo-cms)" $S0CMS 1 "$p"; m=lo-default ;;
                  esac
                  ENVS=("GGZH_SHOWER_TABLE=1" "GGZH_SHOWER_TABLE_DIR=tables")
                  STEP=("$S/run/shower_batch.sh" shower-out "$m" cp5 cms-common -- lhe/*.lhe.gz)
                  OUTP=("tables/table-$s.csv.gz tables/$(sname "${run#table-}")/table-$s.csv.gz") ;;
    # classifier weights (rwt/run.sh clf): CMS-style and default-start source, five inputs or six (leading-jet pT)
    rwt-clf|rwt-clf-default|rwt-clf-jet|rwt-clf-jet-default)
                  res 1 1 8 6; n=${run#rwt-}
                  case $n in *-default) m=lo-cms-default ;; *) m=lo-cms ;; esac
                  t=(); case $n in clf-jet*) t=(--jet-feature) ;; esac
                  need_tables "$m" $S0CMS 100; need_tables nlo-cp5 $S0NLO 200
                  STEP=(rwt/run.sh clf "$n" "$m" nlo-cp5 "${t[@]}")
                  OUTP=("weights-$n.npz rwt/$(sname "weights-$n").npz" "rwt-out/bdt-$n.json rwt/$(sname "bdt-$n").json") ;;
    rwt-summary)  res 1 1 8 4
                  need_tables lo-cms $S0CMS 100; need_tables lo-cms-default $S0CMS 100; need_tables nlo-cp5 $S0NLO 200
                  for n in clf clf-default clf-jet clf-jet-default; do
                    IN+=("rwt-in/weights-$n.npz rwt/$(sname "weights-$n").npz" "rwt-in/bdt-$n.json rwt/$(sname "bdt-$n").json")
                  done
                  IN+=("rwt-in/grids-xsec.json grids/$(sname grids-lo-cms).xsec.json")
                  STEP=(rwt/run.sh summary) ;;
    # the reusable models of the data release: classifier (ONNX) and binned map (correctionlib), nominal and the seven
    # scale variations, trained on all events of the CMS-style source and the target; the ONNX converter wheels
    # (rwt/requirements-export.txt) are staged as products sources/<wheel>
    rwt-model)    res 1 1 8 8
                  need_tables lo-cms $S0CMS 100; need_tables nlo-cp5 $S0NLO 200
                  IN+=("rwt-in/weights-clf.npz rwt/$(sname weights-clf).npz")
                  for x in skl2onnx-1.16.0-py2.py3-none-any.whl onnxconverter_common-1.14.0-py2.py3-none-any.whl; do
                    IN+=("wheels/$x sources/$x")
                  done
                  STEP=(rwt/run.sh model)
                  case $SET in vv*)
                    STEP=(rwt/run.sh model lo-cms nlo-cp5 cms vv)
                    IN+=("rwt-in/bdt-clf.json rwt/$(sname bdt-clf).json"
                         "rwt-in/bdt-clf-jet.json rwt/$(sname bdt-clf-jet).json") ;;
                  esac
                  OUTP=("model.tgz rwt/$(sname model).tgz") ;;
    # ---- the reweighting study with an ATLAS-style source (set atlas13600; rwt/README.md).  LO events of run
    #      events-lo (dynscale 3, PDF4LHC21) showered with A14 as athena sets it (tune a14-athena), the ATLAS base
    #      settings (atlas-common) and the ATLAS POWHEG veto (mode atlas-lo-veto); target: the 13600 production NLO
    #      events with the same A14 settings and the package's NLO veto, PDF4LHC21 nominal.  One row per event.
    table-lo-atlas|table-nlo-a14)
                  case $run in
                    table-lo-atlas) res 100 1 4 6; s=$((S0LO + p)); need_lhe "$(sname lo)" $S0LO 1 "$p"; m=atlas-lo-veto; t=tables/lo-atlas-$E ;;
                    *)              res 200 1 4 6; s=$((S0NLO + p)); need_lhe "nlo-$E" $S0NLO 1 "$p"; m=nlo; t=tables/nlo-a14-$E ;;
                  esac
                  ENVS=("GGZH_SHOWER_TABLE=1" "GGZH_SHOWER_TABLE_DIR=tables")
                  STEP=("$S/run/shower_batch.sh" shower-out "$m" a14-athena atlas-common -- lhe/*.lhe.gz)
                  OUTP=("tables/table-$s.csv.gz $t/table-$s.csv.gz") ;;
    # classifier weights, five inputs or six (leading-jet pT), ATLAS-style source towards the A14 target
    rwt-clf-atlas|rwt-clf-jet-atlas)
                  res 1 1 8 6; n=${run#rwt-}
                  t=(); case $n in clf-jet*) t=(--jet-feature) ;; esac
                  need_tables_as lo-atlas "lo-atlas-$E" $S0LO 100; need_tables_as nlo-a14 "nlo-a14-$E" $S0NLO 200
                  STEP=(rwt/run.sh clf "$n" lo-atlas nlo-a14 "${t[@]}")
                  OUTP=("weights-$n.npz rwt/weights-$n-$E.npz" "rwt-out/bdt-$n.json rwt/bdt-$n-$E.json") ;;
    rwt-summary-atlas)
                  res 1 1 8 4
                  need_tables_as lo-atlas "lo-atlas-$E" $S0LO 100; need_tables_as nlo-a14 "nlo-a14-$E" $S0NLO 200
                  for n in clf-atlas clf-jet-atlas; do
                    IN+=("rwt-in/weights-$n.npz rwt/weights-$n-$E.npz" "rwt-in/bdt-$n.json rwt/bdt-$n-$E.json")
                  done
                  IN+=("rwt-in/grids-xsec.json grids/$(sname grids-lo).xsec.json")
                  STEP=(rwt/run.sh summary-atlas) ;;
    rwt-model-atlas)
                  res 1 1 8 8
                  need_tables_as lo-atlas "lo-atlas-$E" $S0LO 100; need_tables_as nlo-a14 "nlo-a14-$E" $S0NLO 200
                  IN+=("rwt-in/weights-clf.npz rwt/weights-clf-atlas-$E.npz")
                  for n in clf-atlas clf-jet-atlas; do
                    IN+=("rwt-in/bdt-$n.json rwt/bdt-$n-$E.json")
                  done
                  for x in skl2onnx-1.16.0-py2.py3-none-any.whl onnxconverter_common-1.14.0-py2.py3-none-any.whl; do
                    IN+=("wheels/$x sources/$x")
                  done
                  STEP=(rwt/run.sh model lo-atlas nlo-a14 atlas)
                  OUTP=("model.tgz rwt/model-atlas-$E.tgz") ;;
    *) echo "unknown run $run in set $SET" >&2; exit 2 ;;
  esac
  if [ -n "$RWT_AVAILABLE_MANIFEST" ]; then
    ENVS+=("RWT_INPUT_MANIFEST=$RWT_AVAILABLE_MANIFEST")
  fi
  if [[ $SET = parton-zvv-stable-1m-v1 || $SET = parton-zvv-jets-stable-1m-v1 ]]; then
    ENVS+=("RWT_PARTON_CHANNEL=zvv" "RWT_PARTON_RAPIDITY=transverse-mass-rapidity-v1")
    # Keep prepared data and BDT original; require corrected fits for every report.
    IN=("${IN[@]/rwt\/parton-zll-1m-v1/rwt\/parton-zvv-1m-v1}")
    IN=("${IN[@]/rwt\/parton-zvv-1m-v1\/fits/rwt\/parton-zvv-stable-1m-v1\/fits}")
    OUTP=("${OUTP[@]/rwt\/parton-zll-1m-v1/rwt\/$SET}")
    if [[ $SET = parton-zvv-jets-stable-1m-v1 ]]; then
      ENVS+=("RWT_PARTON_CONTEXT=analysis-jets-v1")
      if [[ $run = parton-fit-1m ]]; then res 5 4 24 36; fi
    fi
  elif [[ $SET = parton-zll-jets-1m-v1 || $SET = parton-zvv-jets-1m-v1 ]]; then
    ENVS+=("RWT_PARTON_CONTEXT=analysis-jets-v1")
    if [[ $run = parton-fit-1m ]]; then
      res 5 4 24 36
      OUTP=("${OUTP[@]/rwt\/parton-zll-1m-v1/rwt\/$SET}")
    fi
    if [[ $SET = parton-zvv-jets-1m-v1 ]]; then
      ENVS+=("RWT_PARTON_CHANNEL=zvv")
      IN=("${IN[@]/rwt\/parton-zll-1m-v1/rwt\/parton-zvv-1m-v1}")
    fi
  elif [[ $SET = parton-zvv-1m-v1 ]]; then
    ENVS+=("RWT_PARTON_CHANNEL=zvv")
    IN=("${IN[@]/rwt\/parton-zll-1m-v1/rwt\/parton-zvv-1m-v1}")
    OUTP=("${OUTP[@]/rwt\/parton-zll-1m-v1/rwt\/parton-zvv-1m-v1}")
  elif [[ $SET = parton-zll-1m-recovery-v2 && $run != parton-fit-1m ]]; then
    # Recover missing fit products in the original namespace; publish new reports separately.
    OUTP=("${OUTP[@]/rwt\/parton-zll-1m-v1/rwt\/parton-zll-1m-recovery-v2}")
  fi
  [ -z "$RWT_RETUNED" ] || ENVS+=("RWT_RETUNED=1")
}

known() { runs | tr ' ' '\n' | grep -qx "$1" || { echo "run $1 is not in set $SET" >&2; exit 2; }; }

case $CMD_ in
  list)
    for r in $(runs); do def "$r" 0; echo "$r $JOBS $CPUS $MEM $HOURS"; done ;;
  inputs|outputs|show|run)
    RUN=${3:?run}; P=${4:?job index}; known "$RUN"; def "$RUN" "$P"
    [ "$P" -ge 0 ] && [ "$P" -lt "$JOBS" ] || { echo "job index $P outside 0..$((JOBS - 1))" >&2; exit 2; }
    case $CMD_ in
      inputs)  for x in "${IN[@]}"; do echo "$x"; done ;;
      outputs) for x in "${OUTP[@]}"; do echo "$x"; done ;;
      show)    echo "set $SET run $RUN job $P/$JOBS cpus $CPUS memory ${MEM}GB hours $HOURS"
               for x in "${ENVS[@]}"; do echo "env $x"; done; echo "cmd ${STEP[*]}" ;;
      run)
        echo "production/paper.sh: set $SET run $RUN job $P/$JOBS, package $(head -1 "$PKG/VERSION")"
        for x in "${IN[@]}"; do f=${x%% *}; [ -e "$f" ] || { echo "FATAL: input $f missing (${x#* })" >&2; exit 2; }; done
        export NCPU=$CPUS
        for x in "${ENVS[@]}"; do export "$x"; done
        if [ -f shower-install.tgz ]; then
          mkdir -p shower-install && tar xzf shower-install.tgz -C shower-install || exit 2
          export GGZH_SHOWER_PREFIX=$PWD/shower-install
        fi
        "${STEP[@]}"; rc=$?
        if [ -d shower-out ]; then cp -p shower-out/sum.raw.json "$RUN-$P.raw.json" 2>/dev/null; tar czf "out-$RUN-$P.tgz" shower-out; fi
        for x in "${OUTP[@]}"; do f=${x%% *}; [ "$rc" != 0 ] || [ -e "$f" ] || { echo "FATAL: product $f missing" >&2; rc=3; }; done
        exit $rc ;;
    esac ;;
  *) echo "unknown command $CMD_" >&2; exit 2 ;;
esac
