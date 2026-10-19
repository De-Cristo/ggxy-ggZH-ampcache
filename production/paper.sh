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
# Jobs p = 0 .. <jobs>-1; seeds are fixed functions of p (below), identical in 13600 and 13000.
#
# A job directory contains the unpacked package (VERSION, production/, generator/, shower/) and
# the files listed by "inputs".  Products have storage-independent names: install/<name>.tgz,
# grids/<name>.tgz, lhe/<sample>/events-<seed>.lhe.gz.  Everything else a job writes (out-*.tgz,
# *.json, logs) is its result and stays in the job directory.
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PKG=$(cd "$HERE/.." && pwd)
CMD_=${1:?list|inputs|outputs|show|run}; SET=${2:?set: common|13600|13000|verify|vv13600|vv13000}

TUNES=(monash cp5 a14)
LOMODES=(default scalup scalup-half scalup-quarter power)
MTFAC=(1 2 4)                       # MS-bar top-mass scale mu_t = m_ZH / f

S0NLO=10001                        # first seed of the production NLO events
S0LO=20001; S0HD=10001; S0CMS=40001  # first seeds of the LO, h_damp = inf and CMS-style LO events
VDM=""                             # Z decay setting added to the grids of the CMS-style LO (vv sets: vdecaymode=12)
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
INSTALL=install/ggzh-install-v2; SHOWER_INSTALL=install/shower-install-v2
case $SET in vv*) INSTALL=install/ggzh-install-v3 ;; esac
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
  local k; for k in $(seq 0 $(($3 - 1))); do IN+=("tables/$1/table-$(($2 + k)).csv.gz tables/$(sname "$1")/table-$(($2 + k)).csv.gz"); done
}
res() { JOBS=$1; CPUS=$2; MEM=$3; HOURS=$4; }

# sample names: production lhe/<nlo|lo>-<sqrt(s)>, verification lhe/<...>-nnpdf31-13600, Z -> nu nubar <...>-vv-<sqrt(s)>
sname() { case $SET in verify) echo "$1-nnpdf31-13600" ;; vv*) echo "$1-vv-$E" ;; *) echo "$1-$E" ;; esac; }

def() {
  local run=$1 p=$2 m t
  ENVS=(); STEP=(); IN=(); OUTP=()
  case $run in
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
                  STEP=("$G/condor/run_job.sh" grids 8 3); OUTP=("grids.tgz grids/$(sname grids-lo).tgz") ;;
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
                  OUTP=("model.tgz rwt/$(sname model).tgz") ;;
    *) echo "unknown run $run in set $SET" >&2; exit 2 ;;
  esac
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
