#!/usr/bin/env python3
"""Generate the amplitude-cache patch for POWHEG-BOX-V2/ggxy_ggZH.

Reads the process' pwhg_ggxy.f and Makefile, applies exact string replacements
(each asserted to match exactly once, so a drifted upstream source fails loudly
instead of producing a silently wrong patch), writes the patched files into an
output process directory and a unified diff next to this script.

    make_patched_source.py <process dir> <patched process dir> [patch file]

The Fortran hooks store the alpha_s-stripped squared amplitude keyed on the raw
momenta, the top mass (and mu_R for the two-loop virtual) and the flavour list;
see ampcache.cpp for the rationale.
"""
from __future__ import annotations

import difflib
import shutil
import sys
from pathlib import Path


def replace_once(text: str, old: str, new: str, *, where: str) -> str:
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{where}: expected exactly one match, found {count}:\n{old}")
    return text.replace(old, new)


FORTRAN_EDITS: list[tuple[str, str, str]] = [
    # ---- ggllh1lfunc: declarations and cached evaluation ----------------------
    (
        "ggllh1lfunc declarations",
        "      real(dp),    intent(out) :: mat1l\n"
        "      real(dp)                 :: mt,mut\n",
        "      real(dp),    intent(out) :: mat1l\n"
        "      real(dp)                 :: mt,mut\n"
        "      real(dp)                 :: ckey(4*nlegs+2)\n"
        "      integer                  :: ikey(2),found\n",
    ),
    (
        "ggllh1lfunc body",
        "      ! all other parameters have been initialized in init_couplings\n"
        "      ! mt(mut) is read on the c++ side\n"
        "      mat1l = st_alpha**2*ggllh1l(p)\n"
        "\n"
        "      end subroutine ggllh1lfunc\n",
        "      ! all other parameters have been initialized in init_couplings\n"
        "      ! mt(mut) is read on the c++ side\n"
        "      ! VHTNP amplitude cache (ampcache.cpp): POWHEG's reweighting re-evaluates\n"
        "      ! the same kinematics for every weight; store the alpha_s-stripped value.\n"
        "      ckey(1:4*nlegs) = reshape(p,(/4*nlegs/))\n"
        "      ckey(4*nlegs+1) = mt\n"
        "      ckey(4*nlegs+2) = mut\n"
        "      ikey(1) = 1\n"
        "      ikey(2) = nlegs\n"
        "      call ampcache_lookup(4*nlegs+2,ckey,2,ikey,mat1l,found)\n"
        "      if(found.eq.1) then\n"
        "         mat1l = st_alpha**2*mat1l\n"
        "         return\n"
        "      endif\n"
        "      mat1l = ggllh1l(p)\n"
        "      call ampcache_store(4*nlegs+2,ckey,2,ikey,mat1l)\n"
        "      mat1l = st_alpha**2*mat1l\n"
        "\n"
        "      end subroutine ggllh1lfunc\n",
    ),
    # ---- ggllh2lfunc: the two-loop virtual, keyed additionally on mu_R ---------
    (
        "ggllh2lfunc declarations",
        "      real(dp),    intent(out) :: mat2l\n"
        "      real(dp)                 :: mt,mut\n",
        "      real(dp),    intent(out) :: mat2l\n"
        "      real(dp)                 :: mt,mut\n"
        "      real(dp)                 :: ckey(4*nlegs+3)\n"
        "      integer                  :: ikey(2),found\n",
    ),
    (
        "ggllh2lfunc body",
        "      ! all other parameters have been initialized in init_couplings\n"
        "      ! mt(mut) is read on the c++ side\n"
        "      mat2l = st_alpha**2*ggllh2l(p,dsqrt(st_muren2),mut)\n"
        "\n"
        "      end subroutine ggllh2lfunc\n",
        "      ! all other parameters have been initialized in init_couplings\n"
        "      ! mt(mut) is read on the c++ side\n"
        "      ! VHTNP amplitude cache (ampcache.cpp); the finite two-loop part depends\n"
        "      ! explicitly on mu_R, which is therefore part of the key.\n"
        "      ckey(1:4*nlegs) = reshape(p,(/4*nlegs/))\n"
        "      ckey(4*nlegs+1) = mt\n"
        "      ckey(4*nlegs+2) = mut\n"
        "      ckey(4*nlegs+3) = dsqrt(st_muren2)\n"
        "      ikey(1) = 2\n"
        "      ikey(2) = nlegs\n"
        "      call ampcache_lookup(4*nlegs+3,ckey,2,ikey,mat2l,found)\n"
        "      if(found.eq.1) then\n"
        "         mat2l = st_alpha**2*mat2l\n"
        "         return\n"
        "      endif\n"
        "      mat2l = ggllh2l(p,dsqrt(st_muren2),mut)\n"
        "      call ampcache_store(4*nlegs+3,ckey,2,ikey,mat2l)\n"
        "      mat2l = st_alpha**2*mat2l\n"
        "\n"
        "      end subroutine ggllh2lfunc\n",
    ),
    # ---- calc_matsq_rec: Recola one-loop Born and reals ------------------------
    (
        "calc_matsq_rec declarations",
        "      integer  :: pr,ipowQCD,iorder,prec_flag\n"
        "      real(dp) :: wdp,muR,mt,mut\n",
        "      integer  :: pr,ipowQCD,iorder,prec_flag\n"
        "      real(dp) :: wdp,muR,mt,mut\n"
        "      real(dp) :: ckey(4*nlegs+1)\n"
        "      integer  :: ikey(nlegs+2),found\n",
    ),
    (
        "calc_matsq_rec lookup",
        "      !set scale to some sensible value\n"
        "      muR = mt\n"
        "      call set_alphas_rcl(st_alpha,muR,nf)\n"
        "      call set_mu_uv_rcl(muR)\n"
        "      call set_mu_ir_rcl(muR)\n"
        "\n"
        "      !get pr,ipowQCD and iorder\n"
        "      call get_pr_flav(nlegs,flav,pr,ipowQCD,iorder)\n"
        "\n"
        "      !calculate squared matrix element\n"
        "      call calc_matsq(pr,iorder,0,nlegs,p)\n",
        "      !get pr,ipowQCD and iorder\n"
        "      call get_pr_flav(nlegs,flav,pr,ipowQCD,iorder)\n"
        "\n"
        "      ! VHTNP amplitude cache (ampcache.cpp): at fixed order the squared\n"
        "      ! amplitude is an exact power of alpha_s, so the stored value is\n"
        "      ! alpha_s-free and the current st_alpha is re-applied on a hit.\n"
        "      ! Only the 6-leg real emissions are cached here: setborn caches the Born\n"
        "      ! together with its spin correlations, which Recola computes from the\n"
        "      ! state this call leaves behind (compute_bmunu_rec).\n"
        "      ckey(1:4*nlegs) = reshape(p,(/4*nlegs/))\n"
        "      ckey(4*nlegs+1) = mt\n"
        "      ikey(1:nlegs) = flav(1:nlegs)\n"
        "      ikey(nlegs+1) = nlegs\n"
        "      ikey(nlegs+2) = 3\n"
        "      found = 0\n"
        "      if(nlegs.ge.6) call ampcache_lookup(4*nlegs+1,ckey,nlegs+2,ikey,wlo,found)\n"
        "      if(found.eq.1) then\n"
        "         wlo = wlo*st_alpha**ipowQCD\n"
        "         return\n"
        "      endif\n"
        "\n"
        "      !set scale to some sensible value\n"
        "      muR = mt\n"
        "      call set_alphas_rcl(st_alpha,muR,nf)\n"
        "      call set_mu_uv_rcl(muR)\n"
        "      call set_mu_ir_rcl(muR)\n"
        "\n"
        "      !calculate squared matrix element\n"
        "      call calc_matsq(pr,iorder,0,nlegs,p)\n",
    ),
    (
        "calc_matsq_rec store",
        "            print*,\"Recola + Cuttools(qp) & Oneloop(qp) = \",wlo\n"
        "            print*,\"\"\n"
        "         endif\n"
        "      endif\n"
        "\n"
        "      end subroutine calc_matsq_rec\n",
        "            print*,\"Recola + Cuttools(qp) & Oneloop(qp) = \",wlo\n"
        "            print*,\"\"\n"
        "         endif\n"
        "      endif\n"
        "\n"
        "      if(nlegs.ge.6) call ampcache_store(4*nlegs+1,ckey,nlegs+2,ikey,wlo/st_alpha**ipowQCD)\n"
        "\n"
        "      end subroutine calc_matsq_rec\n",
    ),
]

BORN_EDITS: list[tuple[str, str, str]] = [
    (
        "setborn declarations",
        "      include 'PhysPars.h'\n      integer nlegs\n      parameter (nlegs=nlegborn)\n",
        "      include 'PhysPars.h'\n      include 'pwhg_st.h'\n      integer nlegs\n      parameter (nlegs=nlegborn)\n"
        "      real * 8 ckey(4*nlegs+1),cval(1+16*nlegs)\n      integer ikey(nlegs+1),found\n",
    ),
    (
        "setborn cached Born and bmunu",
        "      if(with_recola) then\n"
        "        call calc_matsq_rec(nlegs,p,bflav,born)\n"
        "        do j=1,nlegs\n"
        "          if(abs(bflav(j)).gt.6) cycle\n"
        "          call compute_bmunu_rec(1,'NLO',2,j,bmunu(0:3,0:3,j))\n",
        "      if(with_recola) then\n"
        "        ! VHTNP amplitude cache: the Born and its spin-correlated tensor are\n"
        "        ! cached together (alpha_s^2 stripped), because compute_bmunu_rec reads\n"
        "        ! Recola's momenta and amplitudes from the last calc_matsq call.\n"
        "        ckey(1:4*nlegs) = reshape(p,(/4*nlegs/))\n"
        "        ckey(4*nlegs+1) = ph_topmass\n"
        "        ikey(1:nlegs) = bflav(1:nlegs)\n"
        "        ikey(nlegs+1) = 4\n"
        "        call ampcache_lookupv(4*nlegs+1,ckey,nlegs+1,ikey,1+16*nlegs,cval,found)\n"
        "        if(found.eq.1) then\n"
        "          born = cval(1)*st_alpha**2\n"
        "          bmunu = reshape(cval(2:1+16*nlegs),(/4,4,nlegs/))*st_alpha**2\n"
        "        else\n"
        "          call calc_matsq_rec(nlegs,p,bflav,born)\n"
        "          do j=1,nlegs\n"
        "            if(abs(bflav(j)).gt.6) cycle\n"
        "            call compute_bmunu_rec(1,'NLO',2,j,bmunu(0:3,0:3,j))\n"
        "          enddo\n"
        "          cval(1) = born/st_alpha**2\n"
        "          cval(2:1+16*nlegs) = reshape(bmunu,(/16*nlegs/))/st_alpha**2\n"
        "          call ampcache_storev(4*nlegs+1,ckey,nlegs+1,ikey,1+16*nlegs,cval)\n"
        "        endif\n"
        "        do j=1,nlegs\n"
        "          if(abs(bflav(j)).gt.6) cycle\n",
    ),
]

MAKEFILE_EDITS: list[tuple[str, str, str]] = [
    (
        "Makefile USER objects",
        "USER=ggxy_interface.o ggxy.o pwhg_ggxy.o init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o \\\n",
        "USER=ampcache.o ggxy_interface.o ggxy.o pwhg_ggxy.o init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o \\\n",
    ),
]


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src = Path(sys.argv[1]).resolve()
    dst = Path(sys.argv[2]).resolve()
    patch_out = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).with_name("0003-vhtnp-amplitude-cache.patch")
    here = Path(__file__).resolve().parent

    if dst.exists():
        raise SystemExit(f"refusing to overwrite {dst}")
    shutil.copytree(src, dst, symlinks=True)

    diff_lines: list[str] = []
    for name, edits in (("pwhg_ggxy.f", FORTRAN_EDITS), ("Born.f", BORN_EDITS), ("Makefile", MAKEFILE_EDITS)):
        original = (src / name).read_text()
        patched = original
        for where, old, new in edits:
            patched = replace_once(patched, old, new, where=where)
        (dst / name).write_text(patched)
        diff_lines += difflib.unified_diff(
            original.splitlines(keepends=True),
            patched.splitlines(keepends=True),
            fromfile=f"a/ggxy_ggZH/{name}",
            tofile=f"b/ggxy_ggZH/{name}",
        )
    shutil.copy2(here / "ampcache.cpp", dst / "ampcache.cpp")
    cpp = (here / "ampcache.cpp").read_text().splitlines(keepends=True)
    diff_lines += difflib.unified_diff([], cpp, fromfile="/dev/null", tofile="b/ggxy_ggZH/ampcache.cpp")
    patch_out.write_text("".join(diff_lines))
    # objects that must be rebuilt: make decides by mtime, and the copies above are new
    for stale in ("pwhg_main",):
        (dst / stale).unlink(missing_ok=True)
    print(f"patched process written to {dst}; patch {patch_out} ({len(diff_lines)} diff lines)")


if __name__ == "__main__":
    main()
