#!/usr/bin/env python3
"""FRESCO benchmarks of physics/transfer (tests/transfer_fresco_test.cpp; docs/upgrade/transfer.md section 9).

Writes, for every case below, a directory fresco/<case>/ with the FRESCO input of the reference run (job.nin, and
fort4.dat when the job reads a table on FRESCO's input file 4) and reference.txt: the numbers the test reads, taken
from FRESCO's output of the reference run and of a refined run (smaller step, larger radii), so that FRESCO's own
convergence is part of the data.  The test never runs FRESCO.

    python3 make_fresco_benchmarks.py --fresco <path to the fresco 3.5 binary> [--out <dir>]

Constants.  FRESCO fixes hbar c = 197.32705 MeV fm and u = 931.49432 MeV in its source; its namelist takes the unit
mass (in u) and 1/alpha.  Every job sets unitmass and finec so that FRESCO's 2u/hbar^2 and e^2 are those of
physics/Constants.h (CODATA 2018): unitmass = (931.49410242/931.49432) (197.32705/197.3269804)^2 and
finec = 197.32705 / (197.3269804/137.035999084).  FRESCO echoes both in its output ("2*amu/hbar^2", "e^2").

Masses are mass numbers (u), so that m_a = m_s + m_x and m_F = m_A + m_x hold exactly and the Jacobi coordinates of
both codes are the same; the binding energies and Q values are inputs.

Lines of reference.txt (whitespace-separated; '#' starts a comment; a case's "drop" leaves out kinds its test does
not read):
  run <name> <settings>                 the FRESCO settings of the following lines (reference or refined)
  constants <2u/hbar^2> <e^2>          as FRESCO prints them
  kin <projectile> <E_cm> <k> <eta>     per partition of the incoming energy (MeV, fm^-1)
  pot <kp> <type> <r> <Re> <Im>         FRESCO's listing of the potential components (fort.34; the Coulomb term per
                                        unit Z1 Z2, the spin-orbit term before 2 l.s), at r = 1 ... 10 fm
  S <l> <2j> <2J> <Re> <Im>             elastic S-matrix elements (fort.7)
  ratio <theta> <sigma/sigma_R>         elastic, c.m. (deg)
  xs <theta> <mb/sr>                    the transfer cross section, c.m. (deg)
  pol <theta> <iT11>                    the projectile's vector analysing power of the transfer (kqmax=1, Madison)
  reaction <mb>                         FRESCO's cumulative reaction cross section
  bound <kn> <scale> <rms> <D0> <ANC>   FRESCO's bound-state listing: depth scale factor (VR), rms radius (fm),
                                        D0 (MeV fm^3/2), ANC (fm^-1/2, at R = rmatch - h)
  ancratio <kn> <r> <value>             u(r)/W(2 kappa r) from fort.46 (the ANC plateau)
  bin <E_min> <E_max> <NK>              a continuum bin's energies (MeV, c.m.) and nodes
  binphase <ik> <k> <E> <delta>         its phase shift (deg) at node ik (every 20th node and the last)
  binnorm <norm> <rms> <D0>             FRESCO's "Norm = ..., RMS = ... & D0 = ..." of the bin (on [0, RSP])
  binwf <N> <u>                         the bin u(r) at r = (N - 1) hcm (every 50th point of FRESCO's listing)

Continuum bins (BE < 0; FRESCO's namelist manual, "CONTINUUM BINS") are built on [0, RSP]; the &FRESCO variable rsp
sets RSP apart from rmatch, so that the bin (the object compared) stays the same while rmatch is converged.  The
bins' widths put a node of the bin's envelope near RSP (docs/upgrade/transfer.md section 15).
"""
import argparse
import datetime
import math
import os
import re
import shutil
import subprocess
import tempfile

HBARC = 197.3269804          # physics/Constants.h
UCONV = 931.49410242
ALPHA = 1.0 / 137.035999084
FRESCO_HBARC = 197.32705     # FRESCO's source (frxx0.f)
FRESCO_AMU = 931.49432
UNITMASS = (UCONV / FRESCO_AMU) * (FRESCO_HBARC / HBARC) ** 2
FINEC = FRESCO_HBARC / (ALPHA * HBARC)
E2 = ALPHA * HBARC

HEAD = """{title}
NAMELIST
 &FRESCO hcm={hcm} rmatch={rmatch} {nonlocal_part}
     jtmin=0.0 jtmax={jtmax} absend=-1 {extra}
     thmin={thmin} thmax=180.0 thinc=5.0
     iter={iter} iblock={iblock} chans=1 smats=2 xstabl=1 treneg=1
     unitmass={unitmass:.10f} finec={finec:.8f}
     elab(1)={elab} /
"""

# The cases: the namelists after &FRESCO, and the reference and refined settings.
CASES = {
    # p + 15N elastic at 5 MeV, spin 1/2, spin-orbit; target spin ignored (0+).
    "p15N_elastic": dict(
        title="p + 15N elastic at E_lab = 5 MeV, spin-orbit", elab=5.0, jtmax=20.5, iter=0, iblock=0, thmin=10.0,
        body="""\
 &PARTITION namep='p' massp=1.0 zp=1 nex=1 namet='15N' masst=15.0 zt=7 qval=0.0 /
 &STATES jp=0.5 ptyp=1 ep=0.0 cpot=1 jt=0.0 ptyt=1 et=0.0 /
 &partition /
 &pot kp=1 type=0 p(1:3)= 15.0 0.0 1.25 /
 &pot kp=1 type=1 p(1:6)= 50.0 1.2 0.65 1.0 1.3 0.5 /
 &pot kp=1 type=2 p(1:6)= 0.0 1.3 0.5 6.0 1.3 0.5 /
 &pot kp=1 type=3 p(1:3)= 6.0 1.1 0.6 /
 &pot /
 &overlap /
 &coupling /
""",
        reference=dict(hcm=0.025, rmatch=25.0), refined=dict(hcm=0.0125, rmatch=30.0)),
    # d + 15N elastic at 10 MeV, spin 1, spin-orbit.
    "d15N_elastic": dict(
        title="d + 15N elastic at E_lab = 10 MeV, spin 1, spin-orbit", elab=10.0, jtmax=30, iter=0, iblock=0,
        thmin=10.0,
        body="""\
 &PARTITION namep='d' massp=2.0 zp=1 nex=1 namet='15N' masst=15.0 zt=7 qval=0.0 /
 &STATES jp=1.0 ptyp=1 ep=0.0 cpot=1 jt=0.0 ptyt=1 et=0.0 /
 &partition /
 &pot kp=1 type=0 p(1:3)= 15.0 0.0 1.3 /
 &pot kp=1 type=1 p(1:6)= 90.0 1.15 0.75 2.0 1.35 0.6 /
 &pot kp=1 type=2 p(1:6)= 0.0 1.35 0.6 10.0 1.35 0.6 /
 &pot kp=1 type=3 p(1:3)= 6.0 1.07 0.66 /
 &pot /
 &overlap /
 &coupling /
""",
        reference=dict(hcm=0.025, rmatch=25.0), refined=dict(hcm=0.0125, rmatch=30.0)),
    # 6Li + 12C elastic at 20 MeV: strong Coulomb and absorption; radii r0 (A_p^1/3 + A_t^1/3).
    "li6c12_elastic": dict(
        title="6Li + 12C elastic at E_lab = 20 MeV", elab=20.0, jtmax=60, iter=0, iblock=0, thmin=5.0,
        body="""\
 &PARTITION namep='6Li' massp=6.0 zp=3 nex=1 namet='12C' masst=12.0 zt=6 qval=0.0 /
 &STATES jp=1.0 ptyp=1 ep=0.0 cpot=1 jt=0.0 ptyt=1 et=0.0 /
 &partition /
 &pot kp=1 type=0 p(1:3)= 12.0 6.0 0.73 /
 &pot kp=1 type=1 p(1:6)= 120.0 0.80 0.85 35.0 0.95 0.85 /
 &pot /
 &overlap /
 &coupling /
""",
        reference=dict(hcm=0.025, rmatch=30.0), refined=dict(hcm=0.0125, rmatch=35.0)),
}

DN_BODY = """\
 &PARTITION namep='d' massp=2.0 zp=1 nex=1 namet='15N' masst=15.0 zt=7 qval=0.0 /
 &STATES jp=1.0 ptyp=1 ep=0.0 cpot=1 jt=0.5 ptyt=-1 et=0.0 /
 &PARTITION namep='n' massp=1.0 zp=0 nex=1 namet='16O' masst=16.0 zt=8 qval=9.9028 /
 &STATES jp=0.5 ptyp=1 ep=0.0 cpot=2 jt=0.0 ptyt=1 et=0.0 /
 &partition /
 &pot kp=1 type=0 p(1:3)= 15.0 0.0 1.3 /
 &pot kp=1 type=1 p(1:6)= 90.0 1.15 0.75 2.0 1.35 0.6 /
 &pot kp=1 type=2 p(1:6)= 0.0 1.35 0.6 10.0 1.35 0.6 /
 &pot kp=2 type=0 p(1:3)= 16.0 0.0 1.3 /
 &pot kp=2 type=1 p(1:6)= 48.0 1.2 0.65 1.0 1.3 0.5 /
 &pot kp=2 type=2 p(1:6)= 0.0 1.3 0.5 6.0 1.3 0.5 /
 &pot kp=3 type=0 p(1:3)= 1.0 0.0 1.0 /
 &pot kp=3 type=1 p(1:3)= 60.0 1.5 0.5 /
 &pot kp=4 type=0 p(1:3)= 15.0 0.0 1.25 /
 &pot kp=4 type=1 p(1:3)= 50.0 1.25 0.65 /
 &pot kp=4 type=3 p(1:3)= 6.0 1.25 0.65 /
 &pot /
 &overlap kn1=1 ic1=1 ic2=2 in=1 kind=0 nn=1 l=0 sn=0.5 j=0.5 kbpot=3 be=2.2246 isc=1 ipc=1 /
 &overlap kn1=2 ic1=1 ic2=2 in=2 kind=0 nn=1 l=1 sn=0.5 j=0.5 kbpot=4 be=12.1274 isc=1 ipc=1 /
 &overlap /
 &coupling icto=2 icfrom=1 {coupling} /
 &cfp in=1 ib=1 ia=1 kn=1 a=1.0 /
 &cfp in=-2 ib=1 ia=1 kn=2 a=1.0 /
 &coupling /
"""

# 12C(6Li,d)16O(g.s.): alpha transfer, alpha + d (2S, one node) and alpha + 12C (0+, two nodes).  {krpot} and {kp5}
# give the post-form interaction its nuclear part only (KRPOT = 5: the alpha-d well at the depth found, with the
# Coulomb term of the binding potential cancelled by a table of -V_C on FRESCO's input file 4).
LI_BODY = """\
 &PARTITION namep='6Li' massp=6.0 zp=3 nex=1 namet='12C' masst=12.0 zt=6 qval=0.0 /
 &STATES jp=1.0 ptyp=1 ep=0.0 cpot=1 jt=0.0 ptyt=1 et=0.0 /
 &PARTITION namep='d' massp=2.0 zp=1 nex=1 namet='16O' masst=16.0 zt=8 qval=5.6883 /
 &STATES jp=1.0 ptyp=1 ep=0.0 cpot=2 jt=0.0 ptyt=1 et=0.0 /
 &partition /
 &pot kp=1 type=0 p(1:3)= 12.0 6.0 0.73 /
 &pot kp=1 type=1 p(1:6)= 120.0 0.80 0.85 35.0 0.95 0.85 /
 &pot kp=2 type=0 p(1:3)= 16.0 0.0 1.3 /
 &pot kp=2 type=1 p(1:6)= 88.0 1.17 0.75 1.0 1.32 0.65 /
 &pot kp=2 type=2 p(1:6)= 0.0 1.32 0.65 12.0 1.32 0.65 /
 &pot kp=3 type=0 p(1:3)= 1.0 0.0 1.9 /
 &pot kp=3 type=1 p(1:3)= 70.0 1.9 0.65 /
 &pot kp=4 type=0 p(1:3)= 12.0 4.0 0.9 /
 &pot kp=4 type=1 p(1:3)= 100.0 0.9 0.65 /
{kp5} &pot /
 &overlap kn1=1 ic1=1 ic2=2 in=1 kind=0 nn=2 l=0 sn=0.0 j=0.0 kbpot=3 {krpot} be=1.4743 isc=1 ipc=1 /
 &overlap kn1=2 ic1=1 ic2=2 in=2 kind=0 nn=3 l=0 sn=0.0 j=0.0 kbpot=4 be=7.1620 isc=1 ipc=1 /
 &overlap /
 &coupling icto=2 icfrom=1 {coupling} /
 &cfp in=1 ib=1 ia=1 kn=1 a=1.0 /
 &cfp in=-2 ib=1 ia=1 kn=2 a=1.0 /
 &coupling /
"""
# The alpha-d depth (MeV) of the nuclear-part vertex: physics/transfer's SolveBoundState for kp=3 (h = 0.05 fm);
# FRESCO's own search for kp=3 prints VR = 1.1114 (77.798 MeV), the same to its four digits.
LI_AD_DEPTH = "77.796542"
LI_KP5 = (" &pot kp=5 type=0 p(1:3)= 1.0 0.0 1.9 /\n"
          f" &pot kp=5 type=1 p(1:3)= {LI_AD_DEPTH} 1.9 0.65 /\n"
          " &pot kp=5 type=1 shape=7 p(1:3)= -1.0 0.0 1.0 /\n")

FR_DN = "kind=7 ip1=0 ip2=0 ip3=0"
CASES.update({
    "dn_finite_range": dict(
        title="15N(d,n)16O(g.s.) at E_lab = 10 MeV, finite range, post, no remnant", elab=10.0, jtmax=30, iter=1,
        iblock=1, thmin=0.0, body=DN_BODY.format(coupling=FR_DN),
        reference=dict(hcm=0.025, rmatch=30.0, rintp=0.05, hnl=0.025, rnl=10.0, centre=0.0, nnu=72),
        refined=dict(hcm=0.02, rmatch=30.0, rintp=0.04, hnl=0.02, rnl=10.0, centre=0.0, nnu=72)),
    "dn_zero_range": dict(
        title="15N(d,n)16O(g.s.) at E_lab = 10 MeV, zero range D0 = -130.31 MeV fm^3/2, recoil (inh=2)", elab=10.0,
        jtmax=30, iter=1, iblock=1, thmin=0.0, extra="inh=2", body=DN_BODY.format(coupling="kind=5 p1=-130.31 p2=0.0"),
        reference=dict(hcm=0.025, rmatch=30.0), refined=dict(hcm=0.0125, rmatch=35.0)),
    "li6d_finite_range": dict(
        title="12C(6Li,d)16O(g.s.) at E_lab = 20 MeV, finite range, post (nuclear part of V_alpha-d), no remnant",
        elab=20.0, jtmax=40, iter=1, iblock=1, thmin=0.0, extra="epc=0.00001 cutl=-6.4", table=True,
        body=LI_BODY.format(kp5=LI_KP5, krpot="krpot=5", coupling=FR_DN),
        reference=dict(hcm=0.0125, rmatch=30.0, rintp=0.05, hnl=0.025, rnl=14.0, centre=-0.1, nnu=72),
        refined=dict(hcm=0.0125, rmatch=30.0, rintp=0.025, hnl=0.0125, rnl=14.0, centre=-0.1, nnu=72),
        also=[dict(hcm=0.025, rmatch=30.0, rintp=0.05, hnl=0.025, rnl=14.0, centre=-0.1, nnu=72, cutl=-3.2),
              dict(hcm=0.05, rmatch=30.0, rintp=0.1, hnl=0.05, rnl=14.0, centre=-0.1, nnu=72, cutl=-1.6)]),
    "li6d_zero_range": dict(
        title="12C(6Li,d)16O(g.s.) at E_lab = 20 MeV, zero range D0 = 68.3056 MeV fm^3/2, recoil (inh=2)", elab=20.0,
        jtmax=40, iter=1, iblock=1, thmin=0.0, extra="inh=2",
        body=LI_BODY.format(kp5="", krpot="", coupling="kind=5 p1=68.3056 p2=0.0"),
        reference=dict(hcm=0.025, rmatch=30.0), refined=dict(hcm=0.0125, rmatch=35.0)),
    # The whole alpha-d binding potential (with its Coulomb tail) as the interaction: FRESCO does not converge in
    # its non-local range here (transfer.md section 9); kept as data for the note, not compared by the test.
    "li6d_finite_range_coulomb": dict(
        title="12C(6Li,d)16O(g.s.) at E_lab = 20 MeV, finite range, post (whole V_alpha-d), no remnant", elab=20.0,
        jtmax=40, iter=1, iblock=1, thmin=0.0, extra="epc=0.00001",
        body=LI_BODY.format(kp5="", krpot="", coupling=FR_DN),
        reference=dict(hcm=0.05, rmatch=30.0, rintp=0.1, hnl=0.05, rnl=22.0, centre=-3.9, nnu=72),
        refined=dict(hcm=0.05, rmatch=30.0, rintp=0.1, hnl=0.05, rnl=26.0, centre=-3.9, nnu=72),
        also=[dict(hcm=0.05, rmatch=30.0, rintp=0.1, hnl=0.05, rnl=16.0, centre=-3.9, nnu=72),
              dict(hcm=0.05, rmatch=30.0, rintp=0.1, hnl=0.05, rnl=12.0, centre=-0.6, nnu=72)]),
})


# 15N(d,n)16O(g.s.) of DN_BODY with a spin-orbit term in both distortions (transfer.md section 19): d + 15N that of
# d15N_elastic, n + 16O a Woods-Saxon derivative of its own; the deuteron's vector analysing power iT11 (kqmax=1, the
# Madison frame: FRESCO's koords=0) is read with the cross section.
DN_SO_BODY = DN_BODY.replace(
    " &pot kp=2 type=0", " &pot kp=1 type=3 p(1:3)= 6.0 1.07 0.66 /\n &pot kp=2 type=0", 1).replace(
    " &pot kp=3 type=0", " &pot kp=2 type=3 p(1:3)= 5.5 1.2 0.65 /\n &pot kp=3 type=0", 1)
CASES.update({
    "dn_finite_range_so": dict(
        title="15N(d,n)16O(g.s.) at E_lab = 10 MeV, spin-orbit in d + 15N and n + 16O, finite range, post, no remnant",
        elab=10.0, jtmax=30, iter=1, iblock=1, thmin=0.0, extra="kqmax=1 pp=0", body=DN_SO_BODY.format(coupling=FR_DN),
        reference=dict(hcm=0.025, rmatch=30.0, rintp=0.05, hnl=0.025, rnl=10.0, centre=0.0, nnu=72),
        refined=dict(hcm=0.02, rmatch=30.0, rintp=0.04, hnl=0.02, rnl=10.0, centre=0.0, nnu=72),
        # larger radii; the radial cut-off cutl (at hcm = 0.0125 without it ERWIN loses 14 digits at high L and the
        # 0 deg cross section moves by 3 %: FRESCO's numerics, transfer.md section 20)
        also=[dict(hcm=0.025, rmatch=35.0, rintp=0.05, hnl=0.025, rnl=12.0, centre=0.0, nnu=96),
              dict(hcm=0.025, rmatch=30.0, rintp=0.05, hnl=0.025, rnl=10.0, centre=0.0, nnu=72, cutl=-3.2),
              dict(hcm=0.0125, rmatch=30.0, rintp=0.025, hnl=0.0125, rnl=10.0, centre=0.0, nnu=72, cutl=-6.4)],
        drop=("S",)),
    "dn_zero_range_so": dict(
        title="15N(d,n)16O(g.s.) at E_lab = 10 MeV, spin-orbit in d + 15N and n + 16O, zero range D0 = -130.31 MeV "
              "fm^3/2, recoil (inh=2)", elab=10.0, jtmax=30, iter=1, iblock=1, thmin=0.0, extra="inh=2 kqmax=1 pp=0",
        body=DN_SO_BODY.format(coupling="kind=5 p1=-130.31 p2=0.0"),
        reference=dict(hcm=0.025, rmatch=30.0), refined=dict(hcm=0.0125, rmatch=35.0), drop=("S",)),
})

# 15N(d,n)[p + 15N] to continuum bins of p + 15N (E_lab = 10 MeV): the projectile overlap of DN_BODY, the bin of
# l, j in the real Woods-Saxon (+ Coulomb) potential kp=5 (no spin-orbit), E_x = 12.1274 MeV + E (the p + 15N
# threshold), the bin weight ISC.  The three bins (transfer.md section 15): non-resonant (l = 1), over a narrow
# d-wave resonance (Gamma < width) and inside a wide s-wave resonance (Gamma > width).
DN_BIN_BODY = """\
 &PARTITION namep='d' massp=2.0 zp=1 nex=1 namet='15N' masst=15.0 zt=7 qval=0.0 /
 &STATES jp=1.0 ptyp=1 ep=0.0 cpot=1 jt=0.5 ptyt=-1 et=0.0 /
 &PARTITION namep='n' massp=1.0 zp=0 nex=1 namet='16O' masst=16.0 zt=8 qval=9.9028 /
 &STATES jp=0.5 ptyp=1 ep=0.0 cpot=2 jt={JF} ptyt={PF} et={EX:.4f} /
 &partition /
 &pot kp=1 type=0 p(1:3)= 15.0 0.0 1.3 /
 &pot kp=1 type=1 p(1:6)= 90.0 1.15 0.75 2.0 1.35 0.6 /
 &pot kp=1 type=2 p(1:6)= 0.0 1.35 0.6 10.0 1.35 0.6 /
 &pot kp=2 type=0 p(1:3)= 16.0 0.0 1.3 /
 &pot kp=2 type=1 p(1:6)= 48.0 1.2 0.65 1.0 1.3 0.5 /
 &pot kp=2 type=2 p(1:6)= 0.0 1.3 0.5 6.0 1.3 0.5 /
 &pot kp=3 type=0 p(1:3)= 1.0 0.0 1.0 /
 &pot kp=3 type=1 p(1:3)= 60.0 1.5 0.5 /
 &pot kp=5 type=0 p(1:3)= 15.0 0.0 1.25 /
 &pot kp=5 type=1 p(1:3)= {V} 1.25 0.65 /
 &pot /
 &overlap kn1=1 ic1=1 ic2=2 in=1 kind=0 nn=1 l=0 sn=0.5 j=0.5 kbpot=3 be=2.2246 isc=1 ipc=1 /
 &overlap kn1=2 ic1=1 ic2=2 in=2 kind=0 nn=1 l={L} sn=0.5 j={J} kbpot=5 be={BE} isc={ISC} nk={{NK}} er={ER} ipc=3 /
 &overlap /
 &coupling icto=2 icfrom=1 {coupling} /
 &cfp in=1 ib=1 ia=1 kn=1 a=1.0 /
 &cfp in=-2 ib=1 ia=1 kn=2 a=1.0 /
 &coupling /
"""
DN_BINS = {
    "nonresonant": dict(L=1, J=0.5, JF=1.0, PF=1, V=50.0, E=2.0, ER=1.89, ISC=2),
    "narrow": dict(L=2, J=2.5, JF=2.0, PF=-1, V=53.0, E=1.9, ER=0.91, ISC=4),
    "wide": dict(L=0, J=0.5, JF=1.0, PF=-1, V=47.0, E=1.58, ER=1.16, ISC=4),
}
# The refined run refines the steps; NK is doubled in a variant of its own, because FRESCO's k rule (equal weights,
# end points included) converges as 1/NK (transfer.md section 15).
BIN_REF = dict(hcm=0.025, rmatch=60.0, rsp=40.0, rintp=0.05, hnl=0.025, rnl=10.0, centre=0.0, nnu=72, cutl=-3.2,
               nk=1600)
BIN_FINE = dict(hcm=0.0125, rmatch=60.0, rsp=40.0, rintp=0.025, hnl=0.0125, rnl=10.0, centre=0.0, nnu=72, cutl=-6.4,
                nk=1600)
for _name, _b in DN_BINS.items():
    _body = dict(JF=_b["JF"], PF=_b["PF"], EX=12.1274 + _b["E"], V=_b["V"], L=_b["L"], J=_b["J"], BE=-_b["E"],
                 ISC=_b["ISC"], ER=-_b["ER"])
    CASES["dn_bin_%s_fr" % _name] = dict(
        title="15N(d,n)[p+15N bin %s, l=%d, E=%g, ER=%g, ISC=%d] at E_lab = 10 MeV, finite range, post" %
              (_name, _b["L"], _b["E"], _b["ER"], _b["ISC"]),
        elab=10.0, jtmax=30, iter=1, iblock=1, thmin=0.0, body=DN_BIN_BODY.format(coupling=FR_DN, **_body),
        reference=BIN_REF, refined=BIN_FINE,
        also=[dict(BIN_REF, rmatch=50.0), dict(BIN_REF, rnl=12.0, nnu=96), dict(BIN_REF, nk=3200)])
    # FRESCO's zero range (KIND=5, INH=2) to a bin is not consistent with its own finite range (transfer.md
    # section 15): kept as data, not compared by the test.
    CASES["dn_bin_%s_zr" % _name] = dict(
        title="15N(d,n)[p+15N bin %s] at E_lab = 10 MeV, zero range D0 = -130.31 MeV fm^3/2, recoil (inh=2)" % _name,
        elab=10.0, jtmax=30, iter=1, iblock=1, thmin=0.0, extra="inh=2",
        body=DN_BIN_BODY.format(coupling="kind=5 p1=-130.31 p2=0.0", **_body),
        reference={k: v for k, v in BIN_REF.items() if k not in ("rintp", "hnl", "rnl", "centre", "nnu")},
        refined={k: v for k, v in BIN_FINE.items() if k not in ("rintp", "hnl", "rnl", "centre", "nnu")})

# 15N(d,n)16O(g.s.) of DN_BODY in the prior form and with the remnant (transfer.md section 17): the core-core
# potential n + 15N (kp=5, FRESCO's IP3 = KPCORE), the full complex remnant (IP2 = -1); FRESCO takes U_sF and U_aA
# from the partitions' optical potentials and, in the prior form, the residual overlap's binding potential (kp=4,
# Coulomb and spin-orbit included) as V_xA.
DN_REM_BODY = DN_BODY.replace(" &pot /\n &overlap kn1=1", """ &pot kp=5 type=0 p(1:3)= 15.0 0.0 1.3 /
 &pot kp=5 type=1 p(1:6)= 48.0 1.2 0.65 1.0 1.3 0.5 /
 &pot kp=5 type=2 p(1:6)= 0.0 1.3 0.5 6.0 1.3 0.5 /
 &pot /
 &overlap kn1=1""", 1)
# The kernel multiplies u_a (not V_sx u_a): the deuteron's tail needs a long non-local range, and FRESCO does not
# converge in it (rnl 20, 26, 32 below; transfer.md section 17).  Data for the note, not compared by the test.
_DN_FR_REF = dict(hcm=0.025, rmatch=30.0, rintp=0.05, hnl=0.025, rnl=20.0, centre=0.0, nnu=96)
_DN_FR_FINE = dict(hcm=0.025, rmatch=30.0, rintp=0.05, hnl=0.025, rnl=26.0, centre=0.0, nnu=96)
for _name, _coupling in (("dn_post_remnant", "kind=7 ip1=0 ip2=-1 ip3=5"),
                         ("dn_prior_remnant", "kind=7 ip1=1 ip2=-1 ip3=5"),
                         ("dn_prior", "kind=7 ip1=1 ip2=0 ip3=0")):
    CASES[_name] = dict(
        title="15N(d,n)16O(g.s.) at E_lab = 10 MeV, finite range, %s" % _coupling, elab=10.0, jtmax=30, iter=1,
        iblock=1, thmin=0.0, extra="epc=0.00001", body=DN_REM_BODY.format(coupling=_coupling), reference=_DN_FR_REF,
        refined=_DN_FR_FINE, also=[dict(_DN_FR_REF, rnl=32.0)])

# 12C(6Li,d)[alpha + 12C] to a bin inside a wide d-wave resonance of alpha + 12C (kp=5, V = 90 MeV): FRESCO loses
# 14 or more digits in the d + 16O* channel at high L (ERWIN) when the bin reaches far; with the bin cut at
# RSP = 12 fm it is stable to its cutl to 0.5 %.  Data for the note (transfer.md section 15), not compared by the test.
LI_BIN_BODY = LI_BODY.replace("et=0.0 /\n &partition /", "et=10.2720 /\n &partition /", 1).replace(
    "jp=1.0 ptyp=1 ep=0.0 cpot=2 jt=0.0 ptyt=1", "jp=1.0 ptyp=1 ep=0.0 cpot=2 jt=2.0 ptyt=1").replace(
    " &pot kp=4 type=0 p(1:3)= 12.0 4.0 0.9 /\n &pot kp=4 type=1 p(1:3)= 100.0 0.9 0.65 /\n",
    " &pot kp=4 type=0 p(1:3)= 12.0 4.0 0.9 /\n &pot kp=4 type=1 p(1:3)= 90.0 0.9 0.65 /\n").replace(
    "nn=3 l=0 sn=0.0 j=0.0 kbpot=4 be=7.1620 isc=1 ipc=1",
    "nn=1 l=2 sn=0.0 j=2.0 kbpot=4 be=-3.11 isc=4 nk={{NK}} er=-0.4 ipc=3")
CASES["li6d_bin_wide_fr"] = dict(
    title="12C(6Li,d)[a+12C bin, l=2, E=3.11, ER=0.4, ISC=4, RSP=12] at E_lab = 20 MeV, finite range, post (nuclear "
          "part of V_alpha-d)", elab=20.0, jtmax=40, iter=1, iblock=1, thmin=0.0, extra="epc=0.00001", table=True,
    body=LI_BIN_BODY.format(kp5=LI_KP5, krpot="krpot=5", coupling=FR_DN),
    reference=dict(hcm=0.025, rmatch=30.0, rsp=12.0, rintp=0.05, hnl=0.025, rnl=14.0, centre=-0.1, nnu=72, cutl=-6.4,
                   nk=400),
    refined=dict(hcm=0.0125, rmatch=30.0, rsp=12.0, rintp=0.025, hnl=0.0125, rnl=14.0, centre=-0.1, nnu=72,
                 cutl=-12.8, nk=400),
    also=[dict(hcm=0.025, rmatch=30.0, rsp=12.0, rintp=0.05, hnl=0.025, rnl=14.0, centre=-0.1, nnu=72, cutl=-9.6,
               nk=400),
          dict(hcm=0.025, rmatch=30.0, rsp=12.0, rintp=0.05, hnl=0.025, rnl=14.0, centre=-0.1, nnu=72, cutl=-3.2,
               nk=400)])


def coulomb_table(z1z2, rc, step, points):
    """-V_C of a uniformly charged sphere (MeV) on r = 0, step, ...: FRESCO's SHAPE=7 input, scaled by P1 = -1."""
    lines = ["V_C of the alpha-d binding potential (MeV), read with p1 = -1 to cancel it (Rc = %g fm)" % rc,
             "%d %g 0.0" % (points, step)]
    vals = []
    for i in range(points):
        r = i * step
        vals.append("%.12e" % (z1z2 * E2 / (2 * rc) * (3 - (r / rc) ** 2) if r < rc else z1z2 * E2 / r))
    for i in range(0, points, 5):
        lines.append(" ".join(vals[i:i + 5]))
    return "\n".join(lines) + "\n"


def table_for(settings):
    """The -V_C table on the job's own grid (FRESCO reads at most as many points as its grid has)."""
    return coulomb_table(2, 1.9, settings["hcm"], int(round(settings["rmatch"] / settings["hcm"])) + 1)


def job_text(case, settings):
    s = dict(settings)
    cutl = s.pop("cutl", None)
    nk = s.pop("nk", None)
    rsp = s.pop("rsp", None)
    nonlocal_keys = ("rintp", "hnl", "rnl", "centre", "nnu")
    nonlocal_part = " ".join("%s=%s" % (k, s[k]) for k in nonlocal_keys if k in s)
    extra = case.get("extra", "")
    if cutl is not None:
        extra = re.sub(r"cutl=\S+", "", extra) + " cutl=%s" % cutl
    if rsp is not None:
        nonlocal_part = ("rsp=%s " % rsp) + nonlocal_part
    body = case["body"].replace("{NK}", str(nk)) if nk is not None else case["body"]
    return HEAD.format(title=case["title"], hcm=s["hcm"], rmatch=s["rmatch"], nonlocal_part=nonlocal_part,
                       jtmax=case["jtmax"], extra=extra, thmin=case["thmin"], iter=case["iter"],
                       iblock=case["iblock"], unitmass=UNITMASS, finec=FINEC, elab=case["elab"]) + body


def run_fresco(fresco, case, settings):
    """Runs FRESCO in a fresh directory; returns (input text, output text, {fort file: text})."""
    text = job_text(case, settings)
    with tempfile.TemporaryDirectory() as d:
        if case.get("table"):
            with open(os.path.join(d, "fort.4"), "w") as f:
                f.write(table_for(settings))
        with open(os.path.join(d, "job.nin"), "w") as f:
            f.write(text)
        with open(os.path.join(d, "job.nin")) as fin, open(os.path.join(d, "job.out"), "w") as fout:
            subprocess.run([fresco], stdin=fin, stdout=fout, stderr=subprocess.STDOUT, cwd=d, check=True)
        out = open(os.path.join(d, "job.out")).read()
        forts = {}
        for n in ("fort.7", "fort.34", "fort.46"):
            p = os.path.join(d, n)
            forts[n] = open(p).read() if os.path.exists(p) else ""
    return text, out, forts


def extract(out, forts):
    lines = []
    m = re.search(r"2\*amu/hbar\^2 =\s*([0-9.]+).*?e\^2 =\s*([0-9.]+)", out)
    if m:
        lines.append("constants %s %s" % (m.group(1), m.group(2)))
    # Kinematics: the first channel line of each partition (columns ... Jtotal E-cm K eta RM*K ...).
    seen = set()
    for ln in out.splitlines():
        f = ln.split()
        if len(f) > 14 and f[0].isdigit() and f[3] == "#" and f[4].endswith(":") and f[1] not in seen:
            seen.add(f[1])  # by projectile; counted from the right, as an incoming channel has an extra "I"
            lines.append("kin %s %s %s %s" % (f[1], f[-10], f[-9], f[-8]))
    # Potentials (fort.34) at r = 1 ... 10 fm.
    kp_type = None
    for ln in forts.get("fort.34", "").splitlines():
        if ln.startswith("#"):
            mm = re.search(r"KP=\s*(\d+) of type\s*(\d+)", ln)
            kp_type = (mm.group(1), mm.group(2)) if mm else None
            continue
        f = ln.split()
        if kp_type and len(f) == 3:
            r = float(f[0])
            if abs(r - round(r)) < 1e-6 and 1 <= round(r) <= 10:
                lines.append("pot %s %s %s %s %s" % (kp_type[0], kp_type[1], f[0], f[1], f[2]))
    for ln in forts.get("fort.7", "").splitlines():
        f = ln.split()
        if len(f) >= 5 and "S(L,J,JT)" in ln:
            lines.append("S %d %d %d %s %s" % (int(f[2]), round(2 * float(f[3])), round(2 * float(f[4])), f[0], f[1]))
    # Cross sections: the blocks "CROSS SECTIONS FOR OUTGOING"; the first is elastic, a second the transfer.
    blocks = out.split("CROSS SECTIONS FOR OUTGOING")[1:]
    for b, block in enumerate(blocks):
        body = block.split("Integrated")[0] if b == 0 else block.split("Finished all xsecs")[0]
        for mm in re.finditer(r"([0-9.]+) deg\.: X-S =\s*([0-9.E+-]+) mb/sr,\s*\n\+\s+(/R|REAC)?\s*=?\s*([0-9.E+-]+)?",
                              body):
            th = mm.group(1)
            if b == 0 and mm.group(3) == "/R":
                lines.append("ratio %s %s" % (th, mm.group(4)))
        if b == 1:
            for mm in re.finditer(r"([0-9.]+) deg\.: X-S =\s*([0-9.E+-]+) mb/sr", body):
                lines.append("xs %s %s" % (mm.group(1), mm.group(2)))
            for mm in re.finditer(r"([0-9.]+) deg\.: X-S =\s*[0-9.E+-]+ mb/sr,\s*& pols =\s*\S+\s+([-0-9.]+)", body):
                lines.append("pol %s %s" % (mm.group(1), mm.group(2)))
    mm = re.search(r"CUMULATIVE REACTION cross section\s*=\s*([0-9.]+)", out)
    if mm:
        lines.append("reaction %s" % mm.group(1))
    for ln in out.splitlines():
        if " VR " in ln and re.match(r"^\s+\d+:", ln):
            t = ln.split(" VR ")[1].split()
            lines.append("bound %s %s %s %s %s" % (ln.split(":")[0].strip(), t[0].rstrip(";"), t[5], t[6], t[8]))
    for mm in re.finditer(r"ANC in ch\s+(\d+) is\s+([-0-9.E+]+)", out):
        lines.append("anc_printed %s %s" % (mm.group(1), mm.group(2)))
    mm = re.search(r"Bin from\s+([0-9.]+) MeV to\s+([0-9.]+) MeV \(cm\) in\s*(\d+) steps", out)
    if mm:
        lines.append("bin %s %s %s" % mm.groups())
        nk = int(mm.group(3))
        for m2 in re.finditer(r"For #\s+(\d+) 1, K,BE =\s+([0-9.]+)\s+([0-9.]+): S =.*?& Del =\s*([-0-9.]+),", out):
            ik = int(m2.group(1))
            if ik % 20 == 1 or ik == nk:
                lines.append("binphase %s %s %s %s" % m2.groups())
    mm = re.search(r"Norm =\s*([0-9.]+) from =\s*[0-9.]+, RMS =\s*([0-9.]+) & D0 =\s*([-0-9.E+]+)", out)
    if mm:
        lines.append("binnorm %s %s %s" % mm.groups())
        listing = out.split("For channel #  1 with L =")[1].split(" No. ")[0]
        for m2 in re.finditer(r"(\d+)\s+([-0-9.]+E[-+]\d+)\s*([-0-9.]+E[-+]\d+)", listing):
            if (int(m2.group(1)) - 1) % 50 == 0:
                lines.append("binwf %s %s" % (m2.group(1), m2.group(2)))
    kn = 0
    for ln in forts.get("fort.46", "").splitlines():
        if ln.startswith("#"):
            kn += 1
            continue
        f = ln.split()
        if len(f) >= 2:
            r = float(f[0])
            if any(abs(r - x) < 1e-6 for x in (8.0, 10.0, 12.0)):
                lines.append("ancratio %d %s %s" % (kn, f[0], f[1]))
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fresco", required=True, help="the fresco binary (version 3.5)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "fresco"))
    ap.add_argument("--case", action="append", help="only these cases")
    args = ap.parse_args()
    for name, case in CASES.items():
        if args.case and name not in args.case:
            continue
        d = os.path.join(args.out, name)
        os.makedirs(d, exist_ok=True)
        runs = [("reference", case["reference"]), ("refined", case["refined"])]
        runs += [("variant%d" % (i + 1), s) for i, s in enumerate(case.get("also", []))]
        lines = ["# %s" % case["title"], "# FRESCO 3.5 (github.com/I-Thompson/fresco c571c9e), %s, by %s" %
                 (datetime.date.today().isoformat(), os.path.basename(__file__))]
        for tag, settings in runs:
            text, out, forts = run_fresco(args.fresco, case, settings)
            if tag == "reference":
                with open(os.path.join(d, "job.nin"), "w") as f:
                    f.write(text)
                if case.get("table"):
                    with open(os.path.join(d, "fort4.dat"), "w") as f:
                        f.write(table_for(settings))
            lines.append("run %s %s" % (tag, " ".join("%s=%s" % kv for kv in settings.items())))
            lines += [l for l in extract(out, forts) if l.split()[0] not in case.get("drop", ())]
        with open(os.path.join(d, "reference.txt"), "w") as f:
            f.write("\n".join(lines) + "\n")
        print(name, len(lines), "lines")


if __name__ == "__main__":
    main()
