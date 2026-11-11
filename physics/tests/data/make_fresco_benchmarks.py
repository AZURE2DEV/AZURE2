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

Lines of reference.txt (whitespace-separated; '#' starts a comment):
  run <name> <settings>                 the FRESCO settings of the following lines (reference or refined)
  constants <2u/hbar^2> <e^2>          as FRESCO prints them
  kin <projectile> <E_cm> <k> <eta>     per partition of the incoming energy (MeV, fm^-1)
  pot <kp> <type> <r> <Re> <Im>         FRESCO's listing of the potential components (fort.34; the Coulomb term per
                                        unit Z1 Z2, the spin-orbit term before 2 l.s), at r = 1 ... 10 fm
  S <l> <2j> <2J> <Re> <Im>             elastic S-matrix elements (fort.7)
  ratio <theta> <sigma/sigma_R>         elastic, c.m. (deg)
  xs <theta> <mb/sr>                    the transfer cross section, c.m. (deg)
  reaction <mb>                         FRESCO's cumulative reaction cross section
  bound <kn> <scale> <rms> <D0> <ANC>   FRESCO's bound-state listing: depth scale factor (VR), rms radius (fm),
                                        D0 (MeV fm^3/2), ANC (fm^-1/2, at R = rmatch - h)
  ancratio <kn> <r> <value>             u(r)/W(2 kappa r) from fort.46 (the ANC plateau)
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
    nonlocal_keys = ("rintp", "hnl", "rnl", "centre", "nnu")
    nonlocal_part = " ".join("%s=%s" % (k, s[k]) for k in nonlocal_keys if k in s)
    extra = case.get("extra", "")
    if cutl is not None:
        extra = re.sub(r"cutl=\S+", "", extra) + " cutl=%s" % cutl
    return HEAD.format(title=case["title"], hcm=s["hcm"], rmatch=s["rmatch"], nonlocal_part=nonlocal_part,
                       jtmax=case["jtmax"], extra=extra, thmin=case["thmin"], iter=case["iter"],
                       iblock=case["iblock"], unitmass=UNITMASS, finec=FINEC, elab=case["elab"]) + case["body"]


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
    mm = re.search(r"CUMULATIVE REACTION cross section\s*=\s*([0-9.]+)", out)
    if mm:
        lines.append("reaction %s" % mm.group(1))
    for ln in out.splitlines():
        if " VR " in ln and re.match(r"^\s+\d+:", ln):
            t = ln.split(" VR ")[1].split()
            lines.append("bound %s %s %s %s %s" % (ln.split(":")[0].strip(), t[0].rstrip(";"), t[5], t[6], t[8]))
    for mm in re.finditer(r"ANC in ch\s+(\d+) is\s+([-0-9.E+]+)", out):
        lines.append("anc_printed %s %s" % (mm.group(1), mm.group(2)))
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
            lines += extract(out, forts)
        with open(os.path.join(d, "reference.txt"), "w") as f:
            f.write("\n".join(lines) + "\n")
        print(name, len(lines), "lines")


if __name__ == "__main__":
    main()
