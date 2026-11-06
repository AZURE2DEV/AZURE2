#!/usr/bin/env python3
"""Park formalism: dS/dE at the level energies against an independent calculation.

Park's overlap J = 1 - sum_c gamma_c^2 dS_c/dE (Phys. Rev. C 104, 064612, Eq. 28)
needs the energy derivative of the shift function at each level energy.  Here
dS_c/dE is read back from the engine as J(gamma_c = 0) - J(gamma_c = 1)
(`park_norms`) for alpha + 12C levels of 16O (a = 5.43 fm) and compared with
mpmath (30 digits; S = rho (F F' + G G')/(F^2 + G^2) above threshold,
S = z W'(z)/W(z), z = 2 kappa a, below it), with and without the Coulomb
table.

Why: until dev 441d14f (2026-10-19) dS/dE was a 1 eV central difference
(gsl_deriv_central) of the exact Coulomb functions, whose rounding noise over
2 eV made dS/dE wrong by +19 % (l = 2) and +54 % (l = 3) at 32.8 MeV above the
alpha threshold, where dS/dE is ~2e-3 MeV^-1; with the gamma^2 ~ 6.6 MeV of a
background pole that moved J by 0.010 (12C+a r3d, Park basis: Marvin elastic
chi2 +27).  Bound levels were off by ~1e-3 (a numerical radial derivative of
the Whittaker function).  The engine now agrees with mpmath to ~2e-6 (table)
and ~2e-5 (exact Coulomb functions).  Both sessions run in one process, so the
check also covers the process-wide derivative memo, which until 2026-10-24 was
keyed without the table (or relativistic kinematics) and handed the second
session the first one's numbers.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/park_dsde_test.py
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
TOL = 1e-4

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)

# (J, parity, E_x, l) -> dS/dE (MeV^-1), mpmath at 30 digits (AZURE2's constants,
# non-relativistic kinematics), alpha + 12C, m = 4.002603254 + 12 u, S_alpha = 7.16192 MeV.
TRUTH = {
    (0.0, 1, 0.0, 0): 0.3263630922827727,          # bound by 7.16 MeV
    (1.0, -1, 7.11685, 1): 0.6114617675897355,     # bound by 45 keV
    (3.0, -1, 11.5, 3): 0.6109096865432108,        # 4.34 MeV above threshold
    (2.0, 1, 40.0, 2): 0.001796701316563697,       # 32.8 MeV above: background poles
    (3.0, -1, 40.0, 3): 0.002883013402134494,
    (5.0, -1, 40.0, 5): 0.007521880384631723,
}

PAIR = "0 1 0 1 0 4.002603254 12 2 6 7.16192 7.16192 0 0 0.0 0 5.43 0 0 0"
LEVELS = [  # J pi E_x (fixed), l, width (eV; the ANC for a bound level), level id
    (0, 1, "0.0", 0, "1.0", 1),
    (1, -1, "7.11685", 1, "1.0", 2),
    (3, -1, "11.5", 3, "100000.0", 3),
    (2, 1, "40.0", 2, "1000000.0", 4),
    (3, -1, "40.0", 3, "1000000.0", 5),
    (5, -1, "40.0", 5, "1000000.0", 6),
]

MODEL = """<config>
true #Perform A-Matrix Calculation
output/ #Full Path to Output Directory
checks/ #Full Path to Checks Directory
none #Compond Nucleus Check
none #Boundary Condition Check
none #Data Check
none #Lo-Matrix and Penetrability Check
none #Legendre Polynomial Check
none #Coulomb Amplitudes Check
none #Reaction Pathway Check
none #Angular Distributions Check
</config>
<levels>
{levels}
</levels>
<segmentsData>
1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs.dat 0 0
</segmentsData>
<segmentsTest>
</segmentsTest>
<targetInt>
</targetInt>
<lastRun>
9997
"output/param.sav"
""
0 0
0 ""
-1. -1. -1.
</lastRun>
"""


def stage(tmp, tag):
    work = os.path.join(tmp, tag)
    for d in ("output", "checks", "data"):
        os.makedirs(os.path.join(work, d))
    # .azr channel lines carry 2s and 2l
    lines = "\n\n".join(f"{J} {pi} {E} 1 1 1 0 {2 * l} {lid} 1 0 {w} {PAIR}" for J, pi, E, l, w, lid in LEVELS)
    with open(os.path.join(work, "park_dsde.azr"), "w") as f:
        f.write(MODEL.format(levels=lines))
    with open(os.path.join(work, "data", "xs.dat"), "w") as f:
        f.write("2.0\t90\t1.0\t0.1\n2.5\t90\t1.0\t0.1\n3.0\t90\t1.0\t0.1\n")
    return work


found = {}
with tempfile.TemporaryDirectory() as tmp:
    for table in (True, False):
        tag = "table" if table else "exact"
        print(f"\nCoulomb table {'on' if table else 'off'}")
        work = stage(tmp, tag)
        M = azure2(os.path.join(work, "park_dsde.azr"), cwd=work, use_brune=False, use_park=True,
                   use_coulomb_table=table)
        P = M.parameters
        levels = [(p.jgroup, p.level, p.J, p.parity, p.value) for p in P if p.kind == "energy"]  # park_norms order
        widths = [p for p in P if p.kind == "width" and p.radiation_type == "P"]
        check(f"{tag}: every width free", all(w.free_index is not None for w in widths))
        x0 = np.asarray(M.params_rwa, float)
        for w in widths:
            x0[w.free_index] = 0.0
        seen = 0
        for li, (jg, lv, J, pi, E) in enumerate(levels):
            for w in (w for w in widths if (w.jgroup, w.level) == (jg, lv)):
                key = (float(J), int(pi), round(float(E), 6), int(w.L))
                if key not in TRUTH:
                    continue
                x = x0.copy()
                j0 = M.park_norms(x)[li]
                x[w.free_index] = 1.0
                d = j0 - M.park_norms(x)[li]
                t = TRUTH[key]
                r = (d - t) / t
                seen += 1
                found[(tag, key)] = d
                check(f"{tag}: J = {J:g}{'+' if pi > 0 else '-'}, E_x = {E:g} MeV, l = {w.L}: dS/dE {d:.10e} "
                      f"(mpmath {t:.10e}, {r:+.1e})", abs(r) < TOL, f"relative error {r:.2e} > {TOL:g}")
        check(f"{tag}: all {len(TRUTH)} levels compared", seen == len(TRUTH), f"compared {seen}")

# The two routes round differently, so identical numbers mean the second session read the first one's
# derivatives (the process-wide memo of CNuc::CalcShiftFunctions was keyed without the table until 2026-10-24).
# Only meaningful where a Coulomb table exists: without one (a fresh CI runner) both sessions take the exact route
# and must agree bit for bit.
table_file = os.environ.get("AZURE2_COULOMB_TABLE") or os.path.expanduser("~/.azure2/coulomb_table.bin")
if os.path.isfile(table_file):
    check("table and exact sessions in one process computed their own dS/dE",
          any(found.get(("table", k)) != found.get(("exact", k)) for k in TRUTH), "bit-identical in every level")
else:
    print(f"  skip  table vs exact in one process: no Coulomb table at {table_file}, both sessions are exact")

print()
if failures:
    print(f"FAILED: {len(failures)} check(s)")
    sys.exit(1)
print("park_dsde: all checks passed")
