#!/usr/bin/env python3
"""Models for the capture analyzing-power coupling-order test (audit item 9, 2026-10-26).

twolevel: n + "15N" (1/2+ on 1/2-: channel spins 0 and 1) into a fictitious compound nucleus with a
0+ ground state, two 1+ levels that decay by M1 to it:  level A (E_cm 1.0 MeV) with a neutron width
in channel spin 0 only, level B (E_cm 1.3 MeV) in channel spin 1 only, both l = 1.  The two levels
share no particle channel, so the level matrix is diagonal and the two capture amplitudes are
T_s = gamma_s gamma_gamma / D_lambda(E) with D = E_lambda - E - gamma^2 [S(E) - S(E_lambda)] - i gamma^2 P(E)
(the Brune form), which verify.py evaluates from the engine's own reduced-width amplitudes
(parameters.out) and the closed-form neutral l = 1 shift and penetrability.  Between the levels
the two amplitudes have different phases, so A_y is not zero and its sign is the coupling order.
The gamma widths are tiny (1 meV) so the photon channel cannot matter anywhere else.

deuteron: the same with a spin-1 projectile and an analyzing-power test segment -- refused.

Run from this directory: python3 make_models.py  (writes models/).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MN, MT, MD = 1.00866, 15.0001, 2.01410
S_N = 2.491
RADIUS = 4.5
EA, EB = S_N + 1.0, S_N + 1.3          # E_x of the two 1+ levels
GA0, GB1, GGAM = 30000.0, -20000.0, 0.001   # eV: neutron widths (s = 0 of A, s = 1 of B), gamma widths

CONFIG = """<config>
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
"""
# pair columns: j1 pi1 j2 pi2 exE m1 m2 z1 z2 sepE(ref) sepE ... pType chRad ...
NPAIR = f"0.5  1  0.5  -1  0  {MN}  {MT}  0  7  {S_N}  {S_N}  0  0  0.0  0  {RADIUS}  0  0  0"
DPAIR = f"1  1  0.5  -1  0  {MD}  {MT}  0  7  {S_N}  {S_N}  0  0  0.0  0  {RADIUS}  0  0  0"
# the last column is the external-capture multipole mask: 0, no external capture
GPAIR = f"1  1  0  1  0  0  16.0061  0  7  {S_N}  0  0  0  0.0  10  0  0  0  0"
DGPAIR = f"1  1  0  1  0  0  17.0061  0  7  {S_N}  0  0  0  0.0  10  0  0  0  0"


def channel(J, pi, E, key, s2, l2, width, tail, lid):
    """One <levels> line: J pi E fixE aa key 2s 2l levelID active fixWidth width <pair columns>."""
    return f"{J}  {pi}  {E:.6g}  1  1  {key}  {s2}  {l2}  {lid}  1  1  {width}  {tail}"


def write(name, levels, tests):
    os.makedirs(os.path.join(HERE, "models"), exist_ok=True)
    txt = (CONFIG + "<levels>\n" + "\n\n".join("\n".join(b) for b in levels) + "\n\n</levels>\n<segmentsData>\n"
           + "</segmentsData>\n<segmentsTest>\n" + "".join(t + "\n" for t in tests)
           + "</segmentsTest>\n<targetInt>\n</targetInt>\n")
    open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)


def elab(ecm, m_p, m_t):
    return ecm * (m_p + m_t) / m_t


def main():
    lvA = [channel(1, 1, EA, 1, 0, 2, GA0, NPAIR, 1), channel(1, 1, EA, 1, 2, 2, 0.0, NPAIR, 1),
           channel(1, 1, EA, 2, 0, 2, GGAM, GPAIR, 1)]
    lvB = [channel(1, 1, EB, 1, 0, 2, 0.0, NPAIR, 2), channel(1, 1, EB, 1, 2, 2, GB1, NPAIR, 2),
           channel(1, 1, EB, 2, 0, 2, GGAM, GPAIR, 2)]
    e0, e1, de = elab(1.05, MN, MT), elab(1.25, MN, MT) + 1e-9, elab(0.05, MN, MT)
    tests = [f"1  1  2  {e0:.6f}  {e1:.6f}  {de:.6f}  0  180  15  7",   # capture A_y
             f"1  1  2  {e0:.6f}  {e1:.6f}  {de:.6f}  0  180  15  5"]   # the c.m. differential capture cross section
    write("twolevel", [lvA, lvB], tests)

    lvD = [channel(1, 1, EA, 1, 1, 2, GA0, DPAIR, 1), channel(1, 1, EA, 2, 0, 2, GGAM, DGPAIR, 1)]
    e0d, e1d = elab(1.05, MD, MT), elab(1.25, MD, MT) + 1e-9
    write("deuteron", [lvD], [f"1  1  2  {e0d:.6f}  {e1d:.6f}  {elab(0.05, MD, MT):.6f}  0  180  15  7"])


if __name__ == "__main__":
    main()
