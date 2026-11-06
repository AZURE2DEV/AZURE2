#!/usr/bin/env python3
"""Models for the lab-frame SUM test (audit item 7, 2026-10-26).

13C + alpha in 17O with two neutron exit pairs of the same masses and different Q values,
n + 16O(6.05, 0+) (key 3) and n + 16O(6.13, 3-) (key 5): the lab-to-c.m. Jacobians of the two
differ at the same lab angle, so a lab-frame differential SUM of the two is the place where the
per-term frame factor (ESegment::TermFrameFactor) matters.  Two levels, 5/2+ and 3/2+, with
widths to all three pairs; no external capture.

  extrap     : test segments (mode 3), lab (isDiff 1) and c.m. (isDiff 5):
               T1 a -> n0+, T2 a -> n3-, T3 = SUM(T1 + T2) lab, T4 = SUM of the same pair twice,
               T5/T6/T7 the c.m. versions of T1/T2/T3, T8 = T1 + 0.5 T2 (scaling token)
  plain_data : two plain lab data segments (a -> n0+, a -> n3-) with a fixed 0.05 MeV energy shift
  sum_data   : the SUM data segment of the two, same shift; check.sh writes data/sum.dat from the
               plain segments' fitted lab values, so its chi2 must be ~0
  sum_free   : sum_data with the widths free, for the AZURE_GRAD_CHECK gradient test

Run from this directory: python3 make_models.py  (writes models/ and models/data/).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MA, M13 = 4.0026, 13.0033
S_ALPHA = 6.35869

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
ALPHA = "0  1  0.5  -1  0  4.0026  13.0033  2  6  4.14313  6.35869  0  0  0.0  0  5.4  0  0  0"
N_0P = "0.5  1  0  1  6.0494  1.00866  15.9949  0  8  4.14313  4.14313  0  0  0.0  0  4.4  0  0  0"
N_3M = "0.5  1  3  -1  6.12989  1.00866  15.9949  0  8  4.14313  4.14313  0  0  0.0  0  4.4  0  0  0"


def channel(J, E, key, s, l, width, tail, free=False, lid=1):
    """One <levels> line: J pi E fixE aa key 2s 2l levelID active fixWidth width <pair columns>."""
    fix = 0 if free else 1
    return f"{J}  1  {E}  1  1  {key}  {int(round(2 * s))}  {2 * l}  {lid}  1  {fix}  {width}  {tail}"


def levels(free):
    lv1 = [channel(2.5, 11.0, 2, 0.5, 3, 3000.0, ALPHA, free), channel(2.5, 11.0, 3, 0.5, 2, 15000.0, N_0P, free),
           channel(2.5, 11.0, 5, 2.5, 1, 20000.0, N_3M, free), channel(2.5, 11.0, 5, 3.5, 1, -10000.0, N_3M, free)]
    lv2 = [channel(1.5, 11.1, 2, 0.5, 1, 5000.0, ALPHA, free, 2), channel(1.5, 11.1, 3, 0.5, 2, -12000.0, N_0P, free, 2),
           channel(1.5, 11.1, 5, 2.5, 1, -15000.0, N_3M, free, 2)]
    return [lv1, lv2]


def write(name, lv, data_segments, test_segments):
    os.makedirs(os.path.join(HERE, "models", "data"), exist_ok=True)
    txt = (CONFIG + "<levels>\n" + "\n\n".join("\n".join(b) for b in lv) + "\n\n</levels>\n<segmentsData>\n"
           + "".join(s + "\n" for s in data_segments) + "</segmentsData>\n<segmentsTest>\n"
           + "".join(s + "\n" for s in test_segments) + "</segmentsTest>\n<targetInt>\n</targetInt>\n")
    open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)


def main():
    # lab energies 5.95-6.30 MeV: E_x 10.9-11.2, both neutron pairs open (thresholds E_x 10.19 and 10.27)
    grid = "5.95  6.30  0.05  0  180  15"
    tests = [f"1  2  3  {grid}  1",                       # T1 a -> n0+, lab
             f"1  2  5  {grid}  1",                       # T2 a -> n3-, lab
             f"1  2  3  {grid}  1  1  0  1  2  5  -999",  # T3 SUM lab
             f"1  2  3  {grid}  1  1  0  1  2  3  -999",  # T4 SUM of the same pair twice
             f"1  2  3  {grid}  5",                       # T5 c.m.
             f"1  2  5  {grid}  5",                       # T6 c.m.
             f"1  2  3  {grid}  5  1  0  1  2  5  -999",  # T7 SUM c.m.
             f"1  2  3  {grid}  1  1  0  -1  1  2  5  -999  0.5"]  # T8 T1 + 0.5 T2, lab
    write("extrap", levels(False), [], tests)

    energies = [5.95 + 0.07 * i for i in range(6)]
    angles = [15.0 + 30.0 * i for i in range(6)]
    rows = [(e, a, 1.0, 0.01) for e in energies for a in angles]
    for fname in ("d1.dat", "d2.dat"):
        with open(os.path.join(HERE, "models", "data", fname), "w") as f:
            for r in rows:
                f.write("  ".join(f"{v:.10g}" for v in r) + "\n")
    seg = "5.9  6.4  0  180  1  1.0  0  0  0.05  0  0"
    write("plain_data", levels(False), [f"1  2  3  {seg}  data/d1.dat", f"1  2  5  {seg}  data/d2.dat"], [])
    write("sum_data", levels(False), [f"1  2  3  {seg}  data/sum.dat  1  0  1  2  5  -999"], [])
    write("sum_free", levels(True), [f"1  2  3  {seg}  data/sum.dat  1  0  1  2  5  -999"], [])


if __name__ == "__main__":
    main()
