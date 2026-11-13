#!/usr/bin/env python3
"""The truth model of p12C_closed: the model with a dense extrapolation of its exit pair 3 (p + 12C[48 MeV]).

    python3 closed_truth.py models/p12C_closed.azr  >  p12C_closed_truth.azr

E_cm from the threshold, 48 MeV, to the rate's upper end, 50 MeV, in steps of 0.125 keV (16000 steps), written as lab
energies (E_lab = E_cm (m1 + m2) / m2) in records of 800 steps.  Mode 3 of the truth model gives sigma(E) without
ReactionRate.cpp; truth.py integrates it.
"""
import sys

M1, M2 = 1.00728, 12.0
E_THRESHOLD, E_TOP, STEPS, PER_RECORD = 48.0, 50.0, 16000, 800


def main():
    text = open(sys.argv[1]).read()
    energies = [E_THRESHOLD + (E_TOP - E_THRESHOLD) * i / STEPS for i in range(STEPS + 1)]
    records = []
    for k in range(0, STEPS, PER_RECORD):
        chunk = energies[k:k + PER_RECORD + 1]
        lab = " ".join("%.12g" % (e * (M1 + M2) / M2) for e in chunk)
        records.append(f"extrapolation x{k // PER_RECORD + 1}\n    reaction       12C(p,p)12C[48]\n"
                       f"    quantity       sigma\n    unit           b\n    energies       list {lab} MeV lab\n")
    sys.stdout.write(text.replace("<extrapolations>\n", "<extrapolations>\n" + "\n".join(records) + "\n", 1))


if __name__ == "__main__":
    main()
