#!/usr/bin/env python3
"""The model of case (c), made from models/p12C.azr when the test runs (a project file kept in tests/ would be a
project of every tree comparison and of tests/upgrade/predictions.txt).

    python3 closed_model.py models/p12C.azr          >  p12C_closed.azr
    python3 closed_model.py models/p12C.azr --truth  >  p12C_closed_truth.azr

p12C_closed is p12C with a third pair, p + 12C in a fictitious state at 48 MeV (an ANC of 0.5 fm^-1/2 on the 1/2+
background level at 20 MeV, s wave), which opens at E_cm = 48 MeV, between the last two backbone nodes of the rate's
grid (47.66 and 50 MeV): one point of the grid has sigma > 0 and the rate takes its per-temperature quadrature.
--truth adds a dense extrapolation of pair 3: E_cm from 48 to 50 MeV in steps of 0.125 keV (16000 steps), written as
lab energies (E_lab = E_cm (m1 + m2) / m2) in records of 800 steps; truth.py integrates its mode-3 output.
"""
import sys

M1, M2 = 1.00728, 12.0
E_THRESHOLD, E_TOP, STEPS, PER_RECORD = 48.0, 50.0, 16000, 800
# the channel of pair 3 on level 4 (1/2+, E_x = 20 MeV): the columns of <levels> (docs/source/reference/data_formats.rst)
PAIR3 = "0.5  1  20  1  1  3  1  0  4  1  1  0.5  0.5  1  0  1  48.0  1.00728  12  1  6  1.94351  1.94351  0  0  0.0  0  3.6  0  0  0"


def closed(text):
    out = []
    for line in text.split("\n"):
        out.append(line)
        t = line.split()
        if len(t) == 31 and t[:3] == ["0.5", "1", "20"] and t[5] == "2" and t[8] == "4":
            out.append(PAIR3)
    if PAIR3 not in out:
        sys.exit("closed_model.py: level 4's pair-2 line not found in the model")
    return "\n".join(out)


def with_truth(text):
    energies = [E_THRESHOLD + (E_TOP - E_THRESHOLD) * i / STEPS for i in range(STEPS + 1)]
    records = []
    for k in range(0, STEPS, PER_RECORD):
        lab = " ".join("%.12g" % (e * (M1 + M2) / M2) for e in energies[k:k + PER_RECORD + 1])
        records.append(f"extrapolation x{k // PER_RECORD + 1}\n    reaction       12C(p,p)12C[48]\n"
                       f"    quantity       sigma\n    unit           b\n    energies       list {lab} MeV lab\n")
    return text.replace("<extrapolations>\n", "<extrapolations>\n" + "\n".join(records) + "\n", 1)


def main():
    text = closed(open(sys.argv[1]).read())
    sys.stdout.write(with_truth(text) if "--truth" in sys.argv[2:] else text)


if __name__ == "__main__":
    main()
