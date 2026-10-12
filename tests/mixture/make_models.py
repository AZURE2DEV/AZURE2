"""Writes the models of tests/mixture/check.sh (run once from this directory; the files are committed).

All derive from tests/convolution_unify/models/fixed.azr (p+12C elastic at 120 deg lab with a 30 keV lab
Gaussian beam spread) with the three 13N levels given their observed widths (see below; the parent's are
sub-keV once the transform reads them as eV) and on the grid 0.35-1.2 MeV lab in 2 keV steps, which holds
the strong 1/2+ s-wave resonance at E_p = 0.46 MeV (Gamma = 32 keV), so the resolution kernels matter at
the 10 % level there.
The mixture block refers to tables that check.sh builds at run time from the model's own bare extrapolation
(bare_wide), so a mixture of "half this model + half a table of this model" must reproduce the plain folded
model.

  gauss        the reference: the Gaussian beam spread on this grid
  identity     gauss + "mixture 1.0 0": the block with no effect
  bare_wide    no experimental effect, 0.2-1.8 MeV lab in 0.1 keV steps (the 20 ns kernel reaches 1.45 MeV
               from the 1.2 MeV end of the grid): the source of table.dat
  mix_gauss    gauss + "mixture 0.5 1 table.dat 0.5"
  mix_two      gauss + two tables of 0.25 each (the same file twice)
  mix_scaled   gauss + table2.dat (= 2 x table.dat) at 0.25
  plain_udr    gauss with the Gaussian replaced by a synthetic time-of-flight kernel (udr.txt): a Gaussian
               delay distribution of 20 ns at a 20 m flight path, i.e. ~4 % in energy (80 keV at 2 MeV,
               clearly distinct from the 30 keV Gaussian)
  mix_udr      plain_udr + "mixture 0.5 1 table.dat 0.5"
  mix_data     the fixed model as a DATA segment (c.m. differential at the c.m. angle of 120 deg lab,
               data.dat written by check.sh from the plain folded values, 1 % errors) with the Gaussian and
               the mixture block: chi2 must vanish and the .out fit column equal mix_gauss's

udr.txt follows tests/reference/udr_reference_test.cpp: blocks of "<E eV>   # eV" then "tau[us] R" rows.
"""
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
base = open(os.path.join(HERE, "..", "convolution_unify", "models", "fixed.azr")).read()
# The parent's widths (0.6, 0.25, 0.2) are read as eV once the parameter transform is on: sub-keV
# resonances that no table could resolve and no kernel affects.  Give the three 13N levels their
# observed widths (1/2+ 2.365: 31.7 keV; 3/2- 3.502: 62 keV; 5/2+ 3.547: 47 keV) so the folds matter.
for ex, w in (("2.365", "31700"), ("3.502", "62000"), ("3.547", "47000")):
    line = [l for l in base.splitlines() if (" " + ex + " ") in l and l.count(" ") > 20]
    assert len(line) == 1, (ex, line)
    tok = line[0].split(); tok[11] = w
    base = base.replace(line[0], " ".join(tok))
TEST = "1 1 1 0.35 1.2 0.002 120 120 1 1"
GAUSS = '1 "1" 50 1 0.030 0 0 "" 0 0 0 0 "" 0 0 0.04 20 50'
UDR = '1 "1" 200 0 0 0 0 "" 0 0 0 0 "" 0 0 0.04 20 50 udr "udr.txt" 20 0 0 1'
assert GAUSS in base


def with_test(text, line):
    return re.sub(r"<segmentsTest>.*?</segmentsTest>", "<segmentsTest>\n" + line + "\n</segmentsTest>", text, flags=re.S)


def with_effect(text, line):
    return re.sub(r"<targetInt>.*?</targetInt>", "<targetInt>\n" + (line + "\n" if line else "") + "</targetInt>", text, flags=re.S)


models = {
    "gauss":     with_test(base, TEST),
    "identity":  with_test(base, TEST + " mixture 1.0 0"),
    "bare_wide": with_effect(with_test(base, "1 1 1 0.2 1.8 0.0001 120 120 1 1"), None),
    "mix_gauss": with_test(base, TEST + ' mixture 0.5 1 "table.dat" 0.5'),
    "mix_two":   with_test(base, TEST + ' mixture 0.5 2 "table.dat" 0.25 "table.dat" 0.25'),
    "mix_scaled": with_test(base, TEST + ' mixture 0.5 1 "table2.dat" 0.25'),
    "plain_udr": with_effect(with_test(base, TEST), UDR),
    "mix_udr":   with_effect(with_test(base, TEST + ' mixture 0.5 1 "table.dat" 0.5'), UDR),
    # data segment: isDiff 4 (c.m. differential), angle window around the c.m. angle, norm 1 fixed
    "mix_data":  with_test(re.sub(r"<segmentsData>\n</segmentsData>",
                                  '<segmentsData>\n1 1 1 0.3 1.3 0 180 4 1 0 0 0 0 0 data.dat 0 0 mixture 0.5 1 "table.dat" 0.5\n</segmentsData>',
                                  base), ""),
}
os.makedirs(os.path.join(HERE, "models"), exist_ok=True)
for name, text in models.items():
    with open(os.path.join(HERE, "models", name + ".azr"), "w") as f:
        f.write(text)

# Synthetic UDR kernel: a Gaussian delay distribution, sigma = 0.020 us, at four energies (same shape).
with open(os.path.join(HERE, "models", "udr.txt"), "w") as f:
    f.write("# synthetic Gaussian time-of-flight delay distribution, sigma 20 ns\n")
    f.write("---------------------------------------------------\n")
    s = 0.020
    for e_ev in (5.0e5, 1.0e6, 2.0e6, 4.0e6):
        f.write(f"{e_ev:.1f}   # eV\n")
        for i in range(0, 801):
            t = (-6.0 + 12.0 * i / 800.0) * s
            f.write(f"{t:.6e} {math.exp(-0.5 * t * t / (s * s)):.8e}\n")
        f.write("\n")
print("written", len(models), "models + udr.txt")
