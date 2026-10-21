#!/usr/bin/env python3
"""Models of tests/electron_screening, built from tests/13N (12C(p,g)13N).

Only the low-energy capture segments stay active: 6 (LUNA HPGe), 7 (LUNA BGO) and 8 (FELS); 13N's own
<targetInt> lines are switched off.  Written next to this script:

  bare.azr        no screening
  scr.azr         screening 1000 eV, fixed, on segments 6 and 7 (8 stays bare)
  conv.azr        the same inside a Gaussian beam spread of 1e-6 MeV (sub-points)
  grad.azr        U_e free with a 200 eV prior about 1000 eV, a beam spread of 2 keV
  bad.azr         a screening block without its sigma (must be refused)

The fit models (fit.azr, fit_prior.azr) need data generated from scr.azr; check.sh writes them.

Run from anywhere: python3 tests/electron_screening/make_models.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
from pyazr.azrfile import AzrModel  # noqa: E402

SRC = os.path.join(ROOT, "tests", "13N", "13N.azr")
KEEP = (6, 7, 8)


def base():
    """13N with only KEEP active and its effects off (data/ is tests/13N/data, linked in by check.sh)."""
    text = open(SRC).read()
    out, block = [], None
    n = 0
    for line in text.splitlines():
        s = line.strip()
        if s in ("<segmentsData>", "<targetInt>", "<segmentsTest>"):
            block = s
        elif s.startswith("</"):
            block = None
        elif block == "<segmentsData>" and s:
            n += 1
            tok = line.split()
            tok[0] = "1" if n in KEEP else "0"
            line = "  ".join(tok)
        elif block == "<targetInt>" and s:
            line = "0" + line.lstrip()[1:]
        elif block == "<segmentsTest>" and s:
            continue
        out.append(line)
    path = os.path.join(HERE, "_base.azr")
    open(path, "w").write("\n".join(out) + "\n")
    m = AzrModel.from_file(path)
    os.remove(path)
    return m


def write(m, name):
    m.write(os.path.join(HERE, name))
    print("wrote", name)


if __name__ == "__main__":
    write(base(), "bare.azr")
    write(base().add_target_effect("6,7", screening=1000.0), "scr.azr")
    write(base().add_target_effect("6,7", gaussian_sigma=1e-6, screening=1000.0), "conv.azr")
    write(base().add_target_effect("6,7", gaussian_sigma=0.002, screening=(1000.0, True, 200.0)), "grad.azr")
    bad = base().add_target_effect("6,7", screening=1000.0)
    bad._suffix = bad._suffix.replace("screening  1000  0  0", "screening  1000  0")
    assert "screening  1000  0\n" in bad._suffix + "\n" or bad._suffix.count("screening") == 1
    write(bad, "bad.azr")
