#!/usr/bin/env python3
"""Models for tests/reaction_rate, derived from tests/13N/13N.azr (p + 12C capture).

Every model has external capture off (EC mask 0 on the gamma lines), no data segments and
every <targetInt> line inactive, so mode 3 (extrapolation) and mode 5 (rate) run in a second
or two on the R-matrix part alone.

  p12C        charged entrance, as 13N:  the Gamow tail of the 1/2+ 2.37 MeV level and the
              others (audit F1: the spline of sigma against E lay 0.77 % low at T9 = 0.01).
  n12C        the same levels with a neutron entrance (z1 = 0, m1 = 1.008665): sigma ~ 1/v
              from the s-wave 1/2+ level (audit F2: sigma = 0 below 1 keV lost 49 % at T9 = 0.01).
  n12C_pert   n12C with useAdaptiveGrid=0, the per-temperature quadrature path (same floor).
  moved       p12C plus a narrow 1/2+ level (Gamma_p = 20 eV s-wave, Gamma_gamma = 2 eV, Gamma = 22 eV) at
              E_x = 2.0935 MeV (E_cm = 0.150 MeV).  moved_azr has it at 2.0985 MeV instead
              (+5 keV = 250 widths); its parameter file is used on moved.azr, so the parameters
              being evaluated put the resonance 5 keV from the .azr value (audit F3: the grid
              clusters followed the .azr value and the resonance fell between backbone nodes).

The truth extrapolation segments (<segmentsTest>, mode 3) are written into the *_truth.azr
copies only: decade segments with a step of 1 % of the decade's low end, and for the narrow
level a 1 eV (c.m.) window of +-20 keV around the resonance (4 eV, 5.5 points per width, was
6e-4 off).

    python3 make_models.py            # writes models/*.azr
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "13N", "13N.azr")
OUT = os.path.join(HERE, "models")

M1_P, M2 = 1.00728, 12.0
M1_N = 1.008665
LAB_P = (M1_P + M2) / M2      # E_lab / E_cm for the proton entrance
LAB_N = (M1_N + M2) / M2

NARROW_EX = 2.0935            # E_x of the narrow 1/2+ level (E_cm = 0.1500 MeV)
NARROW_SHIFT = 0.005          # the parameter file moves it up by this much
NARROW_GP = 20.0              # eV (4 % of the Wigner limit for an s-wave proton at 0.15 MeV)
NARROW_GG = 2.0               # eV


def split_blocks(text):
    """Return the text as (head, levels, mid, segData, segTest, targetInt, tail) pieces."""
    def between(tag):
        a = text.index("<%s>" % tag) + len("<%s>" % tag)
        b = text.index("</%s>" % tag)
        return a, b
    return between


def level_lines(text):
    a = text.index("<levels>") + len("<levels>")
    b = text.index("</levels>")
    return text[:a], [l for l in text[a:b].splitlines() if l.strip()], text[b:]


def set_block(text, tag, body_lines):
    a = text.index("<%s>" % tag) + len("<%s>" % tag)
    b = text.index("</%s>" % tag)
    return text[:a] + "\n" + "\n".join(body_lines) + ("\n" if body_lines else "") + text[b:]


def base_model(text):
    head, lines, tail = level_lines(text)
    out = []
    for l in lines:
        t = l.split()
        t[30] = "0"                       # external capture off
        out.append("  ".join(t))
    text = head + "\n" + "\n".join(out) + "\n" + tail
    text = set_block(text, "segmentsData", [])
    text = set_block(text, "segmentsTest", [])
    a = text.index("<targetInt>") + len("<targetInt>")
    b = text.index("</targetInt>")
    ti = [("0" + l[1:]) for l in text[a:b].splitlines() if l.strip()]
    text = set_block(text, "targetInt", ti)
    return text


def neutron(text):
    head, lines, tail = level_lines(text)
    out = []
    for l in lines:
        t = l.split()
        if t[26] == "0" and t[5] == "1":  # the particle pair 1 lines
            t[17] = "%.6f" % M1_N
            t[19] = "0"
            if t[2] == "0":                # the ground-state ANC: 1.62 fm^-1/2 is too large for a
                t[11] = "0.3"              # neutron bound by 1.94 MeV at 3.6 fm (1 - gamma^2 dS/dE < 0)
        if t[26] == "10":                  # the gamma pair: the compound nucleus is now 13C
            t[18] = "13.00335"
            t[20] = "6"
        out.append("  ".join(t))
    return head + "\n" + "\n".join(out) + "\n" + tail


def with_narrow(text, ex):
    """Add a 1/2+ level at E_x = ex: a copy of the 2.3682 MeV 1/2+ lines (s-wave proton, E1 to
    the ground state) with the narrow widths."""
    head, lines, tail = level_lines(text)
    tmpl = [l for l in lines if l.split()[0] == "0.5" and l.split()[1] == "1" and l.split()[2] == "2.3682"]
    assert len(tmpl) == 2, tmpl
    new = []
    for l in tmpl:
        t = l.split()
        t[2] = "%.4f" % ex
        t[8] = "13"                        # a new level ID
        t[3] = "0"                         # energy free (so the parameter file can move it)
        t[11] = "%g" % (NARROW_GP if t[26] == "0" else NARROW_GG)
        new.append("  ".join(t))
    # keep the J-group together: insert after the last 1/2+ line
    idx = max(i for i, l in enumerate(lines) if l.split()[0] == "0.5" and l.split()[1] == "1")
    lines = lines[:idx + 1] + new + lines[idx + 1:]
    return head + "\n" + "\n".join(lines) + "\n" + tail


def truth_segments(lab, e_low, window=None):
    """Decade segments e_low..50 MeV c.m. (step 1 % of the decade's low end) plus a window."""
    segs = []
    lo = e_low
    while lo < 50.0:
        hi = min(lo * 10.0, 50.0)
        step = lo * 0.01 if hi > lo * 5 else (hi - lo) / 900.0
        segs.append("1 1 2 %.10g %.10g %.10g 0 0 0 0 1 0" % (lo * lab, hi * lab, step * lab))
        lo = hi
    if window:
        e0, half, step = window
        segs.append("1 1 2 %.10g %.10g %.10g 0 0 0 0 1 0" % ((e0 - half) * lab, (e0 + half) * lab, step * lab))
    return segs


def set_potential_key(text, key, value):
    return re.sub(r"(?m)^%s=.*$" % re.escape(key), "%s=%s" % (key, value), text)


def main():
    os.makedirs(OUT, exist_ok=True)
    src = open(SRC).read()
    base = base_model(src)
    models = {}
    models["p12C"] = base
    models["n12C"] = neutron(base)
    models["n12C_pert"] = set_potential_key(models["n12C"], "useAdaptiveGrid", "0")
    models["moved"] = with_narrow(base, NARROW_EX)
    models["moved_azr"] = with_narrow(base, NARROW_EX + NARROW_SHIFT)
    for name, text in models.items():
        open(os.path.join(OUT, name + ".azr"), "w").write(text)
    # truth copies with the dense extrapolation segments
    sepE = 1.94351
    truths = {
        "p12C": (models["p12C"], LAB_P, 1e-4, None),
        "n12C": (models["n12C"], LAB_N, 1e-7, None),
        "moved_azr": (models["moved_azr"], LAB_P, 1e-4,
                      (NARROW_EX + NARROW_SHIFT - sepE, 0.020, 1e-6)),
    }
    for name, (text, lab, e_low, window) in truths.items():
        open(os.path.join(OUT, name + "_truth.azr"), "w").write(
            set_block(text, "segmentsTest", truth_segments(lab, e_low, window)))
    print("wrote", ", ".join(sorted(os.listdir(OUT))))


if __name__ == "__main__":
    main()
