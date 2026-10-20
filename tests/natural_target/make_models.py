#!/usr/bin/env python3
"""Writes the models of the natural-target check (run once; the .azr files are committed).

Two compound nuclei built from the same three-level p+12C skeleton of
tests/convolution_unify/models/fixed.azr:

  A.azr        p + 12C -> 13N (pair key 1), the skeleton as it is
  B.azr        p + 13C -> 14N (pair key 2): the 13C masses and spin, S_p(14N), the levels at
               the same c.m. energies so the resonances sit in the same window, other widths
  M.azr        A's and B's levels in ONE file, <abundances> 0.6 / 0.4 fixed, test segments
                 S1  pair 1 alone          S2  pair 2 alone
                 S3  SUM of (1,1) and the component (2,2), natural 1   = 0.6 S1 + 0.4 S2
                 S4  pair 2 alone, natural 1                             = 0.4 S2
  A_thick.azr, B_thick.azr, M_thick.azr
               the same under a target-integration effect (no beam spread): exact once each
               component's sub-point grid is built with its own pair
  M_samples.azr   M plus a named sample "enr" (0.1 / 0.9, own sum-to-one width): S5 = SUM natural enr
                  = 0.1 S1 + 0.9 S2, S6 = pair 2 natural enr = 0.9 S2
  M_samples_data.azr  two data segments on two samples, four abundances free: the fit must find
                  0.6 / 0.4 and 0.1 / 0.9
  M_data.azr   S3 as a DATA segment reading data.dat (written by check.sh from A and B),
               abundances 0.5 / 0.5 FREE with abundanceSumSigma=0.01: the fit must find 0.6 / 0.4
  bad_level.azr, bad_ratio.azr, bad_key.azr
               a level with channels of both nuclei; "natural 1" on a RATIO; a component naming
               pair key 9 -- each must be refused (ERROR, no output)
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "models")
BASE = os.path.join(HERE, "..", "convolution_unify", "models", "fixed.azr")
TGT = os.path.join(HERE, "..", "convolution_unify", "models", "tgt.azr")
LMAX = 2

# p (1/2+) + 12C (0+): S_p(13N) 1.9435 MeV.  p (1/2+) + 13C (1/2-): S_p(14N) 7.5506 MeV.
PAIR_A = dict(key=1, j1=0.5, pi1=1, j2=0.0, pi2=1, e2=0.0, m1=1.00782, m2=12.0, z1=1, z2=6, sep=1.9435, rad=5.0)
PAIR_B = dict(key=2, j1=0.5, pi1=1, j2=0.5, pi2=-1, e2=0.0, m1=1.00782, m2=13.00335, z1=1, z2=6, sep=7.5506, rad=5.0)
# (J, pi, E_cm, width column): A's are the skeleton's; B's resonances at the same c.m. energies
LEVELS_A = [(0.5, 1, 2.365 - 1.9435, 0.60), (1.5, -1, 3.502 - 1.9435, 0.25), (2.5, 1, 3.547 - 1.9435, 0.20)]
LEVELS_B = [(1.0, 1, 0.70, 0.45), (2.0, -1, 1.25, 0.60), (0.0, 1, 1.85, 0.20)]   # integer J (s = 0, 1); resonances away from 13N's so the two totals differ


def channels(P, j, pi):
    out = []
    s = abs(P["j1"] - P["j2"])
    while s <= P["j1"] + P["j2"] + 1e-9:
        for l in range(LMAX + 1):
            if P["pi1"] * P["pi2"] * (-1) ** l != pi: continue
            if abs(s - l) - 1e-9 <= j <= s + l + 1e-9 and abs((j - abs(s - l)) - round(j - abs(s - l))) < 1e-9:
                out.append((s, l))
        s += 1.0
    return out


def line(P, j, pi, ex, s, l, lid, gamma, efix=0, chfix=0):
    return (f"{j:g} {pi:d} {ex:.6f} {efix} 1 {P['key']} {int(2*s)} {int(2*l)} {lid} 1 {chfix} {gamma:g} "
            f"{P['j1']:g} {P['pi1']} {P['j2']:g} {P['pi2']} {P['e2']:g} {P['m1']} {P['m2']} {P['z1']} {P['z2']} "
            f"{P['sep']} {P['sep']} 0 0 0.0 0 {P['rad']} 0 0 0")


def levels(parts, first_id=1, fixed=False):
    """parts: list of (pair, levels); every level gets every channel its pair allows.
    fixed: energies and widths fixed (the data-mode fit then floats only the abundances)."""
    out, lid = [], first_id
    for P, lv in parts:
        for j, pi, ecm, g in lv:
            chans = channels(P, j, pi)
            assert chans, (P["key"], j, pi)
            for s, l in chans:
                out.append(line(P, j, pi, ecm + P["sep"], s, l, lid, g, efix=int(fixed), chfix=int(fixed)))
            lid += 1
            out.append("")
    return "\n".join(out), lid


def read_block(text, tag):
    a, b = text.index(f"<{tag}>"), text.index(f"</{tag}>")
    return text[a:b + len(tag) + 3]


def model(levels_text, segments_test, abundances=None, config_extra="", target_int="", segments_data=""):
    base = open(BASE).read()
    config = read_block(base, "config")
    if config_extra:
        config = config.replace("</config>", config_extra + "\n</config>")
    tail = base[base.index("<parameterSettings>"):] if "<parameterSettings>" in base else ""
    parts = [config, "<levels>\n" + levels_text + "\n</levels>"]
    if abundances is not None:
        parts.append("<abundances>\n" + "\n".join(abundances) + "\n</abundances>")
    parts.append("<segmentsData>\n" + segments_data + ("\n" if segments_data else "") + "</segmentsData>")
    parts.append("<segmentsTest>\n" + segments_test + "\n</segmentsTest>")
    parts.append("<targetInt>\n" + target_int + ("\n" if target_int else "") + "</targetInt>")
    return "\n".join(parts) + "\n" + tail


def main():
    os.makedirs(OUT, exist_ok=True)
    grid = "1.0 3.0 0.002 0 0 0 0"           # angle-integrated, 1-3 MeV lab, 2 keV (frame-free sums)
    S1 = f"1 1 1 {grid}"
    S2 = f"1 2 2 {grid}"
    S3 = f"1 1 1 {grid} 1 0 1 2 2 -999 natural 1"
    S4 = f"1 2 2 {grid} natural 1"
    # a second, named sample ("enr": 0.1 / 0.9) chosen on the segment line
    S5 = f"1 1 1 {grid} 1 0 1 2 2 -999 sample enr"     # the keyword's name since samples; "natural" stays an alias
    S6 = f"1 2 2 {grid} sample enr"
    lv_a, _ = levels([(PAIR_A, LEVELS_A)])
    lv_b, _ = levels([(PAIR_B, LEVELS_B)])
    lv_m, _ = levels([(PAIR_A, LEVELS_A), (PAIR_B, LEVELS_B)])
    lv_b_from_1, _ = levels([(PAIR_B, LEVELS_B)], first_id=1)
    lv_m_fixed, _ = levels([(PAIR_A, LEVELS_A), (PAIR_B, LEVELS_B)], fixed=True)   # free levels would absorb the weights
    ab = ["7 13 0.6 0 0", "7 14 0.4 0 0"]
    ab_samples = ab + ["enr: 7 13 0.1 0 0", "enr: 7 14 0.9 0 0", "enr: sum 0.01"]
    tgt_line = read_block(open(TGT).read(), "targetInt").splitlines()[1]     # the skeleton's target-integration line
    tgt_m = tgt_line.replace('"1"', '"1-4"', 1)
    files = {
        "A.azr": model(lv_a, S1),
        "B.azr": model(lv_b, S2),
        "M.azr": model(lv_m, "\n".join([S1, S2, S3, S4]), ab),
        "A_thick.azr": model(lv_a, S1, target_int=tgt_line),
        "B_thick.azr": model(lv_b, S2, target_int=tgt_line),
        "M_thick.azr": model(lv_m, "\n".join([S1, S2, S3, S4]), ab, target_int=tgt_m),
        # two samples: the unnamed one (0.6 / 0.4) and "enr" (0.1 / 0.9); S5 = 0.1 S1 + 0.9 S2, S6 = 0.9 S2
        "M_samples.azr": model(lv_m, "\n".join([S1, S2, S3, S4, S5, S6]), ab_samples),
        # data mode with two samples: data.dat (0.6 / 0.4) on the unnamed sample, data_enr.dat (0.1 / 0.9) on "enr",
        # all four abundances free from 0.5, each sample with its own sum-to-one prior
        "M_samples_data.azr": model(lv_m_fixed, "", ["7 13 0.5 1 0", "7 14 0.5 1 0", "enr: 7 13 0.5 1 0", "enr: 7 14 0.5 1 0", "enr: sum 0.01"],
                                    config_extra="abundanceSumSigma=0.01",
                                    segments_data="1 1 1 1.0 3.0 0 0 0 1 0 0 0 0 0 data.dat 1 0 1 2 2 -999 natural 1\n"
                                                  "1 1 1 1.0 3.0 0 0 0 1 0 0 0 0 0 data_enr.dat 1 0 1 2 2 -999 sample enr"),
        # data mode: S3 against data.dat (check.sh writes it), abundances free, sum-to-one prior
        "M_data.azr": model(lv_m_fixed, "", ["7 13 0.5 1 0", "7 14 0.5 1 0"], config_extra="abundanceSumSigma=0.01",
                            segments_data="1 1 1 1.0 3.0 0 0 0 1 0 0 0 0 0 data.dat 1 0 1 2 2 -999 natural 1"),
        # refused: B's levels numbered from 1 again, so level ID 1 has channel lines of both nuclei
        "bad_level.azr": model(lv_a + "\n" + lv_b_from_1, S1, ab),
        # refused: natural on a RATIO composite
        "bad_ratio.azr": model(lv_m, f"1 1 1 {grid} 1 1 1 2 2 -999 natural 1", ab),
        # refused: a component naming a pair key that does not exist
        "bad_key.azr": model(lv_m, f"1 1 1 {grid} 1 0 1 9 9 -999 natural 1", ab),
        # refused: a sample the block does not define
        "bad_sample.azr": model(lv_m, f"1 1 1 {grid} sample nosuch", ab),
    }
    for name, text in files.items():
        with open(os.path.join(OUT, name), "w") as f:
            f.write(text)
        print("wrote", name)


if __name__ == "__main__":
    main()
