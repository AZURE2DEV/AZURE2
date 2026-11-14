#!/usr/bin/env python3
"""pyazr.sources.masses (the AME2020 table) and masses.mass_check / fill_masses on a pyazr.Project.

1. The table reads every nuclide, and the reference rows below, copied from the AME2020 file, come back exactly
   (mass excess in keV, atomic mass in u, the '#' estimated flag).  Mass excess and atomic mass agree to 2e-9 u
   for every measured nuclide, which checks the fixed-width columns are cut where the format says.
2. Names: shorthands, "12C" / "C12" / "c-12", (Z, A); "h" and unknown symbols are refused.
3. Separation energies from mass excesses equal (m1 + m2 - M) times the engine's u (931.49410242 MeV) to 5 eV,
   so a pair filled from the table passes the engine's 1 keV load-time check by five orders of magnitude.
4. masses.fill_masses on tests/13N/13N.azr (a bare-proton 1.00728 u next to atomic masses, S_p 20 eV high):
   the changes are exactly M1 and S of pair 1 and M2 of the gamma pair; nothing else of the pairs and levels
   changes, and every row carries the new S (pair 1's rows, and the entrance-S copy of every row); separation=False keeps S; pairs= restricts; an unknown key is refused.
   A mass-number model (every mass rounded) is filled to the same numbers.

Run from anywhere:  python3 tests/pyazr/masses_test.py
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    import pyazr  # noqa: E402
    from pyazr import Project  # noqa: E402
    from pyazr.sources import masses  # noqa: E402
except ImportError as e:  # as the other pyazr tests: a missing package dependency -> nothing to test
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(77)

FAIL = []


def check(ok, what):
    print(("  ok    " if ok else "  FAIL  ") + what)
    if not ok:
        FAIL.append(what)


# (name, Z, A, mass excess keV, atomic mass u, estimated), from pyazr/data/mass_1.mas20
REFERENCE = [
    ("n", 0, 1, 8071.31806, 1.00866491590, False),
    ("p", 1, 1, 7288.971064, 1.007825031898, False),
    ("d", 1, 2, 13135.722895, 2.014101777844, False),
    ("t", 1, 3, 14949.81090, 3.01604928132, False),
    ("3He", 2, 3, 14931.21888, 3.01602932197, False),
    ("a", 2, 4, 2424.91587, 4.00260325413, False),
    ("12C", 6, 12, 0.0, 12.0, False),
    ("13C", 6, 13, 3125.00933, 13.00335483534, False),
    ("13N", 7, 13, 5345.481, 13.005738609, False),
    ("16O", 8, 16, -4737.00217, 15.99491461926, False),
    ("3Li", 3, 3, 28667.0, 3.030775, True),
]
U_ENGINE = 931.49410242   # include/Constants.h


def table():
    t = masses._load()
    check(len(t) >= 3500, f"table has {len(t)} nuclides (AME2020: 3558)")
    for name, Z, A, dm, mass, est in REFERENCE:
        e = masses.lookup(name)
        check((e.Z, e.A) == (Z, A) and e.mass_excess == dm and e.atomic_mass == mass and e.estimated == est,
              f"{name}: Z={e.Z} A={e.A} excess {e.mass_excess} keV, mass {e.atomic_mass!r} u, estimated {e.estimated}")
    worst = max(abs(e.A + e.mass_excess / 1000 / masses.U_MEV - e.atomic_mass)
                for e in t.values() if not e.estimated)
    check(worst < 2e-9, f"mass excess and atomic mass agree for every measured nuclide (worst {worst:.2e} u)")
    check(masses.AME_VERSION == "AME2020", "AME_VERSION names the evaluation")
    # pyazr's element symbols are the engine's one table (physics/Elements.h, _azure2.element_symbols()); the
    # AME file names every element it lists the same way.
    from pyazr.system.nuclei import SYMBOLS
    odd = sorted({(e.Z, e.symbol) for e in t.values() if e.Z > 0 and SYMBOLS[e.Z] != e.symbol})
    check(len(SYMBOLS) == 119 and SYMBOLS[0] == "n" and not odd,
          f"SYMBOLS (n..Og, {len(SYMBOLS)}) agree with the AME symbols{'; differ: ' + str(odd[:5]) if odd else ''}")


def names():
    cases = {"n": (0, 1), "p": (1, 1), "d": (1, 2), "t": (1, 3), "a": (2, 4), "alpha": (2, 4), "helion": (2, 3),
             "3He": (2, 3), "He3": (2, 3), "12C": (6, 12), "C12": (6, 12), "c-12": (6, 12), "12-C": (6, 12),
             "15N": (7, 15), "N15": (7, 15), "1n": (0, 1), "P31": (15, 31), (8, 16): (8, 16)}
    bad = [k for k, v in cases.items() if masses.parse_nuclide(k) != v]
    check(not bad, f"shorthands, isotope names and (Z, A) parse ({len(cases)} cases){'; wrong: ' + str(bad) if bad else ''}")
    for spec in ("h", "Xx12", "12", ""):
        try:
            masses.parse_nuclide(spec)
            check(False, f"{spec!r} refused")
        except ValueError:
            check(True, f"{spec!r} refused")
    try:
        masses.lookup("12Fe")
        check(False, "a nuclide the table lacks raises KeyError")
    except KeyError:
        check(True, "a nuclide the table lacks raises KeyError")


def separation():
    pairs = [("p", "12C", 1.943490064), ("a", "12C", 7.16191804), ("n", "12C", 4.94630873)]
    for light, heavy, s_ref in pairs:
        s = masses.separation_energy(light, heavy)
        check(s == s_ref, f"S({light} + {heavy}) = {s!r} MeV from the mass excesses")
    worst = 0.0
    for light, heavy in [("p", "12C"), ("a", "12C"), ("n", "13C"), ("a", "13C"), ("p", "7Li"), ("n", "7Be"),
                         ("a", "16O"), ("p", "18O"), ("a", "22Ne"), ("d", "6Li"), ("t", "4He"), ("3He", "4He")]:
        pm = masses.pair_masses(light, heavy)
        comp = masses.lookup(pm.compound)
        s_engine = (pm.M1 + pm.M2 - comp.atomic_mass) * U_ENGINE
        worst = max(worst, abs(pm.S - s_engine))
    check(worst < 5e-6, f"S from mass excesses = (m1 + m2 - M) x 931.49410242 MeV within 5 eV (worst {worst * 1e6:.2f} eV)")


def fill():
    src = os.path.join(ROOT, "tests", "13N", "13N.azr")
    work = tempfile.mkdtemp(prefix="pyazr_masses_")
    before = Project.open(src)
    rows = masses.mass_check(before)
    check([r.pair for r in rows] == [1, 2] and abs(rows[0].dS_kev - 0.0199) < 1e-3 and rows[1].S_table is None,
          f"mass_check on 13N.azr: pair 1 S {rows[0].dS_kev:+.4f} keV from the table, gamma pair not compared")

    m = Project.open(src)
    changes = masses.fill_masses(m)
    want = [(1, "M1", 1.00728, 1.007825031898), (1, "S", 1.94351, 1.943490064), (2, "M2", 13.0057, 13.005738609)]
    check(changes == want, f"fill_masses changes {changes}")
    out = m.save(os.path.join(work, "filled.azr"))
    after = Project.open(out)
    check(all(abs(r.dS_kev) < 1e-9 for r in masses.mass_check(after) if r.S_table is not None)
          and all(r.mass_offset == 0.0 for r in masses.mass_check(after)),
          "the written file agrees with the table exactly")
    # nothing else of the pairs and the levels moved
    import dataclasses
    unchanged = lambda q: dataclasses.replace(q, light_mass=0.0, heavy_mass=0.0, separation_energy=0.0, threshold=0.0)
    check([unchanged(q) for q in before.pairs()] == [unchanged(q) for q in after.pairs()]
          and [(lv.label, [(c.pair, c.l, c.s, c.width_text, c.fixed) for c in lv.channels]) for lv in before.levels()]
          == [(lv.label, [(c.pair, c.l, c.s, c.width_text, c.fixed) for c in lv.channels]) for lv in after.levels()],
          "only the masses and S of the pairs changed")
    rows_after = [after.text.split("\n")[c.line - 1] for lv in after.levels() for c in lv.channels]
    check(all("1.943490064" in r.split() and "1.94351" not in r.split() for r in rows_after),
          "every row carries pair 1's new S (its own column on pair 1's rows, the entrance-S copy on all)")

    m = Project.open(src)
    ch = masses.fill_masses(m, separation=False)
    check(all(c[1] != "S" for c in ch) and masses.mass_check(m)[0].S == 1.94351, "separation=False keeps S")
    m = Project.open(src)
    ch = masses.fill_masses(m, pairs=[2])
    check(ch == [(2, "M2", 13.0057, 13.005738609)], f"pairs=[2] touches pair 2 only ({ch})")
    try:
        masses.fill_masses(Project.open(src), pairs=[7])
        check(False, "an unknown pair key is refused")
    except KeyError:
        check(True, "an unknown pair key is refused")

    d = rows[0].describe()
    check("<-- mass" not in d, "describe: a bare-proton mass (5e-4 u off) is listed, not flagged")
    m = Project.open(src)
    m.set_pair(1, light_mass=1.08664)   # the 17O+a typo for 1.00866
    d = masses.mass_check(m)[0].describe()
    check("<-- mass" in d and abs(masses.mass_check(m)[0].dS_kev - 0.0199) < 1e-3,
          f"describe: a mistyped mass (1.08664 u) is flagged, S still compared through the identified 1H ({d})")

    # a mass-number model fills to the same numbers
    m = Project.open(src)
    for q in m.pairs():
        m.set_pair(q.key, heavy_mass=float(round(q.heavy_mass)),
                   light_mass=float(round(q.light_mass)) if q.light_mass > 0.5 else None)
    masses.fill_masses(m)
    same = [(r.M1, r.M2, r.S) for r in masses.mass_check(m)] == [(r.M1, r.M2, r.S) for r in masses.mass_check(after)]
    check(same, "a model with mass numbers as masses fills to the same masses and S")


def main():
    print("== table");      table()
    print("== names");      names()
    print("== separation"); separation()
    print("== fill");       fill()
    if FAIL:
        print(f"{len(FAIL)} check(s) FAILED")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
