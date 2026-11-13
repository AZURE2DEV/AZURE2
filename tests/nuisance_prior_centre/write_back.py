#!/usr/bin/env python3
"""tests/nuisance_prior_centre/check.sh (g): the write-back of a target's fitted values (the engine's
parameter_value_edits, which pyazr's save_fit and Project use): a screening potential written in keV (the engine
holds eV), a Doppler temperature, an areal density and a sample thickness (the engine holds atoms/b, the record
atoms/b here).  Only the value token changes, in the unit written; the prior stays as written.

Usage: write_back.py WORK   (WORK = check.sh's directory with the runs scr, dop, thk, t2)"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from pyazr.described import _engine  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def target_of(text, needle):
    """The name of the target record holding the line that contains `needle`."""
    name = None
    for line in text.splitlines():
        m = re.match(r"^target\s+(\S+)", line)
        if m:
            name = m.group(1)
        if needle in line:
            return name
    raise ValueError(needle)


def edited(work, run, project, text, values):
    edits, missing, ok = _engine().parameter_value_edits(text.encode(), os.path.join(work, run), values)
    lines = text.splitlines()
    for line, col, old, new in sorted(edits, key=lambda e: (e[0], -e[1])):
        s = lines[line - 1]
        lines[line - 1] = s[:col - 1] + new.decode() + s[col - 1 + len(old):]
    return lines, missing


def main(work):
    cases = [  # run, project, the line before, its parameter's suffix, the value (engine unit), the line written
        ("scr", "grad.azr", ("screening      1000 eV fit prior 1200 200", "screening      1 keV fit prior 1.2 0.2"),
         ".screening", 1100.0, "screening      1.1 keV fit prior 1.2 0.2"),
        ("dop", "dop.azr", None, ".doppler", 310.5, "doppler        free_gas 310.5 K fit prior 320 20"),
        ("t2", "T2.azr", None, ".areal_density", 2.1e17,
         "areal_density  2.1e+17 atoms/cm2 fit prior 2.2e+17 10% grid_max 3e+17"),
        ("thk", "conv.azr", None, ".areal_density", 0.155, "areal_density  0.155 atoms/b fit prior 0.16 0.01"),
    ]
    for run, project, swap, suffix, value, want in cases:
        text = open(os.path.join(work, run, project)).read()
        if swap:
            text = text.replace(*swap)
        needle = want.split()[0] + " "
        needle = {"screening ": "screening      1 keV", "doppler ": "free_gas 300 K fit",
                  "areal_density ": "areal_density  2e+17 atoms/cm2 fit" if run == "t2" else "atoms/b fit prior"}[needle]
        name = target_of(text, needle) + suffix
        lines, missing = edited(work, run, project, text, {name: value})
        got = [s.strip() for s in lines if s.strip().startswith(want.split()[0]) and "prior" in s]
        check(f"{name} = {value:g} written as '{want}'", want in got, (got, missing))
    if failures:
        print(f"{len(failures)} failure(s)")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
