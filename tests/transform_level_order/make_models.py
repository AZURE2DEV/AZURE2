#!/usr/bin/env python3
"""Models for tests/transform_level_order (audit item 13, 2026-10-27).

Derived from tests/capture_ay_coupling/models/twolevel.azr (n + 15N, one 1+ group with two levels,
every width fixed):

  asc.azr       levels in energy order: 3.491 (energy free), 3.791 (energy fixed)
  desc.azr      the same two levels with the 3.791 level listed FIRST
  desc_fit.azr  desc with a total-capture data segment on data/cap.dat (written by check.sh from desc
                with the 3.491 level's formal energy lowered by 0.02 MeV)

Run from anywhere:  python3 tests/transform_level_order/make_models.py
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "capture_ay_coupling", "models", "twolevel.azr"))
OUT = os.path.join(HERE, "models")

TEST = "1 1 2 0.96 1.49 0.0106 0 0 1 4\n"          # total capture, E_cm 0.90-1.40 across both levels
DATA = "1 1 2 0.9 1.6 0 180 3 1.0 0 0 0 0 0 data/cap.dat 0 0\n"


def set_col(block, col, value):
    out = []
    for line in block.splitlines():
        t = line.split()
        t[col] = value
        out.append("  ".join(t))
    return "\n".join(out)


def main():
    os.makedirs(OUT, exist_ok=True)
    s = open(SRC).read()
    m = re.search(r"<levels>\n(.*?)</levels>", s, re.S)
    low, high = m.group(1).strip("\n").split("\n\n")
    low, high = set_col(low, 3, "0"), set_col(high, 3, "1")      # column 4: fix energy
    head, tail = s[:m.start()], s[m.end():]
    tail = re.sub(r"<segmentsTest>\n.*?</segmentsTest>", "<segmentsTest>\n" + TEST + "</segmentsTest>", tail, flags=re.S)

    def write(name, blocks, data=""):
        t = tail.replace("<segmentsData>\n</segmentsData>", "<segmentsData>\n" + data + "</segmentsData>")
        with open(os.path.join(OUT, name), "w") as f:
            f.write(head + "<levels>\n" + "\n\n".join(blocks) + "\n\n</levels>" + t)

    write("asc.azr", [low, high])
    write("desc.azr", [high, low])
    write("desc_fit.azr", [high, low], DATA)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
