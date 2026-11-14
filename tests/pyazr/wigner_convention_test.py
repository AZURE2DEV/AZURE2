#!/usr/bin/env python3
"""AZURE2's Wigner limit is the Teichmann-Wigner sum-rule unit 3 hbar^2 / (2 mu a^2).

`Parameter.wigner_limit` comes from `AChannel::SetWignerLimit`, 1.5 hbar^2 / (mu a^2) in MeV with mu the pair's
reduced mass and a its channel radius; `pyazr.physics.wigner_limit` is the same number (a hbar^2 / (mu a^2) unit is
2/3 of it, and a theta^2 quoted against that unit is 1.5 times larger).

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/wigner_convention_test.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
HBARC = 197.3269804      # MeV fm, Constants.h
UCONV = 931.49410242       # MeV per u, Constants.h

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    from pyazr import Session, physics
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

work = os.path.join(ROOT, "tests", "13N")
m = Session(os.path.join(work, "13N.azr"), cwd=work)
widths = [p for p in m.parameters.of_kind("width") if p.wigner_limit]
pairs_of = {p.number: p for p in m.pairs}
check("the 13N widths carry a Wigner limit", len(widths) >= 3, f"{len(widths)} widths with a limit")
n = 0
for p in widths:
    pair = m.pairs[p.pair - 1]
    if pair.is_photon or pair.channel_radius <= 0:
        continue
    mu = pair.M1 * pair.M2 / (pair.M1 + pair.M2) * UCONV
    expect = 1.5 * HBARC ** 2 / (mu * pair.channel_radius ** 2)
    n += 1
    if n <= 3:
        check(f"{p.name}: wigner_limit = 3 hbar^2 / (2 mu a^2) = {expect:.6f} MeV (pair {pair.key}, a = {pair.channel_radius} fm)",
              abs(p.wigner_limit / expect - 1.0) < 1e-9, f"got {p.wigner_limit}")
    elif abs(p.wigner_limit / expect - 1.0) >= 1e-9:
        check(f"{p.name}: wigner_limit = 3 hbar^2 / (2 mu a^2)", False, f"{p.wigner_limit} vs {expect}")
w = widths[0].wigner_limit
check("pyazr.physics.wigner_limit gives the session's limit", abs(
    float(physics.wigner_limit(pairs_of[widths[0].pair].M1, pairs_of[widths[0].pair].M2,
                               pairs_of[widths[0].pair].channel_radius)) / w - 1.0) < 1e-12)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
