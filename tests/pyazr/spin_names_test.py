#!/usr/bin/env python3
"""pyazr.spin_operator: the operator of a spin observable for its name in the literature (pyazr/spin.py).

- every name of the table, and Aij / Cij for the spins that decide them;
- an operator is returned unchanged;
- Aij where the spins do not decide, and an unknown name, raise pyazr.SpinOperatorError (a pyazr.Error and a
  ValueError) whose message holds its three parts;
- every operator the table gives is one the engine reads: a format-3 project with a spin record of that operator
  validates (p + 3He for spin 1/2 on spin 1/2, d + 4He for a spin-1 beam), and one that the particles' spins do not
  allow (b:yy for a proton beam) is refused by the engine, not by the table.

Run from anywhere:  python3 tests/pyazr/spin_names_test.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

import pyazr  # noqa: E402
from pyazr.spin import LITERATURE_NAMES, spin_operator  # noqa: E402

failures = []


def check(what, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {what}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(what)


expected = {"Ay": "b:y", "Ay0": "b:y", "A0y": "t:y", "Py": "e1:y", "Kyy": "b:y>e1:y", "iT11": "b:iT11",
            "T20": "b:T20", "T21": "b:T21", "T22": "b:T22"}
check("the table of names", LITERATURE_NAMES == expected and pyazr.spin_operator is spin_operator, str(LITERATURE_NAMES))
got = {n: spin_operator(n) for n in expected}
check("each name without spins", got == expected, str(got))
got = {n: spin_operator(n, 1.0, 0.0) for n in expected}
check("each name whatever the spins", got == expected, str(got))

axes = "xyz"
tensor = {f"A{i}{j}": spin_operator(f"A{i}{j}", 1.0, 0.5) for i in axes for j in axes}
check("Aij of a spin-1 beam is the tensor analysing power b:ij",
      all(v == "b:" + k[1:] for k, v in tensor.items()), str(tensor))
for kind in "AC":
    corr = {f"{kind}{i}{j}": spin_operator(f"{kind}{i}{j}", 0.5, 0.5) for i in axes for j in axes}
    check(f"{kind}ij of spin 1/2 on spin 1/2 is the correlation b:i*t:j",
          all(v == f"b:{k[1]}*t:{k[2]}" for k, v in corr.items()), str(corr))
check("Cij without spins is the correlation", spin_operator("Cxz") == "b:x*t:z", spin_operator("Cxz"))
check("an operator is returned unchanged", spin_operator("b:y*t:y") == "b:y*t:y" and spin_operator(" b:T20 ") == "b:T20")


def refused(what, call, *parts):
    try:
        value = call()
    except pyazr.SpinOperatorError as err:
        text = str(err)
        ok = (isinstance(err, pyazr.Error) and isinstance(err, ValueError) and all(p in text for p in parts)
              and all(part and part in text for part in (err.what, err.allowed, err.fix)))
        check(what, ok, text)
    else:
        check(what, False, f"returned {value!r}")


refused("Ayy without spins is refused", lambda: spin_operator("Ayy"), "b:yy", "b:y*t:y", "spin_operator")
refused("Ayy for spin 1/2 on spin 1 is refused", lambda: spin_operator("Ayy", 0.5, 1.0), "b:yy", "b:y*t:y")
refused("Cyy for a spin-1 beam is refused", lambda: spin_operator("Cyy", 1.0, 0.5), "b:yy", "b:y*t:y")
refused("Cyy with the beam's spin alone is refused", lambda: spin_operator("Cyy", 0.5), "target spin None")
refused("an unknown name is refused with the table", lambda: spin_operator("Azzz"), "Kyy", "Aij and Cij")
refused("a lower-case name is refused", lambda: spin_operator("ay"), "Ay0")
check("the error class is described", "SpinOperatorError" in pyazr.kinds()["errors"])

# The engine reads every operator the table gives: the first spin record of a test project, rewritten.
import re  # noqa: E402
import shutil  # noqa: E402
import tempfile  # noqa: E402

work = tempfile.mkdtemp(prefix="spin_names_")
try:
    def errors_with(directory, operator):
        src = os.path.join(ROOT, "tests", directory)
        dst = os.path.join(work, directory)
        if not os.path.isdir(dst):
            shutil.copytree(src, dst)
        path = os.path.join(dst, directory + ".azr")
        text = open(os.path.join(src, directory + ".azr")).read()
        new, n = re.subn(r"(\n\s*operator\s+)\S+", lambda m: m.group(1) + operator, text, count=1)
        assert n == 1, "no operator line in " + directory
        with open(path, "w") as f:
            f.write(new)
        return [d for d in pyazr.validate(path) if d["severity"] == "error"]

    half = [spin_operator(n, 0.5, 0.5) for n in ("Ay", "Ay0", "A0y", "Py", "Kyy")]
    half += [spin_operator(f"{k}{i}{j}", 0.5, 0.5) for k in "AC" for i in axes for j in axes]
    bad = {op: errors_with("spin_observables_p3He", op) for op in sorted(set(half))}
    check(f"p + 3He: the engine reads the {len(bad)} operators of the names for spin 1/2 on spin 1/2",
          not any(bad.values()), str({k: v[0]["message"] for k, v in bad.items() if v}))
    one = [spin_operator(n, 1.0, 1.0) for n in ("Ay", "A0y", "iT11", "T20", "T21", "T22")]
    one += [spin_operator(f"A{i}{j}", 1.0, 1.0) for i in axes for j in axes]
    bad = {op: errors_with("spin_observables_dd", op) for op in sorted(set(one))}
    check(f"d + d: the engine reads the {len(bad)} operators of the names for a spin-1 beam",
          not any(bad.values()), str({k: v[0]["message"] for k, v in bad.items() if v}))
    wrong = errors_with("spin_observables_p3He", "b:yy")
    check("b:yy for a proton beam is refused by the engine (the check bites)",
          len(wrong) >= 1 and any(d.get("field") == "operator" or "operator" in d["message"] for d in wrong),
          str(wrong[:1]))
finally:
    shutil.rmtree(work, ignore_errors=True)

print("FAILED: %d" % len(failures) if failures else "PASS: spin observable names")
sys.exit(1 if failures else 0)
