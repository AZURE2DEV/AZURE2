#!/usr/bin/env python3
"""Composite segments: a zero RATIO denominator, the component scaling factor, and the route pyazr runs.

1. RATIO with a zero denominator (ESegment::CalculateTheoreticalCrossSection).  tests/ratio_frame/models/ratio_data.azr
   (2->3 over a fixed-angle denominator) with the denominator moved to exit pair 5 and every pair-5 width set to 0, so
   the denominator cross section is exactly 0.  Before 2026-11-01 each such point was 0, which a fit compared with the
   datum as if it were a calculated ratio; now it is NaN and a warning is printed once.  Control: the same model with
   the pair-5 widths kept is finite everywhere.

2. Component scaling (SegLine / ExtrapLine, ReadComponentList).  tests/sum_lab_frame/models/sum_data.azr (2->3 plus
   2->5) with the component written in the scaling format, "-1 1 2 5 -999 <s>".  The SUM is linear in s, so
   v(s) - v(0) must equal s * (v(1) - v(0)).  Before 2026-11-01 s went through std::to_string ("%f", six decimals):
   s = 1.5e-7 was read as 0 and 1.2345678e-4 as 1.23e-4 (0.4 % off).  A scaling that is not a number ("abc") was
   read as 1.0 without a word; now the project is refused with a message naming the segment and component.

3. Route.  tests/13N with the R-matrix route asked for (<config> first line false) and ignore_externals=True (the
   pyazr default): the R-matrix route does not implement ignore-externals, so pyazr now switches to the A-matrix with a
   warning (as the CLI does, and as the GUI's dialogs require).  Before, it ran the R-matrix route and the capture
   model was up to 4.5 % off the A-matrix route's.  With ignore_externals=False the R-matrix route runs and equals
   the A-matrix route.  With use_park=True (and ignore_externals=False) the A-matrix is used too: RMatrixFunc has no
   Park matrices and read the Park amplitudes as Brune's (up to a factor 9 off).

Each case runs in its own process (the warnings are printed once per process).
Run from anywhere:  python3 tests/pyazr/composite_terms_test.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
ENERGIES = ["5.95", "6.05", "6.15", "6.25", "6.35"]


def edit(src, dst, fn):
    lines = open(src).read().split("\n")
    open(dst, "w").write("\n".join(fn(lines)))


def worker(azr, cwd, flags, out):
    import numpy as np
    from pyazr import azure2
    try:
        with azure2(azr, cwd=cwd, ignore_externals=flags[0] == "1", use_park=flags[1:2] == "1") as a:
            v = a.calculate_rwa(np.array(a.params_rwa, float))
            res = {"values": [[None if not np.isfinite(x) else float(x) for x in np.ravel(np.asarray(s, float))]
                              for s in v]}
    except Exception as err:
        res = {"error": f"{type(err).__name__}: {err}"}
    json.dump(res, open(out, "w"))


def run(azr, cwd, ignore_externals=True, park=False):
    out = os.path.join(cwd, "res.json")
    if os.path.exists(out):
        os.remove(out)
    flags = ("1" if ignore_externals else "0") + ("1" if park else "0")
    p = subprocess.run([sys.executable, os.path.abspath(__file__), "worker", azr, cwd, flags, out], env=dict(os.environ, OMP_NUM_THREADS="2"), capture_output=True, text=True, timeout=900)
    res = json.load(open(out)) if os.path.exists(out) else {"error": f"exit {p.returncode}: {p.stderr[-600:]}"}
    return res, p.stdout + p.stderr


def project(suite, model, datafile=None):
    work = tempfile.mkdtemp(prefix="pyazr_composite_")
    src = os.path.join(ROOT, "tests", suite, "models") if model else os.path.join(ROOT, "tests", suite)
    name = model or suite
    shutil.copy(os.path.join(src, name + ".azr"), work)
    if os.path.isdir(os.path.join(src, "data")):
        shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
    else:
        os.makedirs(os.path.join(work, "data"))
    if datafile:
        with open(os.path.join(work, "data", datafile), "w") as f:
            for e in ENERGIES:
                f.write(f"{e}  90  1.0  0.01\n")
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    return work, os.path.join(work, name + ".azr")


def main():
    try:
        import numpy as np
        from pyazr import azure2  # noqa: F401
    except Exception as err:  # same skip as the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 0
    failures = []

    def check(name, ok, detail=""):
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
        if not ok:
            failures.append(name)

    works = []
    try:
        # 1. RATIO with a zero denominator
        work, base = project("ratio_frame", "ratio_data", "ratio.dat")
        works.append(work)

        def ratio_lines(zero_widths):
            def fn(lines):
                out, block = [], None
                for l in lines:
                    s = l.split()
                    if l.strip() in ("<levels>", "<segmentsData>"):
                        block = l.strip()
                    elif l.strip() in ("</levels>", "</segmentsData>"):
                        block = None
                    elif block == "<levels>" and len(s) > 12 and s[5] == "5" and zero_widths:
                        s[11] = "0.0"
                        l = "  ".join(s)
                    elif block == "<segmentsData>" and s:
                        assert s[-6:] == ["1", "1", "1", "2", "2", "60"], s
                        l = "  ".join(s[:-6] + ["1", "1", "1", "2", "5", "60"])
                    out.append(l)
                return out
            return fn
        for zero in (True, False):
            azr = os.path.join(work, f"ratio_{int(zero)}.azr")
            edit(base, azr, ratio_lines(zero))
            res, log = run(azr, work)
            if "error" in res:
                check(f"RATIO, pair-5 widths {'0' if zero else 'kept'}: ran", False, res["error"])
                continue
            v = res["values"][0]
            if zero:
                check("RATIO over a zero denominator: every point NaN", all(x is None for x in v), f"{v}")
                check("RATIO over a zero denominator: warning printed",
                      "WARNING: the RATIO denominator of segment" in log and "is 0 at E_lab" in log, log[-600:])
            else:
                check("RATIO control (denominator widths kept): every point finite and non-zero",
                      all(x is not None and x != 0.0 for x in v), f"{v}")

        # 2. Component scaling
        work, base = project("sum_lab_frame", "sum_data", "sum.dat")
        works.append(work)

        def with_scaling(s):
            def fn(lines):
                out = []
                for l in lines:
                    t = l.split()
                    if t[-6:] == ["1", "0", "1", "2", "5", "-999"]:
                        l = "  ".join(t[:-6] + ["1", "0", "-1", "1", "2", "5", "-999", s])
                    out.append(l)
                return out
            return fn
        vals = {}
        for s in ("0", "1", "1.5e-7", "1.2345678e-4"):
            azr = os.path.join(work, f"sum_{s}.azr")
            edit(base, azr, with_scaling(s))
            res, log = run(azr, work)
            vals[s] = np.array(res["values"][0], float) if "error" not in res else None
            if vals[s] is None:
                check(f"SUM, scaling {s}: ran", False, res["error"])
        if all(v is not None for v in vals.values()):
            d1 = vals["1"] - vals["0"]
            check("SUM: the component contributes (v(1) != v(0))", np.all(np.abs(d1) > 0), f"{d1}")
            for s in ("1.5e-7", "1.2345678e-4"):
                want = float(s) * d1
                got = vals[s] - vals["0"]
                rel = np.max(np.abs(got - want) / np.abs(want))
                # v(s) - v(0) is a difference of numbers ~1/s larger: allow for their rounding
                tol = 1e-9 * np.max(np.abs(vals["0"]) / np.abs(want))
                check(f"SUM, scaling {s}: v(s) - v(0) = s (v(1) - v(0)) (worst {rel:.1e}, tol {tol:.0e})", rel <= tol,
                      f"{got} vs {want}")
        azr = os.path.join(work, "sum_bad.azr")
        edit(base, azr, with_scaling("abc"))
        res, log = run(azr, work)
        check("SUM, scaling 'abc': the project is refused", "error" in res, f"{res}")
        check("SUM, scaling 'abc': the engine names segment and component",
              "Data segment #1: component 1 of 1: the scaling factor is not a number" in log, log[-600:])

        # 2b. The two-number component layout of old GUI files ("entrance exit" per component, no angle), read as the
        # GUI reads it.  Before 2026-11-01 the engine took the next component's entrance key as an angle and dropped
        # the last component ("2 2 5 2 3" gave 2->5 at a fixed 2 deg only); the archive has 86 such lines (11B+a ELBE
        # and 13C+a Cierjacks/Fowler neutron totals, 2025).  Too few numbers: a warning, read as a plain segment.
        def with_tail(tail):
            def fn(lines):
                out = []
                for l in lines:
                    t = l.split()
                    if t[-6:] == ["1", "0", "1", "2", "5", "-999"]:
                        l = "  ".join(t[:-6] + tail.split())
                    out.append(l)
                return out
            return fn
        layouts = {}
        for name, tail in (("three", "1 0 2 2 5 -999 2 3 -999"), ("two", "1 0 2 2 5 2 3"),
                           ("short", "1 0 2 2 5"), ("plain", "0")):
            azr = os.path.join(work, f"sum_{name}.azr")
            edit(base, azr, with_tail(tail))
            res, log = run(azr, work)
            layouts[name] = (np.array(res["values"][0], float) if "error" not in res else None, log)
            if layouts[name][0] is None:
                check(f"SUM, {name}-number layout: ran", False, res["error"])
        if all(v[0] is not None for v in layouts.values()):
            d = np.max(np.abs(layouts["two"][0] - layouts["three"][0]) / np.abs(layouts["three"][0]))
            check(f"SUM, two-number layout (2 5, 2 3) = three-number layout (worst {d:.1e})", d <= 1e-12)
            d = np.max(np.abs(layouts["short"][0] - layouts["plain"][0]) / np.abs(layouts["plain"][0]))
            check(f"SUM, too few numbers: read as the plain segment (worst {d:.1e})", d <= 1e-12)
            check("SUM, too few numbers: warning names the segment",
                  "WARNING: Data segment #1: the component list declares 2 component(s) but holds only 2 number(s)"
                  in layouts["short"][1], layouts["short"][1][-600:])

        # 3. Route: R-matrix asked for, with and without ignore_externals
        work, base = project("13N", None)
        works.append(work)

        def route(amatrix):
            def fn(lines):
                c = lines.index("<config>") + 1
                lines[c] = ("true" if amatrix else "false") + "  " + lines[c].split(None, 1)[1]
                return lines
            return fn
        for ign, park in ((True, False), (False, False), (False, True)):
            r = {}
            logs = {}
            for am in (True, False):
                azr = os.path.join(work, f"route_{int(am)}{int(ign)}{int(park)}.azr")
                edit(base, azr, route(am))
                res, logs[am] = run(azr, work, ign, park)
                r[am] = np.concatenate([np.array(x, float) for x in res["values"]]) if "error" not in res else None
                if r[am] is None:
                    check(f"13N route {'A' if am else 'R'}, ignore_externals={ign}: ran", False, res["error"])
            if r[True] is None or r[False] is None:
                continue
            d = np.max(np.abs(r[False] - r[True]) / np.abs(r[True]))
            warned = "WARNING: the R-matrix route (<config> first line false) does not implement" in logs[False]
            if park:
                check("13N, R-matrix asked, use_park=True: warning printed (Park)",
                      warned and "Park's parametrization" in logs[False], logs[False][-600:])
                check(f"13N, R-matrix asked, use_park=True: = A-matrix route (worst {d:.1e})", d <= 1e-12)
            elif ign:
                check("13N, R-matrix asked, ignore_externals=True: warning printed", warned, logs[False][-600:])
                check(f"13N, R-matrix asked, ignore_externals=True: = A-matrix route (worst {d:.1e})", d <= 1e-12)
            else:
                check("13N, R-matrix asked, ignore_externals=False: no warning", not warned)
                check(f"13N, R-matrix route, ignore_externals=False: = A-matrix route (worst {d:.1e})", d <= 1e-10)
    finally:
        for w in works:
            shutil.rmtree(w, ignore_errors=True)
    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(*sys.argv[2:6])
        sys.exit(0)
    sys.exit(main())
