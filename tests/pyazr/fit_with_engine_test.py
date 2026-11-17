#!/usr/bin/env python3
"""Session.fit_with_engine against the program's fit, fitter by fitter (docs/upgrade/fit.md).

Every case runs the program (AZURE2_BIN, one thread) on a copy of a project and fit_with_engine on another copy, from
the same start, and compares what the program writes and prints:
- the end point: the free values of param.sav (printed to 8 significant digits) equal the session's x printed alike;
- the objective: the program's "Total Chi-Squared" line (12 significant digits) equals result.objective printed alike;
- the covariance: with --covariance-band (HESSE after MIGRAD), covariance.dat's band block (7 significant digits)
  equals the session's covariance on those columns; in the program's MINOS run covariance_matrix.out (7 significant
  digits) equals result.covariance; the errors of param.sav equal result.errors (8 significant digits);
- MINOS: param.errors (8 significant digits) equals result.minos.
Cases: tests/13N with migrad (band, call limit 300); T2 (tests/upgrade, THM data, so the program's plain fit is
migrad-gradient) with migrad-gradient (MINOS and band), lm and gsl-lm (band); T2 with a nuisance prior
(<parameterSettings>) on energy_1, 12.4215 +- 0.001 MeV (lm); tests/data_covariance rel.azr, a data covariance with
relvector terms (migrad-gradient
with a call limit of 400, then the re-freeze of the model).  The fitters are deterministic at one thread (the same
code on the same start), so the comparison is exact to the printed digits.  Slow (several minutes at one thread).

Run: python3 tests/pyazr/fit_with_engine_test.py   (AZURE2_BIN: the program; FIT_WITH_ENGINE_CASES=<word>: the
cases whose label has the word, without the example)
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
TESTS = os.path.join(ROOT, "tests")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def copy_project(src, dst, project, extra_dirs=(), edit=None):
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("output", "checks", "expected"))
    for d in extra_dirs:
        shutil.copytree(d, os.path.join(dst, os.path.basename(d)), dirs_exist_ok=True)
    for d in ("output", "checks"):
        os.makedirs(os.path.join(dst, d), exist_ok=True)
    path = os.path.join(dst, project)
    if edit:
        text = open(path).read()
        open(path, "w").write(edit(text))
    return path


def run_program(azr, flags, command):
    proc = subprocess.run([os.environ["AZURE2_BIN"], *command, *flags, os.path.basename(azr)], stdin=subprocess.DEVNULL,
                          text=True, cwd=os.path.dirname(azr), capture_output=True, timeout=3000)
    out = os.path.join(os.path.dirname(azr), "output")
    total = re.findall(r"^Total Chi-Squared: (\S+)", proc.stdout, re.M)
    sav = {}
    with open(os.path.join(out, "param.sav")) as f:
        for line in f:
            p = line.split()
            if len(p) == 3 and not p[0].startswith("#"):
                sav[p[0]] = (float(p[1]), float(p[2]))
    return {"total": float(total[-1]) if total else None, "sav": sav, "out": out, "stdout": proc.stdout}


def sig(v, n):
    return f"{v:.{n - 1}e}"


def same(a, b, n):
    return all(sig(x, n) == sig(y, n) for x, y in zip(a, b)) and len(a) == len(b)


def covariance_matrix_out(path):
    """covariance_matrix.out: the covariance block, as (columns, matrix)."""
    lines = open(path).read().splitlines()
    i = lines.index("Covariance Matrix") + 2
    cols = [int(c) for c in lines[i].split()]
    rows = [[float(v) for v in lines[i + 1 + k].split()[1:]] for k in range(len(cols))]
    return cols, rows


def compare(label, prog, r, s, np, band_file=False, minos=False, hesse=True):
    names = r.names
    check(f"{label}: status {r.status}, objective {r.objective:.12g} = the program's {prog['total']!r}",
          prog["total"] is not None and sig(r.objective, 12) == sig(prog["total"], 12))
    x_prog = [prog["sav"][n][0] for n in names]
    check(f"{label}: the end point equals param.sav ({len(names)} free values, 8 digits)", same(x_prog, r.x, 8),
          f"first difference "
          f"{next(((n, a, b) for n, a, b in zip(names, x_prog, r.x) if sig(a, 8) != sig(b, 8)), None)}")
    # After MIGRAD the program runs HESSE only with --covariance-band; a session always does.
    if r.errors is not None and hesse:
        # The level parameters' errors to the digits of param.sav; those of norms, shifts, abundances and target
        # parameters to 1e-6 relative (3e-8 measured on T2: docs/upgrade/fit.md section 4).
        level = re.compile(r"^(energy|width|rwa|anc)_")
        e_prog = [prog["sav"][n][1] for n in names]
        lv = [k for k, n in enumerate(names) if level.match(n)]
        other = [k for k in range(len(names)) if k not in lv]
        worst = max((abs(r.errors[k] - e_prog[k]) / abs(e_prog[k]) for k in other if e_prog[k]), default=0.0)
        check(f"{label}: the level parameters' errors equal param.sav's (8 digits; {len(lv)} of them), the others to "
              f"{worst:.1e} (at most 1e-6)",
              same([e_prog[k] for k in lv], [r.errors[k] for k in lv], 8) and worst <= 1e-6,
              f"{[(n, a, float(b)) for n, a, b in zip(names, e_prog, r.errors) if sig(a, 8) != sig(b, 8)][:3]}")
    if band_file:
        lines = open(os.path.join(prog["out"], "covariance.dat")).read().splitlines()
        cols = lines[0].split()[1:]
        rows = [[float(v) for v in l.split()] for l in lines if l and not l.startswith("#")]
        idx = [names.index(c) for c in cols]
        mine = [[r.covariance[a, b] for b in idx] for a in idx]
        check(f"{label}: covariance.dat ({len(cols)} band columns) equals the covariance (7 digits)",
              all(same(p, m, 7) for p, m in zip(rows, mine)))
    if minos:
        cols, rows = covariance_matrix_out(os.path.join(prog["out"], "covariance_matrix.out"))
        check(f"{label}: covariance_matrix.out ({len(cols)} x {len(cols)}) equals the covariance (7 digits)",
              len(cols) == len(names) and all(same(p, list(m), 7) for p, m in zip(rows, r.covariance)),
              f"{rows[0][:2]} against {list(r.covariance[0][:2])}")
        errs = {}
        with open(os.path.join(prog["out"], "param.errors")) as f:
            for line in f:
                p = line.split()
                errs[p[0]] = (float(p[2]), float(p[3]))
        lower = [errs[n][0] for n in names]
        upper = [errs[n][1] for n in names]
        check(f"{label}: MINOS equals param.errors (8 digits)",
              same(lower, [abs(r.minos[n][0]) for n in names], 8) and
              same(upper, [abs(r.minos[n][1]) for n in names], 8),
              f"{lower[:2]}/{upper[:2]} against {[r.minos[n] for n in names[:2]]}")


def main():
    try:
        sys.path.insert(0, ROOT)
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import numpy as np
        import pyazr
        pyazr.kinds()
    except Exception as err:  # engine not built
        print(f"skip: engine not available ({type(err).__name__}: {err})")
        return 77
    if not os.environ.get("AZURE2_BIN"):
        print("skip: AZURE2_BIN (the program) is not set")
        return 77
    fitters = pyazr.kinds()["fitters"]  # the engine's table: the program's switch of a method is its row's option

    def switch(method):
        return fitters[method]["option"]

    t2 = os.path.join(TESTS, "upgrade", "T2")
    nuisance_row = "Level 1 Energy (MeV) 12.4215 0 0 0.001 0 1 level 0"

    def with_prior(text):
        return text.replace("<parameterSettings>\n</parameterSettings>",
                            "<parameterSettings>\n" + nuisance_row + "\n</parameterSettings>", 1)

    cases = [
        # label, source dir, project, extra dirs, edit, program flags, command, method, minos, max_calls, band file
        # T2 has THM data: the program's plain fit of it is migrad-gradient; Minuit's numerical gradient on 13N.
        ("13N migrad", os.path.join(TESTS, "13N"), "13N.azr", (), None, ["--covariance-band", "--max-calls=300"],
         ["fit"], "migrad", False, 300, True),
        ("T2 migrad-gradient", t2, "T2.azr", (), None, [switch("migrad-gradient"), "--covariance-band"],
         ["minos", "--delta-chi2=1"], "migrad-gradient", True, None, False),
        ("T2 lm", t2, "T2.azr", (), None, [switch("lm"), "--covariance-band"], ["fit"], "lm", False, None, True),
        ("T2 gsl-lm", t2, "T2.azr", (), None, [switch("gsl-lm"), "--covariance-band"], ["fit"], "gsl-lm",
         False, None, True),
        ("T2 with a nuisance prior on energy_1, lm", t2, "T2.azr", (), with_prior, [switch("lm")], ["fit"], "lm",
         False, None, False),
        ("data covariance with relvector terms, migrad-gradient", os.path.join(TESTS, "data_covariance"), "rel.azr",
         (os.path.join(TESTS, "angle_mapping", "data"),), None, [switch("migrad-gradient"), "--max-calls=400"],
         ["fit"], "migrad-gradient", False, 400, False),
    ]
    # Every fitter of the table has a case (recipe "A fitter" in CLAUDE.md).  The one exception: nlopt, a build option
    # (-DUSE_NLOPT=ON; without it the table has no such row), derivative-free, no covariance: not built where this test
    # was written, so no case was made for it.
    no_case = {"nlopt"}
    check("every fitter of the engine's table has a case here: " + ", ".join(fitters),
          all(m in {c[7] for c in cases} or m in no_case for m in fitters),
          f"no case for {[m for m in fitters if m not in {c[7] for c in cases} and m not in no_case]}")
    with tempfile.TemporaryDirectory() as tmp:
        only = os.environ.get("FIT_WITH_ENGINE_CASES", "")  # a word of the labels: run those cases only
        for k, (label, src, project, extra, edit, flags, command, method, minos, calls, band) in enumerate(cases):
            if only and only not in label:
                continue
            prog = run_program(copy_project(src, os.path.join(tmp, f"p{k}"), project, extra, edit), flags, command)
            with pyazr.Session(copy_project(src, os.path.join(tmp, f"s{k}"), project, extra, edit)) as s:
                r = s.fit_with_engine(method, minos=minos, max_calls=calls)
                compare(label, prog, r, s, np, band_file=band, minos=minos,
                        hesse=not method.startswith("migrad") or "--covariance-band" in flags)
                if "relvector" in label:
                    check(f"{label}: the session re-froze the model as the program did",
                          "relvector model re-frozen" in r.log and "relvector model re-frozen" in prog["stdout"])
        if not only:
            # The example of docs/source/user_guide/pyazr.rst ("The engine's own fitters"), run as written on a
            # copy of T2.
            page = open(os.path.join(ROOT, "docs", "source", "user_guide", "pyazr.rst")).read()
            block = page.split("The engine's own fitters", 1)[1].split(".. code-block:: python", 1)[1]
            lines = []
            for line in block.splitlines()[1:]:
                if line and not line.startswith("   "):
                    break
                lines.append(line[3:])
            code = "\n".join(lines).replace('"T2.azr"', repr(copy_project(t2, os.path.join(tmp, "doc"), "T2.azr")))
            code = code.replace('"T2_fitted.azr"', repr(os.path.join(tmp, "doc", "T2_fitted.azr")))
            scope = {"pyazr": pyazr}
            try:
                exec(code, scope)
                check("the example of the user guide runs: lm, then MIGRAD with a stop request, then save_fit",
                      scope["r"].status in ("converged", "stopped") and
                      os.path.isfile(os.path.join(tmp, "doc", "T2_fitted.azr")), scope["r"].status)
            except Exception as e:  # report, do not stop the other checks
                check("the example of the user guide runs", False, f"{type(e).__name__}: {e}")
            finally:
                if "s" in scope:
                    scope["s"].close()
    print("FAIL" if failures else "OK", f"({len(failures)} failure(s))" if failures else "")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
