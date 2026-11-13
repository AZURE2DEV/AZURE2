#!/usr/bin/env python3
"""Wrong calls raise pyazr's error classes, with a message in three parts (docs/upgrade/pyazr3.md section 2.5).

A table of wrong calls on a copy of tests/13N (and a few without a session).  Each row names the call, the class it
must raise, the builtin class it raised before (so that an existing ``except`` keeps catching it), and three texts the
message must hold: (1) the value given that was wrong, (2) what is allowed, (3) the call that shows or fixes it.
Also checked:

- ``pyazr.Error`` is the base of every class; ``pyazr.AZURE2Error`` is ``pyazr.EngineError`` and
  ``pyazr._azure2.AZURE2Error`` (one class), named "AZURE2Error"; a wrong-length vector is still an ``AZURE2Error``;
- the engine's notice "Could not write to EC Amplitude File." met by a session whose output directory does not exist
  is a ``pyazr.EngineWarning`` that names the file and the consequence, and no longer a line on standard output;
- the documented zeros of ``coulomb_functions`` below threshold are kept (S alone at E <= 0).

Needs the compiled engine; without it exits 77 (ctest: Skipped).  Threads: as the caller's environment sets them.

Run from anywhere:  python3 tests/pyazr/errors_test.py
"""
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({detail})"))
    if not ok:
        failures.append(name)


try:
    import numpy as np
    import pyazr
    from pyazr import _azure2, transform
except Exception as err:  # pragma: no cover
    print(f"skip: the engine module is not available ({type(err).__name__}: {err})")
    sys.exit(77)


def expect(label, call, cls, builtin, parts):
    """Run ``call``; it must raise ``cls`` (also a ``builtin`` and a pyazr.Error) whose text holds every one of
    ``parts`` (what was wrong, what is allowed, which call helps)."""
    try:
        call()
    except Exception as e:  # noqa: BLE001
        text = str(e)
        check(f"{label}: raises {cls.__name__}", isinstance(e, cls), f"{type(e).__name__}: {text[:300]}")
        check(f"{label}: is a {builtin.__name__} and a pyazr.Error",
              isinstance(e, builtin) and isinstance(e, pyazr.Error), type(e).__mro__)
        missing = [p for p in parts if p not in text]
        check(f"{label}: the message has its three parts", not missing, f"missing {missing} in: {text[:400]}")
        check(f"{label}: the message names no repr of an object", "ProjectLevel(" not in text and len(text) < 800,
              text[:200])
        return text
    check(f"{label}: raises {cls.__name__}", False, "no exception")
    return ""


# -- the classes ----------------------------------------------------------------------------------------------------
for name in ("FormatError", "ParameterError", "ParameterKeyError", "ParameterIndexError", "ParameterVectorError",
             "EngineError", "DomainError", "DomainZeroDivision", "SessionClosed", "MissingFile"):
    check(f"pyazr.{name} derives from pyazr.Error", issubclass(getattr(pyazr, name), pyazr.Error))
check("pyazr.EngineWarning is a UserWarning", issubclass(pyazr.EngineWarning, UserWarning))
check("AZURE2Error is EngineError and _azure2.AZURE2Error (one class)",
      pyazr.AZURE2Error is pyazr.EngineError is _azure2.AZURE2Error)
check("the class keeps the name AZURE2Error", pyazr.EngineError.__name__ == "AZURE2Error")
check("_azure2.ParameterVectorError is pyazr.ParameterVectorError",
      _azure2.ParameterVectorError is pyazr.ParameterVectorError)
check("FormatError is pyazr.project.FormatError (one class)", pyazr.project.FormatError is pyazr.FormatError)
for cls, base in ((pyazr.FormatError, ValueError), (pyazr.ParameterError, ValueError),
                  (pyazr.EngineError, RuntimeError), (pyazr.DomainError, ValueError),
                  (pyazr.SessionClosed, RuntimeError), (pyazr.ParameterKeyError, KeyError),
                  (pyazr.ParameterIndexError, IndexError), (pyazr.MissingFile, FileNotFoundError),
                  (pyazr.ParameterVectorError, pyazr.EngineError)):
    check(f"{cls.__name__} is a {base.__name__}", issubclass(cls, base))

# -- wrong calls without a session ----------------------------------------------------------------------------------
tmp = tempfile.mkdtemp(prefix="pyazr_errors_")
try:
    proj_dir = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(ROOT, "tests", "13N"), proj_dir,
                    ignore=shutil.ignore_patterns("output", "checks", "expected"))
    azr = os.path.join(proj_dir, "13N.azr")
    missing = os.path.join(proj_dir, "nope.azr")

    text = expect("azure2 on a missing project file", lambda: pyazr.azure2(missing), pyazr.MissingFile,
                  FileNotFoundError, ["nope.azr", "does not exist", "13N.azr", "pyazr.convert_azr"])
    expect("Project.open on a missing project file", lambda: pyazr.Project.open(missing), pyazr.MissingFile,
           FileNotFoundError, ["nope.azr", "does not exist", "existing format-3 project file"])
    expect("azure2 with a cwd that does not exist", lambda: pyazr.azure2(azr, cwd=os.path.join(tmp, "no_dir")),
           pyazr.MissingFile, FileNotFoundError, ["no_dir", "does not exist", "cwd=None"])
    project = pyazr.Project.open(azr)
    expect("Project[name] of a record it does not have", lambda: project["skowronski_luna_hpgX"],
           pyazr.ParameterKeyError, KeyError,
           ["'skowronski_luna_hpgX'", "did you mean", "Project.records()"])
    expect("penetrability below threshold", lambda: transform.penetrability(0, 1, 6, 0.93, 4.0, -0.1),
           pyazr.DomainError, ValueError, ["energy = -0.1", "> 0", "azure2.coulomb_functions"])
    expect("penetrability with L = -1", lambda: transform.penetrability(-1, 1, 6, 0.93, 4.0, 0.5),
           pyazr.DomainZeroDivision, ZeroDivisionError, ["L = -1", "whole number >= 0", "help(pyazr.penetrability)"])
    expect("whittaker below its domain", lambda: transform.whittaker(0, 1, 6, 0.93, 4.0, 0.0),
           pyazr.DomainZeroDivision, ZeroDivisionError, ["binding_energy = 0.0", "> 0", "pyazr.whittaker"])
    expect("whittaker_wave with rho <= 0", lambda: pyazr.whittaker_wave(0, 1.0, [1.0, -1.0]),
           pyazr.DomainError, ValueError, ["rho = -1.0", "rho finite > 0", "help(pyazr.whittaker_wave)"])

    # -- a session whose output directory does not exist: the engine's notice is a warning ---------------------
    out_dir = os.path.join(proj_dir, "output")
    check("the copy has no output directory", not os.path.exists(out_dir))
    code = ("import warnings, pyazr\n"
            "with warnings.catch_warnings(record=True) as w:\n"
            "    warnings.simplefilter('always')\n"
            f"    s = pyazr.azure2({azr!r})\n"
            "    for x in w:\n"
            "        print('WARNED', x.category.__name__, '|', str(x.message).replace(chr(10), ' '))\n"
            "s.close()\n")
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=600,
                         env=dict(os.environ, PYTHONPATH=ROOT + os.pathsep + os.environ.get("PYTHONPATH", "")))
    warned = [line for line in run.stdout.splitlines() if line.startswith("WARNED EngineWarning")]
    check("a session without its output directory: EngineWarnings", run.returncode == 0 and warned,
          (run.returncode, run.stdout[-600:], run.stderr[-600:]))
    check("each names the file intEC.dat and the consequence",
          warned and all(os.path.join(out_dir, "intEC.dat") in w and "computes them again" in w
                         and "Create the output directory" in w for w in warned), warned)
    printed = "\n".join(line for line in run.stdout.splitlines() if not line.startswith("WARNED"))
    check("the notices are no longer printed on standard output",
          "EC Amplitude File." not in printed, printed[-400:])

    # -- wrong calls on a session ---------------------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", pyazr.EngineWarning)
        s = pyazr.azure2(azr)
    x = np.asarray(s.params_rwa, float)
    n = x.size
    first_free = s.parameters.free[0].name

    expect("parameters.by_name of an unknown name", lambda: s.parameters.by_name("nonexistent"),
           pyazr.ParameterKeyError, KeyError, ["'nonexistent'", "the names are Parameter.name", "table()"])
    expect("pairs.by_number of an unknown pair", lambda: s.pairs.by_number(9), pyazr.ParameterKeyError, KeyError,
           ["no pair 9", "the pair numbers", "azure2.pairs"])
    expect("without_level on a J^pi the model does not have", lambda: s.without_level(x, jpi="9/2+"),
           pyazr.ParameterKeyError, KeyError, ["jpi=9/2+", "the levels:", "azure2.physical_levels()"])
    expect("datasets.by_name of an unknown data set", lambda: s.datasets.by_name("nope"),
           pyazr.ParameterKeyError, KeyError, ["'nope'", "the names:", "azure2.datasets"])
    expect("set_channel_radius of an unknown pair", lambda: s.set_channel_radius(9, 4.0),
           pyazr.ParameterKeyError, KeyError, ["no pair 9", "have [", "azure2.pairs"])
    expect("set_channel_radius with a negative radius", lambda: s.set_channel_radius(1, -1.0),
           pyazr.ParameterError, ValueError, ["got -1.0", "positive length in fm", "channel_radius"])
    expect("objective with one parameter too few", lambda: s.objective(x[:-1]), pyazr.ParameterVectorError,
           pyazr.AZURE2Error, [f"got {n - 1}", f"expected the free parameter vector ({n} values", "params_rwa"])
    expect("calculate_chi2_rwa with one parameter too many", lambda: s.calculate_chi2_rwa(np.append(x, 1.0)),
           pyazr.ParameterVectorError, ValueError, [f"got {n + 1}", f"({n} values", "azure2.parameters.free"])
    bad = x.copy()
    bad[0] = np.nan
    expect("objective with a NaN", lambda: s.objective(bad), pyazr.ParameterVectorError, RuntimeError,
           ["entry 0", first_free, "is nan", "must be a finite number", "params_rwa"])
    bad = x.copy()
    bad[-1] = np.inf
    expect("objective_and_grad with an infinity", lambda: s.objective_and_grad(bad), pyazr.ParameterVectorError,
           ValueError, [f"entry {n - 1}", "is inf", "finite", "params_rwa"])
    expect("ParameterSet.put with the wrong number of values",
           lambda: s.parameters.free.put(x, np.ones(n + 1)),
           pyazr.ParameterError, ValueError, ["value(s) for", "one number for all", "set.free"])
    expect("save_fit with an x of the wrong length", lambda: s.save_fit(os.path.join(tmp, "f.azr"), x=x[:-1]),
           pyazr.ParameterError, ValueError, [f"x has {n - 1} entries", f"{n} free parameters", "params_rwa"])
    expect("update_sav_from_rwa_params without a param.sav", lambda: s.update_sav_from_rwa_params(x),
           pyazr.MissingFile, FileNotFoundError, ["param.sav does not exist", "existing param.sav", "save_fit"])
    E = np.array([0.5, 1.0])
    expect("coulomb_functions on pair 0", lambda: s.coulomb_functions(0, E), pyazr.ParameterIndexError, IndexError,
           ["pair 0", "its pairs are", "azure2.pairs"])
    expect("coulomb_functions on pair 9", lambda: s.coulomb_functions(9, E), pyazr.ParameterIndexError, IndexError,
           ["pair 9", "its pairs are", "azure2.pairs"])
    expect("coulomb_functions with L = -1", lambda: s.coulomb_functions(1, E, L=-1), pyazr.DomainError, ValueError,
           ["L = -1", "whole number >= 0", "azure2.parameters.widths"])
    photon = [p.number for p in s.pairs if p.is_photon][0]
    expect("coulomb_functions on the photon pair", lambda: s.coulomb_functions(photon, E, L=1), pyazr.DomainError,
           ValueError, [f"pair {photon}", "the particle pairs are", "azure2.pairs.particles"])

    # documented zeros are kept: a closed channel gives S alone (and the open energy its numbers)
    cf = s.coulomb_functions(1, np.array([-0.1, 0.5]), L=0)
    check("coulomb_functions below threshold: F, G, P are 0, S is not (documented)",
          cf["F"][0] == 0 and cf["G"][0] == 0 and cf["P"][0] == 0 and cf["S"][0] != 0 and cf["P"][1] > 0, cf)

    s.close()
    expect("objective on a closed session", lambda: s.objective(x), pyazr.SessionClosed, RuntimeError,
           ["azure2.objective", "is closed", "pyazr.azure2("])
    expect("calculate_chi2_rwa on a closed session", lambda: s.calculate_chi2_rwa(x), pyazr.SessionClosed,
           RuntimeError, ["azure2.calculate_chi2_rwa", "cannot be reopened", "is_alive()"])
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"{'FAILED: ' + str(len(failures)) if failures else 'all ok'}")
sys.exit(1 if failures else 0)
