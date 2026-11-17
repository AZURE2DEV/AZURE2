#!/usr/bin/env python3
"""A session on a project with a line window refuses every analytic derivative and keeps every value.

A capture data set with ``line_window`` is integrated over the final state's line shape (src/LineWindow.cpp).  The
adjoint that gives the analytic gradient, the Jacobian of the residual rows and the model sensitivities takes the
T-matrix at the nominal final state instead, so its numbers are wrong (on this model the Jacobian's column of the
final level's energy was 0 where central differences give -2 to -15).  The engine refuses them where every such
derivative starts (GradEnergyTasks.cpp, ``Forward``; LineWindow.h, ``LineWindowDerivativeRefused``), and the binding
raises ``pyazr.EngineError``.  Checked here:

1. values work: ``Objective(x)``, ``Objective.residuals``, ``Objective.log_prob``, ``Session.evaluate``;
2. every derivative raises ``EngineError`` naming ``line_window``: ``Objective.value_and_grad``, ``.grad``,
   ``.jacobian``, ``.residuals_and_jacobian``, ``Session.model_gradients``, and the engine's fitters that read the
   gradient or the rows (``fit_with_engine`` with every row of ``kinds()["fitters"]`` that says so);
3. the same project without ``line_window`` gives its gradient and Jacobian (the refusal is the line window's).

Model: tests/ec_unbound_broad/models/a_w1.azr (12C + p, capture into the 1/2+ state of 13N at E_x = 2.3655 MeV,
Gamma_p = 39 keV) with one windowed data point.  Fails without the refusal: the derivatives return numbers.

Run from anywhere:  python3 tests/pyazr/line_window_derivatives_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import numpy as np  # noqa: E402
    import pyazr  # noqa: E402
    from pyazr import EngineError, Session  # noqa: E402
except Exception as err:
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)

MODEL = os.path.join(ROOT, "tests", "ec_unbound_broad", "models", "a_w1.azr")
DATASET = """dataset lw
    reaction       12C(p,g)13N[2.3655]
    quantity       dsigma_domega
    unit           b/sr
{window}    file           lw.dat
    columns        E_lab[MeV] theta_lab[deg] value error
    norm           1 fixed
"""


def write_project(directory, window):
    """The model with one data point (E_lab = 0.704561 MeV, 90 deg), windowed over +-39 keV or not."""
    os.makedirs(os.path.join(directory, "output"))
    os.makedirs(os.path.join(directory, "checks"))
    with open(MODEL) as f:
        text = f.read()
    block = DATASET.format(window="    line_window    39 keV\n" if window else "")
    edited = text.replace("<data>\n</data>", "<data>\n" + block + "</data>", 1)
    assert edited != text, "the model's empty <data> block was not found"
    path = os.path.join(directory, "lw.azr")
    with open(path, "w") as f:
        f.write(edited)
    with open(os.path.join(directory, "lw.dat"), "w") as f:
        f.write("0.704561 90 1e-6 1e-7\n")
    return path


def refused(name, call):
    try:
        out = call()
    except EngineError as err:
        check(f"{name} raises EngineError naming line_window", "line_window" in str(err), str(err))
        return
    except Exception as err:  # another class is not the house refusal
        check(f"{name} raises EngineError naming line_window", False, f"{type(err).__name__}: {err}")
        return
    check(f"{name} raises EngineError naming line_window", False, f"returned {out!r}"[:200])


work = tempfile.mkdtemp(prefix="line_window_derivatives_")
try:
    s = Session(write_project(os.path.join(work, "windowed"), True))
    f = s.objective()
    x = np.asarray(f.x0, float)

    print("1. values")
    value = f(x)
    check("Objective(x) is a finite number", np.isfinite(value), value)
    rows = f.residuals(x)
    check("Objective.residuals: sum of squares equals the objective", np.isclose(np.sum(rows ** 2), value, rtol=1e-9),
          (np.sum(rows ** 2), value))
    check("Objective.log_prob = -f/2", np.isclose(f.log_prob(x), -0.5 * value, rtol=1e-9), (f.log_prob(x), value))
    ev = s.evaluate(x)
    check("Session.evaluate gives the windowed point's model", len(ev.data) == 1 and np.isfinite(ev.data[0].model[0])
          and ev.data[0].model[0] > 0, ev.data[0].model if ev.data else None)

    print("\n2. every analytic derivative is refused")
    refused("Objective.value_and_grad", lambda: f.value_and_grad(x))
    refused("Objective.grad", lambda: f.grad(x))
    refused("Objective.jacobian", lambda: f.jacobian(x))
    refused("Objective.residuals_and_jacobian", lambda: f.residuals_and_jacobian(x))
    refused("Session.model_gradients", lambda: s.model_gradients(x))
    derivative_fitters = [n for n, r in pyazr.kinds()["fitters"].items() if r["gradient"] or r["residual_rows"]]
    check("the table names the fitters that read the gradient or the rows", "lm" in derivative_fitters
          and "migrad-gradient" in derivative_fitters, derivative_fitters)
    for name in derivative_fitters:
        refused(f"fit_with_engine('{name}')", lambda name=name: s.fit_with_engine(name, max_calls=3))
    check("the values still work after the refusals", np.isclose(f(x), value, rtol=1e-12), (f(x), value))
    s.close()

    print("\n3. without line_window the derivatives are there")
    plain = Session(write_project(os.path.join(work, "plain"), False))
    g = plain.objective()
    y = np.asarray(g.x0, float)
    val, grad = g.value_and_grad(y)
    check("value_and_grad: a finite gradient", np.all(np.isfinite(grad)) and grad.size == y.size, grad)
    check("jacobian: one row per residual", g.jacobian(y).shape == (g.residuals(y).size, y.size))
    plain.close()
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\nFAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("\nPASS: a line window refuses every analytic derivative and keeps the values")
