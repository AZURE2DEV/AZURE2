#!/usr/bin/env python3
"""The analytic derivatives of a project with a line window agree with central differences of its values.

A capture data set with ``line_window`` is integrated over the final state's line shape (src/LineWindow.cpp): its
value is a weighted sum over nodes in the final energy E_f.  The adjoint runs over that sum (GradEnergyTasks.cpp,
``LineUnitAdjoint``): at fixed nodes each node's T-matrix cotangents join the task's bilinear terms, and the line
shape's own dependence on the final level's energy and amplitude -- they place the nodes and set the phase and the
amplitude of the final scattering state -- comes from central differences of the line integration itself
(LineWindow.h, ``LineShapeChain``).  Before this the adjoint took the T-matrix at the nominal final state (the
Jacobian's column of the final level's energy was 0 where central differences give -2 to -15), and every analytic
derivative of such a project was refused.  Checked here, on three models of tests/ec_unbound_broad (12C + p, capture
into the 1/2+ state of 13N at E_x = 2.3655 MeV, Gamma_p = 39 keV):

  a_w1    direct (external) capture alone, differential;
  d_w1    a 5/2+ resonance decaying by E2 into the final state, with external capture, angle-integrated;
  g_w1    the same with external capture into the pair off (the final level found by J^pi and energy);

each with a windowed data set (one point where external capture is on: a moved final-level energy rebuilds that
point's capture integrals, seconds each time; two points on g_w1) and the final level's energy and proton width free:

1. ``Objective.value_and_grad`` and ``.grad`` against central differences of ``Objective(x)``;
2. ``Objective.jacobian`` / ``.residuals_and_jacobian`` against central differences of ``.residuals``;
3. ``Session.model_gradients`` against central differences of the models of ``Session.evaluate``;
4. the columns of the final level's energy and width are not zero (the line shape's share is there);
5. the values after the derivatives are the values before them (the differences leave no trace in the records);
6. the engine's fitters that read the gradient or the rows run (``fit_with_engine``), on a_w1.

Fails without the adjoint over the line shape: the columns of the final level are wrong (and, with the refusal of
before, every derivative raises).

Run from anywhere:  python3 tests/pyazr/line_window_derivatives_test.py
"""
import os
import shutil
import sys
import tempfile
import time

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
    from pyazr import Session  # noqa: E402
except Exception as err:
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)

MODELS = os.path.join(ROOT, "tests", "ec_unbound_broad", "models")
DATASET = """dataset lw
    reaction       12C(p,g)13N[2.3655]
    quantity       {quantity}
    unit           {unit}
    line_window    39 keV
    file           lw.dat
    columns        {columns}
    norm           1 fixed
"""
# model -> (quantity, unit, columns, rows of energy [angle]); the data values are the model's times 1.2 +- 10 %
CASES = {
    "a_w1": ("dsigma_domega", "b/sr", "E_lab[MeV] theta_lab[deg] value error",
             ["0.975546 45"]),
    "d_w1": ("sigma", "b", "E_lab[MeV] _ value error", ["1.755983 0"]),
    "g_w1": ("sigma", "b", "E_lab[MeV] _ value error", ["1.755983 0", "2.16788 0"]),
}


def write_project(directory, model, values=None):
    quantity, unit, columns, rows = CASES[model]
    os.makedirs(os.path.join(directory, "output"), exist_ok=True)
    os.makedirs(os.path.join(directory, "checks"), exist_ok=True)
    with open(os.path.join(MODELS, model + ".azr")) as f:
        text = f.read()
    block = DATASET.format(quantity=quantity, unit=unit, columns=columns)
    edited = text.replace("<data>\n</data>", "<data>\n" + block + "</data>", 1)
    assert edited != text, "the model's empty <data> block was not found"
    path = os.path.join(directory, "lw.azr")
    with open(path, "w") as f:
        f.write(edited)
    with open(os.path.join(directory, "lw.dat"), "w") as f:
        for i, row in enumerate(rows):
            v = 1.0 if values is None else 1.2 * values[i]
            f.write(f"{row} {v:.6e} {0.1 * v:.6e}\n")
    return path


def central(fun, x, i):
    """Central difference of fun (a tuple of arrays) in x[i]: one evaluation of the tuple per side, because a moved
    final-level energy rebuilds the capture integrals of every node of every windowed point."""
    h = 2e-5 * max(abs(x[i]), 1e-3)
    up, dn = x.copy(), x.copy()
    up[i] += h
    dn[i] -= h
    hi, lo = fun(up), fun(dn)
    return [(np.asarray(a, float) - np.asarray(b, float)) / (2 * h) for a, b in zip(hi, lo)]


def close(a, b, rtol):
    """|a - b| <= rtol * the scale of the column (its largest entry), entry by entry."""
    a, b = np.atleast_1d(a), np.atleast_1d(b)
    scale = max(np.max(np.abs(b)), 1e-300)
    return bool(np.all(np.abs(a - b) <= rtol * scale)), float(np.max(np.abs(a - b)) / scale)


def models_of(s, x):
    return np.concatenate([np.asarray(d.model, float) for d in s.evaluate(x).data])


work = tempfile.mkdtemp(prefix="line_window_derivatives_")
try:
    for model in CASES:
        print(f"{model}")
        # the data from the model's own values, so that the residuals are of order one
        s = Session(write_project(os.path.join(work, model), model))
        start = models_of(s, np.asarray(s.objective().x0, float))
        s.close()
        s = Session(write_project(os.path.join(work, model), model, start))
        f = s.objective()
        x = np.asarray(f.x0, float)
        names = [str(n) for n in getattr(f, "names", range(x.size))]
        value = f(x)
        rows = np.asarray(f.residuals(x), float)
        check("the values are finite, the residuals of order one", np.isfinite(value) and np.all(np.isfinite(rows))
              and np.all(np.abs(rows) > 0.1), (value, rows))

        t0 = time.time()
        v2, grad = f.value_and_grad(x)
        jac = np.asarray(f.jacobian(x), float)
        r2, jac2 = f.residuals_and_jacobian(x)
        sens = s.model_gradients(x)
        t_derivatives = time.time() - t0
        t0 = time.time()
        check("value_and_grad gives the value", np.isclose(v2, value, rtol=1e-10), (v2, value))
        check("grad = value_and_grad", np.allclose(f.grad(x), grad, rtol=1e-10, atol=0))
        check("residuals_and_jacobian = residuals, jacobian", np.allclose(r2, rows, rtol=1e-10)
              and np.allclose(jac2, jac, rtol=1e-10, atol=0))
        check("the values after the derivatives are the values before", np.isclose(f(x), value, rtol=1e-12)
              and np.allclose(f.residuals(x), rows, rtol=1e-12, atol=0), (f(x), value))

        fd = [central(lambda y: (f(y), f.residuals(y), models_of(s, y)), x, i) for i in range(x.size)]
        fd_grad = np.array([d[0] for d in fd]).ravel()
        fd_jac = np.array([d[1] for d in fd]).T
        fd_sens = np.array([d[2] for d in fd]).T
        print(f"  times: derivatives {t_derivatives:.1f} s, central differences {time.time() - t0:.1f} s")
        ok, worst = close(grad, fd_grad, 2e-5)
        check("gradient = central differences of the objective", ok, f"worst {worst:.2e}: {grad} vs {fd_grad}")
        try:
            sens = np.concatenate([np.asarray(d, float).reshape(-1, x.size) for d in sens])
        except Exception:
            sens = np.asarray(sens, float).reshape(-1, x.size)
        check("model_gradients has one row per point and one column per parameter", sens.shape == fd_sens.shape,
              (sens.shape, fd_sens.shape))
        for i in range(x.size):
            ok, worst = close(jac[:, i], fd_jac[:, i], 2e-5)
            check(f"Jacobian column {names[i]} = central differences of the residuals", ok,
                  f"worst {worst:.2e}: {jac[:, i]} vs {fd_jac[:, i]}")
            if sens.shape == fd_sens.shape:
                ok, worst = close(sens[:, i], fd_sens[:, i], 2e-5)
                check(f"model gradient {names[i]} = central differences of the models", ok,
                      f"worst {worst:.2e}: {sens[:, i]} vs {fd_sens[:, i]}")
        # the final level: its energy is the parameter whose value is 2.3655, its proton width the one after it
        final = [i for i in range(x.size) if abs(x[i] - 2.3655) < 1e-9]
        check("the final level's energy is free", len(final) == 1, (final, x))
        if len(final) == 1:
            check("the final level's energy column is not zero", np.all(np.abs(jac[:, final[0]]) > 0), jac[:, final[0]])
            check("the final level's width column is not zero", np.all(np.abs(jac[:, final[0] + 1]) > 0),
                  jac[:, final[0] + 1])

        if model == "a_w1":
            fitters = [n for n, r in pyazr.kinds()["fitters"].items() if r["gradient"] or r["residual_rows"]]
            check("the table names the fitters that read the gradient or the rows", "lm" in fitters
                  and "migrad-gradient" in fitters, fitters)
            for name in fitters:
                try:
                    s.fit_with_engine(name, max_calls=3)
                    check(f"fit_with_engine('{name}') runs", True)
                except Exception as err:
                    check(f"fit_with_engine('{name}') runs", False, f"{type(err).__name__}: {err}"[:300])
        s.close()
        print()
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("PASS: the analytic derivatives over a line window agree with central differences")
