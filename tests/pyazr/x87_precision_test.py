#!/usr/bin/env python3
"""A pyazr call computes the same bits in any thread, whatever the thread's x87 precision
(support/include/azure2/support/FloatingPointPrecision.h; api/src/pybind_engine_call.h).

On x86 `long double` is the x87 format, and the mantissa an operation keeps is the precision control of the thread's
control word: 64 bits (0x37F) or 53 (0x27F).  On Windows a Python thread starts with 53 bits and the main thread with
64 (MinGW's start-up), and the same call of T1 gave other bits in a threading.Thread than in the main thread
(session_lock_test.py (b), objective_object_test.py (d); 2026-11-16).  The engine now takes the 64-bit precision for
the length of every call from the binding and in every thread of its OpenMP regions.

(a) What the threads hold (printed for the CI's logs: the precision a Python thread and the OpenMP runtime's threads
    start with on this system): the main thread's control word, a threading.Thread's, and those of the threads of
    an OpenMP region of each, as found and inside the engine's scope.  Checked: inside every scope the precision is
    64 bits, and the calling thread's word is as it was after the call.
(b) T1 at a point: the objective's value and gradient, its rows and Jacobian, and evaluate(x)'s chi^2 and model
    values, computed in the main thread, equal bit for bit those computed
      - in a threading.Thread as it starts,
      - in a threading.Thread whose precision is lowered to 53 bits before the calls (as on Windows; on Linux the
        OpenMP threads that thread starts inherit the lowered word, so the regions' scopes are tested too),
      - in the main thread with its precision lowered to 53 bits (its OpenMP threads, made earlier, keep 64),
    and the lowered word is still in force after the calls.
(c) A representative of each group of the module's functions that compute (repo_consistency check 60 (iii)), called
    in the main thread and in a threading.Thread lowered to 53 bits, bit for bit: rows of the physics table
    (physics_evaluate: Coulomb functions, Whittaker function, penetrability, dS/dE, kinematics), the session's calls
    that release the GIL (transform, channel_functions), those that keep it (Objective(x), log_prob, scattering
    lengths, Park norms), a level door (pyazr.endf.from_session), block physics (collision_matrix,
    optical_scattering), a project's grid (extrapolation_grid, relativistic) and stopping; and the set-up: a Session
    made in a lowered thread gives, called in the main thread, the main session's Objective(x).
Without an x87 (arm64) the control words are 0, (b) compares the threads as they start only and (c) is skipped.

Can fail: with the build before the scope (trial 101), lowering a thread's precision with a ctypes helper moves 27 of
T1's 442 residuals by up to 1.4e-12 relative, its value by 8e-16 and 7 of 13 gradient components (2026-11-16).  (c)
against the build before PrecisionGuard (9d7ec793): the set-up moves Objective(x) by 4.1e-15 relative (2.1e-14
without the Coulomb table); against trial 101 also Objective(x) and log_prob (8.2e-16).  The other representatives
did not move on Linux in either build: their arithmetic is double (SSE) here, or ran in an OpenMP region's scope.

Run: python3 tests/pyazr/x87_precision_test.py
"""
import dataclasses
import os
import shutil
import sys
import tempfile
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
LOWERED = 0x27F  # precision control 10: a 53-bit mantissa, all exceptions masked, round to nearest
PC_BITS, PC_64 = 0x300, 0x300

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import Session
    from pyazr import _azure2
    _azure2._x87_control_words
except Exception as err:  # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def in_thread(fn):
    out = {}

    def run():
        try:
            out["value"] = fn()
        except BaseException as e:  # reported by the caller
            out["error"] = e
    t = threading.Thread(target=run)
    t.start()
    t.join()
    if "error" in out:
        raise out["error"]
    return out["value"]


def words_text(w):
    return (f"calling {w['calling']:#06x}, in the call {w['in_call']:#06x}, after {w['after_call']:#06x}; OpenMP "
            f"threads found {[hex(v) for v in w['omp_found']]}, in their scope {[hex(v) for v in w['omp_in_region']]}")


print("(a) the x87 control words")
main_words = _azure2._x87_control_words()
worker_words = in_thread(_azure2._x87_control_words)
x87 = main_words["calling"] != 0
print(f"  main thread: {words_text(main_words)}")
print(f"  a threading.Thread: {words_text(worker_words)}")
if x87:
    for label, w in (("main thread", main_words), ("a threading.Thread", worker_words)):
        inside = [w["in_call"]] + list(w["omp_in_region"])
        check(f"{label}: 64-bit precision in the call and in every thread of the region",
              all(v & PC_BITS == PC_64 for v in inside), words_text(w))
        check(f"{label}: its control word after the call is the word before it", w["after_call"] == w["calling"],
              words_text(w))
    lowered_words = in_thread(lambda: (_azure2._set_x87_control_word(LOWERED), _azure2._x87_control_words())[1])
    print(f"  a threading.Thread lowered to 53 bits: {words_text(lowered_words)}")
    check("a thread lowered to 53 bits: 64 bits in the call and in every thread of the region, 53 after the call",
          all(v & PC_BITS == PC_64 for v in [lowered_words["in_call"]] + list(lowered_words["omp_in_region"]))
          and lowered_words["calling"] == LOWERED and lowered_words["after_call"] == LOWERED, words_text(lowered_words))
else:
    print("  no x87 on this machine: the words are 0")

print("(b) T1 in the main thread and in other threads, bit for bit")
tmp = tempfile.mkdtemp(prefix="x87_precision_")
try:
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), os.path.join(tmp, "T1"),
                    ignore=shutil.ignore_patterns("output", "checks"))
    os.makedirs(os.path.join(tmp, "T1", "output"), exist_ok=True)
    with Session(os.path.join(tmp, "T1", "T1.azr")) as s:
        f = s.objective()
        x0 = s.parameters.x0
        x = x0 * (1 + 1e-3 * np.cos(np.arange(x0.size)))

        def compute():
            v, g = f.value_and_grad(x)
            r, J = f.residuals_and_jacobian(x)
            ev = s.evaluate(x)
            return {"value": np.array([v]), "gradient": np.asarray(g), "rows": np.asarray(r), "Jacobian": np.asarray(J),
                    "evaluate chi2": np.array([ev.chi2]),
                    "evaluate model": np.concatenate([np.ravel(d.model) for d in ev.data])}

        def lowered():
            _azure2._set_x87_control_word(LOWERED)
            out = compute()
            out["word after"] = _azure2._x87_control_words()["calling"]
            return out

        def differences(a, b):
            """'' when every entry of a equals b's bit for bit, else what differs and by how much."""
            text = []
            for k in b:
                u, v = np.ravel(a[k]).astype(float), np.ravel(b[k]).astype(float)
                if u.shape != v.shape or u.tobytes() != v.tobytes():
                    n = int(np.sum(u != v)) if u.shape == v.shape else -1
                    rel = float(np.max(np.abs(u - v) / np.maximum(np.abs(v), 1e-300))) if n > 0 else 0.0
                    text.append(f"{k}: {n} of {v.size} differ, largest relative {rel:.2e}")
            return "; ".join(text)

        main = compute()
        check("T1 has rows and a Jacobian to compare", main["rows"].size > 100 and main["Jacobian"].size > 1000,
              {k: v.size for k, v in main.items()})
        d = differences(in_thread(compute), main)
        check("a threading.Thread as it starts = the main thread", not d, d)
        if x87:
            w = in_thread(lowered)
            d = differences(w, main)
            check("a threading.Thread lowered to 53 bits = the main thread", not d, d)
            check("  and its word is still the lowered one after the calls", w["word after"] == LOWERED,
                  hex(w["word after"]))
            found = _azure2._x87_control_words()["calling"]
            try:
                m = lowered()
            finally:
                _azure2._set_x87_control_word(found)
            d = differences(m, main)
            check("the main thread lowered to 53 bits = the main thread at 64", not d, d)
finally:
    shutil.rmtree(tmp, ignore_errors=True)


print("(c) a representative of each kind of binding call that computes, in a thread lowered to 53 bits")


def leaves(v, path=""):
    """(path, array) of every number of `v` (dicts by key, sequences, dataclasses, numpy arrays); text as bytes."""
    if dataclasses.is_dataclass(v) and not isinstance(v, type):
        v = {f.name: getattr(v, f.name) for f in dataclasses.fields(v)}
    if isinstance(v, dict):
        for k in sorted(v, key=str):
            yield from leaves(v[k], f"{path}.{k}")
    elif isinstance(v, (list, tuple)) and not all(isinstance(e, (int, float, complex, bool)) for e in v):
        for i, e in enumerate(v):
            yield from leaves(e, f"{path}[{i}]")
    elif isinstance(v, (str, bytes)):
        yield path, np.frombuffer(v.encode() if isinstance(v, str) else v, dtype=np.uint8)
    elif v is None:
        yield path, np.zeros(0)
    else:
        yield path, np.asarray(v)


def leaf_differences(a, b):
    """'' when `a` and `b` hold the same numbers bit for bit, else the first leaves that differ and by how much."""
    la, lb = list(leaves(a)), list(leaves(b))
    if [p for p, _ in la] != [p for p, _ in lb]:
        return "the results have other fields"
    text = []
    for (p, u), (_, v) in zip(la, lb):
        if u.shape != v.shape or u.dtype != v.dtype or u.tobytes() != v.tobytes():
            if u.shape == v.shape and u.dtype.kind in "fc":
                n = int(np.sum(u != v))
                rel = float(np.max(np.abs(u - v) / np.maximum(np.abs(v), 1e-300)))
                text.append(f"{p or 'value'}: {n} of {v.size} differ, largest relative {rel:.1e}")
            else:
                text.append(f"{p or 'value'}: differs")
    return "; ".join(text[:3]) + (f" (and {len(text) - 3} more fields)" if len(text) > 3 else "")


def table_row(name):
    """Row `name` of the physics table on 40 points of its axis argument's sweep, the others at their example."""
    row = next(r for r in _azure2.physics_functions() if r["name"] == name)
    columns = []
    for a in row["arguments"]:
        if a["name"] == row["axis"] and a["sweep"]:
            lo, hi = a["sweep"]["lower"], a["sweep"]["upper"]
            col = np.geomspace(lo, hi, 40) if a["sweep"]["spacing"] == "log" else np.linspace(lo, hi, 40)
        else:
            x = a["example"] if a["example"] is not None else a["default"]
            x = a["choices"].index(x) if isinstance(x, str) else float(x)
            col = np.full(40, x)
        columns.append(np.ascontiguousarray(col, dtype=float))
    values, status = _azure2.physics_evaluate(name, columns)
    return {"values": np.asarray(values), "status": np.asarray(status)}


if not x87:
    print("  no x87 on this machine: nothing to lower")
else:
    from pyazr import physics as ph
    from pyazr import stopping, extrapolation_grid
    from pyazr.endf import model as endf_model

    tmp = tempfile.mkdtemp(prefix="x87_precision_c_")
    try:
        shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), os.path.join(tmp, "T1"),
                        ignore=shutil.ignore_patterns("output", "checks"))
        os.makedirs(os.path.join(tmp, "T1", "output"), exist_ok=True)
        azr = os.path.join(tmp, "T1", "T1.azr")
        with Session(azr) as s:
            f = s.objective()
            x0 = s.parameters.x0
            x = x0 * (1 + 1e-3 * np.cos(np.arange(x0.size)))
            pair = s.pairs[0].key
            e_cm = np.linspace(0.05, 3.0, 25)
            potential = dict(radius_convention="target_mass", target_mass_number=15, coulomb_r0_fm=1.25,
                             real_volume=(50.0, 1.17, 0.75), imaginary_surface=(5.0, 1.32, 0.6),
                             real_spin_orbit=(6.2, 1.01, 0.75))
            levels, gamma = np.array([2.37, 3.5]), np.array([[0.61, 0.12], [-0.83, 0.47]])
            # (group, label, the call): the groups of api/src's computing functions (repo_consistency check 60)
            calls = [
                ("physics table", "coulomb_waves", lambda: table_row("coulomb_waves")),
                ("physics table", "whittaker_wave", lambda: table_row("whittaker_wave")),
                ("physics table", "penetrability", lambda: table_row("penetrability")),
                ("physics table", "shift_derivative", lambda: table_row("shift_derivative")),
                ("physics table", "lab_to_cm_energy (kinematics)", lambda: table_row("lab_to_cm_energy")),
                ("session, GIL released", "Session.transform", lambda: s.transform(x)),
                ("session, GIL released", "Session.channel_functions",
                 lambda: s.channel_functions(pair, e_cm, l=1)),
                ("session, GIL held", "Objective(x)", lambda: f(x)),
                ("session, GIL held", "Objective.log_prob", lambda: f.log_prob(x)),
                ("session, GIL held", "Session.scattering_lengths", lambda: s.scattering_lengths(x)),
                ("session, GIL held", "Session.park_norms", lambda: s.park_norms(x)),
                ("level doors", "pyazr.endf.from_session (formal_levels_of_block)",
                 lambda: endf_model.from_session(s, x)),
                ("block physics", "physics.collision_matrix",
                 lambda: ph.collision_matrix(levels, gamma, 2.35, penetrability=np.array([0.0123, 0.3]),
                                             shift=np.array([-0.71, -0.1]), boundary=np.array([-0.4, 0.05]),
                                             hard_sphere_phase=np.array([-0.002, 0.2]),
                                             coulomb_phase=np.array([0.0, 0.35]))),
                ("block physics", "physics.optical_scattering",
                 lambda: ph.optical_scattering(1, 7, 1.0, 15.0, 4.6875, potential, angles_deg=[30.0, 90.0, 150.0],
                                               projectile_spin=0.5, l_max=16)),
                ("project reading", "extrapolation_grid (c.m. -> lab, relativistic)",
                 lambda: extrapolation_grid("range 0.1 3 MeV cm points 30 log", masses=(1.00782503, 15.0001089),
                                            relativistic=True)),
                ("stopping", "stopping.stopping and range_of",
                 lambda: (stopping.stopping(stopping.Material("Ta2O5"), "p", np.geomspace(0.05, 3.0, 20)),
                          stopping.range_of(stopping.Material("Ta2O5"), "a", np.geomspace(0.1, 5.0, 20)))),
            ]

            def lowered(fn):
                def run():
                    _azure2._set_x87_control_word(LOWERED)
                    return fn()
                return run

            for group, label, fn in calls:
                main = fn()
                d = leaf_differences(in_thread(lowered(fn)), main)
                check(f"{group}: {label}", not d, d)
            # The set-up: a session made in a lowered thread computes, in the main thread, the main session's bits.
            value = f(x)
            made = in_thread(lowered(lambda: Session(azr)))
            try:
                d = leaf_differences(made.objective()(x), value)
            finally:
                made.close()
            check("session set-up: Session(...) made in a lowered thread, Objective(x) in the main thread", not d, d)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok: the same bits in every thread")
