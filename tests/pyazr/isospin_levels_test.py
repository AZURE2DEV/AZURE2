#!/usr/bin/env python3
"""Isospin-labelled levels across compound nuclei (the records of <isospin>, include/IsospinTies.h) in a session.

On tests/isospin_levels/models (README there): iso_mix.azr is the truth, iso_free.azr the model with every reference
parameter, the two shifts and the mix free, with priors on the shifts and the mix.

  1. the parameters: the derived slots of the level block are fixed, the shifts and the mix are parameters of kinds
     isospin_shift and isospin_mix named by their records, with the priors the records write;
  2. the objective is the data chi-squared plus the priors of the shifts and the mix (the objective's term T12), and
     its least-squares rows sum to it;
  3. the ties act: moving a shift moves the chi-squared of the tied nucleus's data set and of no data set of the
     reference nucleus;
  4. the gradient of the objective and the Jacobian of its rows agree with central finite differences for every free
     parameter, the shifts and the mix through the ties;
  5. Project: an isospin_level record removed leaves a level ID across nuclei, which the engine refuses; a second
     reference is refused by the edit (E_ISOSPIN_RECORD).

The data files are written here from a session on the truth model (its extrapolations), as tests/isospin_levels/check.sh
does from a run of the program.

Needs the compiled engine; without it exits 77 (ctest: Skipped).

Run from anywhere:  python3 tests/pyazr/isospin_levels_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def main():
    sys.path.insert(0, ROOT)
    try:
        import numpy as np
        from pyazr import Session, Project, FormatError, AZURE2Error
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"SKIP: engine not importable ({e})")
        return 77

    models = os.path.join(ROOT, "tests", "isospin_levels", "models")
    work = tempfile.mkdtemp(prefix="pyazr_isospin_levels_")
    try:
        for d in ("output", "checks"):
            os.makedirs(os.path.join(work, d))
        shutil.copy(os.path.join(models, "iso_free.azr"), work)
        # Pseudo-data: smooth numbers on each data set's window with 1 % errors (the checks below compare the engine
        # with itself; they need points, not a fit).
        for name, lo, hi, value in (("nt", 2.0, 12.0, 0.5), ("ptn", 1.2, 10.0, 0.05), ("ddp", 0.1, 4.0, 0.02),
                                    ("ddn", 0.1, 4.0, 0.02)):
            with open(os.path.join(work, name + ".dat"), "w") as f:
                for e in np.linspace(lo + 0.05, hi - 0.05, 12):
                    f.write(f"{e:.6f} 0 {value:.6e} {0.01 * value:.6e}\n")
        path = os.path.join(work, "iso_free.azr")

        with Session(path, cwd=work) as s:
            print("1. the parameters")
            params = list(s.parameters)
            by_name = {p.name: p for p in params}
            iso = [p for p in params if p.kind in ("isospin_shift", "isospin_mix")]
            check("two shifts and a mix, named by their records",
                  [p.name for p in iso] == ["shift_4H.shift", "shift_4He.shift", "mix_1.amplitude"], [p.name for p in iso])
            check("they are free", all(not p.fixed for p in iso))
            records = {r["name"]: r for r in s._engine.parameter_records()}
            check("the shift's record: kind, record, nucleus",
                  records["shift_4H.shift"]["kind"] == "isospin_shift" and records["shift_4H.shift"]["isospin_record"] == "shift_4H" and
                  (records["shift_4H.shift"]["nucleus_Z"], records["shift_4H.shift"]["nucleus_A"]) == (1, 4),
                  records["shift_4H.shift"])
            check("the shift's prior is the record's (-0.86, 0.3)",
                  (records["shift_4H.shift"]["prior_centre"], records["shift_4H.shift"]["prior_sigma"]) == (-0.86, 0.3),
                  records["shift_4H.shift"])
            check("the mix's record: the channel it sets (d+d, l = 1, s = 1) and its prior (0, 0.1)",
                  records["mix_1.amplitude"]["kind"] == "isospin_mix" and records["mix_1.amplitude"]["L"] == 1 and
                  records["mix_1.amplitude"]["S"] == 1.0 and records["mix_1.amplitude"]["prior_sigma"] == 0.1,
                  records["mix_1.amplitude"])
            # level 1 (2-, T = 1): energy_1 in 4Li is free, energy_2 (4H) and energy_3 (4He) are derived
            check("the reference member's energy is free, the derived members' energies fixed",
                  not by_name["energy_1"].fixed and by_name["energy_2"].fixed and by_name["energy_3"].fixed,
                  [(n, by_name[n].fixed) for n in ("energy_1", "energy_2", "energy_3")])
            start = {p.name: p.start for p in params}
            check("a derived energy starts on its tie (E_ref + shift)",
                  abs(start["energy_2"] - (start["energy_1"] + start["shift_4H.shift"])) < 1e-12,
                  (start["energy_2"], start["energy_1"], start["shift_4H.shift"]))
            check("the mixed amplitude's slot starts on the mix", start["width_3_3"] == start["mix_1.amplitude"],
                  (start["width_3_3"], start["mix_1.amplitude"]))

            print("\n2. the objective")
            f = s.objective()
            x0 = f.x0
            names = list(f.names)
            at = s.evaluate(x0)
            prior = sum(((x0[names.index(n)] - c) / w) ** 2
                        for n, c, w in (("shift_4H.shift", -0.86, 0.3), ("shift_4He.shift", -0.64, 0.3), ("mix_1.amplitude", 0.0, 0.1)))
            check("objective == data chi2 + the priors of the shifts and the mix",
                  abs(f(x0) - (at.chi2 + prior)) <= 1e-9 * abs(f(x0)), f"{f(x0)} vs {at.chi2} + {prior}")
            rows = f.residuals(x0)
            check("sum of the rows squared == the objective", abs(float(np.sum(rows ** 2)) - f(x0)) <= 1e-9 * abs(f(x0)),
                  f"{float(np.sum(rows ** 2))} vs {f(x0)}")
            check("the last three rows are the priors", np.allclose(sorted(rows[-3:] ** 2), sorted(
                [((x0[names.index(n)] - c) / w) ** 2 for n, c, w in (("shift_4H.shift", -0.86, 0.3), ("shift_4He.shift", -0.64, 0.3),
                                                                       ("mix_1.amplitude", 0.0, 0.1))]), rtol=1e-10, atol=1e-14),
                  rows[-3:])

            print("\n3. the ties act")
            x = x0.copy()
            x[names.index("shift_4H.shift")] += 0.05
            moved = s.evaluate(x)
            chi0 = {d.name: d.chi2 for d in at.data}
            chi1 = {d.name: d.chi2 for d in moved.data}
            check("the 4H shift moves the n+t data set", abs(chi1["nt"] - chi0["nt"]) > 1e-6 * chi0["nt"], (chi0["nt"], chi1["nt"]))
            check("... and no data set of 4He", all(chi1[n] == chi0[n] for n in ("ptn", "ddp", "ddn")),
                  [(chi0[n], chi1[n]) for n in ("ptn", "ddp", "ddn")])

            print("\n4. gradient and Jacobian through the ties")
            value, grad = f.value_and_grad(x0)
            r0, J = f.residuals_and_jacobian(x0)
            worst_g, worst_j = 0.0, 0.0
            for a in range(len(x0)):
                h = 1e-6 * (abs(x0[a]) + 1.0)
                xp, xm = x0.copy(), x0.copy()
                xp[a] += h
                xm[a] -= h
                fd = (f(xp) - f(xm)) / (2 * h)
                worst_g = max(worst_g, abs(grad[a] - fd) / (abs(fd) + abs(grad[a]) + 1e-6 * abs(value) + 1.0))
                col = (f.residuals(xp) - f.residuals(xm)) / (2 * h)
                worst_j = max(worst_j, float(np.max(np.abs(J[:, a] - col) / (np.abs(col) + np.abs(J[:, a]) + 1e-6))))
            check(f"the gradient of every free parameter ({len(x0)}) agrees with finite differences", worst_g < 1e-4, worst_g)
            check("the Jacobian of the rows agrees with finite differences", worst_j < 1e-3, worst_j)
            check("2 J^T r == the gradient", np.allclose(2.0 * J.T @ r0, grad, rtol=1e-6, atol=1e-6 * abs(value)),
                  float(np.max(np.abs(2.0 * J.T @ r0 - grad))))

        print("\n5. Project edits")
        p = Project.open(path)
        try:
            p.add("isospin_reference", "second", nucleus="4He", threshold="p+3H")
            check("a second isospin_reference is refused", False, "no exception")
        except FormatError as e:
            check("a second isospin_reference is refused", "E_ISOSPIN_RECORD" in str(e), str(e))
        p = Project.open(path)
        p.remove("mix_1")
        p.remove("level_1")
        undeclared = p.save(os.path.join(work, "undeclared.azr"))
        try:
            with Session(undeclared, cwd=work):
                pass
            check("a level ID across nuclei without its isospin_level record is refused", False, "the session opened")
        except AZURE2Error as e:
            check("a level ID across nuclei without its isospin_level record is refused",
                  "has channels of two compound nuclei" in str(e), str(e)[-300:])
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
