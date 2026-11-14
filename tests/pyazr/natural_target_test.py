#!/usr/bin/env python3
"""Several compound nuclei in one project, through pyazr.

The models are tests/natural_target/models (format 3): A = p+12C (13N), B = p+13C (14N), M = both in one file with
a target `natural` (isotopes 12C 0.6, 13C 0.4) and extrapolations on it; M_data / M_samples_data the same with data
sets on one or two targets whose abundances are fitted (<target>.isotope.<nucleus>).

Checks:
  1. the file side: Project.levels() gives each level's compound nucleus (13N, 14N); the target's isotopes are
     record fields;
  2. merge_projects(A, B) (the engine's merge of the <levels> rows, experiment::MergeLevelRows) reproduces M through
     the engine: every extrapolation of the merged file equals the standalone project's, a copy of A's on the
     natural target is 0.6 sigma_A; two nuclei with the same J^pi set keep their J-groups apart, Project.levels()
     numbers them as the engine does, and save_fit writes a snapshot of them;
  3. two targets: four abundance parameters, the sum-to-one prior of each target, save_fit writes each target's
     fitted values;
  4. the runtime side: Pair.nucleus, Session.nuclei, the abundance parameters (kind, nucleus, prior sigma),
     terms["abundances"] against the engine's chi-squared;
  5. the abundance column of objective().residuals_and_jacobian against central differences;
  6. save_fit of a free abundance: the target's isotope field carries the fitted value, the .sav carries the
     parameter, the snapshot re-evaluates equal.

Run from anywhere:  python3 tests/pyazr/natural_target_test.py
"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
MODELS = os.path.join(ROOT, "tests", "natural_target", "models")
sys.path.insert(0, ROOT)

try:
    from pyazr import Project, Session
    from pyazr.system import merge_projects
    import pyazr._azure2  # noqa: F401
except ImportError as e:  # as the other pyazr tests: no package or engine module -> nothing to test
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(77)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def main():
    work = tempfile.mkdtemp(prefix="natural_target_")
    try:
        return run(work)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def copy_model(name, work, dest=None):
    path = os.path.join(work, dest or name)
    shutil.copy(os.path.join(MODELS, name), path)
    p = Project.open(path)
    p.set_config(output_dir="output/")
    return p.save()


def grids(path, cwd):
    with Session(path, cwd=cwd, extrapolations=True, use_brune=True) as s:
        x = np.asarray(s.parameters.x0, float)
        return [np.asarray(v, float) for v in [d.model for d in s.evaluate(x).data]]


def isotopes(project, target):
    """{nucleus: abundance} of a target record's isotope fields."""
    f = project[target]["isotope"]
    return {t.split()[0]: float(t.split()[1]) for t in (f if isinstance(f, list) else [f])}


def run(work):
    os.makedirs(os.path.join(work, "output"), exist_ok=True)
    # ---- 1. the file side ------------------------------------------------------
    print("1. file side")
    m = Project.open(os.path.join(MODELS, "M.azr"))
    nuclei = sorted({lv.nucleus for lv in m.levels()})
    check("two nuclei, 13N and 14N", nuclei == [(7, 13), (7, 14)], str(nuclei))
    check("pair keys per nucleus", {lv.nucleus: {c.pair for c in lv.channels} for lv in m.levels()}
          == {(7, 13): {1}, (7, 14): {2}})
    check("the target's isotopes", isotopes(m, "natural") == {"12C": 0.6, "13C": 0.4}, str(isotopes(m, "natural")))

    # ---- 2. merge_projects reproduces M through the engine ----------------------
    print("\n2. merge_projects against the engine")
    a_path, b_path = copy_model("A.azr", work), copy_model("B.azr", work)
    merged = os.path.join(work, "AB.azr")
    rep = merge_projects([a_path, b_path], merged, isotopes={"12C": 0.6, "13C": 0.4})
    check("pair keys offset", rep.pair_map == [{1: 1}, {2: 3}], str(rep.pair_map))
    check("the second x1 renamed", rep.record_map == [{}, {"x1": "x1_2"}], str(rep.record_map))
    ab = Project.open(merged)
    check("merged isotopes", isotopes(ab, "natural") == {"12C": 0.6, "13C": 0.4}, str(isotopes(ab, "natural")))
    check("merged nuclei", sorted({lv.nucleus for lv in ab.levels()}) == [(7, 13), (7, 14)])
    x1 = ab["x1"].fields
    ab.add("extrapolation", "x1_natural", **x1, target="natural")   # a copy of A's x1 on the natural target
    ab.save()
    out_a, out_b, out_m = grids(a_path, work), grids(b_path, work), grids(merged, work)
    check("merged 13N grid == A.azr", np.allclose(out_m[0], out_a[0], rtol=1e-10, atol=0),
          f"{np.max(np.abs(out_m[0] / out_a[0] - 1)):.2e}")
    check("merged 14N grid == B.azr", np.allclose(out_m[1], out_b[0], rtol=1e-10, atol=0),
          f"{np.max(np.abs(out_m[1] / out_b[0] - 1)):.2e}")
    check("plain natural copy of x1 == 0.6 sigma_A", np.allclose(out_m[2], 0.6 * out_a[0], rtol=1e-12, atol=0))

    # ---- 2b. two nuclei with the SAME J^pi set: J-group numbering per nucleus -----
    # The engine numbers J-groups per (J, parity, nucleus); a bookkeeping keyed on (J, parity) alone collapsed them,
    # so save_fit placed nothing (natural neon, 21Ne and 23Ne both half-integer, 2026-10-13).  A twin of A with a
    # 13-mass target has exactly A's J^pi list.
    print("\n2b. same J^pi in two nuclei")
    twin = Project.new(os.path.join(work, "A_twin.azr"), nuclear_part=a_path)   # A's levels, no records
    twin.set_pair(1, heavy_mass=13.00335)
    twin_path = twin.save(os.path.join(work, "A_twin.azr"))
    merged2 = os.path.join(work, "AA.azr")
    merge_projects([a_path, twin_path], merged2)
    aa = Project.open(merged2)
    check("twin forms a second nucleus", sorted({lv.nucleus for lv in aa.levels()}) == [(7, 13), (7, 14)])
    keys = {(lv.group, lv.level): lv for lv in aa.levels() if lv.active}
    with Session(merged2, cwd=work, extrapolations=True, use_brune=True) as s:
        bad = []
        for pp in s.parameters:
            if pp.kind != "energy":
                continue
            lv = keys.get((pp.jgroup, pp.level))
            if lv is None or abs(lv.energy - float(pp.level_energy or 0.0)) > 1e-6:
                bad.append((pp.name, pp.jgroup, pp.level, pp.level_energy, None if lv is None else lv.energy))
        check("Project.levels() numbers J-groups per nucleus like the engine", not bad, str(bad[:4]))
        ngroups = len({pp.jgroup for pp in s.parameters if pp.jgroup is not None and pp.kind in ("energy", "width")})
        check("twice A's J-groups", ngroups == 6, str(ngroups))
        try:
            s.save_fit(os.path.join(work, "AA_snapshot.azr"), np.asarray(s.parameters.x0, float))
            check("save_fit round trip with same-J^pi nuclei", True)
        except Exception as err:   # noqa: BLE001
            check("save_fit round trip with same-J^pi nuclei", False, str(err)[:300])

    # ---- 3. two targets: several abundance sets in one project ----------------
    print("\n3. two targets")
    msd = copy_model("M_samples_data.azr", work)
    e_lab = 1.0 + 0.002 * np.arange(out_a[0].size)
    for f, (wa, wb) in (("data.dat", (0.6, 0.4)), ("data_enr.dat", (0.1, 0.9))):
        y = wa * out_a[0] + wb * out_b[0]
        with open(os.path.join(work, f), "w") as fh:
            for e, v in zip(e_lab, y):
                fh.write(f"{e:.6f} 0 {v:.12e} {0.001 * v:.12e}\n")
    with Session(msd, cwd=work, use_brune=True) as s:
        ab = s.parameters.of_kind("abundance")
        check("four abundance parameters, target by target",
              [p.name for p in ab] == ["natural.isotope.12C", "natural.isotope.13C", "enr.isotope.12C",
                                       "enr.isotope.13C"], str([(p.name, p.sample_name, p.sample_index) for p in ab]))
        x0 = np.asarray(s.parameters.x0, float)
        pen = s.evaluate(x0).terms["abundances"]
        # both targets at 0.5 + 0.5 = 1: no sum penalty, no priors
        check("penalties: two sum-to-one terms at zero", abs(pen) < 1e-12, str(pen))
        x1 = x0.copy()
        x1[ab[2].free_index] = 0.6    # enr sums to 1.1: ((0.1)/0.01)^2 = 100
        pen1 = s.evaluate(x1).terms["abundances"]
        check("penalties: the enr target's own sum width", abs(pen1 - 100.0) < 1e-6, str(pen1))
        azr_path, _ = s.save_fit(os.path.join(work, "M_samples_snapshot.azr"), x1)
        back = Project.open(azr_path)
        check("save_fit writes each target's fitted values",
              abs(isotopes(back, "enr")["12C"] - 0.6) < 1e-12 and isotopes(back, "natural")["12C"] == 0.5,
              (isotopes(back, "enr"), isotopes(back, "natural")))

    # ---- 4. the runtime side ---------------------------------------------------
    print("\n4. runtime")
    mfile = copy_model("M_data.azr", work, "M.azr")
    mm = Project.open(mfile)
    mm.set("natural", isotope=["12C 0.6 fit", "13C 0.4 fit prior 0.4 0.05"])
    mm.save()
    # data = 0.7 sigma_A + 0.3 sigma_B on the 2 keV lab grid of x1 (extrapolation energies are lab)
    y = 0.7 * out_a[0] + 0.3 * out_b[0]
    with open(os.path.join(work, "data.dat"), "w") as f:
        for e, v in zip(e_lab, y):
            f.write(f"{e:.6f} 0 {v:.12e} {0.01 * v:.12e}\n")
    with Session(mfile, cwd=work, use_brune=True) as s:
        pairs = s.pairs
        check("Pair.nucleus", [p.nucleus for p in pairs] == [(7, 13), (7, 14)], str([p.nucleus for p in pairs]))
        check("Session.nuclei", sorted(s.nuclei) == [(7, 13), (7, 14)] and [p.key for p in s.nuclei[(7, 14)]] == [2])
        ab_params = s.parameters.of_kind("abundance")
        check("two abundance parameters after the level block",
              len(ab_params) == 2 and all(p.kind == "abundance" for p in ab_params)
              and ab_params[0].index == sum(1 for p in s.parameters if p.kind in ("energy", "width")) and ab_params[0].nucleus == (7, 13) and ab_params[1].nucleus == (7, 14),
              str(ab_params))
        check("prior sigma and values carried", ab_params[1].prior_sigma == 0.05 and ab_params[0].value == 0.6
              and ab_params[1].value == 0.4 and not ab_params[0].fixed and not ab_params[1].fixed, str(ab_params))
        check("Pair.abundance", [p.abundance for p in pairs] == [0.6, 0.4])
        x = np.asarray(s.parameters.x0, float)
        ia, ib = ab_params[0].free_index, ab_params[1].free_index
        # penalties: at the nominal values the prior is 0 and the sum term 0
        pen = s.evaluate(x)
        check("penalty zero at the nominal abundances", pen.terms["abundances"] == 0.0, str(pen.terms["abundances"]))
        x1 = x.copy()
        x1[ib] = 0.45
        pen1 = s.evaluate(x1).terms["abundances"]
        want = ((0.45 - 0.4) / 0.05) ** 2 + ((0.6 + 0.45 - 1.0) / 0.01) ** 2
        check("penalty formula", abs(pen1 - want) < 1e-9, f"{pen1} vs {want}")
        # the engine's data chi-squared and its abundance term against a pure-data chi-squared formed here from the
        # model values
        ev1 = s.evaluate(x1)
        vals = np.asarray(ev1.data[0].model, float)
        data_chi2 = float(np.sum(((vals - y) / (0.01 * y)) ** 2))
        chi1 = ev1.chi2 + ev1.terms["abundances"]
        check("engine chi2 + its abundance term = data chi2 + abundance penalty", abs(chi1 - (data_chi2 + want)) < 1e-6 * max(1.0, chi1),
              f"{chi1} vs {data_chi2 + want}")
        check("model at the nominal abundances is 0.6 A + 0.4 B",
              np.allclose(np.asarray(s.evaluate(x).data[0].model, float), 0.6 * out_a[0] + 0.4 * out_b[0], rtol=1e-10))

        # ---- 5. the abundance Jacobian column -----------------------------------
        print("\n5. Jacobian")
        n = s.evaluate(x).data[-1].rows.stop   # the data rows; the abundance prior and sum rows follow them
        r, J = s.objective().residuals_and_jacobian(x)
        J = np.asarray(J, float)[:n]
        h = 1e-4
        for col, name in ((ia, "12C"), (ib, "13C")):
            xp = x.copy()
            xp[col] += h
            xm = x.copy()
            xm[col] -= h
            fd = (s.objective().residuals(xp)[:n] - s.objective().residuals(xm)[:n]) / (2 * h)
            an = J[:, col]
            rel = np.max(np.abs(an - fd)) / max(np.max(np.abs(fd)), 1e-300)
            check(f"abundance column of {name} vs central differences", rel < 1e-5, f"rel {rel:.2e}")
        # the data rows carry the model's derivative only (the penalty has rows of its own after them)
        check("abundance column is the nucleus's term over the error",
              np.allclose(J[:, ia], out_a[0] / (0.01 * y), rtol=1e-8),
              f"{np.max(np.abs(J[:, ia] * 0.01 * y / out_a[0] - 1)):.2e}")

        # ---- 6. save_fit round trip ------------------------------------------------
        print("\n6. save_fit")
        xf = x.copy()
        xf[ia] = 0.7
        xf[ib] = 0.3
        # the data part only: the snapshot's isotope fields move to the fitted values, so its prior is 0
        chif = s.evaluate(xf).chi2
        snap, sav = s.save_fit(os.path.join(work, "snap.azr"), xf, verify=False)
        got = isotopes(Project.open(snap), "natural")
        check("the snapshot's isotopes carry the fitted values", got == {"12C": 0.7, "13C": 0.3}, str(got))
        check("param.sav carries the abundances", ab_params[1].name in open(sav).read())
    with Session(snap, cwd=work, use_brune=True) as s2:
        x2 = np.asarray(s2.parameters.x0, float)
        chi2 = s2.evaluate(x2).chi2
        check("snapshot re-evaluates to the same data chi-squared", abs(chi2 - chif) < 1e-6 * max(1.0, chif),
              f"{chi2} vs {chif}")

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
