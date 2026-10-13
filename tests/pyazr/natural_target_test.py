#!/usr/bin/env python3
"""Several compound nuclei in one project, through pyazr.

The models are tests/natural_target/models (make_models.py): A = p+12C (13N),
B = p+13C (14N), M = both in one file with <abundances> 0.6 / 0.4 and a
"natural 1" SUM segment.

Checks:
  1. the file side: AzrModel.nuclei(), abundances, set_abundance, the
     abundanceSumSigma key, add_data_segment(natural=True), set_natural (and its
     RATIO refusal), SegmentSet / TestSegmentSet report the flag;
  2. merge_projects(A, B) reproduces M through the engine: every test segment
     of the merged file equals the standalone project's, and the natural SUM is
     0.6 sigma_A + 0.4 sigma_B;
  3. the runtime side: Pair.nucleus, azure2.nuclei, the abundance parameters
     (kind, nucleus, prior sigma), penalties()["abundance"] against the
     engine's chi-squared;
  4. the abundance column of residual_jacobian against central differences;
  5. save_fit of a free abundance: the <abundances> line carries the fitted
     value, the .sav carries abundance_7_14, the snapshot re-evaluates equal.

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

from pyazr.azrfile import AzrModel                        # noqa: E402
from pyazr.datasets import SegmentSet, TestSegmentSet     # noqa: E402
from pyazr.nuclei import Nucleus, merge_projects          # noqa: E402

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


def run(work):
    # ---- 1. the file side ------------------------------------------------------
    print("1. file side")
    m = AzrModel.from_file(os.path.join(MODELS, "M.azr"))
    nuc = m.nuclei()
    check("two nuclei, 13N and 14N", [n.label for n in nuc] == ["13N", "14N"], str(nuc))
    check("pair keys per nucleus", nuc[Nucleus(7, 13)] == [1] and nuc[Nucleus(7, 14)] == [2], str(nuc))
    check("abundances read", m.abundances[Nucleus(7, 13)] == (0.6, False, 0.0) and m.abundances[Nucleus(7, 14)] == (0.4, False, 0.0),
          str(m.abundances))
    m.set_abundance("14N", 0.35, True, 0.02)
    check("set_abundance edits the line", m.abundances[Nucleus(7, 14)] == (0.35, True, 0.02), str(m.abundances))
    m.abundance_sum_sigma = 0.01
    check("abundanceSumSigma written and read", m.abundance_sum_sigma == 0.01)
    m.abundance_sum_sigma = 0
    check("abundanceSumSigma removed", m.abundance_sum_sigma == 0.0)
    a = AzrModel.from_file(os.path.join(MODELS, "A.azr"))
    check("a single-nucleus project has no block", not a.has_abundance_block and a.abundances[Nucleus(7, 13)] == (1.0, False, 0.0))
    a.set_abundance((7, 13), 0.9, False, 0.0)
    check("set_abundance creates the block after </levels>",
          a.has_abundance_block and a._suffix.splitlines()[:2] == ["</levels>", "<abundances>"], a._suffix[:60])
    try:
        a.set_abundance("23Ne", 0.1)
        check("unknown nucleus refused", False)
    except KeyError:
        check("unknown nucleus refused", True)
    m.add_data_segment("data.dat", 1, 1, "angle-integrated", 1.0, 3.0, natural=True)
    segs = SegmentSet.from_file(m.write(os.path.join(work, "m1.azr")))
    check("add_data_segment(natural=True) writes the token", segs[-1].natural and segs[-1].data_file == "data.dat", str(segs[-1]))
    tsegs = TestSegmentSet.from_file(os.path.join(work, "m1.azr"))
    check("TestSegmentSet reads natural", [t.natural for t in tsegs] == [False, False, True, True], str([t.natural for t in tsegs]))
    check("set_test_natural on and off", m.set_test_natural(1, True) == 1 and
          TestSegmentSet.from_file(m.write(os.path.join(work, "m2.azr")))[0].natural and
          m.set_test_natural(1, False) == 1 and
          not TestSegmentSet.from_file(m.write(os.path.join(work, "m3.azr")))[0].natural)
    r = AzrModel.from_file(os.path.join(MODELS, "bad_ratio.azr"))
    try:
        r.set_test_natural(1, True)
        check("set_natural refuses a RATIO", False)
    except ValueError:
        check("set_natural refuses a RATIO", True)

    try:
        from pyazr import azure2
        import pyazr._azure2  # noqa: F401
    except ImportError as e:
        print(f"pyazr engine module not built ({e}); engine checks skipped")
        return 1 if failures else 0

    # ---- 2. merge_projects reproduces M through the engine ----------------------
    print("\n2. merge_projects against the engine")
    merged = os.path.join(work, "AB.azr")
    rep = merge_projects([os.path.join(MODELS, "A.azr"), os.path.join(MODELS, "B.azr")], merged,
                         {"13N": 0.6, "14N": 0.4})
    check("pair keys offset", rep.pair_map == [{1: 1}, {2: 3}], str(rep.pair_map))
    check("test segments renumbered", rep.test_map == [{1: 1}, {1: 2}], str(rep.test_map))
    text = open(merged).read()
    # the GUI's readers reject an empty <lastRun> or <mcmc> block (natural-neon merge, 2026-10-13)
    check("no empty <lastRun>/<mcmc> block in the merged file", "<lastRun>\n</lastRun>" not in text and "<mcmc>" not in text)
    ab = AzrModel.from_file(merged)
    check("merged abundances", ab.abundances == {Nucleus(7, 13): (0.6, False, 0.0), Nucleus(7, 14): (0.4, False, 0.0)}, str(ab.abundances))
    ab.add_extrapolation(1, 1, 1.0, 3.0, 0.002, "angle-integrated", natural=False)   # a copy of S1, plain
    ab.set_test_natural(3, True)
    ab.write(merged)
    os.makedirs(os.path.join(work, "output"), exist_ok=True)

    def grids(path, cwd):
        with azure2(path, cwd=cwd, data_mode=False, use_brune=True) as s:
            x = np.asarray(s.params_rwa, float)
            return [np.asarray(v, float) for v in s.calculate_rwa(x)], s
    out_a, _ = grids(os.path.join(MODELS, "A.azr"), work)
    out_b, _ = grids(os.path.join(MODELS, "B.azr"), work)
    out_m, _ = grids(merged, work)
    check("merged 13N grid == A.azr", np.allclose(out_m[0], out_a[0], rtol=1e-10, atol=0), f"{np.max(np.abs(out_m[0] / out_a[0] - 1)):.2e}")
    check("merged 14N grid == B.azr", np.allclose(out_m[1], out_b[0], rtol=1e-10, atol=0), f"{np.max(np.abs(out_m[1] / out_b[0] - 1)):.2e}")
    check("plain natural copy of S1 == 0.6 sigma_A", np.allclose(out_m[2], 0.6 * out_a[0], rtol=1e-12, atol=0))

    # ---- 2b. two nuclei with the SAME J^pi set: J-group numbering per nucleus -----
    # The engine numbers J-groups per (J, parity, nucleus); engine_level_keys() keyed on
    # (J, parity) alone collapsed them, so apply_fit/save_fit placed nothing (natural neon,
    # 21Ne and 23Ne both half-integer, 2026-10-13).  A twin of A with a 13-mass target has
    # exactly A's J^pi list.
    print("\n2b. same J^pi in two nuclei")
    twin = os.path.join(work, "A_twin.azr")
    with open(os.path.join(MODELS, "A.azr")) as f, open(twin, "w") as g:
        for ln in f:
            g.write(ln.replace(" 1.00782 12.0 1 6 ", " 1.00782 13.00335 1 6 ") if ln.strip() and ln.split()[0] != "<" and len(ln.split()) > 20 else ln)
    merged2 = os.path.join(work, "AA.azr")
    merge_projects([os.path.join(MODELS, "A.azr"), twin], merged2, {"13N": 0.5, "14N": 0.5})
    aa = AzrModel.from_file(merged2)
    check("twin forms a second nucleus", sorted(n.label for n in aa.nuclei()) == ["13N", "14N"], str(aa.nuclei()))
    lvlmap = aa.engine_level_keys()
    with azure2(merged2, cwd=work, data_mode=False, use_brune=True) as s:
        bad = []
        for pp in s.parameters:
            if pp.kind != "energy":
                continue
            lv = lvlmap.get((pp.jgroup, pp.level))
            if lv is None or abs(lv.energy - float(pp.level_energy or 0.0)) > 1e-6:
                bad.append((pp.name, pp.jgroup, pp.level, pp.level_energy, None if lv is None else lv.energy))
        check("engine_level_keys numbers J-groups per nucleus like the engine", not bad, str(bad[:4]))
        ngroups = len({pp.jgroup for pp in s.parameters if pp.jgroup is not None and pp.kind in ("energy", "width")})
        check("twice A's J-groups", ngroups == 6, str(ngroups))
        try:
            s.save_fit(os.path.join(work, "AA_snapshot.azr"), np.asarray(s.params_rwa, float))
            check("save_fit round trip with same-J^pi nuclei", True)
        except Exception as err:   # noqa: BLE001
            check("save_fit round trip with same-J^pi nuclei", False, str(err)[:300])

    # ---- 3. the runtime side ---------------------------------------------------
    print("\n3. runtime")
    mfile = os.path.join(work, "M.azr")
    mm = AzrModel.from_file(os.path.join(MODELS, "M.azr"))
    mm.set_abundance("14N", 0.4, True, 0.05)
    mm.set_abundance("13N", 0.6, True, 0.0)
    mm.abundance_sum_sigma = 0.01
    mm.add_data_segment("data.dat", 1, 1, "angle-integrated", 1.0, 3.0, natural=True)
    # data = 0.7 sigma_A + 0.3 sigma_B on the 2 keV lab grid of S1 (test energies are lab)
    e_lab = 1.0 + 0.002 * np.arange(out_a[0].size)
    y = 0.7 * out_a[0] + 0.3 * out_b[0]
    with open(os.path.join(work, "data.dat"), "w") as f:
        for e, v in zip(e_lab, y):
            f.write(f"{e:.6f} 0 {v:.12e} {0.01 * v:.12e}\n")
    # the data line must be a SUM over both nuclei: rewrite it with the component
    lines = mm._suffix.splitlines()
    i = max(k for k, ln in enumerate(lines) if "data.dat" in ln)
    lines[i] = lines[i].replace("data.dat  0  0  natural  1", "data.dat 1 0 1 2 2 -999 natural 1")
    mm._suffix = "\n".join(lines) + "\n"
    mm.write(mfile)
    segs = SegmentSet.from_file(mfile)
    check("data segment is a natural SUM", segs[0].natural and segs[0].operation == "sum" and segs[0].components == ((2, 2, -999.0, 1.0),), str(segs[0]))
    with azure2(mfile, cwd=work, use_brune=True) as s:
        pairs = s.pairs
        check("Pair.nucleus", [p.nucleus for p in pairs] == [(7, 13), (7, 14)], str([p.nucleus for p in pairs]))
        check("azure2.nuclei", sorted(s.nuclei) == [(7, 13), (7, 14)] and [p.key for p in s.nuclei[(7, 14)]] == [2])
        ab_params = s.parameters.abundances
        check("two abundance parameters after the level block",
              len(ab_params) == 2 and all(p.kind == "abundance" for p in ab_params)
              and ab_params[0].index == s.n_rmatrix and ab_params[0].nucleus == (7, 13) and ab_params[1].nucleus == (7, 14),
              str(ab_params))
        check("prior sigma and values carried", ab_params[1].prior_sigma == 0.05 and ab_params[0].value == 0.6 and ab_params[1].value == 0.4
              and not ab_params[0].fixed and not ab_params[1].fixed, str(ab_params))
        check("Pair.abundance", [p.abundance for p in pairs] == [0.6, 0.4])
        x = np.asarray(s.params_rwa, float)
        ia, ib = ab_params[0].free_index, ab_params[1].free_index
        # penalties: at the nominal values the prior is 0 and the sum term 0
        pen = s.penalties(x)
        check("penalty zero at the nominal abundances", pen["abundance"] == 0.0, str(pen["abundance"]))
        chi0 = float(s.calculate_chi2_rwa(x)[0])
        x1 = x.copy(); x1[ib] = 0.45
        pen1 = s.penalties(x1)["abundance"]
        want = ((0.45 - 0.4) / 0.05) ** 2 + ((0.6 + 0.45 - 1.0) / 0.01) ** 2
        check("penalty formula", abs(pen1 - want) < 1e-9, f"{pen1} vs {want}")
        # the engine adds the same penalty inside calculate_chi2_rwa: compare the jump against a
        # pure-data evaluation obtained from calculate_rwa
        vals = np.asarray(s.calculate_rwa(x1)[0], float)
        data_chi2 = float(np.sum(((vals - y) / (0.01 * y)) ** 2))
        chi1 = float(s.calculate_chi2_rwa(x1)[0])
        check("engine chi2 = data chi2 + abundance penalty", abs(chi1 - (data_chi2 + want)) < 1e-6 * max(1.0, chi1), f"{chi1} vs {data_chi2 + want}")
        check("model at the nominal abundances is 0.6 A + 0.4 B", np.allclose(np.asarray(s.calculate_rwa(x)[0], float), 0.6 * out_a[0] + 0.4 * out_b[0], rtol=1e-10))

        # ---- 4. the abundance Jacobian column -----------------------------------
        print("\n4. Jacobian")
        r, J = s.residual_jacobian(x)
        J = np.asarray(J, float)
        h = 1e-4
        for col, name in ((ia, "13N"), (ib, "14N")):
            xp = x.copy(); xp[col] += h
            xm = x.copy(); xm[col] -= h
            fd = (np.asarray(s.residual_jacobian(xp)[0], float) - np.asarray(s.residual_jacobian(xm)[0], float)) / (2 * h)
            an = J[:, col]
            rel = np.max(np.abs(an - fd)) / max(np.max(np.abs(fd)), 1e-300)
            check(f"abundance column of {name} vs central differences", rel < 1e-5, f"rel {rel:.2e}")
        # the residual Jacobian must not carry the penalty (a residual row would), only the model's derivative
        check("abundance column is the nucleus's term over the error",
              np.allclose(J[:, ia], out_a[0] / (0.01 * y), rtol=1e-8), f"{np.max(np.abs(J[:, ia] * 0.01 * y / out_a[0] - 1)):.2e}")

        # ---- 5. save_fit round trip ------------------------------------------------
        print("\n5. save_fit")
        xf = x.copy(); xf[ia] = 0.7; xf[ib] = 0.3
        # the data part only: the snapshot's <abundances> nominal moves to the fitted value, so its prior is 0
        chif = float(s.calculate_chi2_rwa(xf)[0]) - s.penalties(xf)["abundance"]
        snap, sav = s.save_fit(os.path.join(work, "snap.azr"), xf, verify=False)
        snapm = AzrModel.from_file(snap)
        check("snapshot <abundances> carries the fitted values",
              snapm.abundances[Nucleus(7, 13)][0] == 0.7 and snapm.abundances[Nucleus(7, 14)][0] == 0.3, str(snapm.abundances))
        check("param.sav carries the abundances", "abundance_7_14" in open(sav).read())
    with azure2(snap, cwd=work, use_brune=True) as s2:
        x2 = np.asarray(s2.params_rwa, float)
        chi2 = float(s2.calculate_chi2_rwa(x2)[0]) - s2.penalties(x2)["abundance"]
        check("snapshot re-evaluates to the same data chi-squared", abs(chi2 - chif) < 1e-6 * max(1.0, chif), f"{chi2} vs {chif}")

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
