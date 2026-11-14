#!/usr/bin/env python3
"""pyazr.endf.write (docs/upgrade/endf.md sections 4, 11, 13.5-13.8): T1, T2 and 13N written as ENDF-6 and read back
(MF=1, MF=2, MF=3, MF=6).

For each project:
(a) the module's own format check of the written file: no error, one warning (LRP=2 in a primary evaluation, p. 42);
(b) MF=3 read back with the module's reader at the data's lab energies and at the midpoints of the written grid
    equals the engine's sigma there to the grid's tolerance (linear interpolation, INT=2);
(c) MF=2 read back (resonances(), to_project) gives the session's levels again (energies and particle widths to
    1e-7), and the MF=3 sections carry the check the writer made (MF=2 reproduces MF=3 to 1e-5);
(d) the list of what is not carried against the file's sections and the project; the MF=1 text names it;
(e) (after the loop) a narrow level is found by the start grid;
(f) MF=6 read back; Material.curves' dsigma/dOmega of the file equals the engine's at the points of the project's
    differential data sets and between the tabulated energies and angles (elastic LAW=5, exits LAW=2, capture final
    states MF=6/MT=102).
With AZURE2_ENDF_UTILITIES and AZURE2_PREPRO set, CHECKR, FIZCON and PSYCHE run on each written file and every
finding is printed (not judged here: the report of the branch reads them).
Exits 77 (skipped) without the compiled engine module.  Run from anywhere:  python3 tests/pyazr/endf_write_test.py
"""
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    import numpy as np
    import pyazr
    from pyazr import _azure2  # noqa: F401
    import pyazr.endf as E
    from pyazr.endf.writer import _pair_keys
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)

warnings.simplefilter("ignore")
PROJECTS = [os.path.join(ROOT, "tests", "upgrade", "T1", "T1.azr"), os.path.join(ROOT, "tests", "upgrade", "T2",
                                                                                  "T2.azr"),
            os.path.join(ROOT, "tests", "13N", "13N.azr")]
TOLERANCE = 1e-3
TMP = tempfile.mkdtemp(prefix="endf_write_")
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:600]}"))
    if not ok:
        failures.append(name)


def engine_dsdo(session, label, e_lab_MeV, angles_cm, folder, pair_keys):
    """The engine's dsigma/dOmega (b/sr, c.m.) of reaction `label` at lab energies x c.m. angles, (n_E, n_angles)."""
    p = pyazr.Project.new(os.path.join(folder, "probe6.azr"), nuclear_part=session.file)
    e_lab_MeV = np.asarray(e_lab_MeV, float)
    np.savetxt(os.path.join(folder, "probe6.dat"), e_lab_MeV, fmt="%.17g")
    p.add_extrapolation("probe", reaction=label, quantity="dsigma_domega", unit="b/sr",
                        energies=f"file {os.path.join(folder, 'probe6.dat')} MeV lab",
                        angles="list " + " ".join(repr(float(a)) for a in angles_cm) + " deg cm", **pair_keys)
    p.save()
    pf = os.path.join(folder, "probe.par")
    P = session.parameters
    with open(pf, "w") as fh:
        fh.write(f"#parametrization {session.basis}\n")
        for q in P:
            fh.write(f"{q.name} {float(q.start if q.fixed else P.x0[q.free_index])!r}\n")
    s = pyazr.Session(p.path, extrapolations=True, parameter_file=pf)
    try:
        d = s.evaluate(s.parameters.x0).data[0]
        return np.asarray(d.model, float).reshape(len(angles_cm), len(e_lab_MeV)).T, np.asarray(d.energies, float)
    finally:
        s.close()


def file_dsdo(m, session, entrance, exit_pair, e_cm, e_lab, angles):
    """The written file's dsigma/dOmega at the points (the module's reader, Material.curves)."""
    by_key = {p.key: p for p in session.pairs}
    req = E.CurveRequest("dsigma_domega", tuple(e_cm), tuple(angles), E.PairSpec.of(by_key[entrance]),
                         E.PairSpec.of(by_key[exit_pair]), tuple(e_lab))
    return np.asarray(m.curves(req).values, float)


def engine_sigma(session, label, e_lab_MeV, folder, pair_keys):
    """The engine's sigma (b) of reaction `label` at lab energies, at the session's start (an extrapolation)."""
    p = pyazr.Project.new(os.path.join(folder, "probe.azr"), nuclear_part=session.file)
    np.savetxt(os.path.join(folder, "probe.dat"), e_lab_MeV, fmt="%.17g")
    p.add_extrapolation("probe", reaction=label, quantity="sigma", unit="b",
                        energies=f"file {os.path.join(folder, 'probe.dat')} MeV lab", **pair_keys)
    p.save()
    pf = os.path.join(folder, "probe.par")
    P = session.parameters
    with open(pf, "w") as fh:
        fh.write(f"#parametrization {session.basis}\n")
        for q in P:
            fh.write(f"{q.name} {float(q.start if q.fixed else P.x0[q.free_index])!r}\n")
    s = pyazr.Session(p.path, extrapolations=True, parameter_file=pf)
    try:
        return np.asarray(s.evaluate(s.parameters.x0).data[0].model, float)
    finally:
        s.close()


for azr in PROJECTS:
    name = os.path.basename(azr)
    print(name)
    d = os.path.join(TMP, name.replace(".azr", ""))
    shutil.copytree(os.path.dirname(azr), d, ignore=shutil.ignore_patterns("output*", "checks*"))
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)
    cwd = os.getcwd()
    os.chdir(d)
    try:
        s = pyazr.Session(name)
        out = os.path.join(d, name.replace(".azr", ".endf"))
        w = E.write(s, out, header={"AUTH": "pyazr tests", "ALAB": "TEST"}, tolerance=TOLERANCE)
        print("\n".join("      " + line for line in str(w).splitlines()))
        # (a)
        found = E.check_format(out).findings
        errors = [f for f in found if f.severity == "error"]
        warns = [f.code for f in found if f.severity == "warning"]
        check(f"{name}: own format check: no error, the one warning LRP=2", not errors and warns == ["LRP"],
              [(f.severity, f.code, f.message) for f in found])
        # (b)
        m = E.read(out)
        project = pyazr.Project.open(azr)
        entrance = next(p.key for p in s.pairs if p.is_entrance)
        lab = np.unique(np.concatenate([np.asarray(x.extra["energy_lab"], float)
                                        for x in s.evaluate(s.parameters.x0).data]))
        for mf, mt, _ in w.sections:
            if mf != 3 or mt == 2:  # MT=2 holds sigma_NI of LAW=5, checked through (f)
                continue
            xs = m.cross_section(mt)
            grid = np.asarray(xs.table.x) / 1e6
            sigma_grid = {mt: np.asarray(xs.table.y)}
            mid = 0.5 * (grid[1:] + grid[:-1])
            exit_key = -1 if mt == 102 else next(
                p.key for p, q in zip([p for p in s.pairs if p.pair_kind == "particle" and p.key != entrance],
                                      [r for r in E.read(out).resonances().pairs[1:]]) if q.mt == mt)
            label = project.reaction(entrance, exit_key)
            keys = _pair_keys(project, entrance, exit_key)
            for what, e in (("the data's lab energies", lab[(lab >= grid[0]) & (lab <= grid[-1])]),
                            ("the midpoints of the grid", mid)):
                ref = engine_sigma(s, label, e, d, keys)
                got = xs(e * 1e6)
                # the grid's rule (pyazr.refine_grid, the engine's adaptive grid): relative to the value, not below
                # 1e-3 of the largest value of the range
                floor = 1e-3 * np.max(np.abs(sigma_grid[mt]))
                dev = np.abs(got - ref) / np.maximum(np.abs(ref), floor)
                big = np.abs(ref) >= floor
                rel = np.abs(got[big] - ref[big]) / np.abs(ref[big])
                print(f"      MT={mt} at {what}: {len(e)} points, max |dev|/max(|sigma|, 1e-3 max) {dev.max():.2e}; "
                      f"relative where sigma >= 1e-3 max: max {rel.max():.2e}, median {np.median(rel):.1e}")
                check(f"{name}: MF=3/MT={mt} at {what} = the engine to {1.5 * TOLERANCE:g} (the grid's rule)",
                      dev.max() < 1.5 * TOLERANCE, dev.max())
        # (c)
        exits = [mt for mf, mt, _ in w.sections if mf == 3 and mt not in (2, 102)]
        check(f"{name}: the writer's check: MF=2 reproduces every particle reaction of MF=3 ({exits}) to 1e-5",
              sorted(w.mf2_check) == sorted(exits) and all(v < 1e-5 for v in w.mf2_check.values()), w.mf2_check)
        back = m.resonances()
        ent = next(p.to_dict() for p in s.pairs if p.is_entrance)
        new, conv = E.to_project(back, os.path.join(d, "back.azr"),
                                 separation_energy=ent["sep_energy"] + ent["excitation"])
        orig = {}
        for lv in project.levels():
            orig.setdefault(lv.jpi, []).append(lv.energy)
        got = {}
        for lv in conv.levels:
            got.setdefault(lv.jpi, []).append(lv.energy_MeV)
        worst = max(abs(a - b) / max(abs(a), 1.0) for k in got for a, b in zip(sorted(orig[k]), sorted(got[k])))
        check(f"{name}: the levels read back from the file: energies to 1e-7 MeV", worst < 1e-7, (orig, got))
        # (d) the list of what is not carried, against the file's own reading: every section a reaction of the
        # project needs (elastic MT=2, each exit, capture MT=102; MF=3 and MF=6) is in the file or named in
        # not_written, never both; the losses name exactly what the project holds (photon pairs, data sets, THM)
        text = " ".join(" ".join(m.header.text).split())
        sections = set(m.sections)
        mts = [2] + [r.mt for r in back.pairs[1:]] + ([102] if any(p.pair_kind == "photon" for p in s.pairs) else [])
        named = {(mf, mt) for mf in (3, 6) for mt in mts if any(t.startswith(f"MF={mf}/MT={mt} ") for t in
                                                                  w.not_written)}
        both = [(mf, mt) for mf in (3, 6) for mt in mts if ((mf, mt) in sections) == ((mf, mt) in named)]
        check(f"{name}: every section the reactions {mts} need is written or named as not written, never both",
              not both and len(w.not_written) == len(named), (both, w.not_written))
        photons = any(p.pair_kind == "photon" for p in s.pairs)
        expect = [photons, True, bool((project.thm() or {}).get("experiments"))]
        got_items = [any("photon channels" in t for t in w.losses),
                     any(t.startswith(f"the {len(s.datasets)} data set(s)") for t in w.losses),
                     any("Trojan-horse" in t for t in w.losses)]
        check(f"{name}: the losses name the photon channels, the {len(s.datasets)} data sets and THM exactly when "
              f"the project has them ({expect})", got_items == expect and len(w.losses) == sum(expect), w.losses)
        check(f"{name}: MF=1 names LRP=2's statement and every loss and section not written",
              "MF=2 (LRF=7, KRM=4, IFG=1) reproduces" in text and "Not carried by this file" in text and
              all(" ".join(("- " + t).split()) in text for t in w.losses + w.not_written), text[-800:])
        check(f"{name}: MF=1: LRP=2, NSUB of the projectile, EMAX the upper energy", m.header.lrp == 2 and
              m.header.nsub == 10 * (1000 * int(ent["Z1"]) + round(ent["M1"])), (m.header.lrp, m.header.nsub))
        # (f) MF=6: every section read back with the module's reader; the file's dsigma/dOmega (Material.curves)
        # against the engine's (i) at the points of the project's differential data sets and (ii) between the
        # tabulated energies and angles: 60 midpoints of the energy grid x 16 c.m. angles 20-170 deg, for the
        # elastic channel (LAW=5), each exit (LAW=2) and each capture final state (MF=6/MT=102 photon)
        for mf, mt, _ in w.sections:
            if mf == 6:
                dist = m.distribution(mt)
                check(f"{name}: MF=6/MT={mt} read back: LCT=2, products {[(p.zap, p.law) for p in dist.products]}",
                      dist.lct == 2 and all(p.law in (2, 4, 5) for p in dist.products))
        targets = []  # (what, exit pair key, label, pair keys)
        targets.append(("elastic MT=2", entrance, project.reaction(entrance, entrance),
                        _pair_keys(project, entrance, entrance)))
        for p in s.pairs:
            if p.key != entrance:
                targets.append((f"{'capture to ' if p.pair_kind == 'photon' else ''}pair {p.key}", p.key,
                                project.reaction(entrance, p.key), _pair_keys(project, entrance, p.key)))
        ev = s.evaluate(s.parameters.x0)
        grid_e = np.asarray(m.cross_section(next(mt for mf, mt, _ in w.sections if mf == 3)).table.x) / 1e6
        mid = 0.5 * (grid_e[1:] + grid_e[:-1])
        mid = mid[np.linspace(0, len(mid) - 1, 60).astype(int)]
        # inside the LAW=5 table: above the angle of its largest cosine the file gives the Coulomb term alone
        angles = np.arange(20.0, 171.0, 10.0)
        for what, key, label, keys in targets:
            # the file holds dsigma/dOmega of LAW=2 as a product of two tables, sigma and p, each interpolated to the
            # grid's rule: 3x the tolerance (measured 1.1e-3 to 3.0e-3); LAW=5 as sigma_c + sigma_NI P_NI with
            # sigma_e - sigma_c interpolated to the rule against the cross section: 5x (measured 2.5e-4 to 3.4e-3)
            bound = (5 if key == entrance else 3) * TOLERANCE
            # (i) the data's points: every dsigma_domega data set of this reaction
            worst_data, n_data = 0.0, 0
            described = {rec.name: rec for rec in s.datasets}
            for dd in ev.data:
                rec = described.get(dd.name)
                if dd.quantity != "dsigma_domega" or rec is None or rec.entrance_key != entrance or \
                        rec.exit_key != key:
                    continue
                lab = np.asarray(dd.extra["energy_lab"], float)
                inside = (lab >= grid_e[0]) & (lab <= grid_e[-1]) & \
                    (np.cos(np.radians(np.asarray(dd.angles, float))) <= E.writer.ELASTIC_MAX_COSINE)
                if not inside.any():
                    continue
                ang = np.asarray(dd.angles, float)
                ref = np.empty(int(inside.sum()))
                for a in np.unique(ang[inside]):  # one engine session per angle of the data set
                    sel = ang[inside] == a
                    ref[sel] = engine_dsdo(s, label, lab[inside][sel], [a], d, keys)[0][:, 0]
                got = file_dsdo(m, s, entrance, key, np.asarray(dd.energies)[inside], lab[inside],
                                np.asarray(dd.angles)[inside])
                worst_data = max(worst_data, float(np.max(np.abs(got - ref) / np.abs(ref))))
                n_data += int(inside.sum())
            if n_data:
                print(f"      {what} at {n_data} data points: largest |file - engine| / engine {worst_data:.2e}")
                check(f"{name}: {what}: dsigma/dOmega of the file at the data's points = the engine to "
                      f"{bound:g}", worst_data < bound, worst_data)
            # (ii) between the tabulated energies and angles
            ref, e_cm = engine_dsdo(s, label, mid, angles, d, keys)
            got = np.array([file_dsdo(m, s, entrance, key, [ec] * len(angles), [el] * len(angles), angles)
                            for ec, el in zip(e_cm[:len(mid)], mid)])
            floor = 1e-3 * np.max(np.abs(ref))
            dev = np.abs(got - ref) / np.maximum(np.abs(ref), floor)
            print(f"      {what} between grid energies and angles ({dev.size} points): largest "
                  f"|file - engine| / max(engine, 1e-3 max) {dev.max():.2e}, median {np.median(dev):.1e}")
            check(f"{name}: {what}: dsigma/dOmega between the tabulated points = the engine to {bound:g}",
                  dev.max() < bound, dev.max())
        if os.environ.get("AZURE2_ENDF_UTILITIES") and os.environ.get("AZURE2_PREPRO"):
            result = E.check(out)
            print("\n".join("      " + line for line in str(result).splitlines()))
        s.close()
    finally:
        os.chdir(cwd)

# (e) a level narrower than the spacing of the start grid is found: the start grid holds the level energies.  A copy
# of T2 with the widths of its 3- level at 13.0734 MeV made 0.1, 0.1 and 0.01 eV (the committed project is not
# changed); the writer's grid of MF=3/MT=800 against the engine on a fine grid across the level.  (A file is not
# written: MF=2 cannot reproduce a level of 0.2 eV to the writer's 1e-5, its ER field holds 1e-3 eV at 1 MeV.)
print("T2 with a narrow level")
d = os.path.join(TMP, "narrow")
shutil.copytree(os.path.dirname(PROJECTS[1]), d, ignore=shutil.ignore_patterns("output*", "checks*"))
for sub in ("output", "checks"):
    os.makedirs(os.path.join(d, sub), exist_ok=True)
cwd = os.getcwd()
os.chdir(d)
try:
    project = pyazr.Project.open("T2.azr")
    level = next(lv for lv in project.levels() if lv.jpi == "3-" and abs(lv.energy - 13.0734) < 1e-6)
    for c in level.channels:
        project.set_channel(level.label, c.pair, c.l, c.s, width=c.width and (0.01 if c.pair == 3 else 0.1))
    project.save()
    s = pyazr.Session("T2.azr")
    ent = next(p.to_dict() for p in s.pairs if p.is_entrance)
    e_lab = float(pyazr.physics.cm_to_lab_energy(ent["M1"], ent["M2"], 13.0734 - ent["sep_energy"]))
    lo, hi = 0.807, 1.207  # the level lies 2.6 keV from an even point and 2.4 keV from a midpoint
    start = E.grid.start_grid(s, None, (lo, hi))
    even = np.linspace(lo, hi, E.grid.START_POINTS)
    check("narrow level: its lab energy is a point of the start grid, > 2 keV from the even points and midpoints",
          min(abs(np.array(start) - e_lab)) < 1e-12 and min(abs(np.concatenate([even, (even[1:] + even[:-1]) / 2])
                                                                - e_lab)) > 2e-3, e_lab)
    work = tempfile.mkdtemp(dir=d)
    sweeps = E.engine_values.Sweeps(s, None, work, 8)
    g = sweeps.energy_grid([E.engine_values.Reaction(800, "exit", project.reaction(1, 2), (), 2)],
                           E.grid.start_grid(s, None, (lo, hi)), TOLERANCE, 4000)
    grid, values = g.energies, {800: g.sigma[(800,)]}
    fine = e_lab + np.linspace(-20e-6, 20e-6, 401)  # +-20 eV lab around the level, 0.1 eV apart
    ref = engine_sigma(s, project.reaction(1, 2), fine, d, _pair_keys(project, 1, 2))
    order = np.argsort(grid)
    got = np.interp(fine, grid[order], values[800][order])
    dev = np.max(np.abs(got - ref)) / np.max(ref)
    print(f"      {len(grid)} grid points; peak of the engine {np.max(ref):.4g} b at {fine[np.argmax(ref)]:.7f} MeV, "
          f"of the grid {np.max(got):.4g} b; largest |grid - engine| / peak over +-20 eV: {dev:.2e}")
    check("narrow level: the grid follows the engine across the level (to 2e-3 of its peak)", dev < 2e-3, dev)
    s.close()
finally:
    os.chdir(cwd)

if os.environ.get("ENDF_TEST_KEEP"):
    print(f"files kept in {TMP}")
else:
    shutil.rmtree(TMP, ignore_errors=True)
print(f"{len(failures)} failure(s)" + (": " + ", ".join(failures) if failures else ""))
sys.exit(1 if failures else 0)
