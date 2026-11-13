#!/usr/bin/env python3
"""Inactive <levels> lines and inactive <segmentsData> lines must not shift
pyazr's bookkeeping off the engine's.

AZURE2 skips every level line with isActive == 0 (CNuc::Fill) and every data
segment with isActive == 0, so neither takes part in the engine's numbering.
Two places in pyazr used to count them anyway:

  * the level bookkeeping numbered levels over all file lines, so save_fit
    wrote a fit one level early in any J-group with an inactive line -- which is how the 13C+alpha archive's 5/2+ block was
    silently corrupted at two successive bakes.

  * the per-data-set priors and chi-squared (Session.evaluate) indexed datasets[i] for i below the
    active count, dropping the normalization penalty of every active segment
    that sits past that position in the file and mislabeling datasets.

Both are exercised here on the 13N model with one inactive level spliced into
the middle of a J-group (Project.add_level, set_level(active=False), its rows
moved in front of the group's second level; Project.levels() numbers the
levels as the engine does, through experiment::LevelEdits) and one inactive segment spliced in front of the
segment whose normalization is freed.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/inactive_lines_test.py
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import Session, Project
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def chi2_line(work):
    out = os.path.join(work, "output", "chiSquared.out")
    line = ""
    if os.path.exists(out):
        for l in open(out):
            if l.startswith("Total-Chi-Squared:"):
                line = l
    return line


with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(ROOT, "tests", "13N"), work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))

    print("1. an inactive level in the middle of a J-group")
    model = Project.open(os.path.join(work, "13N.azr"))
    group = [lv for lv in model.levels() if lv.jpi == "1/2+"]
    check("13N has a two-level 1/2+ group to splice into", len(group) == 2, str(group))
    model.add_level(0.5, +1, 10.0, [dict(pair=c.pair, l=c.l, s=c.s, width=1.0, fixed=c.fixed)
                                    for c in group[0].channels])
    model.set_level("1/2+ 10", active=False)
    spliced = model.save(os.path.join(work, "spliced.azr"))
    # the ghost's rows (the last of the section) moved in front of the 20 MeV pole's first row
    lines = open(spliced).read().split("\n")
    ghost_rows = sorted(c.line - 1 for c in model.level("1/2+ 10").channels)
    rows = [lines[i] for i in ghost_rows]
    for i in reversed(ghost_rows):
        del lines[i]
    at = min(c.line - 1 for c in group[1].channels)
    lines[at:at] = rows
    open(spliced, "w").write("\n".join(lines))
    model = Project.open(spliced)
    ghost = [lv for lv in model.levels() if not lv.active]
    check("the level reports itself inactive", len(ghost) == 1 and abs(ghost[0].energy - 10.0) < 1e-12, str(ghost))
    active = [lv for lv in model.levels() if lv.active]
    check("the engine's numbering skips the inactive level",
          ghost[0].group == 0 and ghost[0].level == 0 and len({(lv.group, lv.level) for lv in active}) == len(active),
          str([(lv.label, lv.group, lv.level) for lv in model.levels()]))
    check("the ghost stands between the two 1/2+ levels",
          min(c.line for c in group[0].channels) < min(c.line for c in ghost[0].channels)
          < min(c.line for c in model.level(("1/2+", 20.0)).channels))

    with Session(spliced, cwd=work) as m:
        x = np.asarray(m.params_rwa, float)
        chi2_ref = m.evaluate(x).chi2
        # every engine (jgroup, level) must land on a file level with the
        # same J^pi and energy
        keys = {(lv.group, lv.level): lv for lv in Project.open(spliced).levels() if lv.active}
        bad = []
        for key, lvls in m.parameters.by_physical_level().items():
            lv = keys.get((key.jgroup, key.level))
            if lv is None:
                bad.append(f"{key}: no file level")
                continue
            p = lvls[0]
            e = 0.0 if p.level_energy is None else float(p.level_energy)
            if lv.jpi != key.jpi or abs(lv.energy - e) > 1e-6:
                bad.append(f"{key}: file has {lv.jpi} at {lv.energy}")
        check("every engine level maps to the matching file level", not bad,
              "; ".join(bad[:4]))

        # perturb the 1/2+ background pole (the level after the ghost) and
        # snapshot: the values must land on the 20 MeV rows, not the ghost
        pole = [p for p in m.parameters.widths
                if p.J is not None and abs(p.J - 0.5) < 1e-9 and p.parity > 0
                and p.level_energy is not None
                and abs(p.level_energy - 20.0) < 1e-6 and not p.fixed]
        check("the 1/2+ pole has a free width to perturb", bool(pole))
        moved = x.copy()
        for p in pole:
            moved[p.free_index] *= 1.2
        chi2_moved = m.evaluate(moved).chi2
        snap, sav = m.save_fit(os.path.join(work, "snap.azr"), moved)
        baked = Project.open(snap)
        g = [lv for lv in baked.levels() if lv.jpi == "1/2+" and not lv.active]
        check("the ghost rows are still inactive and untouched",
              len(g) == 1 and abs(g[0].energy - 10.0) < 1e-9
              and all(c.width == 1.0 for c in g[0].channels))
        check("the pole's width changed on the active 20 MeV rows",
              any(c.width != group[1].channels[0].width
                  for lv in baked.levels() if lv.jpi == "1/2+" and lv.active
                  and abs(lv.energy - 20.0) < 1e-9 for c in lv.channels))

    probe = Project.open(snap)
    probe.set_config(output_dir=os.path.join(work, "output_snap") + "/")
    probe_path = probe.save(os.path.join(work, "snap_probe.azr"))
    os.makedirs(os.path.join(work, "output_snap"), exist_ok=True)
    with Session(probe_path, cwd=work) as m2:
        chi2_back = m2.evaluate(m2.params_rwa).chi2
        check("the snapshot reproduces the perturbed chi-squared",
              np.isclose(chi2_back, chi2_moved, rtol=1e-6),
              f"{chi2_back} vs {chi2_moved} (unperturbed {chi2_ref})")

    purged = Project.open(snap)
    gone = [lv for lv in purged.levels() if not lv.active]
    for lv in gone:
        purged.remove_level(lv)
    check("removing the inactive levels removes exactly the ghost",
          len(gone) == 1 and abs(gone[0].energy - 10.0) < 1e-9
          and not [lv for lv in purged.levels() if not lv.active])
    purged.set_config(output_dir=os.path.join(work, "output_purged") + "/")
    os.makedirs(os.path.join(work, "output_purged"), exist_ok=True)
    purged_path = purged.save(os.path.join(work, "purged.azr"))
    with Session(purged_path, cwd=work) as m3:
        chi2_purged = m3.evaluate(m3.params_rwa).chi2
        check("purging the inactive level changes nothing for the engine",
              np.isclose(chi2_purged, chi2_moved, rtol=1e-9),
              f"{chi2_purged} vs {chi2_moved}")

    print("\n2. an inactive segment in front of a segment with a free norm")
    # the format-3 project: meyer_84's norm free with its 5 % prior (Project.set), an inactive copy of the artemov
    # data set spliced in front of every data set, the outputs to output_seg/
    from pyazr import Project
    proj = Project.open(os.path.join(work, "13N.azr"))
    proj.set("meyer_84", norm="1 fit prior 1 5%")
    segazr = os.path.join(work, "seg.azr")
    proj.save(segazr)
    lines = open(segazr).read().split("\n")
    a0 = lines.index("dataset artemov")
    a1 = lines.index("", a0)
    ghost = ["dataset ghost", "    active         no"] + [l for l in lines[a0 + 1:a1]] + [""]
    lines[a0:a0] = ghost
    c = lines.index("<config>") + 1
    lines[c + 1] = "output_seg/"
    open(segazr, "w").write("\n".join(lines))
    os.makedirs(os.path.join(work, "output_seg"), exist_ok=True)

    with Session(segazr, cwd=work) as m:
        x = np.asarray(m.params_rwa, float)
        active = m.active_datasets
        check("active_datasets has one entry per engine segment",
              len(active) == m.nsegments, f"{len(active)} vs {m.nsegments}")
        check("the inactive line is not among them",
              all(d.active for d in active) and len(m.datasets) == m.nsegments + 1)
        byfile = {d.name: (d.chi2, d.n) for d in m.evaluate(x).data}
        counts = {}
        for i, d in enumerate(active):
            counts[d.name] = counts.get(d.name, 0) + len(m.energies[i])
        check("evaluate().data point counts follow the active segments",
              all(byfile[k][1] == v for k, v in counts.items())
              and set(byfile) == set(counts))
        norm = m.parameters.norms[0]
        seg = m.datasets.by_key(norm.segment_key)
        check("the free norm belongs to meyer_84 by key", seg is not None and "meyer_84" in seg.name)
        moved = x.copy()
        moved[norm.free_index] = 1.10
        pen = m.evaluate(moved)
        check("the penalty lands on meyer_84's active slot",
              np.isclose(np.array([d.norm_chi2 for d in pen.data]).sum(), 4.0, rtol=1e-9)
              and np.isclose(np.array([d.norm_chi2 for d in pen.data])[[i for i, d in enumerate(active)
                                          if d.key == seg.key][0]], 4.0),
              str(np.array([d.norm_chi2 for d in pen.data])))
        chi2_moved = m.evaluate(moved).chi2
        check("objective = chi-squared + penalty",
              np.isclose(m.objective()(moved), chi2_moved + 4.0, rtol=1e-12))
        allrwa = list(np.asarray(m.sess.params_all_rwa(), float))
        it = iter(range(len(moved)))
        sav = os.path.join(work, "moved.sav")
        with open(sav, "w") as fh:
            for i, p in enumerate(m.parameters):
                v = moved[next(it)] if not m.fixed_params[i] else allrwa[i]
                fh.write(f"{p.name:>28s} {float(v): .7e} {0.0: .7e}\n")

    binary = os.environ.get("AZURE2_BIN")  # ctest's: the binary this build made
    for cand in [] if binary else sorted(glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))):
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            binary = cand
            break
    if binary is None:
        print("  skip  no AZURE2 binary to cross-check against")
    else:
        os.makedirs(os.path.join(work, "output_seg"), exist_ok=True)
        subprocess.run([binary, "--no-gui", "--no-readline", "seg.azr"],
                       cwd=work, input="1\nmoved.sav\n\n7\n", text=True,
                       capture_output=True, timeout=900)
        line = ""
        out = os.path.join(work, "output_seg", "chiSquared.out")
        if os.path.exists(out):
            for l in open(out):
                if l.startswith("Total-Chi-Squared:"):
                    line = l
        check("AZURE2 produced a total", bool(line), "no chiSquared.out line")
        if line:
            fields = line.split()
            engine_norm = float(fields[fields.index("Total-Norm-Chi-Squared:") + 1])
            engine_chi2 = float(fields[fields.index("Total-Chi-Squared:") + 1])
            print(f"        engine: {line.strip()}")
            check("pyazr's penalty is the engine's Total-Norm-Chi-Squared",
                  np.isclose(engine_norm, 4.0, rtol=1e-3, atol=1e-6),
                  f"{engine_norm} vs 4.0")
            check("pyazr's chi-squared is the engine's Total-Chi-Squared",
                  np.isclose(engine_chi2, chi2_moved, rtol=1e-4),
                  f"{engine_chi2} vs {chi2_moved}")

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all inactive-line checks passed")
