#!/usr/bin/env python3
"""Inactive <levels> lines and inactive <segmentsData> lines must not shift
pyazr's bookkeeping off the engine's.

AZURE2 skips every level line with isActive == 0 (CNuc::Fill) and every data
segment with isActive == 0, so neither takes part in the engine's numbering.
Two places in pyazr used to count them anyway:

  * AzrModel.engine_level_keys() numbered levels over all file lines, so
    apply_fit / save_fit wrote a fit one level early in any J-group with an
    inactive line -- which is how the 13C+alpha archive's 5/2+ block was
    silently corrupted at two successive bakes.

  * azure2.penalties() / dataset_chi2() indexed datasets[i] for i below the
    active count, dropping the normalization penalty of every active segment
    that sits past that position in the file and mislabeling datasets.

Both are exercised here on the 13N model with one inactive level spliced into
the middle of a J-group and one inactive segment spliced in front of the
segment whose normalization is freed.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/inactive_lines_test.py
"""
import copy
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
    from pyazr import azure2, AzrModel
    from pyazr.azrfile import AzrLevel
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)


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
    model = AzrModel.from_file(os.path.join(work, "13N.azr"))
    group = model.find(jpi="1/2+")
    check("13N has a two-level 1/2+ group to splice into", len(group) == 2,
          str(group))
    ghost = AzrLevel([copy.deepcopy(c) for c in group[0].channels])
    ghost.set_energy(10.0)
    ghost.set_active(False)
    for c in ghost.channels:
        c.gamma = 1.0
    at = model.levels.index(group[1])          # in front of the 20 MeV pole
    model.levels.insert(at, ghost)
    check("the level reports itself inactive", ghost.active is False)
    keys = model.engine_level_keys()
    check("engine_level_keys skips the inactive level",
          ghost not in keys.values() and len(keys) == len(model.active_levels),
          f"{len(keys)} keys for {len(model.active_levels)} active levels")
    spliced = os.path.join(work, "spliced.azr")
    model.write(spliced)

    with azure2(spliced, cwd=work) as m:
        x = np.asarray(m.params_rwa, float)
        chi2_ref = float(np.sum(m.calculate_chi2_rwa(x)))
        # every engine (jgroup, level) must land on a file level with the
        # same J^pi and energy
        fresh = AzrModel.from_file(spliced)
        keys = fresh.engine_level_keys()
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
        # snapshot: the values must land on the 20 MeV line, not the ghost
        pole = [p for p in m.parameters.widths
                if p.J is not None and abs(p.J - 0.5) < 1e-9 and p.parity > 0
                and p.level_energy is not None
                and abs(p.level_energy - 20.0) < 1e-6 and not p.fixed]
        check("the 1/2+ pole has a free width to perturb", bool(pole))
        moved = x.copy()
        for p in pole:
            moved[p.free_index] *= 1.2
        chi2_moved = float(np.sum(m.calculate_chi2_rwa(moved)))
        snap, sav = m.save_fit(os.path.join(work, "snap.azr"), moved)
        baked = AzrModel.from_file(snap)
        g = [lv for lv in baked.levels if lv.jpi == "1/2+" and not lv.active]
        check("the ghost line is still inactive and untouched",
              len(g) == 1 and abs(g[0].energy - 10.0) < 1e-9
              and all(c.gamma == 1.0 for c in g[0].channels))
        check("the pole's width changed on the active 20 MeV line",
              any(c.gamma != group[1].channels[0].gamma
                  for lv in baked.levels if lv.jpi == "1/2+" and lv.active
                  and abs(lv.energy - 20.0) < 1e-9 for c in lv.channels))

    probe = AzrModel.from_file(snap)
    probe.set_output_dir(os.path.join(work, "output_snap"))
    probe_path = probe.write(os.path.join(work, "snap_probe.azr"))
    with azure2(probe_path, cwd=work) as m2:
        chi2_back = float(np.sum(m2.calculate_chi2_rwa(m2.params_rwa)))
        check("the snapshot reproduces the perturbed chi-squared",
              np.isclose(chi2_back, chi2_moved, rtol=1e-6),
              f"{chi2_back} vs {chi2_moved} (unperturbed {chi2_ref})")

    purged = AzrModel.from_file(snap)
    gone = purged.purge_inactive_levels()
    check("purge_inactive_levels removes exactly the ghost",
          len(gone) == 1 and abs(gone[0].energy - 10.0) < 1e-9
          and not [lv for lv in purged.levels if not lv.active])
    purged.set_output_dir(os.path.join(work, "output_purged"))
    purged_path = purged.write(os.path.join(work, "purged.azr"))
    with azure2(purged_path, cwd=work) as m3:
        chi2_purged = float(np.sum(m3.calculate_chi2_rwa(m3.params_rwa)))
        check("purging the inactive level changes nothing for the engine",
              np.isclose(chi2_purged, chi2_moved, rtol=1e-9),
              f"{chi2_purged} vs {chi2_moved}")

    print("\n2. an inactive segment in front of a segment with a free norm")
    model = AzrModel.from_file(os.path.join(work, "13N.azr"))
    model.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    # splice an inactive copy of the artemov line at the top of the block
    lines = model._suffix.splitlines()
    top = lines.index("<segmentsData>") + 1
    ghost_seg = lines[top].split()
    ghost_seg[0] = "0"
    lines.insert(top, " ".join(ghost_seg))
    model._suffix = "\n".join(lines)
    model.set_output_dir(os.path.join(work, "output_seg"))
    segazr = model.write(os.path.join(work, "seg.azr"))

    with azure2(segazr, cwd=work) as m:
        x = np.asarray(m.params_rwa, float)
        active = m.active_datasets
        check("active_datasets has one entry per engine segment",
              len(active) == m.nsegments, f"{len(active)} vs {m.nsegments}")
        check("the inactive line is not among them",
              all(d.active for d in active) and len(m.datasets) == m.nsegments + 1)
        byfile = m.dataset_chi2(x)
        counts = {}
        for i, d in enumerate(active):
            counts[d.name] = counts.get(d.name, 0) + len(m.energies[i])
        check("dataset_chi2 point counts follow the active segments",
              all(byfile[k][1] == v for k, v in counts.items())
              and set(byfile) == set(counts))
        norm = m.parameters.norms[0]
        seg = m.datasets.by_key(norm.segment_key)
        check("the free norm belongs to meyer_84 by key", "meyer_84" in seg.name)
        moved = x.copy()
        moved[norm.free_index] = 1.10
        pen = m.penalties(moved)
        check("the penalty lands on meyer_84's active slot",
              np.isclose(pen["norm"].sum(), 4.0, rtol=1e-9)
              and np.isclose(pen["norm"][[i for i, d in enumerate(active)
                                          if d.key == seg.key][0]], 4.0),
              str(pen["norm"]))
        chi2_moved = float(np.sum(m.calculate_chi2_rwa(moved)))
        check("objective = chi-squared + penalty",
              np.isclose(m.objective(moved), chi2_moved + 4.0, rtol=1e-12))
        allrwa = list(np.asarray(m.sess.params_all_rwa(), float))
        it = iter(range(len(moved)))
        sav = os.path.join(work, "moved.sav")
        with open(sav, "w") as fh:
            for i, p in enumerate(m.parameters):
                v = moved[next(it)] if not m.fixed_params[i] else allrwa[i]
                fh.write(f"{p.name:>28s} {float(v): .7e} {0.0: .7e}\n")

    binary = None
    for cand in sorted(glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))):
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
