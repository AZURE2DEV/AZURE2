#!/usr/bin/env python3
"""The parameter table and the pairs, built from the engine's records (pyazr.Parameters, pyazr.Pair).

A session gives each fit parameter and each particle pair as a record with named fields; the fields are one table
of the engine (api/src/pybind_records.cpp: _azure2.parameter_record_fields(), pair_record_fields()).  Checked:

  1. the table: every kind of fit parameter of the engine (_azure2.parameter_kinds()) is described by a field shown
     as its owner, so that a new kind cannot be added without saying what its parameters belong to; every field
     names known kinds, has a doc and a unit entry; pyazr.kinds() carries the fields and each kind's owner fields;
  2. on tests/upgrade/T2, tests/13N, the gamma_decay projects of tests/gamma_decay/models (partial.azr, record.azr)
     and tests/18O_p_a_thm with a linear cbackground=, against the program run on the same project (menu 1):
     the names and start values of Parameters, row by row, are those of the program's param.par; Parameters.table()
     has one line per parameter with its name and start; the fixed flags are the project's (through
     Project.levels(): a level energy fixed as written or below every particle threshold of its nucleus, a channel
     row fixed as written or written as 0, a channel without a row fixed at zero; a data set's nuisance through the
     engine's description of the project; a THM background term through the program's thm_experiments.out);
  3. on the same projects: x0 is params_rwa (a copy), free / of_kind / of_level / by_name select rows, put and take
     write and read x by row, write and read a parameter file give x back, a fixed row is refused by put; on 13N
     pyazr.snapshot.write_values puts x into a Project in memory (the file untouched) and the project saved
     from it starts at x;
  4. the pair records against the project's pairs (Project.pairs(): key, masses, separation energy, excitations,
     channel radius), and the session's pairs in order of their numbers.

Needs the compiled engine and the program (AZURE2_BIN); without them exits 77 (ctest: Skipped).

Run from anywhere:  AZURE2_BIN=<build>/src/AZURE2 python3 tests/pyazr/parameter_records_test.py
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
os.environ.setdefault("OMP_NUM_THREADS", "2")

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({str(detail)[:600]})"))
    if not ok:
        failures.append(name)


try:
    import numpy as np
    import pyazr
    from pyazr import _azure2, Session, Project
    import thm_records
    from pyazr.described import describe_project
    from pyazr.snapshot import write_values
except Exception as err:  # pragma: no cover
    print(f"skip: the engine module is not available ({type(err).__name__}: {err})")
    sys.exit(77)

BINARY = os.environ.get("AZURE2_BIN")
if not BINARY or not os.access(BINARY, os.X_OK):
    print("skip: AZURE2_BIN names no program to compare against")
    sys.exit(77)


# -- 1. the table -------------------------------------------------------------------------------------------------
print("1. the table of record fields")
kinds = [k["key"] for k in _azure2.parameter_kinds()]
fields = _azure2.parameter_record_fields()
undescribed = [k for k in kinds if not any(f["owner"] and k in f["kinds"] for f in fields)]
check("every kind of parameter has an owner field in the record table", not undescribed, undescribed)
check("every field names known kinds", all(set(f["kinds"]) <= set(kinds) for f in fields),
      [f["name"] for f in fields if not set(f["kinds"]) <= set(kinds)])
check("every field has a doc", all(f["doc"] for f in fields), [f["name"] for f in fields if not f["doc"]])
names = [f["name"] for f in fields]
check("field names are unique", len(names) == len(set(names)))
check("the common fields describe every kind",
      all(f["every_kind"] and set(f["kinds"]) == set(kinds) for f in fields
          if f["name"] in ("index", "name", "kind", "unit", "start", "fixed", "free_index", "bounds")))
k = pyazr.kinds()
check("pyazr.kinds(): parameter_fields and pair_fields are the engine's",
      list(k["parameter_fields"]) == names and list(k["pair_fields"]) == [f["name"] for f in _azure2.pair_record_fields()])
check("pyazr.kinds(): each kind's owner fields are the fields that describe it",
      all(k["parameters"][key]["owner"] == [f["name"] for f in fields if not f["every_kind"] and key in f["kinds"]]
          for key in kinds))
check("every pair field has a doc", all(f["doc"] for f in _azure2.pair_record_fields()))


# -- 2. to 4. on projects ------------------------------------------------------------------------------------------
def copy_project(src_dir, dst):
    shutil.copytree(src_dir, dst, ignore=shutil.ignore_patterns("output", "checks", "run.log", "expected"))
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(dst, sub), exist_ok=True)
    return dst


def run_program(path, cwd):
    """The program's menu 1 on PATH: its output directory's param.par as [(name, value)], and the directory."""
    proc = subprocess.run([BINARY, "--no-gui", "--no-readline", os.path.basename(path)], cwd=cwd,
                          input="1\n\n\n7\n", text=True, capture_output=True, timeout=3600)
    out = os.path.join(cwd, "output")
    par = os.path.join(out, "param.par")
    if not os.path.isfile(par):
        print(proc.stdout[-2000:])
        raise RuntimeError(f"the program wrote no param.par for {path}")
    rows = []
    for line in open(par):
        t = line.split()
        if len(t) >= 2 and not t[0].startswith("#"):
            rows.append((t[0], float(t[1])))
    return rows, out


def expected_fixed(project, path, P, outdir):
    """{name: fixed} from the project (levels, the engine's description of its data sets) and the program's
    thm_experiments.out; names no source states are left out."""
    want = {}
    levels = {(lv.group, lv.level): lv for lv in project.levels() if lv.active}
    pairs = project.pairs()
    for p in P.of_kind("energy", "width"):
        lv = levels[(p.jgroup, p.level)]
        if p.kind == "energy":
            # fixed as written, and a level below every particle threshold of its nucleus (format 3's rule)
            open_pairs = [q for q in pairs if q.kind == 0 and (q.light[0] + q.heavy[0], q.light[1] + q.heavy[1])
                          == tuple(lv.nucleus) and lv.energy > q.threshold]
            want[p.name] = lv.fixed or not open_pairs
        else:
            # fixed as written, and a width written as 0 (no row: the engine's zero input)
            rows = [c for c in lv.channels if c.channel == p.channel]
            want[p.name] = rows[0].fixed or rows[0].width == 0.0 if rows else True
    for d in describe_project(path, os.path.dirname(path)).get("datasets", []):
        for q in d.get("parameters", []):
            if q.get("status") in ("fit", "fixed"):
                want.setdefault(q["name"], q["status"] == "fixed")
    thm = os.path.join(outdir, "thm_experiments.out")
    if os.path.isfile(thm) and P.of_kind("cbkg"):
        flags = [line.rstrip().endswith("fixed") for line in open(thm)
                 if line.strip() and not line.startswith("#") and "cbkg" in line]
        for p, f in zip(P.of_kind("cbkg"), flags):
            want[p.name] = f
    return want


def check_project(label, path, cwd, write_back=False):
    print(label)
    par, outdir = run_program(path, cwd)
    project = Project.open(path)
    with Session(path, cwd=cwd) as s:
        P = s.parameters
        n = len(P)
        check(f"{label}: one row per parameter of param.par ({n})", n == len(par), (n, len(par)))
        check(f"{label}: the names are param.par's, in its order", P.names == [r[0] for r in par],
              [(a, b[0]) for a, b in zip(P.names, par) if a != b[0]][:5])
        start = np.array([p.start for p in P])
        ref = np.array([r[1] for r in par])
        bad = [(P[i].name, start[i], ref[i]) for i in range(min(n, len(par)))
               if abs(start[i] - ref[i]) > 1e-6 * max(abs(ref[i]), 1e-12)]
        check(f"{label}: the start values are param.par's (rel 1e-6, its 7 digits)", not bad, bad[:5])
        want = expected_fixed(project, path, P, outdir)
        wrong = [(name, f) for name, f in want.items() if P.by_name(name).fixed != f]
        check(f"{label}: the fixed flags are the project's ({len(want)} of {n} stated)", not wrong and want, wrong[:5])
        text = P.table().splitlines()
        check(f"{label}: table(): a header and one line per parameter, with its name and start",
              len(text) == n + 1 and all(line.split()[:2] == [str(p.index), p.name] and f"{p.start:.6g}" in line
                                          for line, p in zip(text[1:], P)), text[:3])
        kinds_seen = {p.kind for p in P}
        check(f"{label}: table(): every row says what it belongs to (an owner field of its kind)",
              all(p.owner_text() for p in P), [p.name for p in P if not p.owner_text()][:5])
        x0 = P.x0
        check(f"{label}: x0 is params_rwa, and a copy",
              np.array_equal(x0, np.asarray(s.params_rwa, float)) and P.x0 is not P.x0 and not np.shares_memory(x0, P.x0))
        free = P.free
        check(f"{label}: free is the rows of x in order", [p.free_index for p in free] == list(range(len(x0))))
        check(f"{label}: names, kinds, units, free_mask, bounds: one per row",
              len(P.names) == len(P.kinds) == len(P.units) == len(P.free_mask) == len(P.bounds) == n
              and int(P.free_mask.sum()) == len(x0))
        for kind in sorted(kinds_seen):
            sel = P.of_kind(kind)
            assert all(p.kind == kind for p in sel)
        check(f"{label}: of_kind selects every kind ({sorted(kinds_seen)}) and is again a Parameters",
              sum(len(P.of_kind(kd)) for kd in kinds_seen) == n and isinstance(P.of_kind(*kinds_seen), pyazr.Parameters))
        lab = next((p.level_label for p in P if p.level_label), None)
        if lab:
            lv = P.of_level(lab)
            check(f"{label}: of_level('{lab}') is the level's energy and channels",
                  len({(p.jgroup, p.level) for p in lv}) == 1 and lv[0].kind == "energy"
                  and lab in [l.label for l in project.levels()])
        target = free[len(free) // 2] if len(free) else None
        if target is not None:
            x1 = P.put(x0, target, 2.5)
            x2 = P.put(x0, target.name, P.take(x0, target) + 1.0)
            check(f"{label}: put and take by row and by name",
                  x1[target.free_index] == 2.5 and np.count_nonzero(x1 != x0) <= 1
                  and x2[target.free_index] == x0[target.free_index] + 1.0 and P.take(x1, [target])[0] == 2.5)
            path_sav = os.path.join(cwd, "records.sav")
            P.write(path_sav, x1)
            back = P.read(path_sav)
            check(f"{label}: write then read gives x back", np.allclose(back, x1, rtol=1e-7, atol=0))
        if target is not None and write_back:
            before = open(path).read()
            edited = write_values(s, Project.open(path), x1)
            check(f"{label}: write_values edits the project in memory, the file untouched",
                  open(path).read() == before and edited.text != before)
            edited.save(os.path.join(cwd, "written.azr"))
            with Session(os.path.join(cwd, "written.azr"), cwd=cwd) as back:
                check(f"{label}: the project written with x1 starts at x1 (rel 1e-6)",
                      np.allclose(back.parameters.x0, x1, rtol=1e-6, atol=1e-12),
                      (back.parameters.x0[target.free_index], x1[target.free_index]))
        fixed = next((p for p in P if p.fixed), None)
        if fixed is not None:
            try:
                P.put(x0, fixed, 1.0)
                check(f"{label}: put refuses a fixed row", False)
            except pyazr.ParameterError as err:
                check(f"{label}: put refuses a fixed row", "fixed" in str(err) and fixed.name in str(err), err)
        pairs = s.pairs
        ppairs = {q.key: q for q in project.pairs()}
        ok = [q.number for q in pairs] == list(range(1, len(pairs) + 1)) and len(pairs) == len(s.sess.pair_records())
        for q in pairs:
            r = ppairs[q.key]
            ok &= (abs(q.M1 - r.light_mass) < 1e-9 and abs(q.M2 - r.heavy_mass) < 1e-9
                   and abs(q.sep_energy - r.separation_energy) < 1e-9 and abs(q.channel_radius - r.radius) < 1e-9
                   and abs(q.light_excitation - r.light_excitation) < 1e-9
                   and abs(q.heavy_excitation - r.heavy_excitation) < 1e-9
                   and q.ptype == r.kind and q.pair_kind == {0: "particle", 10: "photon", 20: "beta"}[r.kind])
        check(f"{label}: the pair records are the project's pairs ({len(pairs)})", ok, [q.to_dict() for q in pairs][:2])
        return kinds_seen


seen = set()
with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
    d = copy_project(os.path.join(ROOT, "tests", "upgrade", "T2"), os.path.join(tmp, "T2"))
    seen |= check_project("tests/upgrade/T2", os.path.join(d, "T2.azr"), d)
    d = copy_project(os.path.join(ROOT, "tests", "13N"), os.path.join(tmp, "13N"))
    seen |= check_project("tests/13N", os.path.join(d, "13N.azr"), d, write_back=True)
    for model in ("partial", "record"):
        d = copy_project(os.path.join(ROOT, "tests", "gamma_decay", "models"), os.path.join(tmp, "gd_" + model))
        seen |= check_project(f"tests/gamma_decay/models/{model}.azr", os.path.join(d, model + ".azr"), d)
    d = copy_project(os.path.join(ROOT, "tests", "18O_p_a_thm"), os.path.join(tmp, "thm"))
    text = open(os.path.join(d, "18O_p_a_thm.azr")).read()
    text = thm_records.thm_set(text, "experiment[E] cbackground=1/2+:2:0.5,0,0.5,1:linear=0.25,-0.5f,0,1e-3",
                               thm_records.datasets(text))
    open(os.path.join(d, "cb.azr"), "w").write(text)
    seen |= check_project("tests/18O_p_a_thm with a linear cbackground=", os.path.join(d, "cb.azr"), d)
check("the projects hold levels, radiative numbers, nuisances and a THM background",
      {"energy", "width", "radiative", "norm", "cbkg"} <= seen, sorted(seen))

print(f"\n{len(failures)} failure(s)" if failures else "\nall passed")
sys.exit(1 if failures else 0)
