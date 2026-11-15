#!/usr/bin/env python3
"""Data files outside data/ (model.data_dirs / link_data_dirs), 2026-11-07 (the evaluator: "Yes, fix it in rmfit").

Before: every working directory rmfit writes a model copy into linked data/ only, so segments reading e.g.
akdata/ or thermal/ were lost there (4 of 75 segments in one fit; pyazr refused the segment count).
Fixture: 13N with two of its data files moved to other/ and sub/deeper/ (segment paths rewritten), plus a shift
function whose expression holds a "/" (must not count as a path).
Part A (no engine): data_dirs, link_data_dirs (idempotent, broken links left alone), shards, single_variant.
Part B (engine): the moved model gives the same objective as the original, single session and Sharded (3).
Run: PYTHONPATH=<dir holding rmfit> python3 -I tests/test_data_dirs.py [--no-engine]
"""
import os, re, shutil, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def make_fixture(tmp):
    from rmfit import AZURE2_ROOT
    fix = next(d for d in (AZURE2_ROOT, os.path.expanduser("~/AZURE2")) if os.path.isdir(os.path.join(d, "tests", "13N")))
    src = os.path.join(fix, "tests", "13N")
    orig = os.path.join(tmp, "orig"); moved = os.path.join(tmp, "moved")
    for d in (orig, moved):
        shutil.copytree(os.path.join(src, "data"), os.path.join(d, "data"))
        os.makedirs(os.path.join(d, "output"))
    text = open(os.path.join(src, "13N.azr"), encoding="latin-1").read()
    open(os.path.join(orig, "13N.azr"), "w", encoding="latin-1").write(text)
    os.makedirs(os.path.join(moved, "other")); os.makedirs(os.path.join(moved, "sub", "deeper"))
    shutil.move(os.path.join(moved, "data", "meyer_84.dat"), os.path.join(moved, "other", "meyer_84.dat"))
    shutil.move(os.path.join(moved, "data", "kettner_55.dat"), os.path.join(moved, "sub", "deeper", "kettner_55.dat"))
    t2 = text.replace("data/meyer_84.dat", "other/meyer_84.dat").replace("data/kettner_55.dat", "sub/deeper/kettner_55.dat")
    assert t2 != text
    open(os.path.join(moved, "13N.azr"), "w", encoding="latin-1").write(t2)
    return orig, moved


def part_a(orig, moved, tmp):
    from rmfit.model import data_dirs, link_data_dirs, build_shards, single_variant
    print("Part A: data_dirs / link_data_dirs")
    check("original model: data only", data_dirs(os.path.join(orig, "13N.azr")) == ["data"], data_dirs(os.path.join(orig, "13N.azr")))
    dd = data_dirs(os.path.join(moved, "13N.azr"))
    check("moved model: data, other, sub", dd == ["data", "other", "sub"], dd)
    # a "/" inside an expression or a non-existent dir is not a data directory
    p = os.path.join(moved, "expr.azr")
    t = open(os.path.join(moved, "13N.azr"), encoding="latin-1").read()
    m = re.search(r"<segmentsData>\n(.*?)\n", t)
    line = m.group(1)
    t = t.replace(line, line + '  shiftfn "a0/x+b/c"  1 0 0.02 1  nothere/x.dat', 1)
    open(p, "w", encoding="latin-1").write(t)
    check("expressions and missing dirs ignored", data_dirs(p) == ["data", "other", "sub"], data_dirs(p))
    w = os.path.join(tmp, "w1")
    got = link_data_dirs(os.path.join(moved, "13N.azr"), w)
    check("links made", got == ["data", "other", "sub"] and all(os.path.islink(os.path.join(w, s)) for s in got), got)
    check("deep file reachable", os.path.isfile(os.path.join(w, "sub", "deeper", "kettner_55.dat")))
    check("idempotent", link_data_dirs(os.path.join(moved, "13N.azr"), w) == got)
    os.remove(os.path.join(w, "other")); os.symlink("/nonexistent/x", os.path.join(w, "other"))
    try:
        link_data_dirs(os.path.join(moved, "13N.azr"), w)
        check("a broken existing link is left alone (no FileExistsError)", os.readlink(os.path.join(w, "other")) == "/nonexistent/x")
    except FileExistsError as e:
        check("a broken existing link is left alone (no FileExistsError)", False, str(e))
    sv = single_variant(os.path.join(moved, "13N.azr"), "sv", os.path.join(tmp, "sv"))
    check("single_variant dir has other/ and sub/",
          all(os.path.isdir(os.path.join(tmp, "sv", s)) for s in ("data", "other", "sub")))
    shards = build_shards(os.path.join(moved, "13N.azr"), 3, os.path.join(tmp, "shards_a"))
    check("shard dir has other/ and sub/",
          all(os.path.isdir(os.path.join(tmp, "shards_a", s)) for s in ("data", "other", "sub")), os.listdir(os.path.join(tmp, "shards_a")))
    return sv


def part_b(orig, moved, tmp):
    from rmfit.engine import SingleSession, Sharded
    from rmfit.model import single_variant
    print("Part B: objectives (13N)")
    vals = {}
    for tag, d in (("orig", orig), ("moved", moved)):
        sv = single_variant(os.path.join(d, "13N.azr"), f"b_{tag}", os.path.join(tmp, f"b_{tag}"))
        with SingleSession(sv, cwd=os.path.join(tmp, f"b_{tag}"), name=f"b_{tag}") as ev:
            x0 = np.array(ev.x0, float)
            vals[tag] = (ev.objective(x0), ev.npts, x0, list(ev.keys))
    check("same point count", vals["orig"][1] == vals["moved"][1], f"{vals['orig'][1]} vs {vals['moved'][1]}")
    check("single session: same objective", np.isclose(vals["orig"][0], vals["moved"][0], rtol=1e-10),
          f"{vals['orig'][0]} vs {vals['moved'][0]}")
    sh = Sharded(os.path.join(moved, "13N.azr"), os.path.join(tmp, "shards_b"), nshards=3)
    try:
        from rmfit.keys import take
        xs = take(vals["orig"][2], vals["orig"][3], sh.keys, sh.x0)
        o = sh.objective(np.array(xs, float))
        check("Sharded (3) on the moved model = single original", np.isclose(o, vals["orig"][0], rtol=1e-8),
              f"{o} vs {vals['orig'][0]}")
        check("Sharded point count", sh.npts == vals["orig"][1], f"{sh.npts}")
    finally:
        sh.close()


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_datadirs_")
    try:
        orig, moved = make_fixture(tmp)
        part_a(orig, moved, tmp)
        if "--no-engine" not in sys.argv:
            part_b(orig, moved, tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("ALL OK" if not failures else f"FAILURES: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
