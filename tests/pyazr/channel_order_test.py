#!/usr/bin/env python3
"""Channel order inside a level is presentation, not physics.

A J-group's channels can be written in any order -- the group owns the set, and
every level of the group carries a width for each of them. Listing them
differently must not change a result, and must certainly not crash.

It did both. CompoundSystem::TransformIn and CompoundSystem::TransformOut each padded their
per-channel arrays for a non-particle channel by copying entry [0], to keep the
array indexed by channel number. When channel 1 was *itself* a photon or beta
channel that array was still empty, so the copy read off the end of an empty
vector: a segfault on a file whose <levels> simply happens to list its capture
channel first. Every reader of those entries is guarded by RadType == 'P', so
the padding value is never used and zero does as well.

TransformIn is on the default path. TransformOut's copy is reached only with
the Brune formalism off, which is why neither had been seen.

Reversing every level keeps the levels of a group in step. One level listing
its channels in a different order from the rest of its group, or leaving out a
channel the group carries, is the harder case: until 2026-10-17 CompoundSystem::Fill
stored each width in the position of its line, so the first put widths on the
wrong channels without a message and the second left the level's width vector
one short (TransformIn then wrote past its end). Widths are now matched to
their channel by (pair, l, s) and an unlisted channel gets a zero, fixed width;
both must give the as-written chi-squared, and save_fit must write a snapshot
of the short model that reads back.

pyazr.Project offers no reordering of a level's rows (the order is not physics); the variants are made here by
moving and deleting whole rows of the file, at the lines Project.levels() gives for each channel.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/channel_order_test.py
"""
import importlib.util
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


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import Session
    from pyazr import Project
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)

# 13N is the fixture with a photon channel to move: reversing each level's
# channels puts it first, which is the condition that used to crash.
PROJECT = os.path.join(ROOT, "tests", "13N")

with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "13N")
    shutil.copytree(PROJECT, work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))

    plain = os.path.join(work, "13N.azr")
    levels = Project.open(plain).levels()
    pairs = {q.key: q for q in Project.open(plain).pairs()}
    text = open(plain).read().split("\n")

    def variant_of(order):
        """The file with each level's rows replaced by `order(level)` (a list of its channels, or None to keep)."""
        lines = list(text)
        for lv in levels:
            new = order(lv)
            if new is None:
                continue
            at = sorted(c.line - 1 for c in lv.channels)
            rows = [text[c.line - 1] for c in new]
            for i in reversed(at):
                del lines[i]
            lines[at[0]:at[0]] = rows
            # the following levels' lines move with the length of this one
            if len(rows) != len(at):
                raise AssertionError("one level shortened at most, the last edit")
        return "\n".join(lines)

    def write(name, body):
        path = os.path.join(work, name)
        open(path, "w").write(body)
        return path

    reversed_path = write("reversed.azr", variant_of(lambda lv: list(reversed(lv.channels))))
    photon_first = sum(1 for lv in levels if pairs[lv.channels[-1].pair].kind == 10)
    check("the reversed file leads with a photon channel", photon_first > 0,
          "nothing to test: no level starts with one")

    # One level out of step with its group: the second level of a J-group
    # (whose channel order the group took from its first level) reversed.
    seen, odd = set(), []
    for lv in levels:
        if (lv.J, lv.parity) in seen and len(lv.channels) > 1 and not odd:
            odd.append(lv.label)
        seen.add((lv.J, lv.parity))
    one_path = write("one_level.azr", variant_of(lambda lv: list(reversed(lv.channels)) if lv.label in odd else None))
    check("one level reversed", len(odd) == 1)

    # A short level: a later level of its group without the row of a zero-width
    # channel that is neither its first nor its last (rows deleted from the end of the file backwards, so the
    # line numbers of the levels before it hold).
    seen, dropped = set(), None
    for lv in levels:
        if (lv.J, lv.parity) in seen and dropped is None and len(lv.channels) > 2:
            mid = [c for c in lv.channels[1:-1] if c.width == 0.0]
            if mid:
                dropped = (lv.jpi, lv.energy, mid[0].pair, mid[0].l, mid[0].line)
        seen.add((lv.J, lv.parity))
    lines = list(text)
    if dropped is not None:
        del lines[dropped[4] - 1]
    short_path = write("short.azr", "\n".join(lines))
    check("one level short", dropped is not None)

    # Every variant gets its own output directory. AZURE2 caches its
    # external-capture integrals in output/intEC.dat and reuses them when the
    # file is there, so runs sharing one directory are not independent -- on
    # this model that alone moves chi-squared by 0.6, which would swamp what
    # this test is trying to see.
    results = {}
    for name, source in (("as written", plain), ("reversed", reversed_path),
                         ("one level", one_path), ("short", short_path)):
        for brune in (True, False):
            tag = f"{name.replace(' ', '_')}_{brune}"
            outdir = os.path.join(work, "out_" + tag)
            os.makedirs(outdir, exist_ok=True)
            variant = Project.open(source)
            variant.set_config(output_dir=outdir + "/")
            path = variant.save(os.path.join(work, tag + ".azr"))
            try:
                with Session(path, cwd=work, use_brune=brune) as m:
                    x = np.asarray(m.parameters.x0, float)
                    results[(name, brune)] = m.evaluate(x).chi2
            except Exception as err:
                results[(name, brune)] = f"{type(err).__name__}: {err}"

    for key, value in sorted(results.items(), key=lambda kv: str(kv[0])):
        print(f"        {key[0]:11} brune={str(key[1]):5} -> {value}")

    for brune in (True, False):
        a = results[("as written", brune)]
        for name in ("reversed", "one level", "short"):
            b = results[(name, brune)]
            check(f"{name} loads (brune={brune})", isinstance(b, float), str(b))
            if isinstance(a, float) and isinstance(b, float):
                check(f"and gives the same chi-squared (brune={brune})",
                      abs(a - b) <= 1e-9 * max(abs(a), 1.0), f"{a} vs {b}")

    # save_fit writes every parameter, the short level's unlisted zero width
    # included: it has no line to go on, and must not make the snapshot refuse.
    outdir = os.path.join(work, "out_short_save")
    os.makedirs(outdir, exist_ok=True)
    variant = Project.open(short_path)
    variant.set_config(output_dir=outdir + "/")
    variant.save(os.path.join(work, "short_save.azr"))
    try:
        with Session(os.path.join(work, "short_save.azr"), cwd=work) as m:
            snap, _ = m.save_fit(os.path.join(work, "short_snapshot.azr"))
        # read back from an empty output directory, not the save run's (intEC cache, above)
        fresh = os.path.join(work, "out_short_back")
        os.makedirs(fresh, exist_ok=True)
        back_model = Project.open(snap)
        back_model.set_config(output_dir=fresh + "/")
        back_model.save(snap)
        with Session(snap, cwd=work) as back:
            c = back.evaluate(np.asarray(back.parameters.x0, float)).chi2
        a = results[("as written", True)]
        check("save_fit of the short model reads back",
              isinstance(a, float) and abs(c - a) <= 1e-6 * max(abs(a), 1.0), f"{c} vs {a}")
    except Exception as err:
        check("save_fit of the short model", False, f"{type(err).__name__}: {err}")

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all channel-order checks passed")
