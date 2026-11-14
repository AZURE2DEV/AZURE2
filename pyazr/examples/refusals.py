#!/usr/bin/env python3
"""What a refusal looks like, and how a script (or an agent) reads it instead of parsing the message.

Every deliberate error of pyazr is a ``pyazr.Error`` that also derives from the builtin the call raised before
(``ValueError``, ``KeyError``, ...).  It carries three parts beside its message: ``what`` was wrong, what is
``allowed``, which call shows or ``fix``es it (``None`` where the error does not give one).  A ``FormatError`` of the
project's grammar carries the engine's diagnostics as records: ``codes`` and ``diagnostics`` (each a dict with
``severity``, ``code``, ``record``, ``field``, ``found``, ``allowed``, ``message``, the ``file`` and ``line``), and the
project is left unchanged.  ``Project.validate()`` gives the same records for the project as it stands (an empty list:
nothing to refuse).

The script makes three refusals on a project -- a data set the grammar refuses, a pair key that does not exist, a
fit vector of the wrong length -- and prints what a caller reads of each.  Nothing is written.

Usage (from the repository's root):

    python3 pyazr/examples/refusals.py tests/upgrade/T3/T3.azr
"""
import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np  # noqa: E402

import pyazr  # noqa: E402


def show(title, call):
    """Run ``call``; print the class, the builtin it also is, and the three parts of the refusal."""
    try:
        call()
    except pyazr.Error as e:
        builtin = next(c.__name__ for c in type(e).__mro__ if c.__module__ == "builtins")
        print(f"\n{title}\n  {type(e).__name__} (a pyazr.Error and a {builtin})")
        print(f"  message  {str(e).splitlines()[0][:150]}")
        for part in ("what", "allowed", "fix"):
            print(f"  {part:<8} {getattr(e, part)}")
        if isinstance(e, pyazr.FormatError):
            print(f"  codes    {e.codes}")
            for d in e.diagnostics:
                print(f"  record {d['record']!r}, field {d['field']!r}: {d['code']} (found {d['found']!r})")
        return e
    raise SystemExit(f"{title}: not refused")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr)")
    args = ap.parse_args()

    p = pyazr.Project.open(args.project)
    pair = p.pairs()[0]
    data = next(r for r in p.records().values() if r.record_kind == "dataset")
    show(f"1. a data set of the angle-integrated elastic cross section of {pair.label} (Project.add_dataset)",
         lambda: p.add_dataset("refused", reaction=p.reaction(pair.key, pair.key), quantity="sigma", unit="b/sr",
                               file=data["file"], columns="E_lab[MeV] _ value error", norm="1 fixed"))
    print(f"  the project is unchanged: {'refused' not in p.records()}; Project.validate() -> {p.validate()}")

    with pyazr.Session(args.project) as s:
        keys = [q.key for q in s.pairs]
        show("2. channel functions of a pair key the project does not have (Session.channel_functions)",
             lambda: s.channel_functions(max(keys) + 1, [1.0]))
        show("3. an evaluation at a fit vector of the wrong length (Session.evaluate)",
             lambda: s.evaluate(np.zeros(len(s.parameters.x0) + 1)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
