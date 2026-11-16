#!/usr/bin/env python3
"""Fetch a data set from EXFOR and put it into a project as a new data set, with the paper behind it.

``pyazr.sources.nds`` end to end:

1. search EXFOR for a reaction (``search_exfor``; network only),
2. fetch a data set and write it as an AZURE2 data file (``fetch_exfor``, ``ExforData.to_azr``: lab energy in MeV,
   angle in deg, cross section in b, error in b),
3. add it to a copy of the project as a ``dataset`` record (``Project.add_dataset``) and calculate it,
4. the level scheme of the compound nucleus from ENSDF (``fetch_levels``; network only),
5. the paper behind the data set (``reference``; ``resolve_doi``, network only).

The IAEA services are on the public internet.  ``--local DIR`` reads a local copy of the EXFOR master files instead
(``fetch_exfor(..., source="local")``, ``$AZURE2_EXFOR_MASTER``) and skips the steps that need the network; the
repository holds a synthetic entry for that, ``tests/pyazr/data/exfor_master`` (T0001, made up: not physics).  The
data file and the project copy go to ``--out`` (default: a new temporary directory, printed); the copy's session reads
the project's own data files from the project's directory (``Session(..., cwd=)``).

Usage (from the repository's root):

    python3 pyazr/examples/exfor_fetch.py tests/upgrade/T1/T1.azr --entry T0001003.1 \\
        --local tests/pyazr/data/exfor_master [--out DIR]
    python3 pyazr/examples/exfor_fetch.py tests/upgrade/T1/T1.azr --entry <EXFOR data set> --target C-12 --reaction p,g
"""
import argparse
import os
import sys
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "2")

from pyazr import Project, Session  # noqa: E402
from pyazr.sources import nds as N  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr) the data set is added to")
    ap.add_argument("--entry", required=True, help="the EXFOR data set (e.g. O2599004; a pointer as O2599004.1)")
    ap.add_argument("--entrance", type=int, default=1, help="the entrance pair's key")
    ap.add_argument("--exit", dest="exit_pair", type=int, default=2, help="the exit pair's key")
    ap.add_argument("--target", default="C-12", help="EXFOR target for the search (network only)")
    ap.add_argument("--reaction", default="p,g", help="EXFOR reaction for the search (network only)")
    ap.add_argument("--local", default=None, help="a local copy of the EXFOR master files: no network")
    ap.add_argument("--out", default=None, help="the directory of the data file and the copy (default: temporary)")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(args.project))
    out_dir = os.path.abspath(args.out or tempfile.mkdtemp(prefix="exfor_fetch."))
    os.makedirs(out_dir, exist_ok=True)
    source = "auto"
    if args.local:
        os.environ["AZURE2_EXFOR_MASTER"] = os.path.abspath(args.local)
        source = "local"

    if not args.local:                                                     # 1. search
        print(f"== search_exfor({args.target}, {args.reaction}, SIG) ==")
        for hit in N.search_exfor(target=args.target, reaction=args.reaction, quantity="SIG", limit=5):
            print(" ", hit)

    data = N.fetch_exfor(args.entry, source=source)                        # 2. fetch and convert
    print(f"\n== fetch_exfor({args.entry!r}) ==\n  reaction: {data.reaction}\n"
          f"  columns: {[(c.name, c.unit) for c in data.columns]}")
    kw = data.to_azr(out_dir, entrance=args.entrance, exit=args.exit_pair)
    print(f"  written: {kw['data_file']} (lab E in MeV, angle in deg, sigma and error in b)")
    with open(kw["data_file"], encoding="utf-8") as f:
        for line in f.readlines()[:3]:
            print("   ", line.rstrip())

    p = Project.open(args.project)                                         # 3. a data set of a copy
    name = "exfor_" + args.entry.replace(".", "_").lower()
    p.add_dataset(name, reaction=p.reaction(args.entrance, args.exit_pair), quantity="sigma", unit="b",
                  norm="1 fixed", file=kw["data_file"], columns="E_lab[MeV] _ value error")
    path = p.save(os.path.join(out_dir, os.path.basename(args.project)))
    with Session(path, cwd=here) as s:
        d = next(d for d in s.evaluate().data if d.name == name)
    print(f"\n== the copy {path} ==\n  {name}: {d.n} points, chi2 = {d.chi2:.4g} at the project's parameters")

    if not args.local:                                                     # 4. levels
        print("\n== fetch_levels of the compound nucleus ==")
        nucleus = data.reaction.split(")")[-1].split(",")[0].replace("-", "").lower()
        for lv in N.fetch_levels(nucleus)[:6]:
            print(f"  {lv.energy:8.3f} keV  {lv.jp:<5}  {lv.half_life}")

    ref = N.reference(args.entry.split(".")[0], source=source)             # 5. the paper
    print(f"\n== reference ==\n  entry: {ref.entry}\n  code: {ref.x4ref}\n  year: {ref.year}")
    if not args.local:
        doi = N.resolve_doi(ref)
        print(f"  DOI: {doi} -> https://doi.org/{doi}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
