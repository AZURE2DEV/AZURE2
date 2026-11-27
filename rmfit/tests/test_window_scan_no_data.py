#!/usr/bin/env python3
"""A level candidate whose window holds no data point is untestable, not a crash (2026-11-27).

levelscan --refit on a candidate above the top of the data (a refit half a MeV above the highest data point) built a
window project with every data set deactivated; Sharded found no shard block and raised IndexError in _assemble.
scan_jpi_on_window now returns one "crashed" result saying so (run_level_add records it as scan_failed, which the drive
counts as untestable, not as a marginal add), without opening an engine session.
On the 13N fixture (AZURE2_ROOT/tests/13N, data well below Ex 40 MeV): a candidate at Ex 40 MeV.
Shared by both backends; no engine session is opened.
Run:  AZURE2_ROOT=<tree with tests/13N> python3 rmfit/tests/test_window_scan_no_data.py     (seconds)"""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit import AZURE2_ROOT  # noqa: E402
from rmfit.structure import scan_jpi_on_window, LevelCandidate  # noqa: E402

FAILS = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""), flush=True)
    if not ok:
        FAILS.append(name)


def main():
    src = os.path.join(AZURE2_ROOT, "tests", "13N")
    if not os.path.isfile(os.path.join(src, "13N.azr")):
        print(f"SKIP: no tests/13N under AZURE2_ROOT={AZURE2_ROOT}")
        sys.exit(77)
    tmp = tempfile.mkdtemp(prefix="rmfit_nodata_")
    try:
        work = os.path.join(tmp, "13N")
        shutil.copytree(src, work)
        lines = []
        cand = LevelCandidate(energy=40.0, jpi="1/2+", channels=[], source="test", note="far above the data")
        res = scan_jpi_on_window(os.path.join(work, "13N.azr"), os.path.join(tmp, "scan"), cand,
                                 {"use_brune": True}, jpis=["1/2+"], log=lines.append, nshards=2)
        check("one result, no exception", len(res) == 1, str(res))
        r = res[0] if res else None
        check("status crashed (-> scan_failed in run_level_add), no vector", r is not None and r.status == "crashed" and r.x is None, str(r))
        check("the note says the window has no data", r is not None and "no data point" in r.note, str(r))
        check("the log says so too", any("no data point" in ln for ln in lines), "\n".join(lines))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if FAILS:
        print(f"{len(FAILS)} FAILED"); sys.exit(1)
    print("all window-without-data checks passed")


if __name__ == "__main__":
    main()
