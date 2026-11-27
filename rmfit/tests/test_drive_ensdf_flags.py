#!/usr/bin/env python3
"""rmfit drive: ENSDF candidate widths and J^pi, and review flags of a level that moved (2026-11-21).

In a benchmark, pyazr.fetch_levels left width_kev None for every 16O level above S_alpha (ENSDF quotes the
width in the half-life field; half_life_sec = hbar ln 2 / Gamma), so the broad-first sort fell back to energy order, the
level_tests_max cut dropped the top of the window, and the queue opened with two unnatural-parity levels that collapsed
and ended it ("2 marginal adds in a row") before any level near the deleted ones was tried; "4(-)" / "4(+)" were read
as no J^pi.  The REVIEW step then waited on the base model's own broad-0+ flag because the pole had moved 220 keV.
  (a) CANDIDATES: a width from half_life_sec when width_kev is None; tentative J^pi "3/2(+)", "(7/2)+" read; a J list
      ("(1/2,3/2)+", NNDC's "(1 2)+") stays None; broad first, so the budget keeps the broadest levels.
  (b) review.same_flag_subject: the same J^pi level within max(50 keV, half the larger width) is the same flag.
  (c) REVIEW: a HIGH flag on a base-model level that moved within its width is inherited (no decision); a flag on
      another level still waits.
  (d) bench.new_review_flags uses the same matching.
Engine-free, seconds.  Run:  python3 rmfit/tests/test_drive_ensdf_flags.py"""
import json, os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
os.environ.pop("JOB_ID", None)
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402  (shared by both backends since 2026-11-26)

FAILS = []
HL = 6.582119569e-19 * 0.6931471805599453        # hbar ln 2, keV s


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""))
    sys.stdout.flush()
    if not ok:
        FAILS.append(name)


class Lv:
    def __init__(self, e, jp, width_kev=None, gamma_kev=None):
        self.energy_mev, self.jp, self.width_kev = e, jp, width_kev
        self.half_life_sec = HL / gamma_kev if gamma_kev else None
        self.half_life = None


def main():
    import pyazr
    if _PYAZR3:
        import test_p3_drive_dry as tdd
    else:
        import test_drive_dry as tdd
    from rmfit import drive as dv
    from rmfit import review as rv
    from rmfit.spec import Finding
    tmp = tempfile.mkdtemp(prefix="rmfit_ensdf_flags_")
    try:
        base = tdd.make_dir(tmp)                   # 13N, fit window Ex 1.9-4.5, level_tests_max 3 (set per call below)
        d = dv.Driver(tmp, base=base, local=True)
        # (a) CANDIDATES from a fake ENSDF table (no level of the 13N model within 30 keV of these)
        fake = [Lv(2.60, "3/2(+)", gamma_kev=500.0), Lv(2.90, "(7/2)+", width_kev=30.0), Lv(3.10, "1/2-", gamma_kev=1e-3),
                Lv(4.20, "(1/2,3/2)+", gamma_kev=1000.0), Lv(4.40, "5/2-"), Lv(4.45, "(1 2)+", gamma_kev=900.0)]
        orig = getattr(pyazr, "fetch_levels", None)
        cpath = os.path.join(tmp, "candidates.json")

        def build(nmax):
            if os.path.exists(cpath):
                os.remove(cpath)
            d.budget["level_tests_max"] = nmax
            pyazr.fetch_levels = lambda cn: fake
            try:
                d._s_candidates()
            finally:
                if orig is not None:
                    pyazr.fetch_levels = orig
            return json.load(open(cpath))["candidates"]
        bank = build(10)                           # the whole table
        by_e = {c["energy"]: c for c in bank}
        check("width from half_life_sec when width_kev is None (2.6 MeV: 500 keV)", 2.6 in by_e and abs((by_e[2.6]["width_kev"] or 0) - 500.0) < 0.5,
              str(bank))
        check("tentative J^pi read: '3/2(+)' -> 3/2+, '(7/2)+' -> 7/2+", by_e.get(2.6, {}).get("jpi") == "3/2+" and by_e.get(2.9, {}).get("jpi") == "7/2+",
              str([(c["energy"], c["jpi"]) for c in bank]))
        check("J lists stay without a J^pi: '(1/2,3/2)+', and NNDC's '(1 2)+' is not 12+",
              4.2 in by_e and by_e[4.2]["jpi"] is None and 4.45 in by_e and by_e[4.45]["jpi"] is None, str([by_e.get(4.2), by_e.get(4.45)]))
        check("broad first over the whole table: 1000, 900, 500, 30 keV, 1 eV, unknown width",
              [c["energy"] for c in bank] == [4.2, 4.45, 2.6, 2.9, 3.1, 4.4], str([(c["energy"], c["width_kev"]) for c in bank]))
        bank = build(3)
        check("the budget (3) keeps the three broadest (1000, 900, 500 keV)", [c["energy"] for c in bank] == [4.2, 4.45, 2.6],
              str([(c["energy"], c["width_kev"]) for c in bank]))
        # (b) the matcher
        S = rv.same_flag_subject
        check("same 0+ moved 220 keV, G ~2 MeV -> same flag", S("0+ 14.3323 G=2.09 MeV", "0+ 14.5548 G=1.71 MeV"))
        check("other J^pi at the same energy -> different", not S("2+ 14.3323 G=2.09 MeV", "0+ 14.3323 G=2.09 MeV"))
        check("3- moved 570 keV, G 763 keV -> different", not S("3- 11.4336 G=763 keV", "3- 12.0036 G=763 keV"))
        check("narrow 6+ moved 61 keV -> different (50 keV floor)", not S("6+ 11.4389 G=28 eV", "6+ 11.5000 G=28 eV"))
        check("non-level subjects: equal text only", S("norm SETB25", "norm SETB25") and not S("norm SETB25", "norm ND5U26"))
        # (c) REVIEW with a stub review
        d._campaign = tdd.StubCampaign(tmp, d.decisions)
        d.state["baseline_review"] = [["broad-level-averages-peaks", "0+ 14.5548 G=1.71 MeV"]]
        findings = [Finding("HIGH", "broad-level-averages-peaks", "0+ 14.3323 G=2.09 MeV: a level 2.09 MeV wide spans 24 visible unfitted peaks")]
        orig_rc = rv.review_campaign
        rv.review_campaign = lambda c, **kw: dict(findings=list(findings), counts={"HIGH": len(findings)})
        try:
            n0 = len(d.decisions.all())
            r1 = d._s_review()
            check("REVIEW: the base model's 0+ flag after the pole moved is inherited -> next, no decision",
                  r1[0] == "next" and len(d.decisions.all()) == n0 and "inherited" in r1[1], str(r1))
            findings.append(Finding("HIGH", "broad-level-averages-peaks", "2+ 12.9500 G=1.4 MeV: a level 1.4 MeV wide spans 7 visible unfitted peaks"))
            r2 = d._s_review()
            check("REVIEW: a flag on another level still waits", r2[0] == "stay" and len(d.decisions.all()) == n0 + 1, str(r2))
        finally:
            rv.review_campaign = orig_rc
        # (d) bench
        sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "bench")))
        from rmfit.bench import bench as bb
        run = {("HIGH", "broad-level-averages-peaks", "0+ 14.3323 G=2.09 MeV"), ("MEDIUM", "broad-level-averages-peaks", "3- 11.4336 G=763 keV")}
        ref = {("HIGH", "broad-level-averages-peaks", "0+ 14.5547 G=1.7 MeV"), ("HIGH", "broad-level-averages-peaks", "0+ 13.4734 G=2.48 MeV")}
        new = bb.new_review_flags(run, ref)
        check("bench: only the 3- flag is new", new == [("MEDIUM", "broad-level-averages-peaks", "3- 11.4336 G=763 keV")], str(new))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if FAILS:
        print(f"{len(FAILS)} FAILED"); sys.exit(1)
    print("all ENSDF-candidate / review-flag checks passed")


if __name__ == "__main__":
    main()
