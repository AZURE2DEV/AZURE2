"""review.py (2026-11-07): the engine-free review of a finished CLI run, on the AZURE2 13N test output.  Seconds."""
import os, sys, tempfile, shutil
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3 (its pyazr 3 version is test_p3_review.py)")
    sys.exit(77)
from rmfit import AZURE2_ROOT
from rmfit import review as rv

SRC = os.path.join(AZURE2_ROOT, "tests", "13N")
FAILS = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""))
    if not ok:
        FAILS.append(name)


def main():
    azr = os.path.join(SRC, "13N.azr"); odir = os.path.join(SRC, "output")
    # 1. the point table reproduces chiSquared.out (parsing, block assignment and the n*d / n*e convention)
    tab = rv.table_from_output(azr, odir)
    chisq = rv.read_chisq_out(os.path.join(odir, "chiSquared.out"))
    ntot = sum(d["n"] for d in chisq.values()); ctot = sum(d["chi2"] for d in chisq.values())
    check("point count = chiSquared.out Total-N", tab.n == ntot, f"{tab.n} vs {ntot}")
    check("plain chi2 = chiSquared.out total (0.1 %)", abs(tab.chi2.sum() - ctot) / ctot < 1e-3, f"{tab.chi2.sum():.1f} vs {ctot:.1f}")
    per = {k: float(tab.chi2[tab.key == k].sum()) for k in tab.keys()}
    worst = max(abs(per[k] - chisq[k]["chi2"]) / max(chisq[k]["chi2"], 1) for k in per)
    check("per-segment chi2 matches (1e-3)", worst < 1e-3, f"worst rel {worst:.2e}")
    check("angles carried for differential segments", np.isfinite(tab.ang[tab.key == 3]).all() and np.isnan(tab.ang[tab.key == 1]).all())
    # 1b. a c.m. differential segment (isDiff 4, "differential-cm") carries its angles too (missing from ANGLE_OBS until
    #     2026-11-15): the same output read through a copy of the model with segment 3 marked c.m.
    tmpc = tempfile.mkdtemp(prefix="rmfit_review_cm_")
    try:
        lines = open(azr).read().split("\n")
        i3 = [i for i, ln in enumerate(lines) if "meyer_84.dat" in ln][0]
        f = lines[i3].split()
        f[7] = "4"
        lines[i3] = " ".join(f)
        azr_cm = os.path.join(tmpc, "13N.azr")
        open(azr_cm, "w").write("\n".join(lines))
        os.symlink(os.path.join(SRC, "data"), os.path.join(tmpc, "data"))
        tabc = rv.table_from_output(azr_cm, odir)
        from rmfit.model import base_segments
        segc = [s for s in base_segments(azr_cm) if s.key == 3][0]
        check("differential-cm segment carries its angles", segc.observable == "differential-cm" and
              np.isfinite(tabc.ang[tabc.key == 3]).all() and np.allclose(tabc.ang[tabc.key == 3], tab.ang[tab.key == 3]),
              segc.observable)
    finally:
        shutil.rmtree(tmpc, ignore_errors=True)
    # 2. levels from parameters.out: the 1/2+ 2.3682 has G = 34 keV + 480 meV
    lv = rv.levels_from_parameters_out(os.path.join(odir, "parameters.out"))
    l1 = [l for l in lv if l.jpi == "1/2+" and abs(l.energy - 2.3682) < 1e-3]
    check("parameters.out levels parsed", len(lv) >= 6 and len(l1) == 1, str(len(lv)))
    check("1/2+ 2.3682 total width 34.00048 keV", l1 and abs(l1[0].gamma_mev - (34.0e-3 + 480e-9)) < 1e-9, l1 and str(l1[0].gamma_mev))
    # 3. levels from the .azr's physical gamma fields (open channels only)
    la = rv.levels_from_azr(azr)
    check(".azr levels carry widths", la and all(l.gamma_mev is not None for l in la), str(len(la)))
    # 4. the whole review runs and writes REVIEW.md + figures
    tmp = tempfile.mkdtemp(prefix="rmfit_review_")
    try:
        r = rv.review_output(azr, odir, out_dir=tmp, figures=2)
        txt = open(r["path"]).read()
        check("REVIEW.md written with findings", os.path.exists(r["path"]) and sum(r["counts"].values()) > 0 and "## Per-point summary" in txt)
        check("sections present", all(h in txt for h in ("## Normalizations", "## Residual runs", "## Average versus structure",
                                                          "## Zigzag in the calculated curve", "## Levels inside the data range", "## Figures")))
        figs = [f for f in os.listdir(os.path.join(tmp, "figures")) if f.endswith(".png")]
        check("figures written (residual map + 2 data sets)", len(figs) == 3, str(figs))
        codes = {f.code for f in r["findings"]}
        check("the seed's misfit shows as systematic-misfit + residual-run", {"systematic-misfit", "residual-run"} <= codes, str(codes))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    # 5. synthetic norm-opposite case (norm 1.32 while the fit is 0.64 x the data)
    n = 40
    ex = np.linspace(8.0, 9.0, n); d = np.full(n, 10.0); e = np.full(n, 0.3); c = 0.64 * d
    tab2 = rv.PointTable(key=np.full(n, 1), name=np.array(["SETB25"] * n, dtype=object), ecm=ex, ex=ex, ang=np.full(n, np.nan),
                         fit=c, data=d, err=e, norm=np.full(n, 1.32))

    class Seg:
        key = 1; observable = "angle-integrated"; vary_norm = True; norm_error = 9.7; norm_group = "SETB25"; covariance = None; data_file = "SETB25.dat"
        entrance_key = 1; exit_key = 2

    class Segs(list):
        def by_key(self, k): return self[0]

    rev = rv.Review(tab2, [], Segs([Seg()])).run(None)
    opp = [f for f in rev.findings if f.code == "norm-opposite"]
    check("norm-opposite fires at norm 1.32 vs n* 0.64", len(opp) == 1 and "0.640" in opp[0].message, [f.line() for f in rev.findings][:3])
    # 5b. width-tied close pair (gap 1 keV, widths ~220 keV): flagged at the default 10 % scale, not when
    #     the gap exceeds it, not for another J^pi, not outside the data range, and never with the check switched off
    exc = np.linspace(10.0, 11.0, 30)
    tabc = rv.PointTable(key=np.full(30, 1), name=np.array(["X"] * 30, dtype=object), ecm=exc, ex=exc, ang=np.full(30, np.nan),
                         fit=np.full(30, 5.0), data=np.full(30, 5.0), err=np.full(30, 0.5), norm=np.ones(30))
    lvc = [rv.Level("7/2+", 10.4198, 0.220, [(1, 4, 0.5, 0.2), (2, 3, 0.5, 0.02)]), rv.Level("7/2+", 10.4208, 0.300, [(1, 4, 0.5, 0.3)]),
           rv.Level("7/2+", 10.60, 0.300), rv.Level("5/2-", 10.4200, 0.250), rv.Level("5/2+", 10.70, 0.100), rv.Level("5/2+", 10.72, 0.100),
           rv.Level("1/2+", 15.0, 5.0), rv.Level("1/2+", 15.1, 5.0), rv.Level("3/2+", 10.5, None), rv.Level("3/2+", 10.5001, 0.2)]
    revc = rv.Review(tabc, lvc, None); revc.close_pairs()
    cl = [f for f in revc.findings if f.code == "level-close-width"]
    check("level-close-width flags only the 7/2+ pair at 10.42", len(cl) == 1 and cl[0].severity == "MEDIUM" and "7/2+" in cl[0].message
          and "10.41980" in cl[0].message and "10.42080" in cl[0].message, [f.line() for f in revc.findings])
    check("its section is written", any(h.startswith("Close same-J^pi pairs") for h, _ in revc.sections))
    rev20 = rv.Review(tabc, lvc, None, close_fraction=0.25); rev20.close_pairs()
    check("a larger scale also flags the 5/2+ pair (20 keV of 100 keV)", len(rev20.findings) == 2, [f.line() for f in rev20.findings])
    rev0 = rv.Review(tabc, lvc, None, close_fraction=0); rev0.close_pairs()
    check("close_fraction 0 switches the check off", not rev0.findings and not rev0.sections)
    # 6. zigzag: a saw-tooth fit over finely spaced points
    ex = np.linspace(11.40, 11.47, 36); c = 5.0 * (1 + 0.03 * np.where(np.arange(36) % 2, 1, -1)); d = np.full(36, 5.0); e = np.full(36, 0.5)
    tab3 = rv.PointTable(key=np.full(36, 89), name=np.array(["ND_an_t1"] * 36, dtype=object), ecm=ex, ex=ex, ang=np.full(36, np.nan),
                         fit=c, data=d, err=e, norm=np.ones(36))
    rev3 = rv.Review(tab3, [rv.Level("7/2-", 11.435, 606e-6)], None).run(None)
    zz = [f for f in rev3.findings if f.code == "fit-zigzag"]
    check("fit-zigzag fires on a 3 % saw-tooth and names the narrow level", len(zz) == 1 and "7/2- 11.4350" in zz[0].message, [f.line() for f in rev3.findings][:3])
    print("\nFAILED: " + ", ".join(FAILS) if FAILS else "\nPASS")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
