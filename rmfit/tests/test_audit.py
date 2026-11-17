"""audit.py (2026-11-07): engine-free, seconds.  A copy of the 13N test model with injected problems must raise the
expected codes; the clean model must raise none of them."""
import os, re, shutil, sys, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit import AZURE2_ROOT
from rmfit.audit import audit_model, write_audit
from rmfit.azr3 import Model as AzrModel

SRC = os.path.join(AZURE2_ROOT, "tests", "13N")


def codes(fs, sev=None):
    return sorted({f.code for f in fs if sev is None or f.severity == sev})


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_audit_")
    try:
        shutil.copytree(os.path.join(SRC, "data"), os.path.join(tmp, "data"))
        os.makedirs(os.path.join(tmp, "output"))
        m = AzrModel.from_file(os.path.join(SRC, "13N.azr"))
        m.set_output_dir(os.path.join(tmp, "output"))
        clean = m.write(os.path.join(tmp, "clean.azr"))
        f0 = audit_model(clean, run_target_energies=False)
        c0 = codes(f0)
        print("clean model codes:", c0)
        # the 13N test model is not "clean" by every check (placeholder zero widths, integer masses, test data with
        # round-percent errors): assert only on the codes it genuinely lacks, and on the HIGH ones
        for bad in ("norm-prior-tight", "level-merge", "data-file", "output-path", "extrapolation-folded", "legacy-sav", "norm-nominal", "model-setup"):
            assert bad not in c0, (bad, [f.message for f in f0 if f.code == bad])
        assert "model" in c0
        print("clean model HIGH:", [(f.code, f.message[:60]) for f in f0 if f.severity == "HIGH"])
        # --- inject problems ---------------------------------------------------------------------------
        # (a) nominal norm 0.9, free, with a prior 0.03 % wide (a fraction written where per cent was meant), on the
        # first active data set
        import pyazr
        bad = os.path.join(tmp, "bad.azr")
        mb = AzrModel.from_file(clean)
        first = mb.datasets()[0]
        mb.project.set(first, norm=pyazr.Nuisance(0.9, "fit", prior=(0.9, "0.03%")))
        data_file = str(mb.project[first].get("file"))
        # (b) second level of the same J^pi within 1 keV, and a free channel width exactly 0
        lv = mb.active_levels[0]
        chans = [dict(pair=c.pair, L=c.L, S=c.S, gamma=0.0, fixed=False) for c in lv.channels[:1]]
        mb.add_level(lv.J, lv.parity, float(lv.energy), chans, level_fixed=False)           # read as one level: level-merge
        lv2 = mb.active_levels[1]
        mb.add_level(lv2.J, lv2.parity, float(lv2.energy) + 0.0005, [dict(pair=c.pair, L=c.L, S=c.S, gamma=0.0, fixed=False) for c in lv2.channels[:1]],
                     level_fixed=False)                                                     # other ID, 0.5 keV: level-close
        mb.set_output_dir(os.path.join(tempfile.mkdtemp(prefix="rmfit_elsewhere_"), "output"))  # elsewhere: another campaign's output
        mb.write(bad)
        # (c) a data file whose errors are an integer-percent FLAG column, used by the first segment
        dfile = os.path.join(tmp, "data", os.path.basename(data_file))
        d = np.loadtxt(dfile)
        d = d[:40] if len(d) >= 40 else d
        d[:, 3] = d[:, 2] * np.array([1, 2, 3, 4, 5] * (len(d) // 5 + 1))[:len(d)] / 100.0
        d[3, 3] = d[3, 2] * 1e-5                     # one implausibly tiny error
        d = np.vstack([d, d[:2]])                    # duplicate rows
        np.savetxt(dfile, d)
        # (d) legacy parameter names in a .sav next to the model
        open(os.path.join(tmp, "old.sav"), "w").write("j=1_la=0_energy 2.37\nj=1_la=0_ch=0_rwa 0.1\n")
        f1 = audit_model(bad, run_target_energies=False)
        c1 = codes(f1)
        print("bad model codes:", c1)
        for want in ("norm-nominal", "norm-prior-tight", "level-merge", "zero-width", "data-error-pattern", "data-error-tiny",
                     "data-duplicates", "legacy-sav", "output-path"):
            assert want in c1, (want, c1)
        assert "norm-prior-tight" in codes(f1, "MEDIUM") and "level-merge" in codes(f1, "HIGH") and "level-close" in codes(f1, "MEDIUM") and "legacy-sav" in codes(f1, "HIGH")
        out = os.path.join(tmp, "AUDIT.md")
        counts = write_audit(f1, out, bad)
        txt = open(out).read()
        assert txt.startswith("# Pre-flight audit") and "norm-prior-tight" in txt and counts["HIGH"] >= 3, counts
        print("test_audit: PASS")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
