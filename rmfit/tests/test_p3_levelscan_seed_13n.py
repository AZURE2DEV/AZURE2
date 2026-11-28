#!/usr/bin/env python3
"""levelscan dummy seeds must be physical: no failed Brune transform (2026-11-20).

template_variant used to write every particle channel of every dummy level at a flat 1 eV.  A
channel just above its threshold cannot hold 1 eV (Gamma_W far below it), so AZURE2's transform of
the dummy from physical widths to Brune reduced widths had no real solution: "Denominator less than
zero while transforming ... The transformation may not have been successful for this level" (hundreds
of times in one scan).  The seed is now min(gamma_ref_ev, 1e-2 * Gamma_W), and a J^pi with
no channel of Gamma_W >= 0.1 eV is not written (a physical seed there is so narrow that AZURE2's target-
integration grid hangs at init; an exactly-zero seed drops the parameter key -- both found on this fixture).

On 13N (12C+p, S_p = 1.944 MeV):
  A. Ex 2.1 MeV, the OLD behaviour (seed_theta2 huge -> a flat 1 eV): the transform warning (minimal input);
  B. Ex 2.1 MeV, default seeds: the J^pi that cannot carry 0.1 eV are skipped, the rest load (well within the time
     limit) with no transform warning, every seed has theta^2 <= 1e-2 against the ENGINE's Gamma_W, and the engine's
     rwa of each seed equals the Wigner-based amplitude of that width (ratio ~1: the transform succeeded);
  D. Ex 1.97 MeV (26 keV above threshold): every J^pi skipped, the variant loads;
  C. Ex 9.0 MeV (1 eV was harmless there): the same parameter keys as the 1 eV seeds and the same objective with the
     dummy amplitudes zeroed.
Run: python3 rmfit/tests/test_p3_levelscan_seed_13n.py (pyazr 3; the pyazr 2 version is test_levelscan_seed_13n.py)
"""
import os, shutil, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if not _PYAZR3:
    print("SKIP: a pyazr 3 test; this process runs pyazr 2 (the pyazr 2 version is test_levelscan_seed_13n.py)")
    sys.exit(77)
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.engine import Sharded
from rmfit.levelscan import template_variant, dummy_labels
from rmfit.structure import resolve_jpis
from rmfit.azr3 import Model as AzrModel
FLAGS = {"use_brune": True}
WARN = "Denominator less than zero while transforming"

failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok: failures.append(name)


def load_capturing(path, shard_dir):
    """Sharded evaluator of ``path`` with the engine's stdout AND stderr captured into a file (pyazr 3 prints the
    transform warnings on stderr, in the workers it starts; pyazr 2 printed them on stdout -- capturing stdout alone made
    checks B and D pass without looking, 2026-11-26); returns (ev, text)."""
    cap = os.path.join(shard_dir + "_stdout.txt")
    os.makedirs(os.path.dirname(cap), exist_ok=True)
    sys.stdout.flush(); sys.stderr.flush()
    saved1, saved2 = os.dup(1), os.dup(2)
    fd = os.open(cap, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    os.dup2(fd, 1); os.dup2(fd, 2)
    try:
        ev = Sharded(path, shard_dir, nshards=2, flags=FLAGS)
        ev.objective(ev.x0)
    finally:
        sys.stdout.flush(); sys.stderr.flush()
        os.dup2(saved1, 1); os.dup2(saved2, 2); os.close(saved1); os.close(saved2); os.close(fd)
    return ev, open(cap).read()


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_seed_")
    work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True); os.makedirs(os.path.join(work, "output"))
    base = os.path.join(work, "13N.azr")
    jpis = resolve_jpis(AzrModel.from_file(base), None)
    print(f"J^pi groups: {jpis}")
    EX = 2.1

    print(f"A. flat 1 eV dummies at Ex {EX} MeV (the behaviour before 2026-11-20)")
    var_old = template_variant(base, os.path.join(work, "old.azr"), EX, jpis=jpis, seed_theta2=1e300, min_open_ev=0.0,
                               output_dir=os.path.join(work, "out_old"))
    seeds_old = [c["gamma"] for lv in var_old["levels"].values() for c, ph in zip(lv["channels"], lv["photon"]) if not ph]
    check("A: the old behaviour writes 1 eV in every particle channel of every J^pi",
          len(var_old["levels"]) == len(jpis) and all(g == 1.0 for g in seeds_old), (len(var_old["levels"]), seeds_old))
    ev, out = load_capturing(var_old["path"], os.path.join(tmp, "sh_old"))
    ev.close()
    print(f"   transform warnings: {out.count(WARN)}")
    check("A: the flat 1 eV seed reproduces the transform failure", out.count(WARN) > 0, out[-400:])

    print(f"B. default seeds at Ex {EX} MeV")
    var = template_variant(base, os.path.join(work, "new.azr"), EX, jpis=jpis, output_dir=os.path.join(work, "out_new"))
    print(f"   written {sorted(var['levels'])}; skipped {var['skipped']}")
    check("B: some J^pi skipped and some written", var["skipped"] and var["levels"], var["skipped"])
    check("B: written + skipped = every J^pi", sorted(list(var["levels"]) + list(var["skipped"])) == sorted(jpis))
    t0 = time.time()
    ev, out = load_capturing(var["path"], os.path.join(tmp, "sh_new"))
    print(f"   load + one evaluation {time.time() - t0:.0f} s; transform warnings {out.count(WARN)}")
    check("B: no transform warning with the new seeds", out.count(WARN) == 0, out[-400:])
    labs = dummy_labels(ev, var)
    check("B: every written dummy found in the engine's key table", set(labs) == set(var["levels"]), f"{labs} vs {list(var['levels'])}")
    gW_amp = ev.wigner_amplitudes(); gW_wid = ev.wigner_widths(ev.x0)
    ratios, t2s, nz = [], [], 0
    for jpi, info in var["levels"].items():
        lab = labs.get(jpi)
        for j, k in enumerate(ev.keys):
            if k.kind != "width" or k.level != lab:
                continue
            spec = [c for c, ph in zip(info["channels"], info["photon"])
                    if int(c["pair"]) == k.pair and int(c["L"]) == k.L and abs(float(c["S"]) - k.S) < 1e-6 and not ph]
            if not spec:
                continue
            g = float(spec[0]["gamma"])
            nz += int(ev.x0[j] != 0.0)
            if j in gW_wid and gW_wid[j] > 0 and gW_amp.get(j, 0) > 0:
                t2s.append(g / gW_wid[j])
                ratios.append(abs(float(ev.x0[j])) / np.sqrt(g) / (gW_amp[j] / np.sqrt(gW_wid[j])))
    ev.close()
    print(f"   {len(ratios)} particle channels; seed theta^2 (engine Gamma_W) {min(t2s):.3g}..{max(t2s):.3g}; "
          f"x0 / Wigner-based amplitude {min(ratios):.4f}..{max(ratios):.4f}; nonzero dummy amplitudes {nz}")
    check("B: every particle seed has theta^2 <= 1.01 * seed_theta2 (engine Gamma_W)", ratios and max(t2s) <= 1.01 * var["seed_theta2"], max(t2s))
    check("B: the engine's rwa of each seed equals the Wigner-based amplitude (transform succeeded)",
          ratios and 0.98 < min(ratios) and max(ratios) < 1.02, f"{min(ratios)}..{max(ratios)}")
    check("B: dummy amplitudes are nonzero (keys survive)", nz == len(ratios) and nz > 0, (nz, len(ratios)))

    print("D. Ex 1.97 MeV (26 keV above the proton threshold)")
    vd = template_variant(base, os.path.join(work, "thr.azr"), 1.97, jpis=jpis, output_dir=os.path.join(work, "out_thr"))
    check("D: every J^pi skipped", not vd["levels"] and sorted(vd["skipped"]) == sorted(jpis), vd["skipped"])
    ev, out = load_capturing(vd["path"], os.path.join(tmp, "sh_thr"))
    ev.close()
    check("D: the variant loads, no transform warning", out.count(WARN) == 0, out[-300:])

    print("C. Ex 9.0 MeV: default seeds vs the old flat 1 eV")
    vo = template_variant(base, os.path.join(work, "old9.azr"), 9.0, jpis=jpis, seed_theta2=1e300, min_open_ev=0.0, output_dir=os.path.join(work, "o9"))
    vn = template_variant(base, os.path.join(work, "new9.azr"), 9.0, jpis=jpis, output_dir=os.path.join(work, "n9"))
    eo, out_o = load_capturing(vo["path"], os.path.join(tmp, "sh_o9"))
    en, out_n = load_capturing(vn["path"], os.path.join(tmp, "sh_n9"))
    check("C: same parameter keys as with the 1 eV seeds", [str(k) for k in eo.keys] == [str(k) for k in en.keys])
    check("C: no transform warning with the new seeds at 9 MeV", out_n.count(WARN) == 0, out_n[-300:])
    lo, ln = dummy_labels(eo, vo), dummy_labels(en, vn)
    xo, xn = eo.x0.copy(), en.x0.copy()
    for j, k in enumerate(eo.keys):
        if k.kind == "width" and k.level in set(lo.values()): xo[j] = 0.0
    for j, k in enumerate(en.keys):
        if k.kind == "width" and k.level in set(ln.values()): xn[j] = 0.0
    oo, on = eo.objective(xo), en.objective(xn)
    eo.close(); en.close()
    print(f"   zero-amplitude dummies: objective {oo:.6f} (1 eV variant) vs {on:.6f} (new seeds)")
    check("C: the incumbent objective on the variant does not depend on the seeds", abs(oo - on) <= 1e-6 * max(1.0, abs(oo)), f"{oo} vs {on}")

    shutil.rmtree(tmp, ignore_errors=True)
    print("FAILED: " + ", ".join(failures) if failures else "ALL OK")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
