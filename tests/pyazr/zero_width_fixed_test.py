#!/usr/bin/env python3
"""A width given as 0 is not fitted, in every basis (docs/upgrade/external_capture.md section 9).

The engine fixes a width parameter whose input is zero.  It used to test the transformed amplitude instead, so in the
formal basis a zero input came out as 0 or as rounding (1e-15 ... 1e-40) by the order of the sums, and rounding decided
which parameters a fit varied: tests/beam_profile_kernel freed 35 such parameters on 5055bec0 and 48 on a18fe89e.

Checks:
  (a) beam_profile_kernel: the free parameters of the formal basis (--no-brune) are those of the Brune basis, by name;
  (b) tests/13N with the photon width of the 20 MeV 1/2+ level set to 0 and --no-ignore-externals (where the
      amplitude is gamma^int = -Re gamma^ext, not 0): that parameter is fixed, in the Brune and formal bases;
  (c) apply_fit in strict mode, fixed parameters included, places the formal-basis vector save_fit writes for
      beam_profile_kernel.  This guards the round trip only: no test project has a fixed width on a channel its level
      does not list (scan of the 75 projects that load in the formal basis, 2026-11-08), so the azrfile.py branch
      that skips such a width whatever its value is not exercised here, and (c) passes on a18fe89e too; that branch
      is tested in apply_fit_test.py ("fixed 5, no line").

Run from anywhere:  python3 tests/pyazr/zero_width_fixed_test.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)


def prepare(work, project):
    src = os.path.join(ROOT, "tests", project)
    shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    if project == "13N":
        from pyazr import AzrModel
        m = AzrModel.from_file(os.path.join(src, "13N.azr"))
        level = [lv for lv in m.levels if lv.J == 0.5 and lv.parity == +1 and lv.energy == 20.0][0]
        [c for c in level.channels if c.pair == 2][0].gamma = 0.0
        path = os.path.join(work, "m.azr")
        m.write(path)
        return path
    shutil.copy(os.path.join(src, project + ".azr"), work)
    return os.path.join(work, project + ".azr")


def worker(project, brune, ignore, out):
    import numpy as np
    from pyazr import AzrModel, azure2
    work = tempfile.mkdtemp(prefix="pyazr_zero_fixed_")
    try:
        azr = prepare(work, project)
        res = {}
        with azure2(azr, cwd=work, use_brune=brune == "1", ignore_externals=ignore == "1") as a:
            res["free"] = sorted(p.name for p in a.parameters
                                 if p.kind in ("energy", "width") and not p.fixed)
            res["fixed"] = sorted(p.name for p in a.parameters if p.kind == "width" and p.fixed)
            if project == "beam_profile_kernel" and brune == "0":
                # the vector save_fit writes: every parameter, physical (apply_fit refuses include_fixed with transform=)
                want = np.asarray(a.transform_all_rwa(a._all_rwa(np.array(a.params_rwa, float)), include_fixed=True), float)
                try:
                    AzrModel.from_file(azr).apply_fit(a.parameters, want, physical=True, pairs=a.pairs, strict=True,
                                                      include_fixed=True)
                    res["apply_fit"] = "ok"
                except Exception as err:
                    res["apply_fit"] = f"{type(err).__name__}: {str(err)[:300]}"
        json.dump(res, open(out, "w"))
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    try:
        import numpy  # noqa: F401
        from pyazr import azure2  # noqa: F401
    except Exception as err:  # same skip as the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 0
    failures = []

    def check(name, ok, detail=""):
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
        if not ok:
            failures.append(name)

    def run(project, brune, ignore=True):
        d = tempfile.mkdtemp(prefix="pyazr_zero_fixed_out_")
        out = os.path.join(d, "res.json")
        try:
            p = subprocess.run([sys.executable, os.path.abspath(__file__), "worker", project, str(int(brune)),
                                str(int(ignore)), out], capture_output=True, text=True, timeout=900)
            if p.returncode != 0 or not os.path.exists(out):
                print(p.stderr[-1500:])
                return None
            return json.load(open(out))
        finally:
            shutil.rmtree(d, ignore_errors=True)

    bpk_brune, bpk_formal = run("beam_profile_kernel", True), run("beam_profile_kernel", False)
    n_brune, n_formal = run("13N", True, ignore=False), run("13N", False, ignore=False)
    if None in (bpk_brune, bpk_formal, n_brune, n_formal):
        check("every configuration ran", False, "see the traceback above")
    else:
        extra = sorted(set(bpk_formal["free"]) - set(bpk_brune["free"]))
        missing = sorted(set(bpk_brune["free"]) - set(bpk_formal["free"]))
        check(f"(a) beam_profile_kernel: formal free set = Brune free set ({len(bpk_formal['free'])} and "
              f"{len(bpk_brune['free'])})", not extra and not missing,
              f"formal only: {extra[:6]}; Brune only: {missing[:6]}")
        # width_4_2: the photon channel of the 20 MeV 1/2+ level (level 4 in the file's order, channel 2).
        for label, r in (("Brune", n_brune), ("formal", n_formal)):
            check(f"(b) 13N, photon width 0, --no-ignore-externals, {label}: width_4_2 fixed",
                  "width_4_2" in r["fixed"] and "width_4_2" not in r["free"])
        check(f"(c) apply_fit strict, formal basis, fixed included: {bpk_formal.get('apply_fit')}",
              bpk_formal.get("apply_fit") == "ok")
    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(*sys.argv[2:6])
        sys.exit(0)
    sys.exit(main())
