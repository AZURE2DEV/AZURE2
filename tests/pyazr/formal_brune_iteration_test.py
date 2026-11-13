#!/usr/bin/env python3
"""The formal -> Brune iteration inside an evaluation does not jump (docs/upgrade/level_parameters.md section 9).

Model: tests/beam_profile_kernel in the formal basis (--no-brune), where external capture goes to final states that
share their J^pi blocks and the --ignore-externals masks are projectors, so every evaluation runs brune_from_formal.
Raising width_9_2 (a 1- particle width, 1730) by 1e-6 relative once made the 1- level with Brune energy 20.0 MeV stop
at steps of 7.7e-12 MeV against a 1e-12 MeV tolerance; its formal values were then kept, the capture masks changed
and chi2 jumped 40.44 -> 42.93 (external_capture.md section 8).

Check: chi2 at the start and at +-1e-6 relative on width_9_2 is smooth -- the second difference is below 1e-2 of the
first (1e-4 on 5055bec0, which had no projectors; 1.0 on the builds with the jump) -- and finite.

Run from anywhere:  python3 tests/pyazr/formal_brune_iteration_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
PROJECT = os.path.join(ROOT, "tests", "beam_profile_kernel")


def main():
    try:
        import numpy as np
        from pyazr import Session
    except Exception as err:  # same skip as the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 77
    failures = []

    def check(name, ok, detail=""):
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
        if not ok:
            failures.append(name)

    work = tempfile.mkdtemp(prefix="pyazr_brune_iteration_")
    try:
        shutil.copytree(os.path.join(PROJECT, "data"), os.path.join(work, "data"))
        shutil.copy(os.path.join(PROJECT, "beam_profile_kernel.azr"), work)
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        with Session(os.path.join(work, "beam_profile_kernel.azr"), cwd=work, use_brune=False) as a:
            x = np.array(a.params_rwa, float)
            target = [p for p in a.parameters if str(p).find("width_9_2,") >= 0 and not p.fixed]
            check("width_9_2 is a free parameter", len(target) == 1)
            if target:
                i = target[0].free_index
                h = 1e-6 * (abs(x[i]) + 1.0)          # the step of pyazr's finite-difference fallback
                c = []
                for d in (-h, 0.0, h):
                    v = x.copy()
                    v[i] += d
                    c.append(a.calculate_chi2_rwa(v)[0])
                first, second = c[2] - c[0], c[2] - 2.0 * c[1] + c[0]
                check(f"chi2 finite at the start and +-h ({c[0]:.10g}, {c[1]:.10g}, {c[2]:.10g})",
                      all(np.isfinite(c)))
                ratio = abs(second) / max(abs(first), 1e-300)
                check(f"chi2 smooth in width_9_2: |second difference| / |first| = {ratio:.1e}", ratio < 1e-2)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
