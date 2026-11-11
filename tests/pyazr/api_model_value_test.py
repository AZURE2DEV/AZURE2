#!/usr/bin/env python3
"""The API reports every segment's model value, also for an E1-only or E2-only segment.

The value a data point is compared with is its segment's model value (Results::model, ModelPass.h): for an
angle-integrated-E1 / -E2 segment (isDiff 5 / 6) the E1 or E2 part of the capture cross section.  The binding's
calculated_segments (update_segments, update_segments_rwa, update_segments_all_rwa) returned the base point's leaf
cross section -- the total -- for such a segment, while chi^2 used the part; pyazr's calculate* swapped the part in
on its own side, calculate_sfactor* did not.  Since the plan branch (director's decision, HANDOFF) the engine returns
the model value for every segment.

Model: tests/channel_space/channel_space.azr plus E1-only and E2-only copies of its capture data set (`multipole E1`,
`E2`).  Checks, for each
of the three entry points: calculated_segments of the E1 (E2) segment equals its calculated_segments_e1 (_e2), bit
for bit; the plain capture segment's total differs from its E1 part (control); pyazr's calculate_rwa and
calculate_sfactor_rwa give the part (times the S-factor conversion).  The start-head engine fails the first checks.

Needs the compiled engine; skips cleanly without it.  Run from anywhere:  python3 tests/pyazr/api_model_value_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
failures = 0


def check(name, ok, detail=""):
    global failures
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}{'' if ok or not detail else '  -- ' + detail}")
    if not ok:
        failures += 1


def main():
    sys.path.insert(0, ROOT)
    try:
        import numpy as np
        from pyazr import azure2
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"skip: engine not importable ({e})")
        return 77

    work = tempfile.mkdtemp(prefix="pyazr_api_model_value_")
    try:
        src = os.path.join(ROOT, "tests", "channel_space")
        shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
        for d in ("output", "checks"):
            os.makedirs(os.path.join(work, d), exist_ok=True)
        azr = os.path.join(work, "parts.azr")
        text = open(os.path.join(src, "channel_space.azr")).read()
        extra = "".join(
            f"dataset capture_{m}\n    reaction       12C(p,g)13N\n    quantity       sigma\n    unit           b\n"
            f"    multipole      {m}\n    file           data/capture.dat\n    columns        E_lab[MeV] _ value error\n"
            f"    energy_window  0.5 7 MeV\n    norm           1 fixed\n\n" for m in ("E1", "E2"))
        assert text.count("</data>\n") == 1
        text = text.replace("</data>\n", extra + "</data>\n", 1)
        with open(azr, "w") as f:
            f.write(text)

        E1, E2 = "angle-integrated cross section, E1 part", "angle-integrated cross section, E2 part"
        with azure2(azr, cwd=work) as m:
            s = m.sess
            obs = [d.observable for d in m.active_datasets]
            check("the model has a total, an E1-only and an E2-only capture segment",
                  E1 in obs and E2 in obs, f"{obs}")
            k1, k2 = obs.index(E1), obs.index(E2)
            kt = k1 - 1  # the plain capture segment, written just before the two copies
            x = np.asarray(m.params_rwa, float)
            entry = {
                "update_segments_rwa": lambda: s.update_segments_rwa(x),
                "update_segments": lambda: s.update_segments(np.asarray(m.params, float)),
                "update_segments_all_rwa": lambda: s.update_segments_all_rwa(np.asarray(s.params_all_rwa(), float)),
            }
            for name, call in entry.items():
                n = int(call())
                check(f"{name}: {n} segments", n == len(obs))
                for k, part, getter in ((k1, "E1", s.calculated_segments_e1), (k2, "E2", s.calculated_segments_e2)):
                    v, p = np.asarray(s.calculated_segments(k)), np.asarray(getter(k))
                    same = v.shape == p.shape and v.tobytes() == p.tobytes()
                    check(f"{name}: the {part}-only segment reports its {part} part", same,
                          f"max |value - part| = {np.max(np.abs(v - p)):.3g} b" if v.shape == p.shape else "")
                tot, e1 = np.asarray(s.calculated_segments(kt)), np.asarray(s.calculated_segments_e1(kt))
                check(f"{name}: control -- the plain capture segment's total is not its E1 part",
                      not np.array_equal(tot, e1), f"max |total/E1 - 1| = {np.max(np.abs(tot / e1 - 1)):.3g}")
            v = m.calculate_rwa(x)
            check("calculate_rwa: the E2-only segment is the E2 part",
                  np.array_equal(np.asarray(v[k2]), np.asarray(s.calculated_segments_e2(k2))))
            sf = m.calculate_sfactor_rwa(x)
            want = np.asarray(s.calculated_segments_e2(k2)) * np.asarray(s.calculated_conv(k2))
            check("calculate_sfactor_rwa: the E2-only segment is the E2 part times the conversion",
                  np.array_equal(np.asarray(sf[k2]), want))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print(f"api_model_value: {'FAIL' if failures else 'all checks passed'} ({failures} failure{'' if failures == 1 else 's'})")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
