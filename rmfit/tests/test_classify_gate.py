#!/usr/bin/env python3
"""The 'unphysical' gate in structure.classify (2026-09-16).

In one high-energy window two levels were accepted with large gains, with
channels sitting exactly ON the theta^2 cap and the rest collapsed to ~1e-35; their
Brune transformation failed and the model could not be baked.  The blowup test
fires only ABOVE the cap.  A level whose released polish leaves any of its width
amplitudes on a bound, or whose theta^2 is within 0.1% of the cap, must now come
back "unphysical" -- and the ordinary outcomes must be unchanged.
Pure function, no engine.  Run:  python3 rmfit/tests/test_classify_gate.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.structure import classify  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def main():
    cap, tau = 0.5, 1000.0
    before = 930_000.0

    # 1. the round-3 case: huge gain, rho fine, but three channels ON the cap
    th = {"g[7/2-#11 p2 L3 S0.5]": 0.29, "g[7/2-#11 p5 L2 S2.5]": 0.5,
          "g[7/2-#11 p6 L1 S1.5]": 0.5, "g[7/2-#11 p8 L1 S2.5]": 0.5, "g[7/2-#11 p3 L3 S0.5]": 1e-35}
    out = classify(before, before - 320_000, th, 0.0, cap, tau, 0.08)
    check("theta2 exactly at the cap -> unphysical (theta2 path)", out == "unphysical", out)

    # 2. same, told through the polish's own at-bound list, theta2 just under the cap
    th2 = dict(th, **{"g[7/2-#11 p5 L2 S2.5]": 0.49, "g[7/2-#11 p6 L1 S1.5]": 0.49, "g[7/2-#11 p8 L1 S2.5]": 0.49})
    out = classify(before, before - 320_000, th2, 0.0, cap, tau, 0.08, at_bound=["g[7/2-#11 p8 L1 S2.5]"])
    check("width on a bound (bounds_hit) -> unphysical", out == "unphysical", out)

    # 3. theta2 within 0.1% of the cap -> unphysical; 1% below -> not
    out = classify(before, before - 50_000, {"a": cap * 0.9995}, 0.0, cap, tau, 0.1)
    check("theta2 within 0.1% of cap -> unphysical", out == "unphysical", out)
    out = classify(before, before - 50_000, {"a": cap * 0.99}, 0.0, cap, tau, 0.1)
    check("theta2 1% below cap, clean gain -> accept", out == "accept", out)

    # 4. the ordinary ladder is unchanged
    good = {"a": 0.39, "b": 0.093, "c": 0.024}          # the accepted 5/2+ 12.626, cap 1.0
    check("clean level -> accept", classify(before, before - 135_000, good, 0.0, 1.0, 5300, 0.07) == "accept")
    check("above cap -> blowup", classify(before, before - 1e5, {"a": 1.2}, 0.0, 1.0, tau, 0.1) == "blowup")
    check("all dead -> collapse", classify(before, before - 100, {"a": 1e-6}, 0.0, cap, tau, 0.1) == "collapse")
    check("moved > 50 keV -> relocated", classify(before, before - 1e5, good, 0.06, 1.0, tau, 0.1) == "relocated")
    check("gain below tau -> marginal", classify(before, before - 100, good, 0.0, 1.0, tau, 0.1) == "marginal")
    check("worse -> live-harmful", classify(before, before + 7715, good, 0.0, 1.0, tau, 1.62) == "live-harmful")
    check("rho > 0.5 -> redistributes", classify(before, before - 1e5, good, 0.0, 1.0, tau, 0.7) == "redistributes")
    check("empty at_bound tuple is inert", classify(before, before - 1e5, good, 0.0, 1.0, tau, 0.1, at_bound=()) == "accept")

    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all classify-gate checks passed")


if __name__ == "__main__":
    main()
