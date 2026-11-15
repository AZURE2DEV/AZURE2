#!/usr/bin/env python3
"""structure.classify with the level's physical width (2026-09-20): "collapse" means vanished,
not narrow; "narrow" = passes everything but is below the campaign's width floor."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.structure import classify

ok = True
def check(name, got, want):
    global ok
    good = got == want; print(f"  {'ok  ' if good else 'FAIL'}  {name}" + ("" if good else f"  -- got {got!r}, want {want!r}")); ok &= good

def main():
    tiny = {"n2 L1": 5.6e-5, "a0 L3": 4.8e-6}            # a narrow level: ~35 eV, Gamma_W 0.3-1 MeV
    base = dict(cap=1.0, tau=610.0, rho=0.21)
    check("no width given -> the old theta^2-only gate (backward compatible)", classify(305048, 304370, tiny, -0.0014, **base), "collapse")
    check("35 eV wide with tiny theta^2 is NOT a collapse", classify(305048, 304370, tiny, -0.0014, level_gamma_ev=35.0, **base), "accept")
    check("... and is 'narrow' under a 10 keV campaign floor", classify(305048, 304370, tiny, -0.0014, level_gamma_ev=35.0, narrow_width_ev=1e4, **base), "narrow")
    check("a genuinely vanished level is still a collapse", classify(305048, 304370, tiny, 0.0, level_gamma_ev=0.02, narrow_width_ev=1e4, **base), "collapse")
    broad = {"n1 L0": 0.155}
    check("a broad level above the floor is accepted", classify(306742, 303917, broad, -0.004, level_gamma_ev=1.06e5, narrow_width_ev=1e4, cap=1.0, tau=732.0, rho=0.40), "accept")
    check("narrow does not mask a marginal gain", classify(305048, 304900, tiny, 0.0, level_gamma_ev=35.0, narrow_width_ev=1e4, **base), "marginal")
    check("narrow does not mask redistribution", classify(305048, 304370, tiny, 0.0, level_gamma_ev=35.0, narrow_width_ev=1e4, cap=1.0, tau=610.0, rho=0.7), "redistributes")
    check("unphysical still comes before everything after blowup", classify(305048, 304370, {"a": 0.9995}, 0.0, level_gamma_ev=5e4, narrow_width_ev=1e4, **base), "unphysical")
    check("floor 0 / None = off", classify(305048, 304370, tiny, 0.0, level_gamma_ev=35.0, narrow_width_ev=0.0, **base), "accept")
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
