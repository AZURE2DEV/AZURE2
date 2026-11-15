#!/usr/bin/env python3
"""classify(..., broad_width_ev): the 'broad' verdict (2026-09-22, the evaluator) -- mirror of 'narrow'.
A level that passes every test but ends wider than BoundsPolicy.max_level_width_ev is "broad" and
not adopted.  Off by default; ordering of every other verdict unchanged.  Pure.
Run:  python3 rmfit/tests/test_classify_broad.py"""
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.normpath(os.path.join(HERE,"..","..")))
from rmfit.structure import classify
from rmfit.local import BoundsPolicy
failures=[]
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok: failures.append(name)
def main():
    good={"a":0.02,"b":0.01}; B=930_000.0; A=B-50_000; tau=1000.0
    check("policy field exists, default off", getattr(BoundsPolicy(), "max_level_width_ev", None) == 0.0)
    check("default: no broad verdict", classify(B,A,good,0.0,0.5,tau,0.1,level_gamma_ev=2e6) == "accept")
    check("wider than the ceiling -> broad", classify(B,A,good,0.0,0.5,tau,0.1,level_gamma_ev=2e6,broad_width_ev=5e5) == "broad")
    check("narrower than the ceiling -> accept", classify(B,A,good,0.0,0.5,tau,0.1,level_gamma_ev=3e5,broad_width_ev=5e5) == "accept")
    check("width unknown -> accept (no verdict without a width)", classify(B,A,good,0.0,0.5,tau,0.1,level_gamma_ev=None,broad_width_ev=5e5) == "accept")
    check("marginal still wins over broad", classify(B,B-100,good,0.0,0.5,tau,0.1,level_gamma_ev=2e6,broad_width_ev=5e5) == "marginal")
    check("unphysical still wins over broad", classify(B,A,{"a":0.5},0.0,0.5,tau,0.1,level_gamma_ev=2e6,broad_width_ev=5e5) == "unphysical")
    check("splits still wins over broad", classify(B,A,good,0.0,0.5,tau,0.1,level_gamma_ev=2e6,broad_width_ev=5e5,splits=True) == "splits")
    check("narrow still wins over broad when both set (impossible in practice)", classify(B,A,good,0.0,0.5,tau,0.1,level_gamma_ev=100.0,narrow_width_ev=1e3,broad_width_ev=5e5) == "narrow")
    print(); 
    if failures: print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all classify-broad checks passed")
if __name__=="__main__": main()
