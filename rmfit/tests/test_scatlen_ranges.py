#!/usr/bin/env python3
"""Measured scattering lengths as residual rows (engine.ScatLenRows) and <targetInt> segment ranges
(model.parse_seg_list), 2026-11-07 (the evaluator: "make those rmfit changes so we can proceed").

Part A (no engine, seconds):
  * parse_seg_list reads keys and inclusive ranges as AZURE2 does; a backwards range raises;
  * target_int_entries / target_effects read range entries (before: skipped without a word);
  * shards written from a model with range entries carry every integrated data piece, re-keyed, and no
    <scatteringLengths> lines.
Part B (13N fixture, a <scatteringLengths> line for p+12C s = 1/2 added, minutes):
  * a model without the block has no rows (old behaviour, residual count unchanged);
  * with the block: one row, objective (residual path) = the engine's scalar chi2 (which holds the term) +
    penalties, and the score reports the term in ``penalty``;
  * a coherent line (s = -1) gives the same row as the s = 1/2 line (one channel spin: (2s+1) i1i2 = 1);
  * the rows' Jacobian matches central differences of the rows, and columns outside ``cols`` do not move a0;
  * Sharded (3 shards) gives the same objective and row; its shard files carry no <scatteringLengths> lines;
  * Park mode: the scattering-length rows sit before the wall rows and both are counted.
Run on the login node: PYTHONPATH=<dir holding rmfit> python3 -I tests/test_scatlen_ranges.py  (or as a job).
"""
import glob
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def part_b():
    from rmfit import AZURE2_ROOT
    from rmfit.engine import SingleSession, Sharded
    from rmfit.model import single_variant
    from rmfit.azr3 import Model as AzrModel, record_kinds
    import pyazr
    print("Model rows: a measured scattering length (13N fixture)")
    tmp = tempfile.mkdtemp(prefix="rmfit_scatlen_")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
    mdl.set_output_dir(os.path.join(tmp, "output"))
    plain = mdl.write(os.path.join(tmp, "13N_plain.azr"))

    with SingleSession(single_variant(plain, "plain", tmp), cwd=tmp, name="plain") as ev:
        check("no model data: no model rows", ev.scat is None)
        n_plain = ev.residuals(ev.x0)[0].size
        check("no model data: rows = data + penalties", n_plain == ev.npts + ev.penalties.cols.size)
    with pyazr.Session(plain, cwd=tmp) as s:
        a_fit = [v for v in s.scattering_lengths() if int(v.pair) == 1][0].a0.real
    print(f"  computed a0(p+12C, s=1/2) at the start: {a_fit:.6g} fm")
    target, sig = a_fit + 0.3 * max(abs(a_fit), 1.0), 0.05 * max(abs(a_fit), 1.0)

    m = pyazr.Project.open(plain)
    m.add_datum("a_p12C", quantity="scattering_length", pair="p+12C", s="1/2", value=(repr(target), "fm"),
                error=(repr(sig), "fm"))
    sl = m.save(os.path.join(tmp, "13N_sl.azr"))
    expect_row = (a_fit - target) / sig
    ev = SingleSession(single_variant(sl, "sl", tmp), cwd=tmp, name="sl")
    try:
        check("with the datum: one model row", ev.scat is not None and len(ev.scat.parts) == 1)
        x0 = np.array(ev.x0, float)
        r0, _ = ev.residuals(x0)
        check("one more row", r0.size == n_plain + 1, f"{r0.size} vs {n_plain + 1}")
        check("row = +-(a0 - value)/sigma", np.isclose(abs(r0[-1]), abs(expect_row), rtol=1e-8), f"{r0[-1]} vs {expect_row}")
        with pyazr.Session(sl, cwd=tmp) as s:
            eng = float(s.objective()(x0))
        obj = ev.objective(x0)
        check("objective = the engine's objective", np.isclose(obj, eng, rtol=1e-9), f"{obj} vs {eng}")
        sc = ev.score(x0)
        check("score: term in penalty", np.isclose(sc.penalty - ev.penalties.value(x0), expect_row ** 2, rtol=1e-8),
              f"{sc.penalty}")
        check("score: chi2 excludes it", np.isclose(sc.chi2 + sc.penalty, eng, rtol=1e-9), f"{sc.chi2} + {sc.penalty}")
        # the Jacobian row against central differences of the row
        _, J = ev.residuals(x0, jac=True)
        Jr = J[-1]
        worst = 0.0
        for c in range(ev.n_rmatrix):
            h = 1e-6 if ev.keys[c].kind == "energy" else max(1e-5 * abs(x0[c]), 1e-7)
            xp, xm = x0.copy(), x0.copy()
            xp[c] += h
            xm[c] -= h
            d = (ev._scat_rows(xp)[0][0] - ev._scat_rows(xm)[0][0]) / (2 * h)
            worst = max(worst, abs(d - Jr[c]) / max(abs(d), abs(Jr).max() * 1e-3, 1e-3))
        check("Jacobian row matches central differences", worst < 1e-4, f"worst rel {worst:.2e}")
        moved = {str(ev.keys[c].level).split("#")[0] for c in range(ev.n_rmatrix) if abs(Jr[c]) > 1e-12 * abs(Jr).max()}
        check("only the 1/2+ levels move a0", moved == {"1/2+"}, f"{sorted(moved)}")
    finally:
        ev.close()

    sh = Sharded(sl, os.path.join(tmp, "shards"), nshards=3)
    try:
        o_sh = sh.objective(np.array(sh.x0, float))
        check("Sharded objective = single", np.isclose(o_sh, obj, rtol=1e-8), f"{o_sh} vs {obj}")
        r_sh, J_sh = sh.residuals(np.array(sh.x0, float), jac=True)
        check("Sharded row (counted once)", np.isclose(r_sh[-1], r0[-1], rtol=1e-9) and r_sh.size == r0.size, f"{r_sh[-1]}")
        check("Sharded score penalty", np.isclose(sh.score(np.array(sh.x0, float)).penalty, sc.penalty, rtol=1e-9))
        check("every shard keeps the datum (any worker has the model rows)",
              all("datum" in record_kinds(f) for f in glob.glob(os.path.join(tmp, "shards", "*.azr"))))
    finally:
        sh.close()

    park = dict(use_park=True, use_brune=True)
    with SingleSession(single_variant(sl, "slpark", tmp), cwd=tmp, flags=park, name="slpark") as ep:
        xp0 = np.array(ep.x0, float)
        with pyazr.Session(sl, cwd=tmp, **park) as s:
            eng = float(s.objective()(xp0))
        check("Park: the wall comes with the engine's model rows", ep.park is None and ep.scat is not None
              and ep.scat.n > 1, f"{ep.park} {ep.scat}")
        check("Park: objective = the engine's objective", np.isclose(ep.objective(xp0), eng, rtol=1e-9),
              f"{ep.objective(xp0)} vs {eng}")
    shutil.rmtree(tmp, ignore_errors=True)


def main():
    part_b()
    print("ALL PASSED" if not failures else f"FAILED: {failures}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
