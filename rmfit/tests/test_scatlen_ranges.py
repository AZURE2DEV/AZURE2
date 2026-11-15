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


def part_a():
    from rmfit.model import parse_seg_list, target_int_entries, target_effects, parse_target_effect
    print("Part A: <targetInt> segment ranges")
    check("plain list", parse_seg_list('"8,9"') == [8, 9])
    check("ranges", parse_seg_list('"41,47-50,54-56"') == [41, 47, 48, 49, 50, 54, 55, 56])
    check("spaces", parse_seg_list(' "1, 3 - 4" ') == [1, 3, 4])
    try:
        parse_seg_list('"5-3"')
        check("backwards range raises", False)
    except ValueError:
        check("backwards range raises", True)
    # a minimal <targetInt> block in an archive model's format
    line_r = '1  "41,47-50,54-56"  30  0  0  1  4.6e18  "a0*exp(-a1*x)"  2  1  1  0  0  1  ""  0'
    line_p = '1  "40"  30  0  0  1  4.6e18  "a0*exp(-a1*x)"  2  1  1  0  0  1  ""  0'
    text = f"<targetInt>\n{line_r}\n{line_p}\n</targetInt>\n"
    ent = target_int_entries(text)
    check("target_int_entries reads the range entry", len(ent) == 2 and ent[0][1] == [41, 47, 48, 49, 50, 54, 55, 56],
          f"{ent}")
    check("target_effects reads the range entry", [e.keys for e in target_effects(text)][:1] == [[41, 47, 48, 49, 50, 54, 55, 56]],
          f"{[e.keys for e in target_effects(text)]}")
    check("parse_target_effect keys", parse_target_effect(line_r).keys[-1] == 56)

    # shards of a real archive model (no engine: plan + write only)
    base = os.path.join(os.environ.get("RMFIT_ARCHIVE", os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))), "an archived campaign")  # archive model; RMFIT_ARCHIVE
    if not os.path.exists(base):
        print("  (archive model not found, shard check skipped)")
        return
    from rmfit.model import build_shards
    tmp = tempfile.mkdtemp(prefix="rmfit_ranges_")
    try:
        text = open(base, encoding="latin-1").read()
        integ = set()
        for en, keys, line in target_int_entries(text):
            if en:
                integ.update(keys)
        n_range = sum(1 for en, keys, line in target_int_entries(text)
                      if en and "-" in line.split('"')[1])
        check("archive model: the 9 enabled range entries read", n_range == 9, f"{n_range}")
        shards = build_shards(base, 8, tmp)
        covered, sl_lines = set(), 0
        for s in shards:
            st = open(s.azr, encoding="latin-1").read()
            for en, keys, line in target_int_entries(st):
                for k in keys:
                    p = s.pieces[k - 1]
                    covered.add(p.base_key)
            if "<scatteringLengths>" in st:
                a, b = st.index("<scatteringLengths>"), st.index("</scatteringLengths>")
                sl_lines += sum(1 for l in st[a:b].splitlines()[1:] if l.strip() and not l.strip().startswith("#"))
        data_keys = {p.base_key for s in shards for p in s.pieces if p.kind == "data"}
        missing = sorted((integ & data_keys) - covered)
        check("every integrated data piece keeps its target effect in its shard", not missing, f"missing {missing[:20]}")
        check("no shard carries <scatteringLengths> lines", sl_lines == 0, f"{sl_lines}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def part_b():
    from rmfit import AZURE2_ROOT
    from rmfit.engine import SingleSession, Sharded
    from rmfit.model import single_variant
    from pyazr import AzrModel
    print("Part B: scattering-length rows (13N fixture)")
    fix = next(d for d in (AZURE2_ROOT, os.path.expanduser("~/AZURE2")) if os.path.isdir(os.path.join(d, "tests", "13N")))
    tmp = tempfile.mkdtemp(prefix="rmfit_scatlen_")
    shutil.copytree(os.path.join(fix, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(fix, "tests", "13N", "13N.azr"))
    mdl.set_output_dir(os.path.join(tmp, "output"))
    plain = mdl.write(os.path.join(tmp, "13N_plain.azr"))

    with SingleSession(single_variant(plain, "plain", tmp), cwd=tmp, name="plain") as ev:
        check("no block: no rows", ev.scat is None)
        n_plain = ev.residuals(ev.x0)[0].size
        check("no block: rows = data + penalties", n_plain == ev.npts + ev.penalties.cols.size)
        a_fit = [r for r in ev.worker.call("scatlen", ev.x0) if int(r[0]) == 1 and abs(r[1] - 0.5) < 1e-6][0][3]
    print(f"  computed a0(p+12C, s=1/2) at the fit: {a_fit:.6g} fm")
    target, sig = a_fit + 0.3 * max(abs(a_fit), 1.0), 0.05 * max(abs(a_fit), 1.0)

    def with_line(s, name):
        m = AzrModel.from_file(plain)
        m.add_scattering_length(1, s, target, sig)
        return m.write(os.path.join(tmp, f"13N_{name}.azr"))

    sl = with_line(0.5, "sl")
    coh = with_line(-1, "coh")
    expect_row = (a_fit - target) / sig
    ev = SingleSession(single_variant(sl, "sl", tmp), cwd=tmp, name="sl")
    try:
        check("with block: rows built", ev.scat is not None and len(ev.scat.parts) == 1)
        x0 = np.array(ev.x0, float)
        r0, _ = ev.residuals(x0)
        check("one more row", r0.size == n_plain + 1, f"{r0.size} vs {n_plain + 1}")
        check("row = (a0 - value)/sigma", np.isclose(r0[-1], expect_row, rtol=1e-9), f"{r0[-1]} vs {expect_row}")
        eng = float(ev.worker.call("eval", x0)) + ev.penalties.value(x0)
        obj = ev.objective(x0)
        check("objective = engine scalar chi2 (holds the term) + penalties", np.isclose(obj, eng, rtol=1e-9),
              f"{obj} vs {eng}")
        sc = ev.score(x0)
        check("score: term in penalty", np.isclose(sc.penalty - ev.penalties.value(x0), expect_row ** 2, rtol=1e-9),
              f"{sc.penalty}")
        check("score: chi2 excludes it", np.isclose(sc.chi2 + sc.penalty, eng, rtol=1e-9), f"{sc.chi2} + {sc.penalty}")
        # Jacobian rows vs central differences of the row; columns outside cols must not move a0
        _, J = ev.residuals(x0, jac=True)
        Jr = J[-1]
        worst, outside = 0.0, 0.0
        cols = set(ev.scat.cols)
        for c in range(ev.n_rmatrix):
            h = 1e-6 if ev.keys[c].kind == "energy" else max(1e-5 * abs(x0[c]), 1e-7)
            xp, xm = x0.copy(), x0.copy()
            xp[c] += h
            xm[c] -= h
            d = (ev._scat_rows(xp)[0][0] - ev._scat_rows(xm)[0][0]) / (2 * h)
            if c in cols:
                worst = max(worst, abs(d - Jr[c]) / max(abs(d), 1e-3))
            else:
                outside = max(outside, abs(d))
        check("Jacobian row matches central differences", worst < 1e-4, f"worst rel {worst:.2e}")
        check("columns outside cols do not move a0", outside < 1e-10, f"max |d| {outside:.2e}")
        check("cols = the 1/2+ levels only", all(str(ev.keys[c].level).startswith("1/2+") for c in cols) and cols,
              f"{sorted({str(ev.keys[c].level) for c in cols})}")
    finally:
        ev.close()

    with SingleSession(single_variant(coh, "coh", tmp), cwd=tmp, name="coh") as ec:
        rc, _ = ec.residuals(ec.x0)
        check("coherent line (s=-1) = s=1/2 line for one channel spin", np.isclose(rc[-1], expect_row, rtol=1e-9),
              f"{rc[-1]} vs {expect_row}")

    sh = Sharded(sl, os.path.join(tmp, "shards"), nshards=3)
    try:
        o_sh = sh.objective(np.array(sh.x0, float))
        check("Sharded objective = single", np.isclose(o_sh, obj, rtol=1e-8), f"{o_sh} vs {obj}")
        r_sh, J_sh = sh.residuals(np.array(sh.x0, float), jac=True)
        check("Sharded row", np.isclose(r_sh[-1], expect_row, rtol=1e-9), f"{r_sh[-1]}")
        check("Sharded score penalty", np.isclose(sh.score(np.array(sh.x0, float)).penalty, sc.penalty, rtol=1e-9))
        lines = 0
        for f in glob.glob(os.path.join(tmp, "shards", "*.azr")):
            lines += len([e for e in AzrModel.from_file(f).scattering_lengths if e.active])
        check("shard files carry no scattering-length lines", lines == 0, f"{lines}")
    finally:
        sh.close()

    park = dict(use_park=True, use_brune=True)
    with SingleSession(single_variant(sl, "slpark", tmp), cwd=tmp, flags=park, name="slpark") as ep:
        rp, _ = ep.residuals(ep.x0)
        L = len(ep.park.labels)
        check("Park: scattering-length row just before the wall rows", rp.size == n_plain + 1 + L, f"{rp.size}")
        eng = float(ep.worker.call("eval", np.array(ep.x0, float))) + ep.penalties.value(np.array(ep.x0, float))
        check("Park: objective = engine chi2 (wall + term) + penalties",
              np.isclose(ep.objective(np.array(ep.x0, float)), eng, rtol=1e-9),
              f"{ep.objective(np.array(ep.x0, float))} vs {eng}")
        check("Park: scat_value reads the right row", np.isclose(ep.scat_value(rp), rp[-L - 1] ** 2, rtol=1e-12))
    shutil.rmtree(tmp, ignore_errors=True)


def main():
    part_a()
    if "--no-engine" not in sys.argv:
        part_b()
    print("ALL OK" if not failures else f"FAILURES: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
