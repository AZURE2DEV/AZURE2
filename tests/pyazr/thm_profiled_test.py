#!/usr/bin/env python3
"""pyazr treats a THM segment's free norm the way the CLI does: profiled.

A THM (HOES) excitation function has an arbitrary overall scale, so when its
norm is free the engine does not fit it: every evaluation sets it to the
chi-squared optimum n* = S_mm/S_md and counts S_dd - S_md^2/S_mm, with no norm
penalty, and the norm is not a Minuit parameter (ESegment::IsProfiledNorm).
pyazr used to evaluate such a model at the norm written in the file (212924
instead of 2137.83 on tests/7Li_p_a), list the profiled norm in
parameter_info() (so ``m.parameters`` failed), refuse residual_jacobian, and
leave the THM points out of chi2_and_grad altogether.

Checked here on tests/7Li_p_a (one THM segment, free norm):

  1. calculate_chi2_rwa equals a fresh CLI run (rel 1e-6), and the CLI still
     lands on the recorded pin; parameters / parameter_info / params_rwa
     agree in count and the THM norm is not among them; residuals and
     segment_chi2 add up to chi2; segment_norms reports n*;
     write_output_files reproduces the CLI's chiSquared.out and AZUREOut_*.
  2. On a copy with all but one level fixed (so the THM finite differences
     are cheap): residual_jacobian against central differences of the
     residuals (rel 1e-4), chi2_and_grad against 2 J^T r and against
     differences of chi2, model_gradients against differences of the model.

Needs the compiled engine and an AZURE2 binary; skips cleanly without them.

Run from anywhere:  python3 tests/pyazr/thm_profiled_test.py
"""
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "7Li_p_a")
PIN = 2138.48          # tests/7Li_p_a/expected/chiSquared.out (re-pinned after the dS/dE step fix)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2, AzrModel
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)

binary = os.environ.get("AZURE2_BIN")
if not binary:
    cands = [c for c in glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))
             if os.path.isfile(c) and os.access(c, os.X_OK)]
    binary = max(cands, key=os.path.getmtime) if cands else None
if binary is None:
    print("skip: no AZURE2 binary to compare against")
    sys.exit(77)


def fresh_copy(dst):
    shutil.copytree(SOURCE, dst)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(dst, junk), ignore_errors=True)
        os.makedirs(os.path.join(dst, junk))
    for f in glob.glob(os.path.join(dst, "run.log")):
        os.remove(f)


def output_files(work):
    out = os.path.join(work, "output")
    return {os.path.basename(f): open(f).read()
            for f in sorted(glob.glob(os.path.join(out, "*")))
            if os.path.isfile(f) and not os.path.basename(f).startswith("intEC")}


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


with tempfile.TemporaryDirectory() as tmp:
    cli_dir = os.path.join(tmp, "cli")
    py_dir = os.path.join(tmp, "py")
    fresh_copy(cli_dir)
    fresh_copy(py_dir)

    # -- the CLI, fresh -------------------------------------------------------
    proc = subprocess.run([binary, "--no-gui", "--no-readline", "7Li_p_a.azr"],
                          cwd=cli_dir, input="1\n\n\n7\n", text=True,
                          capture_output=True, timeout=1800)
    log = proc.stdout.replace("\r", "\n")
    hit = re.findall(r"Total Chi-Squared:\s*([0-9.eE+-]+)", log)
    from_cli = output_files(cli_dir)
    if not hit or "chiSquared.out" not in from_cli:
        print("FAIL: the CLI run produced no chi-squared")
        print(log[-2000:])
        sys.exit(1)
    cli_chi2 = float(hit[-1])
    cli_norm = float(from_cli["chiSquared.out"].splitlines()[1].split(",")[3])
    print(f"CLI: chi2 {cli_chi2!r}, norm {cli_norm}")

    print("1. the profiled THM norm, as the CLI has it")
    # run_tests.sh's tolerance: the pin is the regression suite's, and moves
    # within it when e.g. the shift-function derivative gets more accurate.
    check("the CLI still lands on the pin (rel 1e-3)", rel(cli_chi2, PIN) < 1e-3,
          f"{cli_chi2} vs {PIN}")
    with azure2(os.path.join(py_dir, "7Li_p_a.azr"), cwd=py_dir) as m:
        x = np.asarray(m.params_rwa, float)
        chi2 = float(np.sum(m.calculate_chi2_rwa(x)))
        print(f"        pyazr chi2 {chi2!r}")
        check("calculate_chi2_rwa == CLI (rel 1e-6)", rel(chi2, cli_chi2) < 1e-6,
              f"{chi2} vs {cli_chi2}")

        n_info = np.asarray(m.sess.parameter_info()).size // 16
        check("parameter_info has one record per parameter",
              n_info == len(m.fixed_params), f"{n_info} vs {len(m.fixed_params)}")
        try:
            ps = m.parameters
            check("m.parameters builds", True)
            check("free parameters == params_rwa",
                  len(ps.free) == x.size, f"{len(ps.free)} vs {x.size}")
            check("the profiled THM norm is not a parameter",
                  not any(p.segment_key == 1 for p in ps.norms),
                  [p.name for p in ps.norms])
        except Exception as err:
            check("m.parameters builds", False, f"{type(err).__name__}: {err}")
        check("norm_indices is empty (the only norm is profiled)",
              len(m.norm_indices()) == 0, m.norm_indices())

        r = m.residuals(x)
        check("one residual per point", r.size == sum(len(e) for e in m.energies),
              f"{r.size}")
        check("sum r^2 == chi2", rel(float(np.sum(r ** 2)), chi2) < 1e-9,
              f"{np.sum(r ** 2)} vs {chi2}")
        seg = m.segment_chi2(x)
        check("segment_chi2 sums to chi2", rel(float(seg.sum()), chi2) < 1e-9,
              f"{seg.sum()} vs {chi2}")
        pen = m.penalties(x)
        check("no norm penalty on a profiled norm", float(pen["norm"].sum()) == 0.0,
              pen["norm"])
        check("objective == chi2 (no penalties here)",
              rel(m.objective(x), chi2) < 1e-12, f"{m.objective(x)}")
        norms = m.segment_norms(x)
        check("segment_norms gives n* (the CLI's norm, 6 digits)",
              rel(float(norms[0]), cli_norm) < 1e-5, f"{norms[0]} vs {cli_norm}")
        d = np.asarray(m.cross[0]); e = np.asarray(m.cross_err[0])
        mod = np.asarray(m.calculate_rwa(x)[0])
        nstar = np.sum(mod * mod / e ** 2) / np.sum(mod * d / e ** 2)
        check("n* = S_mm/S_md", rel(float(norms[0]), nstar) < 1e-10,
              f"{norms[0]} vs {nstar}")

        m.write_output_files(x)
    from_py = output_files(py_dir)
    owned = {"chiSquared.out"} | {n for n in from_cli if n.startswith("AZUREOut_")}
    check("write_output_files wrote the CLI's data files", owned <= set(from_py),
          f"missing {sorted(owned - set(from_py))}")
    for name in sorted(owned & set(from_py)):
        check(f"{name} is byte-identical to the CLI's", from_py[name] == from_cli[name])

    # -- derivatives, on a copy with only one level free -----------------------
    print("\n2. derivatives through the profiled norm (2+ 20.1 MeV level free)")
    small = os.path.join(tmp, "small")
    fresh_copy(small)
    model = AzrModel.from_file(os.path.join(small, "7Li_p_a.azr"))
    for lv in model.levels:
        keep = lv.jpi == "2+" and abs(lv.energy - 20.1) < 1e-6
        lv.set_fixed(not keep)
        for c in lv.channels:
            c.channel_fixed = not keep
    model.write(os.path.join(small, "small.azr"))

    with azure2(os.path.join(small, "small.azr"), cwd=small) as m:
        x = np.asarray(m.params_rwa, float)
        print(f"        {x.size} free parameters: "
              + ", ".join(p.name for p in m.parameters.free))
        check("a handful of free parameters", 2 <= x.size <= 12, x.size)
        chi2 = float(np.sum(m.calculate_chi2_rwa(x)))
        check("same chi2 as the full model (same values)", rel(chi2, cli_chi2) < 1e-6,
              f"{chi2}")

        r, J = m.residual_jacobian(x)
        check("residual_jacobian no longer refuses THM", J.shape == (r.size, x.size),
              J.shape)
        check("its residuals are the profiled ones", rel(float(np.sum(r ** 2)), chi2) < 1e-9,
              f"{np.sum(r ** 2)}")

        energy = [p.free_index for p in m.parameters.free if p.kind == "energy"]
        widths = [p.free_index for p in m.parameters.free if p.kind == "width"]
        cols = energy[:1] + widths[:2]
        for c in cols:
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (m.residuals(xp) - m.residuals(xm)) / (2 * h)
            err = np.linalg.norm(J[:, c] - fd) / max(np.linalg.norm(fd), 1e-300)
            check(f"J[:, {c}] ({m.parameters.free[c].name}) vs central differences "
                  f"(rel {err:.1e})", err < 1e-4)

        c2, g = m.chi2_and_grad(x)
        check("chi2_and_grad value == chi2", rel(c2, chi2) < 1e-9, f"{c2} vs {chi2}")
        gj = 2.0 * J.T @ r
        check("gradient == 2 J^T r", np.allclose(g, gj, rtol=1e-8,
                                                 atol=1e-10 * np.max(np.abs(gj))),
              f"max diff {np.max(np.abs(g - gj))}")
        for c in cols:
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (m.calculate_chi2_rwa(xp)[0] - m.calculate_chi2_rwa(xm)[0]) / (2 * h)
            check(f"d chi2/d x[{c}] vs central differences (rel {rel(g[c], fd):.1e})",
                  rel(g[c], fd) < 1e-4, f"{g[c]} vs {fd}")

        G = m.model_gradients(x)[0]
        check("model_gradients covers the THM points", G.shape == (r.size, x.size),
              G.shape)
        for c in cols[:2]:
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (np.asarray(m.calculate_rwa(xp)[0]) - np.asarray(m.calculate_rwa(xm)[0])) / (2 * h)
            err = np.linalg.norm(G[:, c] - fd) / max(np.linalg.norm(fd), 1e-300)
            check(f"d model/d x[{c}] vs central differences (rel {err:.1e})", err < 1e-4)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM-profiling checks passed")
