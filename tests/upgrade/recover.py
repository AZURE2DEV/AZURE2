#!/usr/bin/env python3
"""Recovery check: a fit with the present engine returns the true parameters of T1 and T3.

    python3 tests/upgrade/recover.py [--binary AZURE2] [--threads N] [--max-pull 4] [T1 T3]

For each project the true model (T<n>.azr) is calculated once (menu 1) to obtain the true reduced-width vector
(param.par: level energies E_x in MeV, reduced width amplitudes in MeV^1/2, ANCs as reduced widths of the bound
channels), and the start model (T<n>_start.azr: every free parameter moved, energies by +4 keV (T1) or +10 keV (T3),
widths and ANCs by a factor 1.15) is fitted with MIGRAD on the analytic gradient (menu 2, --use-gradient).  The fit
passes when

  (a) MIGRAD reports a valid minimum;
  (b) chi^2 at the fitted parameters <= chi^2 at the truth + 1e-3 (the fit is no worse than the truth on the same
      noisy data, so it is not stuck in a secondary minimum);
  (c) every free parameter satisfies |p_fit - p_true| <= max_pull * sigma_p, sigma_p the MIGRAD uncertainty in
      param.sav (Brune basis, the engine's own parameter vector), max_pull = 4 by default;
  (d) the same fit to noise-free data (the committed data files with each value replaced by the truth, the
      uncertainties kept, and the measured scattering lengths of T3 (scattering_length data of <modelData>) replaced
      by the truth's) returns every free parameter to |p_fit - p_true| <= 0.05 sigma_p.

(c) tests the fit on the committed noisy data; with 13 (T1) and 8 (T3) parameters a 3 sigma excursion is expected in a
few per cent of noise realizations, hence 4.  An ensemble of 8 noise seeds for T1 (2026-11-06) gave pull means within
+-0.7 and standard deviations 0.7-1.5, so the MIGRAD uncertainties are the right scale.  (d) has no noise and is the
sharp test: a wrong gradient, a parameter the start model does not reach or a secondary minimum fails it.

The binary is --binary, $AZURE2_BIN, $AZURE2_BINARY, else the newest build*/src/AZURE2.  Exit status 1 on failure.
"""
import argparse
import importlib.util
import os
import re
import shutil
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))
import upgrade_project as up  # noqa: E402


def load_generator(name):
    spec = importlib.util.spec_from_file_location(f"{name}_generate", os.path.join(HERE, name, "generate.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def copy_project(name, work):
    src = os.path.join(HERE, name)
    for f in (f"{name}.azr", f"{name}_start.azr", "temperatures.dat"):
        shutil.copy(os.path.join(src, f), work)
    shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))


def write_noise_free(work, project):
    """Replace every data value by the model at the current output (the truth calculation), in the data frame:
    AZUREOut column 4 is the c.m. model, column 6 the c.m. datum, so lab model = col4 * datum_lab / col6."""
    rows = up.segment_outputs(os.path.join(work, "output"), project.data)
    for seg, r in zip(project.data, rows):
        path = os.path.join(work, "data", seg.name + ".dat")
        lines = [ln.split() for ln in open(path) if ln.strip() and not ln.startswith("#")]
        if len(lines) != len(r):
            raise RuntimeError(f"{seg.name}: {len(lines)} data rows, {len(r)} output rows")
        with open(path, "w") as f:
            for x, rr in zip(lines, r):
                f.write(f"{x[0]} {x[1]} {rr[3] * float(x[2]) / rr[5]:.12e} {x[3]}\n")


def noise_free_scattering_lengths(work, azr):
    """Set the scattering-length data of <modelData> in azr (in work) to the truth: the lengths in
    output/scatteringLengths.out of the truth calculation (s-wave, J = s), matched by pair label and channel spin."""
    a0 = {}
    for ln in open(os.path.join(work, "output", "scatteringLengths.out")):
        p = ln.split()
        if len(p) >= 6 and p[0].isdigit() and p[2][-1] in "+-":
            a0[(p[1], p[2][:-1])] = (p[4], p[5])
    path = os.path.join(work, azr)
    text = open(path).read()
    if "<modelData>" not in text:
        return

    def fix(m):
        records = re.split(r"(?m)^(?=datum )", m.group(1))
        out = []
        for rec in records:
            field = dict(re.findall(r"(?m)^ +(\w+) +(\S+)", rec))
            if field.get("quantity") == "scattering_length":
                re_, im_ = a0[(field["pair"], field["s"])]   # integer s in T3
                rec = re.sub(r"(?m)^( +value +)\S+", lambda v: v.group(1) + re_, rec)
                rec = re.sub(r"(?m)^( +imag +)\S+", lambda v: v.group(1) + im_, rec)
            out.append(rec)
        return "<modelData>" + "".join(out) + "</modelData>"
    text = re.sub(r"<modelData>(.*?)</modelData>", fix, text, flags=re.S)
    open(path, "w").write(text)


def check_pulls(truth, fit, free, max_pull):
    failures = 0
    print(f"        {'parameter':12s} {'truth':>14s} {'fit':>14s} {'sigma':>11s} {'pull':>9s}")
    for n in free:
        t, (f, s) = truth[n][0], fit[n]
        pull = (f - t) / s if s > 0 else float("inf")
        good = abs(pull) <= max_pull
        failures += not good
        print(f"  {'ok  ' if good else 'FAIL'}  {n:12s} {t:14.7e} {f:14.7e} {s:11.3e} {pull:+9.3f}")
    return failures


def recover(binary, name, threads, max_pull):
    gen = load_generator(name)
    free = [n for n, ref in gen.model().param_names() if ref is not None]
    project = gen.model()
    work = tempfile.mkdtemp(prefix=f"{name}_recover_")
    try:
        copy_project(name, work)
        up.run(binary, work, f"{name}.azr", up.CALC, threads=threads)
        truth = up.read_param_file(os.path.join(work, "output", "param.par"))
        chi_truth, _ = up.read_chi2(os.path.join(work, "output", "chiSquared.out"))
        t0 = time.time()
        log = up.run(binary, work, f"{name}_start.azr", up.FIT, flags=["--use-gradient"], threads=threads)
        seconds = time.time() - t0
        fit = up.read_param_file(os.path.join(work, "output", "param.sav"))
        chi_fit, _ = up.read_chi2(os.path.join(work, "output", "chiSquared.out"))
        # (d) noise-free: values from the truth calculation, then the same fit
        up.run(binary, work, f"{name}.azr", up.CALC, threads=threads)
        write_noise_free(work, project)
        noise_free_scattering_lengths(work, f"{name}_start.azr")
        log0 = up.run(binary, work, f"{name}_start.azr", up.FIT, flags=["--use-gradient"], threads=threads)
        fit0 = up.read_param_file(os.path.join(work, "output", "param.sav"))
        chi_fit0, _ = up.read_chi2(os.path.join(work, "output", "chiSquared.out"))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    failures = 0
    for label, lg in (("noisy data", log), ("noise-free data", log0)):
        m = re.search(r"MIGRAD converged: (valid|INVALID)[^\n]*", lg)
        valid = bool(m) and m.group(1) == "valid"
        failures += not valid
        print(f"{name} {label}: {'ok  ' if valid else 'FAIL'} {m.group(0) if m else 'no MIGRAD summary in the log'}")
    print(f"  fit time {seconds:.1f} s on {threads} threads")
    ok = chi_fit <= chi_truth + 1e-3
    failures += not ok
    print(f"  {'ok  ' if ok else 'FAIL'}  (b) chi2 at the fit {chi_fit:.3f} <= chi2 at the truth {chi_truth:.3f}")
    print(f"  (c) noisy data, |pull| <= {max_pull:g}")
    failures += check_pulls(truth, fit, free, max_pull)
    print(f"  (d) noise-free data (chi2 at the fit {chi_fit0:.3e}), |pull| <= 0.05")
    failures += check_pulls(truth, fit0, free, 0.05)
    return failures


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("projects", nargs="*", default=["T1", "T3"])
    ap.add_argument("--binary")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--max-pull", type=float, default=4.0)
    args = ap.parse_args()
    binary = up.find_binary(args.binary)
    failures = sum(recover(binary, n, args.threads, args.max_pull) for n in args.projects)
    print(f"recovery: {'FAIL' if failures else 'ok'} ({failures} problems)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
