#!/usr/bin/env python3
"""T2 recovery check: a MIGRAD fit of T2.azr must recover the parameters the data were made with.

Runs, in a temporary directory:
  1. a calculation of T2_true.azr (the truth), to read the true level parameters in the units of
     param.sav (Brune basis: E in MeV, reduced-width amplitudes in MeV^1/2);
  2. the fit of T2.azr (AZURE2 fit) with --use-lm, started off the truth.  T2's THM segments have no
     analytic residual Jacobian in today's engine, so --use-lm falls back to MIGRAD with Minuit's
     numerical gradient (about 1100 calls: 585 s at 8 threads with exact Coulomb functions, 72 s with a
     Coulomb table, the same parameters to the digits of param.sav).  Without a flag the CLI switches on the analytic-gradient
     MIGRAD for a project with THM segments; on T2 that path spends minutes per gradient call and
     did not converge in 10 min at 8 threads (2026-11-06), so it is not used here;
and compares every free parameter p with its true value t through the pull (p - t)/dp, dp the
fit's own uncertainty in param.sav.

Tolerances:
  noisy data (the committed data/):    |pull| <= 3 for every free parameter, chi2/N <= 1.5;
  --asimov (noiseless data from make_t2.py --asimov):  |pull| <= 0.5.  The Asimov optimum is not
     exactly the truth, because the priors are centred on nominal values, not on the truth; the
     shift this causes is (offset/sigma_prior)(dp/sigma_prior) in units of dp, below 0.5 here.

    AZURE2_BINARY=/path/to/AZURE2 python3 tests/upgrade/T2/check_recovery.py [--asimov] [--threads 8]
          [--keep DIR] [--json FILE]

Exit status 0 if every check passes.  numpy only (and make_t2.py for --asimov).
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = ["T2.azr", "T2_true.azr", "cov_diff.txt", "prior_cov.f3.txt", "data"]

# Free parameters: param.sav name (format 3; the old name in the comment) -> how the truth is found ("level" = from
# the truth run's param.par, else the key of truth.json's "nuisance" table).
FREE = {
    "energy_1": "level", "width_1_3": "level",
    "width_2_1": "level", "width_2_4": "level",
    "energy_5": "level", "width_5_3": "level",
    "s01_int_thin.norm": "s1_norm",                     # segment_1_norm
    "s01_int_thin.shift.offset": "s1_shift",            # segment_1_energy_shift
    "D.value": "D_norm",                                # segment_2_norm (normalization group D)
    "G.offset": "G_shift",                              # segment_2_energy_shift (shift group G)
    "G.a0": "G_fn_a0",                                  # segment_2_energy_shift_fn0
    "s04_thick.norm": "s4_norm",                        # segment_4_norm
    "s04_thick.shift.sqrt": "s4_sqrt",                  # segment_4_energy_shift_sqrt
    "tgt_s05_implanted.areal_density": "s5_density",    # density_4
    "s06_gas.norm": "s6_norm",                          # segment_6_norm
    "s07_profile.norm": "s7_norm",                      # segment_7_norm
    "s08_sum.norm": "s8_norm",                          # segment_8_norm
    "natN.isotope.15N": "abundance",                    # abundance_natN_8_16
    "s13_thm_oes.norm": "s13_norm",                     # segment_13_norm
}


def read_params(path):
    out = {}
    for line in open(path):
        f = line.split()
        if len(f) >= 3 and not f[0].startswith("#"):
            try:
                out[f[0]] = (float(f[1]), float(f[2]))
            except ValueError:
                pass
    return out


def chi2_total(path):
    for line in open(path):
        if line.startswith("Total-Chi-Squared:"):
            f = line.split()
            return float(f[1]), int(f[f.index("Total-N:") + 1])
    raise RuntimeError("no Total-Chi-Squared in " + path)


def run(binary, workdir, azr, command, threads, extra=()):
    for d in ("output", "checks"):
        os.makedirs(os.path.join(workdir, d), exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), LC_ALL="C")
    t0 = time.time()
    r = subprocess.run([binary, command, "--no-gui", *extra, azr], stdin=subprocess.DEVNULL, cwd=workdir, env=env,
                       capture_output=True, text=True, timeout=3600)
    open(os.path.join(workdir, azr + ".log"), "w").write(r.stdout + r.stderr)
    return time.time() - t0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--binary", default=os.environ.get("AZURE2_BINARY", ""))
    ap.add_argument("--asimov", action="store_true", help="fit noiseless data (made by make_t2.py --asimov)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--keep", help="copy the work directory here")
    ap.add_argument("--engine-arg", action="append", default=[], help="extra AZURE2 option for the fit (repeatable)")
    ap.add_argument("--json", help="write the numbers to this file")
    args = ap.parse_args()
    found = shutil.which(args.binary) if args.binary else None  # which(): on Windows the path may lack its .exe
    if not found:
        sys.exit("set AZURE2_BINARY or pass --binary path/to/AZURE2")
    binary = os.path.abspath(found)
    work = tempfile.mkdtemp(prefix="t2rec.")
    fail = 0

    def check(name, ok, detail):
        nonlocal fail
        print(f"  {'ok  ' if ok else 'FAIL'} {name}: {detail}")
        fail += 0 if ok else 1

    try:
        for f in FILES:
            src = os.path.join(HERE, f)
            (shutil.copytree if os.path.isdir(src) else shutil.copy)(src, os.path.join(work, f))
        if args.asimov:      # noiseless data, covariance and prior mean at the truth; same T2.azr
            subprocess.run([sys.executable, os.path.join(HERE, "make_t2.py"), "--binary", binary,
                            "--asimov", work], check=True, capture_output=True, text=True)
        truthdir, fitdir = os.path.join(work, "truth"), os.path.join(work, "fit")
        for d in (truthdir, fitdir):
            os.makedirs(d)
            for f in FILES:
                src = os.path.join(work, f)
                (shutil.copytree if os.path.isdir(src) else shutil.copy)(src, os.path.join(d, f))
        run(binary, truthdir, "T2_true.azr", "calculate", args.threads)
        t_fit = run(binary, fitdir, "T2.azr", "fit", args.threads, ["--use-lm"] + args.engine_arg)
        truth_lv = read_params(os.path.join(truthdir, "output", "param.par"))
        nuis = json.load(open(os.path.join(HERE, "truth.json")))["nuisance"]
        sav = os.path.join(fitdir, "output", "param.sav")
        if not os.path.isfile(sav):
            print(open(os.path.join(fitdir, "T2.azr.log")).read()[-3000:])
            sys.exit("FAIL: the fit wrote no param.sav")
        start = read_params(os.path.join(fitdir, "output", "param.par"))
        fit = read_params(sav)
        chi2, npts = chi2_total(os.path.join(fitdir, "output", "chiSquared.out"))
        tol = 0.5 if args.asimov else 3.0
        print(f"T2 recovery ({'Asimov' if args.asimov else 'noisy data, seed in truth.json'}), "
              f"{args.threads} threads, fit {t_fit:.1f} s, chi2 {chi2:.2f} for N = {npts}")
        rows = []
        for name, src in FREE.items():
            t = truth_lv[name][0] if src == "level" else nuis[src]["true"]
            p, dp = fit[name]
            pull = (p - t) / dp if dp > 0 else float("inf")
            rows.append(dict(name=name, start=start[name][0], truth=t, fit=p, error=dp, pull=pull))
            check(name, abs(pull) <= tol, f"start {start[name][0]:.6g}, truth {t:.6g}, fit {p:.6g} +- {dp:.3g}, "
                                          f"pull {pull:+.2f} (|pull| <= {tol})")
        if not args.asimov:
            check("chi2/N", chi2 / npts <= 1.5, f"{chi2:.2f}/{npts} = {chi2 / npts:.3f} (<= 1.5)")
        if args.json:
            json.dump(dict(asimov=args.asimov, threads=args.threads, fit_seconds=t_fit, chi2=chi2, n=npts,
                           params=rows), open(args.json, "w"), indent=1)
    finally:
        if args.keep:
            shutil.copytree(work, args.keep, dirs_exist_ok=True)
        shutil.rmtree(work, ignore_errors=True)
    print("PASS" if fail == 0 else f"FAIL ({fail})")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
