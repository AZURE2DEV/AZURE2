#!/usr/bin/env python3
"""Plain-shell job templates (site setting "scheduler": "shell", rmfit/site.py): no scheduler, the script is the job.

  (a) scheduler() accepts "shell"; submit_cmd() is bash; submit_hint() names `nohup ./script &` (with RMFIT_CMD) here and
      qsub / sbatch on the other kinds; RMFIT_SHELL_JOB counts as a job for site, cli and drive (the login-node refusal);
  (b) the three writers write bash scripts with no scheduler directive and no unfilled field (`bash -n` accepts them),
      the site's setup lines, the job marker, a lock directory and a log file; the drive script loops legs (no
      resubmission);
  (c) the generated scripts run under bash with a stub python3 put on PATH by the site's setup lines:
      the campaign script logs to rmfit_shell.log and passes RMFIT_CMD with the marker set; a second run is refused while
      the lock is held; the lock is gone afterwards; a wrong engine md5 stops it before python3 -m rmfit.cli (lock
      released); the drive script runs legs while the state is CONTINUE and stops on DONE / WAITING; SIGTERM during a
      step releases the lock.
Engine-free, seconds.  Run:  python3 rmfit/tests/test_site_shell.py"""
import json, os, re, shutil, signal, subprocess, sys, tempfile, time

# The generated job scripts are POSIX shell for a Linux cluster or login node, and this test runs them.  On Windows
# "bash" from a native Python resolves to the WSL launcher in System32 ("Windows Subsystem for Linux has no installed
# distributions", the CI runner 2026-11-17), and signals, locks and md5sum differ: skipped there.
if os.name == "nt":
    print("skip: rmfit's job scripts are POSIX shell; not run on Windows")
    sys.exit(77)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from rmfit import site                                     # noqa: E402

FAILS = []
REAL_PY = sys.executable


def check(name, ok, info=""):
    print(("ok   " if ok else "FAIL ") + name + ("" if ok else f"  [{info}]"))
    if not ok:
        FAILS.append(name)


STUB = r"""#!/bin/bash
# stub python3: '-m rmfit.cli . <cmd>' is faked (recorded, drive state advanced); everything else goes to the real Python
if [ "$1" = "-m" ] && [ "$2" = "rmfit.cli" ]; then
  shift 3
  echo "STUB rmfit.cli $* marker=${RMFIT_SHELL_JOB:-unset}" >> "$STUB_LOG"
  if [ "$1" = "sleep" ]; then sleep "$2"; exit 0; fi
  if [ "$1" = "drive" ]; then
    n=$(( $(cat drive_legs 2>/dev/null || echo 0) + 1 )); echo $n > drive_legs
    st=CONTINUE; [ $n -ge ${STUB_LEGS:-3} ] && st=${STUB_FINAL:-DONE}
    echo "{\"status\": \"$st\"}" > drive_state.json
  fi
  exit 0
fi
exec REALPY "$@"
"""


def writers(tmp, cfg=None):
    from rmfit.campaign import Campaign
    from rmfit import drive
    from rmfit.spec import FitSpec
    from rmfit.bench import bench
    cdir = os.path.join(tmp, "X+p", "camp"); os.makedirs(cdir, exist_ok=True)
    c = object.__new__(Campaign); c.dir = cdir; c.cfg = dict({"base_azr": "x.azr", "nshards": 6}, **(cfg or {}))
    paths = {"campaign": c.write_job_script(overwrite=True),
             "drive": drive.write_job_script(cdir, FitSpec({"reaction": "X+p", "nshards": 4}), dict(c.cfg), overwrite=True)}
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    case = bench.load_case(os.path.join(here, "bench", "cases", "13N_tier0.json")); case["_base"] = tmp
    paths["bench"] = bench.write_job(case, next(iter(case["runs"])))
    return cdir, paths


def run(script, env_extra=None, timeout=60):
    env = {k: v for k, v in os.environ.items() if k not in ("JOB_ID", "SLURM_JOB_ID", "RMFIT_SHELL_JOB", "RMFIT_CMD")}
    env.update(env_extra or {})
    return subprocess.run(["bash", script], capture_output=True, text=True, env=env, timeout=timeout)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_site_shell_")
    keys = ("RMFIT_MAIL", "RMFIT_SITE", "RMFIT_SCHEDULER", "SLURM_JOB_ID", "JOB_ID", "RMFIT_SHELL_JOB", "RMFIT_CMD")
    saved = {k: os.environ.pop(k, None) for k in keys}
    try:
        sf = os.path.join(tmp, "site.json")
        stubbin = os.path.join(tmp, "stubbin"); os.makedirs(stubbin)
        open(os.path.join(stubbin, "python3"), "w").write(STUB.replace("REALPY", REAL_PY)); os.chmod(os.path.join(stubbin, "python3"), 0o755)
        json.dump({"scheduler": "shell", "mail": "you@example.org", "sge_queue": "q", "setup": [f'export PATH="{stubbin}:$PATH"']}, open(sf, "w"))
        os.environ["RMFIT_SITE"] = sf

        # (a)
        from rmfit import cli, drive
        check("(a) scheduler shell, submit_cmd bash", site.scheduler() == "shell" and site.submit_cmd() == "bash")
        check("(a) submit_hint (shell)", site.submit_hint("run_crc_drive") == "nohup ./run_crc_drive >/dev/null 2>&1 &"
              and site.submit_hint("run_crc_rmfit", "uq ...") == 'RMFIT_CMD="uq ..." nohup ./run_crc_rmfit >/dev/null 2>&1 &',
              site.submit_hint("run_crc_rmfit", "uq ..."))
        os.environ["RMFIT_SCHEDULER"] = "sge"
        h_sge = (site.submit_hint("run_crc_drive"), site.submit_hint("run_crc_rmfit", "uq ..."))
        os.environ["RMFIT_SCHEDULER"] = "slurm"
        h_slurm = site.submit_hint("run_crc_rmfit", "uq ...")
        del os.environ["RMFIT_SCHEDULER"]
        check("(a) submit_hint (sge, slurm)", h_sge == ("qsub run_crc_drive", 'qsub -v RMFIT_CMD="uq ..." run_crc_rmfit')
              and h_slurm == 'sbatch --export=ALL,RMFIT_CMD="uq ..." run_crc_rmfit', str((h_sge, h_slurm)))
        check("(a) no marker: not in a job", not site.in_batch_job() and not cli.in_sge_job_env() and not drive.in_sge_job())
        os.environ["RMFIT_SHELL_JOB"] = "4242"
        check("(a) RMFIT_SHELL_JOB: in a job (site, cli, drive)", site.in_batch_job() and cli.in_sge_job_env() and drive.in_sge_job())
        del os.environ["RMFIT_SHELL_JOB"]

        # (b)
        cdir, paths = writers(tmp)
        for name, p in paths.items():
            text = open(p).read(); L = text.splitlines()
            check(f"(b) {name}: bash, no scheduler directive, no mail", L[0] == "#!/bin/bash" and not any(l.startswith(("#$", "#SBATCH")) for l in L)
                  and "you@example.org" not in text)
            check(f"(b) {name}: no unfilled field", not re.search(r"(?<!\$)\{[a-z_]+\}", text), re.findall(r"(?<!\$)\{[a-z_]+\}", text))
            check(f"(b) {name}: setup line, marker, lock, log", f'export PATH="{stubbin}:$PATH"' in L and "export RMFIT_SHELL_JOB=$$" in text
                  and 'mkdir "$LOCK"' in text and "exec >> " in text)
            r = subprocess.run(["bash", "-n", p], capture_output=True, text=True)
            check(f"(b) {name}: bash -n", r.returncode == 0, r.stderr)
        d = open(paths["drive"]).read()
        check("(b) drive: legs loop in the script, no resubmission", "while true; do" in d and "qsub" not in d and "sbatch" not in d)
        check("(b) campaign: lock in the campaign directory, shard count named",
              f'LOCK="{os.path.join(cdir, site.SHELL_LOCK)}"' in open(paths["campaign"]).read() and "with 6 shard processes" in open(paths["campaign"]).read())

        # (c) behaviour under bash, with the stub python3
        stublog = os.path.join(tmp, "stub.log")
        env = {"STUB_LOG": stublog}
        lock = os.path.join(cdir, site.SHELL_LOCK)
        r = run(paths["campaign"], dict(env, RMFIT_CMD="polish --nfev 3"))
        log = open(os.path.join(cdir, "rmfit_shell.log")).read() if os.path.exists(os.path.join(cdir, "rmfit_shell.log")) else ""
        calls = open(stublog).read() if os.path.exists(stublog) else ""
        check("(c) campaign: exit 0, RMFIT_CMD passed with the marker set", r.returncode == 0 and re.search(r"STUB rmfit.cli polish --nfev 3 marker=\d+", calls) is not None,
              r.stdout + r.stderr + calls)
        check("(c) campaign: terminal names the log; log has the import guard and done", "output in" in r.stdout and "import guard ok" in log and "done (exit 0)" in log,
              r.stdout + log[-500:])
        check("(c) campaign: lock released", not os.path.exists(lock))
        os.makedirs(lock); open(os.path.join(lock, "pid"), "w").write("999999\n")
        r = run(paths["campaign"], dict(env, RMFIT_CMD="status"))
        check("(c) second run refused while the lock is held (pid named), lock kept", r.returncode == 1 and "999999" in r.stderr and os.path.isdir(lock), r.stderr)
        shutil.rmtree(lock)

        # wrong engine md5: stops before rmfit.cli, releases the lock
        so = os.path.join(tmp, "eng", "pyazr", "_azure2.cpython-312-x86_64-linux-gnu.so"); os.makedirs(os.path.dirname(so)); open(so, "wb").write(b"x")
        _, p2 = writers(os.path.join(tmp, "pin"), {"engine": {"azure2_root": os.path.join(tmp, "eng"), "so_md5": "f" * 32}})
        n0 = len(open(stublog).read().splitlines())
        r = run(p2["campaign"], dict(env, RMFIT_CMD="status"))
        pdir = os.path.dirname(p2["campaign"])
        plog = open(os.path.join(pdir, "rmfit_shell.log")).read()
        check("(c) wrong engine md5: exit 1 before rmfit.cli, lock released", r.returncode == 1 and "ENGINE PIN FAILED" in plog
              and len(open(stublog).read().splitlines()) == n0 and not os.path.exists(os.path.join(pdir, site.SHELL_LOCK)), plog[-300:])

        # drive: legs while CONTINUE
        for final, want in (("DONE", 3), ("WAITING", 2)):
            for f in ("drive_legs", "drive_state.json", "drive_shell.log"):
                if os.path.exists(os.path.join(cdir, f)):
                    os.remove(os.path.join(cdir, f))
            r = run(paths["drive"], dict(env, STUB_LEGS=str(want), STUB_FINAL=final, DRIVE_MAX_MIN="5"))
            dlog = open(os.path.join(cdir, "drive_shell.log")).read()
            legs = int(open(os.path.join(cdir, "drive_legs")).read())
            check(f"(c) drive: {want} legs, the last {final}, then the script ends", r.returncode == 0 and legs == want
                  and f"leg {want} done (exit 0, status {final})" in dlog and dlog.count("next leg") == want - 1
                  and not os.path.exists(lock), dlog[-400:])
        check("(c) drive: --max-minutes from DRIVE_MAX_MIN", "drive --max-minutes 5" in open(stublog).read())

        # SIGTERM during a step: the lock goes
        env2 = dict({k: v for k, v in os.environ.items() if k not in keys or k == "RMFIT_SITE"}, STUB_LOG=stublog, RMFIT_CMD="sleep 2")
        pr = subprocess.Popen(["bash", paths["campaign"]], env=env2, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        t0 = time.time()
        while not os.path.isdir(lock) and time.time() - t0 < 20:
            time.sleep(0.1)
        time.sleep(0.5)
        held = os.path.isdir(lock)
        pr.send_signal(signal.SIGTERM); rc = pr.wait(timeout=30)
        check("(c) SIGTERM: script ends (143) and releases the lock", held and rc == 143 and not os.path.exists(lock), f"held {held} rc {rc}")
    finally:
        for k, v in saved.items():
            os.environ.pop(k, None)
            if v is not None:
                os.environ[k] = v
        shutil.rmtree(tmp, ignore_errors=True)
    if FAILS:
        print(f"{len(FAILS)} FAILED"); sys.exit(1)
    print("all site shell checks passed")


if __name__ == "__main__":
    main()
