#!/usr/bin/env python3
"""Slurm job templates (site setting "scheduler": "slurm", rmfit/site.py).

  (a) scheduler(): default sge, RMFIT_SCHEDULER and the site file select it, anything else raises; submit_cmd();
      in_batch_job() / cli.in_sge_job_env() / drive.in_sge_job() see SLURM_JOB_ID;
  (b) slurm_minutes / slurm_time_str; drive_leg_minutes caps the leg at slurm_time - 30 min;
  (c) the three writers (Campaign.write_job_script, drive.write_job_script, bench.write_job) write bash scripts with the
      #SBATCH header from the site file (partition, account, time, mail, extra, cpus = shards), the setup lines, no SGE
      directive and no unfilled field; `bash -n` accepts them; the drive script resubmits with sbatch;
  (d) enginepin.guard_lines_sh run by bash: passes on the right .so md5, exits 1 on a wrong one.
Engine-free, seconds; no Slurm needed.  The SGE scripts are covered by test_site_mail.py.
Run:  python3 rmfit/tests/test_site_slurm.py"""
import hashlib, json, os, re, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from rmfit import site, enginepin                         # noqa: E402

FAILS = []


def check(name, ok, info=""):
    print(("ok   " if ok else "FAIL ") + name + ("" if ok else f"  [{info}]"))
    if not ok:
        FAILS.append(name)


def writers(tmp):
    from rmfit.campaign import Campaign
    from rmfit import drive
    from rmfit.spec import FitSpec
    from rmfit.bench import bench
    cdir = os.path.join(tmp, "X+p", "camp"); os.makedirs(cdir, exist_ok=True)
    c = object.__new__(Campaign); c.dir = cdir; c.cfg = {"base_azr": "x.azr", "nshards": 6}
    out = {"campaign": open(c.write_job_script(overwrite=True)).read(),
           "drive": open(drive.write_job_script(cdir, FitSpec({"reaction": "X+p", "nshards": 4}), {}, overwrite=True)).read()}
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    case = bench.load_case(os.path.join(here, "bench", "cases", "13N_tier0.json")); case["_base"] = tmp
    out["bench"] = open(bench.write_job(case, next(iter(case["runs"])))).read()
    return out


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_site_slurm_")
    saved = {k: os.environ.pop(k, None) for k in ("RMFIT_MAIL", "RMFIT_SITE", "RMFIT_SCHEDULER", "SLURM_JOB_ID", "JOB_ID")}
    try:
        sf = os.path.join(tmp, "site.json")
        os.environ["RMFIT_SITE"] = os.path.join(tmp, "absent.json")
        check("(a) default scheduler sge, qsub", site.scheduler() == "sge" and site.submit_cmd() == "qsub")
        os.environ["RMFIT_SCHEDULER"] = "Slurm"
        check("(a) RMFIT_SCHEDULER selects slurm, sbatch", site.scheduler() == "slurm" and site.submit_cmd() == "sbatch")
        os.environ["RMFIT_SCHEDULER"] = "pbs"
        try:
            site.scheduler(); raised = False
        except ValueError:
            raised = True
        check("(a) an unknown scheduler raises", raised)
        del os.environ["RMFIT_SCHEDULER"]
        from rmfit import cli, drive
        check("(a) no job env: not in a job", not site.in_batch_job() and not cli.in_sge_job_env() and not drive.in_sge_job())
        os.environ["SLURM_JOB_ID"] = "123"
        check("(a) SLURM_JOB_ID: in a job (site, cli, drive)", site.in_batch_job() and cli.in_sge_job_env() and drive.in_sge_job())
        del os.environ["SLURM_JOB_ID"]

        check("(b) slurm_minutes", [site.slurm_minutes(s) for s in ("30", "12:00:00", "2-00:00:00", "1-12", "x")] == [30, 720, 2880, 2160, None])
        check("(b) slurm_time_str", site.slurm_time_str(630) == "0-10:30:00")
        json.dump({"scheduler": "slurm", "slurm_time": "8:00:00"}, open(sf, "w")); os.environ["RMFIT_SITE"] = sf
        check("(b) drive leg capped by slurm_time", site.drive_leg_minutes(600) == (450, 480) and site.drive_leg_minutes(120) == (120, 150))

        json.dump({"scheduler": "slurm", "mail": "you@example.org", "slurm_partition": "compute", "slurm_account": "proj1",
                   "slurm_time": "2-00:00:00", "slurm_extra": ["--mem=0"], "setup": ["module load python/3.12"]}, open(sf, "w"))
        for name, text in writers(tmp).items():
            L = text.splitlines()
            sb = [l for l in L if l.startswith("#SBATCH")]
            check(f"(c) {name}: bash, #SBATCH header from the site file",
                  L[0] == "#!/bin/bash" and {"#SBATCH --partition=compute", "#SBATCH --account=proj1", "#SBATCH --mail-user=you@example.org",
                                             "#SBATCH --mem=0", "#SBATCH --nodes=1"} <= set(sb), sb)
            first_cmd = next(i for i, l in enumerate(L[1:], 1) if l.strip() and not l.startswith("#"))
            check(f"(c) {name}: every #SBATCH line before the first command", all(L.index(s) < first_cmd for s in sb))
            check(f"(c) {name}: setup line, no SGE directive, no unfilled field",
                  "module load python/3.12" in L and not any(l.startswith("#$") for l in L) and not re.search(r"(?<!\$)\{[a-z_]+\}", text),
                  re.findall(r"(?<!\$)\{[a-z_]+\}", text))
            p = os.path.join(tmp, f"{name}.sh"); open(p, "w").write(text)
            r = subprocess.run(["bash", "-n", p], capture_output=True, text=True)
            check(f"(c) {name}: bash -n", r.returncode == 0, r.stderr)
        w = writers(tmp)
        check("(c) campaign: cpus = the campaign's shards", "#SBATCH --cpus-per-task=6" in w["campaign"].splitlines())
        check("(c) campaign: mail on begin, end, fail", "#SBATCH --mail-type=BEGIN,END,FAIL" in w["campaign"].splitlines())
        check("(c) drive: --time = leg + 30 min, sbatch resubmission",
              "#SBATCH --time=0-10:30:00" in w["drive"].splitlines() and 'DRIVE_MAX_MIN="${DRIVE_MAX_MIN:-600}"' in w["drive"]
              and re.search(r"if ! sbatch \S+/run_crc_drive; then", w["drive"]) is not None)

        so = os.path.join(tmp, "eng", "pyazr", "_azure2.cpython-312-x86_64-linux-gnu.so"); os.makedirs(os.path.dirname(so))
        open(so, "wb").write(b"fake engine")
        md5 = hashlib.md5(b"fake engine").hexdigest()
        for label, rec_md5, want in (("right md5 passes", md5, 0), ("wrong md5 exits 1", "f" * 32, 1)):
            script = enginepin.guard_lines_sh({"azure2_root": os.path.join(tmp, "eng"), "so_md5": rec_md5}) + 'echo "guard passed $AZURE2_ROOT"\n'
            r = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
            check(f"(d) sh guard: {label}", r.returncode == want and (("guard passed" in r.stdout) == (want == 0)), r.stdout + r.stderr)
        r = subprocess.run(["bash", "-c", enginepin.guard_lines_sh(None) + 'echo "$AZURE2_ROOT"'], capture_output=True, text=True,
                           env=dict(os.environ, HOME="/h", AZURE2_ROOT=""))
        check("(d) sh guard without a pin: $HOME/AZURE2", r.stdout.strip() == "/h/AZURE2", r.stdout)
    finally:
        for k, v in saved.items():
            os.environ.pop(k, None)
            if v is not None:
                os.environ[k] = v
        shutil.rmtree(tmp, ignore_errors=True)
    if FAILS:
        print(f"{len(FAILS)} FAILED"); sys.exit(1)
    print("all site slurm checks passed")


if __name__ == "__main__":
    main()
