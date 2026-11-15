"""Site settings: values that belong to a cluster site or a user, not to a campaign.

Precedence: environment variable > site file > default.  The site file is ``$RMFIT_SITE`` or ``~/.rmfit/site.json``,
for example::

    {"mail": "you@example.org"}

Keys:
  mail             job mail address (env ``RMFIT_MAIL``).  SGE: the ``-M`` line; Slurm: ``--mail-user``.  Unset: no address
                   line, so the scheduler mails the submitting user (its own default).
  scheduler        "sge" (default), "slurm" or "shell" (env ``RMFIT_SCHEDULER``): which job templates rmfit writes and
                   how its messages say to start them.  SGE templates: ``templates/run_crc_*.csh`` (csh); Slurm templates:
                   ``templates/run_slurm_*.sh`` (bash); "shell": ``templates/run_shell_*.sh``, bash scripts with no
                   scheduler that run where they are started (a workstation, or a node already allocated), log to a file
                   in the campaign directory, hold a lock directory so that two runs cannot share a campaign, and run the
                   drive's legs in a loop instead of resubmitting.  No mail on "shell".
  sge_queue        SGE queue for the ``-q`` line (env ``RMFIT_SGE_QUEUE``), e.g. "long".  Unset: no ``-q`` line, so SGE
                   picks the queue from its defaults and the job's requests.
  sge_pe           SGE parallel environment for the ``-pe <pe> <slots>`` line (env ``RMFIT_SGE_PE``; default "smp", a
                   common name for a one-node shared-memory PE: ``qconf -spl`` lists a site's PEs).  Slurm needs neither:
                   ``--cpus-per-task`` and ``slurm_partition`` play these parts.
  slurm_partition  ``--partition`` (omitted when unset: the cluster's default partition)
  slurm_account    ``--account``
  slurm_time       ``--time`` for campaign steps and bench jobs (Slurm formats MM, HH:MM:SS, D-HH:MM:SS; omitted when
                   unset).  A drive leg asks for DRIVE_MAX_MIN + 30 minutes, capped by slurm_time with DRIVE_MAX_MIN
                   shortened to fit.
  slurm_extra      list of further ``#SBATCH`` options, e.g. ["--mem=0", "--qos=long"]
  setup            list of shell lines run before Python in every job script rmfit writes (modules, a Python with numpy
                   and scipy on PATH, ...), in that script's shell: csh for the SGE templates, bash for the Slurm and
                   shell ones.  Unset: a comment, and the job uses the environment it starts with.

Example for an SGE site::

    {"mail": "you@example.org", "sge_queue": "long", "sge_pe": "smp", "setup": ["module load python/3.12 gsl"]}

Example for a workstation::

    {"scheduler": "shell", "setup": ["source ~/venvs/rmfit/bin/activate"]}

Example for a Slurm site::

    {"scheduler": "slurm", "mail": "you@example.org", "slurm_partition": "compute", "slurm_time": "2-00:00:00",
     "slurm_extra": ["--mem=0"], "setup": ["module load python/3.12 gsl"]}
"""
import json
import os

SITE_FILE_DEFAULT = os.path.join("~", ".rmfit", "site.json")


def site_file() -> str:
    return os.path.expanduser(os.environ.get("RMFIT_SITE") or SITE_FILE_DEFAULT)


def load() -> dict:
    """The site file's settings ({} when it is absent); a malformed file raises with its path."""
    path = site_file()
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as f:
            data = json.load(f)
    except ValueError as e:
        raise ValueError(f"rmfit site file {path}: {e}") from None
    if not isinstance(data, dict):
        raise ValueError(f"rmfit site file {path}: expected a JSON object")
    return data


def job_mail() -> str:
    """The job mail address: ``RMFIT_MAIL``, else the site file's ``mail``, else '' (none)."""
    return (os.environ.get("RMFIT_MAIL") or str(load().get("mail") or "")).strip()


def mail_line() -> str:
    """The SGE ``-M`` line for a job script, or a comment saying how to set one (SGE templates)."""
    addr = job_mail()
    if addr:
        return f"#$ -M {addr}"
    return "# no -M line: set RMFIT_MAIL or \"mail\" in ~/.rmfit/site.json (rmfit/site.py); SGE mails the submitting user"


def sge_queue() -> str:
    """The SGE queue: ``RMFIT_SGE_QUEUE``, else the site file's ``sge_queue``, else '' (no -q line)."""
    return (os.environ.get("RMFIT_SGE_QUEUE") or str(load().get("sge_queue") or "")).strip()


def queue_line() -> str:
    """The SGE ``-q`` line for a job script, or a comment saying how to set one (SGE templates)."""
    q = sge_queue()
    if q:
        return f"#$ -q {q}"
    return "# no -q line: set RMFIT_SGE_QUEUE or \"sge_queue\" in ~/.rmfit/site.json (rmfit/site.py); SGE picks the queue"


SGE_PE_DEFAULT = "smp"


def sge_pe() -> str:
    """The SGE parallel environment: ``RMFIT_SGE_PE``, else the site file's ``sge_pe``, else 'smp'."""
    return (os.environ.get("RMFIT_SGE_PE") or str(load().get("sge_pe") or SGE_PE_DEFAULT)).strip()


def pe_line(slots: int) -> str:
    """The SGE ``-pe <pe> <slots>`` line (SGE templates)."""
    return f"#$ -pe {sge_pe()} {int(slots)}"


SCHEDULERS = ("sge", "slurm", "shell")


def scheduler() -> str:
    """'sge' or 'slurm': ``RMFIT_SCHEDULER``, else the site file's ``scheduler``, else 'sge'."""
    s = (os.environ.get("RMFIT_SCHEDULER") or str(load().get("scheduler") or "sge")).strip().lower()
    if s not in SCHEDULERS:
        raise ValueError(f"rmfit scheduler {s!r}: expected one of {SCHEDULERS} (RMFIT_SCHEDULER or 'scheduler' in {site_file()})")
    return s


def submit_cmd() -> str:
    """The command that starts a job script (``bench run --submit``): qsub, sbatch, or bash for "shell"."""
    return {"sge": "qsub", "slurm": "sbatch", "shell": "bash"}[scheduler()]


def submit_hint(script: str, cmd: str = "") -> str:
    """How to start ``script`` (with ``RMFIT_CMD`` = cmd if given), for messages."""
    s = scheduler()
    if s == "sge":
        return f'qsub -v RMFIT_CMD="{cmd}" {script}' if cmd else f"qsub {script}"
    if s == "slurm":
        return f'sbatch --export=ALL,RMFIT_CMD="{cmd}" {script}' if cmd else f"sbatch {script}"
    return (f'RMFIT_CMD="{cmd}" ' if cmd else "") + f"nohup ./{script} >/dev/null 2>&1 &"


def in_batch_job() -> bool:
    """True inside a job: SGE ``JOB_ID``, Slurm ``SLURM_JOB_ID``, or ``RMFIT_SHELL_JOB`` (set by the "shell" scripts,
    which are the job on a site without a scheduler), whatever the site setting."""
    return bool(os.environ.get("JOB_ID") or os.environ.get("SLURM_JOB_ID") or os.environ.get("RMFIT_SHELL_JOB"))


def slurm_minutes(spec: str):
    """Minutes of a Slurm time string (MM, MM:SS, HH:MM:SS, D-HH, D-HH:MM, D-HH:MM:SS); None if it cannot be read."""
    s = str(spec).strip()
    try:
        days = 0
        if "-" in s:
            d, s = s.split("-", 1)
            days = int(d)
            parts = [int(x) for x in s.split(":")]
            parts += [0] * (3 - len(parts))                  # D-HH, D-HH:MM
            h, m, sec = parts
        else:
            parts = [int(x) for x in s.split(":")]
            if len(parts) == 1:
                h, m, sec = 0, parts[0], 0
            elif len(parts) == 2:
                h, m, sec = 0, parts[0], parts[1]
            else:
                h, m, sec = parts
        return days * 1440 + h * 60 + m + (1 if sec else 0)
    except ValueError:
        return None


def slurm_time_str(minutes: int) -> str:
    d, rem = divmod(int(minutes), 1440)
    h, m = divmod(rem, 60)
    return f"{d}-{h:02d}:{m:02d}:00"


def slurm_header(jobname: str, slots: int, log: str, mail_events: str = "abe", time_minutes=None) -> str:
    """``#SBATCH`` lines for a one-node job.  mail_events in SGE letters (a = FAIL, b = BEGIN, e = END);
    time_minutes overrides slurm_time (the drive's computed leg length)."""
    cfg = load()
    lines = [f"#SBATCH --job-name={jobname}", "#SBATCH --nodes=1", "#SBATCH --ntasks=1",
             f"#SBATCH --cpus-per-task={int(slots)}", f"#SBATCH --output={log}", "#SBATCH --open-mode=append"]
    for key, opt in (("slurm_partition", "--partition"), ("slurm_account", "--account")):
        if cfg.get(key):
            lines.append(f"#SBATCH {opt}={cfg[key]}")
    if time_minutes is not None:
        lines.append(f"#SBATCH --time={slurm_time_str(time_minutes)}")
    elif cfg.get("slurm_time"):
        lines.append(f"#SBATCH --time={cfg['slurm_time']}")
    addr = job_mail()
    if addr:
        types = [t for c, t in (("b", "BEGIN"), ("e", "END"), ("a", "FAIL")) if c in mail_events]
        lines.append(f"#SBATCH --mail-user={addr}")
        if types:
            lines.append(f"#SBATCH --mail-type={','.join(types)}")
    else:
        lines.append("# no --mail-user: set RMFIT_MAIL or \"mail\" in ~/.rmfit/site.json (rmfit/site.py)")
    for opt in cfg.get("slurm_extra") or []:
        lines.append(f"#SBATCH {opt}")
    return "\n".join(lines)


def drive_leg_minutes(max_minutes: int):
    """(DRIVE_MAX_MIN, Slurm --time minutes) for a drive leg: the leg plus 30 minutes for the last state to finish and
    the job to resubmit, shortened to fit slurm_time when that is smaller."""
    cap = slurm_minutes(load().get("slurm_time")) if load().get("slurm_time") else None
    leg = int(max_minutes)
    if cap is not None and leg + 30 > cap:
        leg = max(cap - 30, 10)
    return leg, (leg + 30 if cap is None else min(leg + 30, cap))


SHELL_LOCK = ".rmfit_shell.lock"


def shell_lock_lines(lockdir: str) -> str:
    """bash lines for a "shell" script: take the lock directory (mkdir is atomic) or stop, remove it on exit."""
    return (f'LOCK="{lockdir}"\n'
            'if ! mkdir "$LOCK" 2>/dev/null; then\n'
            '  echo "$(date) another rmfit shell run holds $LOCK (pid $(cat "$LOCK/pid" 2>/dev/null || echo ?)); '
            'remove the directory if that process is gone" >&2\n'
            '  exit 1\n'
            'fi\n'
            'echo $$ > "$LOCK/pid"\n'
            'trap \'rm -rf "$LOCK"\' EXIT\n'
            "trap 'exit 143' TERM INT HUP")


def setup_lines() -> str:
    """The site's setup lines for a job script (or a comment saying where to put them); csh for SGE, bash otherwise."""
    lines = load().get("setup") or []
    if not lines:
        return ("# site setup (modules, a Python with numpy/scipy on PATH): \"setup\" in ~/.rmfit/site.json (rmfit/site.py),\n"
                "# in this script's shell (csh for SGE, bash for Slurm and shell)")
    return "\n".join(str(x) for x in lines)
