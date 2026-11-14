# Running rmfit jobs on a shared cluster

rmfit steps run as batch jobs (Sun Grid Engine `qsub` in the examples) on shared nodes and a shared group file system.
This reference lists what a job script must contain and the process, quota and authentication traps of such a
cluster. Nothing here depends on one machine; replace the scheduler commands with your site's.

## 1. The job script

- **Prove the environment on the compute node.** A cluster's default Python can change, and a module can fail to load
  silently, leaving the system Python. The first command of every job imports what the run needs and exits otherwise:

  ```csh
  python3 -c "import sqlite3, sys; sys.path.insert(0, '$AZURE2_ROOT'); import pyazr; pyazr.describe()" || exit 1
  ```

  A login-node dry run proves nothing: login shells inherit another interpreter and other libraries (a compute node's
  `sqlite3` once lacked `sqlite3_deserialize`). Pin the interpreter and, if needed, a newer SQLite library on
  `LD_LIBRARY_PATH` in the script.
- **Pin the build.** rmfit reads `AZURE2_ROOT` and, when unset, falls back to a default checkout which it puts first on
  `sys.path`. With only `PYTHONPATH` set, the parent imports one pyazr and the shard workers another -- the shared copy
  that installs replace under running jobs. Set both, and assert before any work:

  ```csh
  setenv AZURE2_ROOT <frozen build directory>
  setenv PYTHONPATH  ${AZURE2_ROOT}:${RMFIT_HOME}
  ```
  ```python
  assert pyazr.__file__.startswith(rmfit.AZURE2_ROOT)
  ```

  Shared binaries change under you: when a result stops reproducing, check the modification time of the engine
  module (`pyazr.describe()["engine"]["module"]`, with its git version) and of the CLI binary first.
- **`if __name__ == "__main__": main()`** in every script that opens a `SingleSession`, `Sharded` or `Campaign`. The
  engine starts workers with multiprocessing's spawn method, which re-imports the script in each child; module-level
  code runs again there and the parent dies ("An attempt has been made to start a new process before the current
  process has finished its bootstrapping phase", or EOFError). Heredoc scripts fail for the same reason.
- **csh quoting.** Never put JSON in a shell variable or backquotes (csh brace-expands `{...}`: "Missing '}'"); pass
  literal single-quoted JSON or `@file.json`. A value starting with a minus is read by argparse as an option: use
  `--opt=value`. csh `eval` of a step string drops everything after `#` (as in `9/2+#2`) and splits quoted values:
  write steps as explicit quoted lines and dry-run them into Python's `argv` before submitting.
- **Run the job's exact command line once on the login node for a few seconds** before submitting it.
- **Queue behind a running job** with `qsub -hold_jid <id>`; a held step starts at once if the job it waits for
  fails, so it checks its inputs (the previous export, an install marker) and stops if they are missing.
- **Confirm the submission.** Check the scheduler's reply for the job id ("Your job N") before writing it anywhere;
  a failed `qsub` can print only "Exiting.", and a script must never log "submitted" unconditionally.
- **Mail on start, end and abort** (`-m abe`).

## 2. Authentication

On clusters where submission needs a Kerberos/AFS token:
- A detached (`nohup`) submitter has no token ("job does not provide an AFS token"): submit from a live session's shell.
- A restarted session can inherit a `KRB5CCNAME` that points at a vanished cache (`klist`: "No credentials cache
  found"). List the caches newest first, take the one whose `klist` shows a ticket-granting ticket valid now, export
  `KRB5CCNAME=FILE:<path>` for the submission, or ask the user to `kinit`. A reboot removes the caches; never hard-code a
  cache path across days.
- Long login-node computations launched from a tool die with the session: use `nohup`, a done-file, and watch the
  done-file.

## 3. Processes and quotas

- **`pkill -f` / `pgrep -f` patterns are regular expressions, and reaction names contain `+`.** `pkill -f "13C+a/..."`
  matches nothing (`C+` = one or more C) and looks exactly like "no such process". A pattern that appears in the
  calling shell's own command line kills that shell. Find the pid with a fixed-string match
  (`ps -eo pid,cmd | grep -F "<name>"`) or escape the plus, kill by pid, confirm with `kill -0 <pid>` failing, and only
  then remove its log.
- **Launch background programs with `< /dev/null` and a size cap** (`2>&1 | head -c 200000 > log`). A program that
  reads from an open stdin can wait forever (AZURE2's commands read nothing from it); a log deleted under a live process keeps growing as a hidden
  `.nfs...` file (one reached 748 GB and filled a 5 TB group quota, killing every job writing there). `ls -la` shows
  the `.nfs` file's size; `ls -l /proc/<pid>/fd` finds the holder.
- **When a batch of jobs dies at once, check the group quota first** (the group's own quota, not the size of the whole
  file system that `df` may show). Truncated logs ending after the module lines are the sign.

## 4. Watching jobs

- `campaign.log` and the ledger are the truth; the scheduler's stdout file lags. Do not kill a job for silence in the
  stdout file.
- Poll, do not `tail -F`, on NFS.
- Check the queue for other users' and sessions' jobs on the same campaign before every submission.
- After any rmfit install, read the first new job's log past its evaluator start-up; running drivers keep the modules
  they loaded, so a bad install shows up only in new jobs.
