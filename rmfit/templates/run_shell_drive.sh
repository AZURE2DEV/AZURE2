#!/bin/bash
# rmfit drive with no scheduler (written by `rmfit <dir> drive` when the site scheduler is "shell"; rmfit/site.py).
# The script runs on the machine where it is started, with {slots} shard processes: start it where that many cores are free
# (a workstation, or a node already allocated to you).  Usage:
#   nohup ./run_crc_drive >/dev/null 2>&1 &                (resumes at the state in drive_state.json; the first run needs
#                                                           `rmfit <dir> drive --base <azr> [--sav ..] --dry-run` once)
#   DRIVE_MAX_MIN=240 nohup ./run_crc_drive >/dev/null 2>&1 &   (shorter legs)
# Legs run one after another in this script (no resubmission): a leg ends with DRIVE_STATUS CONTINUE (wall cap reached:
# the next leg starts), WAITING (a DECISIONS.md entry needs the evaluator: answer with `rmfit <dir> decide`, then start
# this script again), DONE, STOP (STOP_DRIVE file, budget, or a refused check) or ERROR.  `touch <dir>/STOP_DRIVE` stops
# after the current state.  DRIVE.md is rewritten after every state; output goes to {campaign_dir}/drive_shell.log.
# One run at a time per campaign: the lock directory below refuses a second one.
# Engine pin: AZURE2_ROOT and the _azure2 .so md5 below come from campaign.json["engine"] (or spec.engine before init).

DRIVE_MAX_MIN="${DRIVE_MAX_MIN:-{max_minutes}}"
cd {campaign_dir} || exit 1
{lock_lines}
echo "rmfit drive: output in {campaign_dir}/drive_shell.log"
exec >> drive_shell.log 2>&1
echo "$(date) rmfit drive on $(hostname) (shell, pid $$)"
{setup_lines}
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export RMFIT_SHELL_JOB=$$                # this script is the job here: the drive refuses engine states outside one
{engine_guard}
export PYTHONPATH={pythonpath}:${AZURE2_ROOT}
if ! python3 -c "import sqlite3, rmfit, rmfit.engine, rmfit.ledger, rmfit.drive, pyazr._azure2, pyazr; assert pyazr.__file__.startswith(rmfit.AZURE2_ROOT), pyazr.__file__; print('import guard ok', pyazr.__file__, rmfit.__file__)"; then
  echo "$(date) IMPORT GUARD FAILED -- aborting"
  exit 1
fi
leg=0
while true; do
  leg=$((leg + 1))
  python3 -m rmfit.cli . drive --max-minutes "${DRIVE_MAX_MIN}"
  rc=$?
  st=$(python3 -c "import json; print(json.load(open('drive_state.json')).get('status', 'ERROR'))" 2>/dev/null || echo ERROR)
  echo "$(date) leg $leg done (exit $rc, status $st)"
  [ "$st" = "CONTINUE" ] || break
  echo "$(date) next leg"
done
exit $rc
