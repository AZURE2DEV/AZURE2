#!/bin/bash
# rmfit campaign step with no scheduler (written by `rmfit init` when the site scheduler is "shell"; rmfit/site.py).
# The script runs on the machine where it is started, with {nshards} shard processes: start it where that many cores are
# free (a workstation, or a node already allocated to you).  Usage:
#   nohup ./run_crc_rmfit >/dev/null 2>&1 &                               (runs the command in RMFIT_CMD; default: status)
#   RMFIT_CMD="round --rounds 2" nohup ./run_crc_rmfit >/dev/null 2>&1 &
#   RMFIT_CMD="polish --nfev 150 --tighten" ./run_crc_rmfit                (in the foreground)
# Output goes to {campaign_dir}/rmfit_shell.log.  Every command is resumable from the ledger; STATUS.md is rewritten at the
# end.  One run at a time per campaign: the lock directory below refuses a second one.  To stop: touch STOP (sign-round
# loops), or kill the python3 process; the lock goes when the script exits.
# Engine pin: AZURE2_ROOT and the _azure2 .so md5 below come from campaign.json["engine"]; the run stops before Python if
# the pinned build is not the one on disk.

RMFIT_CMD="${RMFIT_CMD:-status}"
cd {campaign_dir} || exit 1
{lock_lines}
echo "rmfit ${RMFIT_CMD}: output in {campaign_dir}/rmfit_shell.log"
exec >> rmfit_shell.log 2>&1
echo "$(date) rmfit ${RMFIT_CMD} on $(hostname) (shell, pid $$)"
{setup_lines}
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export RMFIT_SHELL_JOB=$$                # this script is the job here: rmfit's engine commands refuse to run outside one
{engine_guard}
export PYTHONPATH={pythonpath}:${AZURE2_ROOT}
if ! python3 -c "import sqlite3, rmfit, rmfit.engine, rmfit.ledger, pyazr._azure2, pyazr; assert pyazr.__file__.startswith(rmfit.AZURE2_ROOT), pyazr.__file__; print('import guard ok', pyazr.__file__, rmfit.__file__)"; then
  echo "$(date) IMPORT GUARD FAILED -- aborting"
  exit 1
fi
set -f                                   # RMFIT_CMD is split into words like the csh template, but never globbed
python3 -m rmfit.cli . ${RMFIT_CMD}
rc=$?
echo "$(date) done (exit $rc)"
exit $rc
