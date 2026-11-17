#!/bin/bash
# rmfit campaign step on one node, Slurm (written by `rmfit init` when the site scheduler is slurm; rmfit/site.py).  Usage:
#   sbatch run_crc_rmfit                                      (runs the command in RMFIT_CMD below; default: status)
#   sbatch --export=ALL,RMFIT_CMD="round --rounds 2" run_crc_rmfit
#   sbatch --export=ALL,RMFIT_CMD="polish --nfev 150 --tighten" run_crc_rmfit
#   sbatch --dependency=afterany:<jobid> --export=ALL,RMFIT_CMD="..." run_crc_rmfit   (queue behind another step)
# --export splits on commas: a command containing a comma goes in a file (rmfit's JSON options accept @file).
# Every command is resumable from the ledger; STATUS.md is rewritten at the end.
# One campaign command at a time per directory; the login node is for status/report/log.
# Engine pin: AZURE2_ROOT and the _azure2 .so md5 below come from campaign.json["engine"]; the job stops before Python if
# the pinned build is not the one on disk.

{slurm_header}

{setup_lines}
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
RMFIT_CMD="${RMFIT_CMD:-status}"
echo "$(date) rmfit ${RMFIT_CMD} on $(hostname) (job ${SLURM_JOB_ID:-none})"
cd {campaign_dir} || exit 1
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
