#!/bin/bash
# rmfit drive: one leg of the closed loop on one node, Slurm (written by `rmfit <dir> drive` when the site scheduler is
# slurm; rmfit/site.py).  Usage:
#   sbatch run_crc_drive                                  (resumes at the state in drive_state.json; the first leg needs
#                                                          `rmfit <dir> drive --base <azr> [--sav ..] --dry-run` once)
#   sbatch --export=ALL,DRIVE_MAX_MIN=240 run_crc_drive   (shorter leg; keep it below the job's --time minus 30 min)
#   sbatch --dependency=afterany:<jobid> run_crc_drive    (queue behind another step)
# The leg ends with DRIVE_STATUS CONTINUE (wall cap reached: this script resubmits itself), WAITING (a DECISIONS.md entry
# needs the evaluator: answer with `rmfit <dir> decide`, then sbatch again), DONE, STOP (STOP_DRIVE file, budget, or a
# refused check) or ERROR.  `touch <dir>/STOP_DRIVE` stops after the current state.  DRIVE.md is rewritten after every state.
# Engine pin: AZURE2_ROOT and the _azure2 .so md5 below come from campaign.json["engine"] (or spec.engine before init).
{leg_note}
{slurm_header}

echo "$(date) rmfit drive leg on $(hostname) (job ${SLURM_JOB_ID:-none})"
{setup_lines}
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
DRIVE_MAX_MIN="${DRIVE_MAX_MIN:-{max_minutes}}"
cd {campaign_dir} || exit 1
{engine_guard}
export PYTHONPATH={pythonpath}:${AZURE2_ROOT}
if ! python3 -c "import sqlite3, rmfit, rmfit.engine, rmfit.ledger, rmfit.drive, pyazr._azure2, pyazr; assert pyazr.__file__.startswith(rmfit.AZURE2_ROOT), pyazr.__file__; print('import guard ok', pyazr.__file__, rmfit.__file__)"; then
  echo "$(date) IMPORT GUARD FAILED -- aborting"
  exit 1
fi
python3 -m rmfit.cli . drive --max-minutes "${DRIVE_MAX_MIN}"
rc=$?
st=$(python3 -c "import json; print(json.load(open('drive_state.json')).get('status', 'ERROR'))" 2>/dev/null || echo ERROR)
echo "$(date) leg done (exit $rc, status $st)"
if [ "$st" = "CONTINUE" ]; then
  echo "$(date) resubmitting for the next leg"
  if ! sbatch {campaign_dir}/run_crc_drive; then
    echo "$(date) RESUBMIT FAILED (sbatch not allowed from compute nodes here?) -- sbatch run_crc_drive from a login node"
  fi
fi
exit $rc
