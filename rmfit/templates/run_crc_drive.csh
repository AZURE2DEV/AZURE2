#!/bin/csh
# rmfit drive: one leg of the closed loop on one node (written by `rmfit <dir> drive`).  Usage:
#   qsub run_crc_drive                       (resumes at the state in drive_state.json; the first leg needs
#                                             `rmfit <dir> drive --base <azr> [--sav ..] --dry-run` once to record the base)
#   qsub -v DRIVE_MAX_MIN=240 run_crc_drive  (shorter leg)
#   qsub -hold_jid <jobid> run_crc_drive     (queue behind another step)
# The leg ends with DRIVE_STATUS CONTINUE (wall cap reached: this script resubmits itself), WAITING (a DECISIONS.md entry
# needs the evaluator: answer with `rmfit <dir> decide`, then qsub again), DONE, STOP (STOP_DRIVE file, budget, or a refused
# check) or ERROR.  `touch <dir>/STOP_DRIVE` stops after the current state.  DRIVE.md is rewritten after every state.
# Engine pin (2026-11-07): AZURE2_ROOT and the _azure2 .so md5 below come from campaign.json["engine"] (or spec.engine before
# init); the job stops before Python if the pinned build is not the one on disk.

{mail_line}
#$ -m a
{pe_line}
{queue_line}
#$ -N {jobname}
#$ -cwd
#$ -j y
#$ -o drive_sge.log

echo "`date` rmfit drive leg on `hostname` (job $JOB_ID)"
{setup_lines}
setenv OMP_NUM_THREADS 1
setenv PYTHONUNBUFFERED 1
if (! $?DRIVE_MAX_MIN) setenv DRIVE_MAX_MIN {max_minutes}
cd {campaign_dir}
{engine_guard}
setenv PYTHONPATH {pythonpath}:${AZURE2_ROOT}
python3 -c "import sqlite3, rmfit, rmfit.engine, rmfit.ledger, rmfit.drive, pyazr._azure2, pyazr; assert pyazr.__file__.startswith(rmfit.AZURE2_ROOT), pyazr.__file__; print('import guard ok', pyazr.__file__, rmfit.__file__)"
if ($status != 0) then
  echo "`date` IMPORT GUARD FAILED -- aborting"
  exit 1
endif
python3 -m rmfit.cli . drive --max-minutes ${DRIVE_MAX_MIN}
set rc = $status
set st = `python3 -c "import json; print(json.load(open('drive_state.json')).get('status', 'ERROR'))"`
echo "`date` leg done (exit $rc, status $st)"
if ("$st" == "CONTINUE") then
  echo "`date` resubmitting for the next leg"
  qsub {campaign_dir}/run_crc_drive
  if ($status != 0) echo "`date` RESUBMIT FAILED (no AFS/Kerberos token on the node?) -- qsub run_crc_drive from a login node"
endif
exit $rc
