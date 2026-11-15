#!/bin/csh
# rmfit campaign step on one node (written by `rmfit init`).  Usage:
#   qsub run_crc_rmfit                            (runs the command in RMFIT_CMD below)
#   qsub -v RMFIT_CMD="round --rounds 2" run_crc_rmfit
#   qsub -v RMFIT_CMD="polish --nfev 150 --tighten" run_crc_rmfit
#   qsub -hold_jid <jobid> -v RMFIT_CMD="..." run_crc_rmfit   (queue behind another step)
# Every command is resumable from the ledger; STATUS.md is rewritten at the end.
# One campaign command at a time per directory; the login node is for status/report/log.
# Engine pin (2026-11-07): AZURE2_ROOT and the _azure2 .so md5 below come from campaign.json["engine"];
# the job stops before Python if the pinned build is not the one on disk.

{mail_line}
#$ -m abe
{pe_line}
{queue_line}
#$ -N {jobname}
#$ -cwd
#$ -j y
#$ -o rmfit_sge.log

{setup_lines}
setenv OMP_NUM_THREADS 1
setenv PYTHONUNBUFFERED 1
if (! $?RMFIT_CMD) setenv RMFIT_CMD "status"
echo "`date` rmfit ${RMFIT_CMD} on `hostname`"
cd {campaign_dir}
{engine_guard}
setenv PYTHONPATH {pythonpath}:${AZURE2_ROOT}
python3 -c "import sqlite3, rmfit, rmfit.engine, rmfit.ledger, pyazr._azure2, pyazr; assert pyazr.__file__.startswith(rmfit.AZURE2_ROOT), pyazr.__file__; print('import guard ok', pyazr.__file__, rmfit.__file__)"
if ($status != 0) then
  echo "`date` IMPORT GUARD FAILED -- aborting"
  exit 1
endif
python3 -m rmfit.cli . ${RMFIT_CMD}
echo "`date` done (exit $status)"
