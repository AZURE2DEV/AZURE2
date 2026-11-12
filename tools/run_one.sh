#!/bin/bash
#
# One job of tools/runall.sh: run every .azr of one folder with one engine in one basis, menu 1 (calculate with data),
# one OpenMP thread, in a private copy of the folder's top-level test folder (so nested folders such as X and
# X/models never share a work tree).  Each run is cut off after 900 s and its log after 2 MB.
#
#   tools/run_one.sh SRC_TESTS CMP TAG ENGINE MODE DIR
#
#   SRC_TESTS  the tests/ directory of the source tree
#   CMP        the results root ($WORK/cmp in runall.sh)
#   TAG        the run's name
#   ENGINE     the AZURE2 binary (absolute path)
#   MODE       the basis flag: "" (Brune), --no-brune or --use-park
#   DIR        the folder, relative to SRC_TESTS
#
# Output: CMP/<basis>/TAG/.res-<DIR with "/" as "__">/<file>.azr/{output,checks,run.log} for each .azr of DIR, where
# <basis> is brune, no-brune or use-park.  The work copy CMP/<basis>/TAG/job-<...> is removed at the end.
# Environment: AZURE2_COULOMB_TABLE, the Coulomb table file (required).
# Exit 2 on a usage error or without the table, 1 when DIR is not found in the copy.

[ $# -eq 6 ] || { echo "usage: run_one.sh SRC_TESTS CMP TAG ENGINE MODE DIR" >&2; exit 2; }
# The reference outputs were made with the universal Coulomb table: without it (unset, or the default file of another
# format) the engine uses the exact Coulomb functions and every charged-channel number moves at about 1e-7.
[ -n "${AZURE2_COULOMB_TABLE:-}" ] && [ -f "$AZURE2_COULOMB_TABLE" ] ||
  { echo "run_one.sh: set AZURE2_COULOMB_TABLE to the Coulomb table file the reference runs used" >&2; exit 2; }
SRC=$1; WORK=$2; TAG=$3; BIN=$4; MODE=$5; DIR=$6
M=${MODE:-brune}; M=${M#--}
FLAT=${DIR//\//__}; TOP=${DIR%%/*}
J=$WORK/$M/$TAG/job-$FLAT
rm -rf $J; mkdir -p $J; cp -r $SRC/$TOP $J/
cd $J/$DIR || exit 1
R=$WORK/$M/$TAG/.res-$FLAT
for azr in *.azr; do
  [ -f "$azr" ] || continue
  rm -rf output checks; mkdir -p output checks
  # menu 1, then the prompts of the run answered with their defaults, then 7 (exit); MODE unquoted: "" is no flag
  printf '1\n\n\n7\n' | OMP_NUM_THREADS=1 timeout 900 $BIN --no-gui --no-readline $MODE $azr 2>&1 | head -c 2000000 > run.log
  mkdir -p $R/$azr
  mv output checks run.log $R/$azr/ 2>/dev/null
done
rm -rf $J
