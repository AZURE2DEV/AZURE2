#!/bin/bash
#
# Run every project of a test tree with one engine, in each basis, for a file-by-file comparison of two engines
# (tools/cls.sh, tools/classify.py).  Every folder of $SRC_TESTS that holds an .azr (output/ folders left out) is one
# job of tools/run_one.sh: menu 1 (calculate with data), one OpenMP thread, in a private copy of the folder's top-level
# test folder.
#
#   tools/runall.sh ENGINE TAG [MODE ...]
#
#   ENGINE  the AZURE2 binary, by an absolute path (each job runs it from its own copy of the tree)
#   TAG     the name of this run; every run gets a fresh tag (the results of an earlier run with it are not removed)
#   MODE    the bases, as engine flags: default the three, "" (Brune), --no-brune (formal), --use-park
#
# Environment:
#   SRC_TESTS  the tests/ directory of the source tree whose projects are run (required)
#   WORK       scratch (required): the results land in $WORK/cmp/<basis>/<TAG>/.res-<folder>/<file>.azr/, which
#              holds output/, checks/ and run.log; <basis> is brune, no-brune or use-park, <folder> the project's
#              folder under tests/ with "/" written "__"
#   JOBS       runs at once, one thread each (default 1: the caller states its thread budget)
#
# Exit 2 on a usage error; otherwise 0 (a failed run shows as missing or different files in the comparison).

usage() { echo "usage: SRC_TESTS=<tests dir> WORK=<scratch> [JOBS=n] runall.sh ENGINE TAG [MODE ...]" >&2; exit 2; }
[ $# -ge 2 ] || usage
[ -n "${SRC_TESTS:-}" ] || { echo "runall.sh: set SRC_TESTS to the tests/ directory of the source tree" >&2; usage; }
[ -n "${WORK:-}" ] || { echo "runall.sh: set WORK to a scratch directory" >&2; usage; }
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$1; TAG=$2; shift 2
MODES=("$@"); [ ${#MODES[@]} -eq 0 ] && MODES=("" "--no-brune" "--use-park")
DIRS=$(cd "$SRC_TESTS" && find . -name '*.azr' -not -path '*/output/*' | sed 's|/[^/]*$||; s|^\./||' | sort -u)
for m in "${MODES[@]}"; do
  for d in $DIRS; do echo "$SRC_TESTS $WORK/cmp $TAG $BIN '$m' $d"; done
done | xargs -P "${JOBS:-1}" -L 1 bash "$HERE/run_one.sh"
echo runs done $TAG
