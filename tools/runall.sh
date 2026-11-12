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
#   TAG     the name of this run; every run gets a fresh tag: a tag whose results exist in one of the bases is refused
#           (exit 2), unless REPLACE_TAG=1, which removes that tag's results in the bases of this run first
#   MODE    the bases, as engine flags: default the three, "" (Brune), --no-brune (formal), --use-park
#
# Environment:
#   SRC_TESTS  the tests/ directory of the source tree whose projects are run (required)
#   WORK       scratch (required): the results land in $WORK/cmp/<basis>/<TAG>/.res-<folder>/<file>.azr/, which
#              holds output/, checks/ and run.log; <basis> is brune, no-brune or use-park, <folder> the project's
#              folder under tests/ with "/" written "__"
#   JOBS       runs at once, one thread each (default 1: the caller states its thread budget)
#   AZURE2_COULOMB_TABLE  the Coulomb table file (required, as by the reference runs; checked by run_one.sh)
#
#   REPLACE_TAG=1  run under a tag that exists, after removing its results (else refused: a folder of a project
#              that was renamed since would survive and be compared)
#
# Exit 2 on a usage error or an existing tag; 1 when a job failed (its message printed by run_one.sh; a failed run
# of the engine itself shows as missing or different files in the comparison); otherwise 0.

usage() { echo "usage: SRC_TESTS=<tests dir> WORK=<scratch> [JOBS=n] runall.sh ENGINE TAG [MODE ...]" >&2; exit 2; }
[ $# -ge 2 ] || usage
[ -n "${SRC_TESTS:-}" ] || { echo "runall.sh: set SRC_TESTS to the tests/ directory of the source tree" >&2; usage; }
[ -n "${WORK:-}" ] || { echo "runall.sh: set WORK to a scratch directory" >&2; usage; }
[ -f "${AZURE2_COULOMB_TABLE:-}" ] || { echo "runall.sh: set AZURE2_COULOMB_TABLE to the Coulomb table file" >&2; usage; }
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$1; TAG=$2; shift 2
MODES=("$@"); [ ${#MODES[@]} -eq 0 ] && MODES=("" "--no-brune" "--use-park")
# tests/convert/ holds the originals of converted projects (fixtures of the converter, without data): not run
for m in "${MODES[@]}"; do
  B=${m:-brune}; B=$WORK/cmp/${B#--}/$TAG
  if [ -e "$B" ]; then
    [ "${REPLACE_TAG:-0}" = 1 ] || { echo "runall.sh: the tag $TAG exists ($B): take a fresh tag, or REPLACE_TAG=1" >&2; exit 2; }
    rm -rf "$B"
  fi
done
DIRS=$(cd "$SRC_TESTS" && find . -name '*.azr' -not -path '*/output/*' -not -path './convert/*' | sed 's|/[^/]*$||; s|^\./||' | sort -u)
for m in "${MODES[@]}"; do
  for d in $DIRS; do echo "$SRC_TESTS $WORK/cmp $TAG $BIN '$m' $d"; done
done | xargs -P "${JOBS:-1}" -L 1 bash "$HERE/run_one.sh"
rc=$?
if [ $rc -ne 0 ]; then
  echo "runs done $TAG: FAILED (a job ended with an error, above; xargs status $rc)" >&2
  exit 1
fi
echo runs done $TAG
