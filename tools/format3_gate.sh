#!/usr/bin/env bash
#
# The converted-tree gate of project format 3 (docs/upgrade/project_format3.md section 8.5): the tests/ tree of a
# source tree is converted to format 3 with --exact-grids (tools/convert_tree.sh, without the old engine's pass), made
# runnable (each original replaced by its conversion, tools/format3_gate.py tree), run in the three bases (menu 1,
# runall.sh) and every output file compared with those of the original tree's run under the reference tag
# (tools/format3_gate.py compare): identical, or a difference of a class listed in format3_gate.py's CLASSES.
#
#   tools/format3_gate.sh <AZURE2 engine> <source tree (repository root)> <tag>
#
# Environment (as runall.sh):
#   WORK     scratch (required): the runs land in $WORK/cmp/<basis>/tab-<tag>, the converted tree in $WORK/gate-<tag>
#   JOBS     parallel runs and conversions, one thread each (default 1)
#   REF      the reference tag in $WORK/cmp, the original tree's run by tools/runall.sh with the same engine (required)
#   RUNALL   the runall.sh to run the converted tree with (default: tools/runall.sh, beside this script)
#
# Prints the counts per basis; the full list in $WORK/gate-<tag>/gate.txt.  Exit 0 when no file differs and none is
# missing, 1 otherwise, 2 on a usage or set-up error.

usage() { echo "usage: WORK=<scratch> REF=<reference tag> [JOBS=n] format3_gate.sh ENGINE SOURCE_TREE TAG${1:+ ($1)}" >&2; exit 2; }
[ $# -eq 3 ] || usage
[ -n "${WORK:-}" ] || usage "set WORK to the scratch directory of runall.sh"
[ -n "${REF:-}" ] || usage "set REF to the reference tag"
ENGINE=$1; SRC=$2; TAG=$3
HERE=$(cd "$(dirname "$0")" && pwd)
RUNALL="${RUNALL:-$HERE/runall.sh}"
JOBS="${JOBS:-1}"
ENGINE="$(cd "$(dirname "$ENGINE")" && pwd)/$(basename "$ENGINE")"
SRC="$(cd "$SRC" && pwd)"
WORK="$(cd "$WORK" && pwd)"
if [ ! -x "$ENGINE" ] || [ ! -d "$SRC/tests" ] || [ ! -f "$RUNALL" ]; then
  echo "format3_gate.sh: need an engine ($ENGINE), a source tree with tests/ ($SRC) and runall.sh ($RUNALL)" >&2
  exit 2
fi
for b in brune no-brune use-park; do
  if [ ! -d "$WORK/cmp/$b/$REF" ]; then
    echo "format3_gate.sh: no reference run $WORK/cmp/$b/$REF" >&2
    exit 2
  fi
done

G="$WORK/gate-$TAG"
rm -rf "$G"
mkdir -p "$G"
JOBS=$JOBS OLD_ENGINE=0 CONVERT_OPTIONS="--exact-grids" "$HERE/convert_tree.sh" "$ENGINE" "$SRC" "$G/converted" \
  > "$G/convert.log" 2>&1
python3 "$HERE/format3_gate.py" tree "$G/converted" "$G/tree" > "$G/tree.txt" || exit 2
for b in brune no-brune use-park; do rm -rf "$WORK/cmp/$b/tab-$TAG"; done
SRC_TESTS="$G/tree/tests" WORK="$WORK" JOBS="$JOBS" bash "$RUNALL" "$ENGINE" "tab-$TAG" > "$G/runall.log" 2>&1
python3 "$HERE/format3_gate.py" compare "$WORK/cmp" "$REF" "tab-$TAG" "$G/converted" > "$G/gate.txt"
status=$?
grep '^##' "$G/gate.txt"
exit $status
