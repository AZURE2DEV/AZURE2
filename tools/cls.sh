#!/bin/bash
#
# Compare the outputs of two runs of tools/runall.sh file by file (tools/classify.py), in the three bases.
#
#   tools/cls.sh REFTAG NEWTAG
#
# Environment: WORK, the scratch directory runall.sh wrote to (required; the runs are in $WORK/cmp/<basis>/<TAG>).
# Output: the table in $WORK/cls-NEWTAG.tsv (one line per file that differs, then the "## " summary lines per basis);
# its last five lines are printed.  Exit 2 on a usage error, else 0 (the table says what differs).

[ $# -eq 2 ] || { echo "usage: WORK=<scratch> cls.sh REFTAG NEWTAG" >&2; exit 2; }
[ -n "${WORK:-}" ] || { echo "cls.sh: set WORK to the scratch directory used by runall.sh" >&2; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd)
python3 "$HERE/classify.py" "$WORK/cmp" "$1" "$2" > "$WORK/cls-$2.tsv"
tail -5 "$WORK/cls-$2.tsv"
