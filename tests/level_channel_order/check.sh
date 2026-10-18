#!/usr/bin/env bash
#
# A level's width lines are matched to its J-group's channels by (pair, l, s), not
# by their position in the file.  Two edits of the 13N reference model must give
# exactly the reference chi-squared (CLI defaults; tests/pyazr/channel_order_test.py
# covers both parametrizations through pyazr):
#   reordered  the three channel lines of the 3/2- level at 20 MeV (level ID 6, the
#              second level of its group, which takes its channel order from ID 5)
#              listed last-first;
#   short      the 5/2+ level at 20 MeV (ID 10) without the line of its zero-width
#              p+12C(gamma) l=4 channel (the middle channel of its group).
# Before the fix the widths were stored in line order: "reordered" put
# the 5 MeV proton width on a capture channel (chi2 312,264 instead of 336,521), and
# "short" put the l=6 width in the l=4 slot and left the width vector one short
# (out-of-bounds write in TransformIn).  Both of those widths are zero here, so the
# old binary's "short" chi2 is still right; there only the warning tells them apart
# (valgrind sees the write).  Each edit must print its warning.
#
#   ./tests/level_channel_order/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# macOS has no timeout(1); run without a limit there
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/level_channel_order.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# The edited models, built from the reference at run time.
awk '/<levels>/ { f = 1; print; next } /<\/levels>/ { f = 0 }
     f && $9 == 6 { held[++n] = $0; if (n == 3) for (i = 3; i >= 1; i--) print held[i]; next }
     { print }' "$REF/13N.azr" > "$WORK/reordered.azr"
awk '/<levels>/ { f = 1; print; next } /<\/levels>/ { f = 0 }
     f && $9 == 10 && $6 == 2 && $8 == 4 { next }
     { print }' "$REF/13N.azr" > "$WORK/short.azr"
# reversing three lines leaves the middle one in place: two lines differ
[ "$(diff "$REF/13N.azr" "$WORK/reordered.azr" | grep -c '^>')" = 2 ] || { bad "reordered model not built"; exit 1; }
[ "$(diff "$REF/13N.azr" "$WORK/short.azr" | grep -c '^<')" = 1 ] || { bad "short model not built"; exit 1; }

# run NAME MODEL -> $WORK/NAME/output/chiSquared.out, $WORK/NAME/run.log
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$REF/data" "$dir/data"
  cp "$2" "$dir/t.azr"
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -s "$dir/output/chiSquared.out" ]; then
    bad "$1: no chiSquared.out written"; tail -15 "$dir/run.log"; return 1
  fi
}
total() { grep -o 'Total-Chi-Squared: [^ ]*' "$1/output/chiSquared.out" | awk '{ print $2 }'; }

run ref       "$REF/13N.azr"        || exit 1
run reordered "$WORK/reordered.azr" || exit 1
run short     "$WORK/short.azr"     || exit 1
ref="$(total "$WORK/ref")"
for m in reordered short; do
  if cmp -s "$WORK/ref/output/chiSquared.out" "$WORK/$m/output/chiSquared.out"; then
    ok "$m: chiSquared.out identical to the reference (total $ref)"
  else
    bad "$m: total chi2 $(total "$WORK/$m"), reference $ref"
  fi
done
if grep -q 'lists its channels in a different order from its J-group' "$WORK/reordered/run.log"; then
  ok "reordered: warning printed"; else bad "reordered: no reordering warning"; fi
if grep -q 'has no line for channel 2 of its J-group' "$WORK/short/run.log"; then
  ok "short: warning printed"; else bad "short: no missing-channel warning"; fi
if grep -q 'WARNING: level' "$WORK/ref/run.log"; then bad "reference: unexpected level warning"; fi

if [ "$failures" -eq 0 ]; then echo "PASS: level channel order"; else echo "FAIL: $failures check(s)"; exit 1; fi
