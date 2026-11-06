#!/usr/bin/env bash
#
# Total capture is the exit key -1 together with the total-capture data type (audit 2026-10-25, item 32; 2026-10-27).
#
# The reader expands a total-capture segment (exit key -1 AND data type 3 on a data line / 4 on a test line) into one
# term per capture pair.  The key -1 with another data type left a segment whose exit pair does not exist:
# GetPair(GetPairNumFromKey(-1)) read pairs_[-1], and in mode 3 the external-capture lookup for that point crashed
# (SIGSEGV on tests/13N; a model without external capture gave meaningless numbers instead).  Both halves are now
# checked against each other and the mismatch is refused.
#
#   (a) tests/13N + a TEST line "1 1 -1 0.5 2.0 0.5 0 0 0 0" (exit -1, type 0), mode 3: ERROR, no output;
#   (b) a DATA line with exit -1 and type 0 (angle-integrated), mode 1: ERROR, no chiSquared.out;
#   (c) the proper test line "1 1 -1 0.5 2.0 0.5 0 0 0 4" beside "1 1 2 ... 0": AZUREOut_aa=1_TOTAL_CAPTURE.extrap
#       equals AZUREOut_aa=1_R=2.extrap (13N has one capture pair) point by point;
#   (d) the type with a real exit key, "1 1 2 0.5 2.0 0.5 0 0 0 4": runs as total capture with a WARNING that the key
#       is ignored (as every build did silently), and equals (c).
# dev 510e6d1: (a) segmentation fault (exit 139) with no output; (b) the same read of pairs_[-1].
#
#   ./tests/total_capture_code/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/total_capture_code.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
# model NAME MODE LINE...: tests/13N with the given lines as its only test lines (mode 3) or an extra data line (mode 1);
# every <targetInt> line inactive (they fold test segments by number).
model() {
  local name="$1" mode="$2"; shift 2
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"; cp -r "$BASE/data" "$WORK/$name/"
  # the lines go through the environment: BSD awk (macOS) refuses a newline in -v ("newline in string")
  export LINES="$(printf '%s\n' "$@")"
  if [ "$mode" = 3 ]; then
    awk '/<segmentsTest>/{print; print ENVIRON["LINES"]; skip=1; next} /<\/segmentsTest>/{skip=0} !skip{print}' "$BASE/13N.azr"
  else
    awk '{print} /<segmentsData>/{print ENVIRON["LINES"]}' "$BASE/13N.azr"
  fi | awk '/<targetInt>/{p=1; print; next} /<\/targetInt>/{p=0} p{ $1 = 0 } {print}' > "$WORK/$name/m.azr"
  (cd "$WORK/$name" && printf '%s\n\n\n7\n' "$mode" | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1; echo "exit $?" >> run.log)
}
errors() { grep -c "total-capture data type" "$WORK/$1/run.log"; }
status() { tail -1 "$WORK/$1/run.log"; }

model a 3 "1 1 -1 0.5 2.0 0.5 0 0 0 0"
if [ "$(errors a)" -ge 1 ] && [ -z "$(ls "$WORK/a/output" | grep AZUREOut)" ]; then ok "(a) test line with exit -1 and data type 0: refused, no output"
else bad "(a) exit -1 with type 0 on a test line: errors $(errors a), $(status a), files $(ls "$WORK/a/output" | tr '\n' ' ')"; fi

model b 1 "1 1 -1 0 5 0 180 0 1 0 0 0 0 0 data/artemov.dat 0 0"
if [ "$(errors b)" -ge 1 ] && [ ! -s "$WORK/b/output/chiSquared.out" ]; then ok "(b) data line with exit -1 and data type 0: refused, no chiSquared.out"
else bad "(b) exit -1 with type 0 on a data line: errors $(errors b), $(status b), total $(grep '^Total' "$WORK/b/output/chiSquared.out" 2> /dev/null)"; fi

model c 3 "1 1 -1 0.5 2.0 0.5 0 0 0 4" "1 1 2 0.5 2.0 0.5 0 0 0 0"
ft="$WORK/c/output/AZUREOut_aa=1_TOTAL_CAPTURE.extrap"; fr="$WORK/c/output/AZUREOut_aa=1_R=2.extrap"
if [ -s "$ft" ] && [ -s "$fr" ] && paste <(awk 'NF{print $4}' "$ft") <(awk 'NF{print $4}' "$fr") | awk '{ d = $1 - $2; if (d < 0) d = -d; if (d > 1e-9 * ($2 > 0 ? $2 : 1)) bad = 1 } END { exit bad }'; then
  ok "(c) total capture with data type 4 equals capture to the single gamma pair ($(awk 'NF' "$ft" | wc -l | tr -d ' ') energies)"
else bad "(c) total capture (type 4) vs capture to pair 2: $(status c); $(head -2 "$ft" 2> /dev/null | awk '{print $4}' | tr '\n' ' ') vs $(head -2 "$fr" 2> /dev/null | awk '{print $4}' | tr '\n' ' ')"; fi

model d 3 "1 1 2 0.5 2.0 0.5 0 0 0 4"
fd="$(ls "$WORK/d/output"/AZUREOut_aa=1_TOTAL_CAPTURE.extrap "$WORK/d/output"/AZUREOut_aa=1_R=2.extrap 2> /dev/null | head -1)"
if [ -s "$fd" ] && [ "$(grep -c 'the exit key is ignored' "$WORK/d/run.log")" -ge 1 ] && paste <(awk 'NF{print $4}' "$ft") <(awk 'NF{print $4}' "$fd") | awk '{ d = $1 - $2; if (d < 0) d = -d; if (d > 1e-9 * ($2 > 0 ? $2 : 1)) bad = 1 } END { exit bad }'; then
  ok "(d) the total-capture type with exit pair 2: runs as total capture with a warning, equals (c)"
else bad "(d) type 4 with exit 2: $(status d); warnings $(grep -c 'the exit key is ignored' "$WORK/d/run.log"); file $(ls -la "$fd" 2> /dev/null | awk '{print $5}') bytes"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
