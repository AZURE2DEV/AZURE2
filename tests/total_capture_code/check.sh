#!/usr/bin/env bash
#
# Total capture is the exit key -1 together with the total-capture data type (audit 2026-10-25, item 32; 2026-10-27).
#
# The reader expands a total-capture segment (exit key -1 AND data type 3 on a data line / 4 on a test line) into one
# term per capture pair.  The key -1 with another data type left a segment whose exit pair does not exist:
# GetPair(GetPairNumFromKey(-1)) read pairs_[-1], and in mode 3 the external-capture lookup for that point crashed
# (SIGSEGV on tests/13N).  In format 3 the two halves are one token: the residual state [total] of the reaction
# (12C(p,g)13N[total]), so the mismatch cannot be written; the old lines that had it are the converter's
# (docs/upgrade/format3_tree.md section 4).
#
#   (c) tests/13N with two extrapolations, 12C(p,g)13N[total] and 12C(p,g)13N, mode 3: AZUREOut_aa=1_TOTAL_CAPTURE.extrap
#       equals AZUREOut_aa=1_R=2.extrap (13N has one capture pair) point by point;
#   (d) the old test line "1 1 2 0.5 2.0 0.5 0 0 0 4" (the total-capture type with a real exit key, which the old engine
#       ran as total capture with a warning), converted with AZURE2 --convert ... --to 3: a [total] extrapolation that
#       equals (c);
#   (a), (b) the old lines "1 1 -1 0.5 2.0 0.5 0 0 0 0" (test line) and a data line with exit -1 and type 0, which the
#       old engine refused: AZURE2 --convert refuses them with the old engine's reason and a non-zero exit, no file.
# dev 510e6d1: the old lines of (c) with exit -1 and type 0 gave a segmentation fault (exit 139) with no output.
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
# model NAME RECORDS: tests/13N with the given extrapolation records (format 3) as its only extrapolations, mode 3
model() {
  local name="$1"
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"; cp -r "$BASE/data" "$WORK/$name/"
  # the records go through the environment: BSD awk (macOS) refuses a newline in -v ("newline in string")
  export RECORDS="$2"
  awk '{ sub(/\r$/, "") } /^<\/extrapolations>$/ { print ENVIRON["RECORDS"] } { print }' "$BASE/13N.azr" > "$WORK/$name/m.azr"
  run "$name"
}
run() { (cd "$WORK/$1" && printf '3\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1; echo "exit $?" >> run.log); }
status() { tail -1 "$WORK/$1/run.log"; }
same() { [ -s "$1" ] && [ -s "$2" ] && paste <(awk 'NF{print $4}' "$1") <(awk 'NF{print $4}' "$2") | awk '{ d = $1 - $2; if (d < 0) d = -d; if (d > 1e-9 * ($2 > 0 ? $2 : 1)) bad = 1 } END { exit bad }'; }
grid="    quantity       sigma
    unit           b
    energies       range 0.5 2.0 MeV lab step 0.5"

model c "extrapolation total
    reaction       12C(p,g)13N[total]
$grid
extrapolation to_pair_2
    reaction       12C(p,g)13N
$grid"
ft="$WORK/c/output/AZUREOut_aa=1_TOTAL_CAPTURE.extrap"; fr="$WORK/c/output/AZUREOut_aa=1_R=2.extrap"
if same "$ft" "$fr"; then
  ok "(c) total capture ([total]) equals capture to the single gamma pair ($(awk 'NF' "$ft" | wc -l | tr -d ' ') energies)"
else bad "(c) total capture vs capture to pair 2: $(status c); $(head -2 "$ft" 2> /dev/null | awk '{print $4}' | tr '\n' ' ') vs $(head -2 "$fr" 2> /dev/null | awk '{print $4}' | tr '\n' ' ')"; fi

# (d) the old line, through the converter: its original project (tests/convert/fixtures) with the line as its only
# test line and its <targetInt> lines inactive (they fold test segments by number)
OLD="$HERE/../convert/fixtures/tests/13N/13N.azr"
mkdir -p "$WORK/d/output" "$WORK/d/checks"; cp -r "$BASE/data" "$WORK/d/"
awk '{ sub(/\r$/, "") } /<segmentsTest>/ { print; print "1 1 2 0.5 2.0 0.5 0 0 0 4"; skip = 1; next } /<\/segmentsTest>/ { skip = 0 }
     /<targetInt>/ { p = 1; print; next } /<\/targetInt>/ { p = 0 } p { $1 = 0 } !skip { print }' "$OLD" > "$WORK/d/old.azr"
(cd "$WORK/d" && run_limited 300 "$AZURE2_BIN" --convert old.azr m.azr --to 3 > convert.log 2>&1)
fd="$WORK/d/output/AZUREOut_aa=1_TOTAL_CAPTURE.extrap"
if grep -q '^    reaction       12C(p,g)13N\[total\]$' "$WORK/d/m.azr" 2> /dev/null && run d && same "$ft" "$fd"; then
  ok "(d) the old total-capture type with exit pair 2 converts to [total] and equals (c)"
else bad "(d) type 4 with exit 2: $(tail -1 "$WORK/d/convert.log" 2> /dev/null); $(status d); file $(ls -la "$fd" 2> /dev/null | awk '{print $5}') bytes"; fi

# (a), (b) the old lines with exit -1 and a type that is not total capture, which the old engine refused: the converter
# refuses them too, with the old engine's reason and a non-zero exit, and writes no file
refused() {  # NAME SECTION LINE: the original with LINE as the first line of SECTION (segmentsTest or segmentsData)
  local d="$WORK/$1"
  mkdir -p "$d"; cp -r "$BASE/data" "$d/"
  awk -v sec="<$2>" -v line="$3" '{ sub(/\r$/, "") } { print } $0 == sec { print line }' "$OLD" > "$d/old.azr"
  (cd "$d" && run_limited 300 "$AZURE2_BIN" --convert old.azr m.azr --to 3 > convert.log 2>&1; echo "exit $?" >> convert.log)
  [ ! -f "$d/m.azr" ] && ! grep -q '^exit 0$' "$d/convert.log" && grep -q "the exit key -1 (total capture) needs the total-capture data type" "$d/convert.log"
}
if refused a segmentsTest "1 1 -1 0.5 2.0 0.5 0 0 0 0"; then ok "(a) test line with exit -1 and type 0: the converter refuses it with the old engine's reason"
else bad "(a) exit -1 with type 0 on a test line: $(tail -2 "$WORK/a/convert.log" | tr '\n' ' ')"; fi
if refused b segmentsData "1 1 -1 0 5 0 180 0 1 0 0 0 0 0 data/artemov.dat 0 0"; then ok "(b) data line with exit -1 and type 0: the converter refuses it with the old engine's reason"
else bad "(b) exit -1 with type 0 on a data line: $(tail -2 "$WORK/b/convert.log" | tr '\n' ' ')"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
