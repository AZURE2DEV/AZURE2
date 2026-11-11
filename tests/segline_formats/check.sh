#!/usr/bin/env bash
#
# Segment-line layouts the reader must accept (audit 2026-10-25, 06_parameters_outputs_api C7/C8, checklist item 29;
# 2026-10-27).
#
# A <segmentsData> line's three shift fields (shift shiftError vary) are optional: older files carry one shift value or
# none before the data file.  SegLine.h decided the layout by whether ">> double" succeeded, so (i) a file whose name
# starts with digits in the no-shift layout ("... normErr 12C_p.dat") was read as a 12 MeV shift with the file "C_p.dat",
# and (ii) the one-shift layout ("... normErr 0.001 data/foo.dat") left the stream failed and the file empty -- in both
# cases EData printed "Could Not Fill Segment" and the fit went on without the data.  Now the layout is decided by whole
# tokens.  And a <segmentsTest> line's keyword reader stopped at a data-only block (sqrtshift, shiftgroup, ...), losing a
# "sample" or "total" written after it; those blocks are skipped now.
#
#   (a) tests/13N with segment 3 in the no-shift layout and its file 84deg_meyer.dat in the model directory: the total equals
#       expected/chiSquared.out (the segment is read, no shift);
#   (b) segment 3 in the one-shift layout "0.001 data/meyer_84.dat": the total equals the three-field layout
#       "0.001 0 0 data/meyer_84.dat" and differs from (a) (the 1 keV shift is applied, so the check has teeth);
#   (c) tests/natural_neon in mode 3 with an n + 20Ne test segment written "sqrtshift 0 0 0 shiftgroup G sample 1":
#       the extrapolation equals the one written "sample 1" (weighted by the abundance), which itself differs from the
#       bare line.
# dev e7d72bd: (a) and (b) drop the segment, (c) loses the sample weighting.
#
#   ./tests/segline_formats/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../13N"
NEON="$HERE/../natural_neon"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/segline_formats.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
near() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(a != "" && b != "" && d <= t) }'; }

SEG3='1              1              1              0              5              84.3           84.3           1              1              0              5              0              0.005          0               data/meyer_84.dat 0 0'
grep -qF "$SEG3" "$BASE/13N.azr" || { echo "  FAIL  tests/13N/13N.azr no longer holds segment 3 as '$SEG3'"; exit 1; }
# model NAME NEW-SEGMENT-3-LINE [RENAMED-FILE]
model() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$BASE/data" "$WORK/$1/"
  [ -n "${3:-}" ] && cp "$BASE/data/meyer_84.dat" "$WORK/$1/$3"   # in the model directory: the name itself starts with digits
  awk -v old="$SEG3" -v new="$2" '{ sub(/\r$/, "") } index($0, old) == 1 { print new; next } { print }' "$BASE/13N.azr" > "$WORK/$1/t.azr"
}
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
total() { awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2> /dev/null; }
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$BASE/expected/chiSquared.out")"

model noshift "1 1 1 0 5 84.3 84.3 1 1 0 5 84deg_meyer.dat" 84deg_meyer.dat
model oneshift "1 1 1 0 5 84.3 84.3 1 1 0 5 0.001 data/meyer_84.dat"
model threefield "1 1 1 0 5 84.3 84.3 1 1 0 5 0.001 0 0 data/meyer_84.dat 0 0"
run noshift; run oneshift; run threefield
a="$(total noshift)"; b="$(total oneshift)"; c="$(total threefield)"
if near "$a" "$want" 1e-5 && ! grep -q "Could Not Fill" "$WORK/noshift/run.log"; then ok "(a) no-shift layout with a file name starting with digits: segment read, total $a = $want"; else bad "(a) no-shift layout with '84deg_meyer.dat': total '$a' vs $want (filled: $(grep -c 'Could Not Fill' "$WORK/noshift/run.log") 'Could Not Fill')"; fi
if near "$b" "$c" 1e-9 && ! near "$b" "$want" 1e-5 && ! grep -q "Could Not Fill" "$WORK/oneshift/run.log"; then ok "(b) one-shift layout: total $b = the three-field layout $c, and the 1 keV shift moved it from $want"; else bad "(b) one-shift layout: total '$b', three-field '$c', unshifted $want"; fi

# ---- (c) a test line with data-only keyword blocks before "sample 1" -------------------------------------------------
# (natural_neon has an <abundances> block: "sample 1" on a plain n + 20Ne test line weights it by the 20Ne abundance)
tmodel() {  # NAME TEST-LINE
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$NEON/data" "$WORK/$1/"
  awk -v line="$2" '{ sub(/\r$/, "") } /^<segmentsTest>$/ { print; print line; skip = 1; next } /^<\/segmentsTest>$/ { skip = 0 } !skip { print }' "$NEON/natural_neon.azr" > "$WORK/$1/t.azr"
}
trun() { (cd "$WORK/$1" && printf '3\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
tmodel tbare "1 2 2 0.5 1.0 0.1 0 0 0 0"
tmodel tplain "1 2 2 0.5 1.0 0.1 0 0 0 0 sample 1"
tmodel tkeys "1 2 2 0.5 1.0 0.1 0 0 0 0 sqrtshift 0 0 0 shiftgroup G sample 1"
trun tbare; trun tplain; trun tkeys
ext() { awk '{ sub(/\r$/, "") } NF >= 4 { printf "%.6f %.10e\n", $1, $4 }' "$WORK/$1/output/AZUREOut_aa=2_R=2.extrap" 2> /dev/null; }
ext tbare > "$WORK/tbare.x"; ext tplain > "$WORK/tplain.x"; ext tkeys > "$WORK/tkeys.x"
n="$(wc -l < "$WORK/tplain.x" | tr -d ' ')"
maxrel() { paste "$1" "$2" | awk 'BEGIN { m = 0 } { r = ($2 - $4) / $2; if (r < 0) r = -r; if (r > m) m = r } END { printf "%.3e", m }'; }
w="$(maxrel "$WORK/tbare.x" "$WORK/tplain.x")"; d="$(maxrel "$WORK/tplain.x" "$WORK/tkeys.x")"
if [ "$n" -ge 5 ] && awk -v w="$w" 'BEGIN { exit !(w > 0.01) }'; then ok "(c) 'sample 1' weights the plain n + 20Ne test line by the abundance (moves it by $w)"; else bad "(c) 'sample 1' does not change the test line ($w): the check has no teeth"; fi
if [ "$n" -ge 5 ] && awk -v d="$d" 'BEGIN { exit !(d <= 1e-9) }'; then ok "(c) with sqrtshift and shiftgroup before 'sample 1' the line equals the plain 'sample 1' line ($n energies, $d)"; else bad "(c) the test line with data-only keywords differs from the plain one by $d ($n energies): 'sample 1' was lost"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
