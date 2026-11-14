#!/usr/bin/env bash
#
# AZURE2 --validate predicts what a run refuses (docs/upgrade/project_format3.md section 6.2): each variant of two small
# format-3 projects is given to --validate and then run (menu 1); validation must report it with its code where the
# run refuses it (or, for a row the run drops, prints the run's message too).
#
# Models (format 3, *.f3 so that tools/runall.sh does not take them as projects; copied to t.azr in a scratch dir):
#   rth.f3     : alpha + 12C ratio to Rutherford on the hard sphere (tests/ratio_to_rutherford), data/rth.dat
#   capture.f3 : 16O(p,g)17F(0.495 MeV) with E1 external capture (tests/ec_final_state_listed), data/ecm.dat
# Variants (made here):
#   (a) rth.f3 as it is: valid, runs
#   (b) two rows of data/rth.dat that do not parse: E_DATA_ROW twice (each row, with line and column); the run stops
#       reading the file at the first ("Cannot parse line 3") and drops the data set, its only one: nothing is computed
#   (c) a covariance file that does not exist: E_COVARIANCE_FILE; the run is refused
#   (d) a covariance file with a variance 0 on a row inside the window: E_COVARIANCE_FILE; the run is refused
#   (e) a row with uncertainty 0: W_DATA_ERROR (a warning, the file is valid); the run skips the point in the chi^2
#   (f) capture.f3 without its final-state level: E_FINAL_STATE_NOT_LEVEL; the run is refused
#   (g) a second level of the final state's J^pi 0.2 keV away: E_FINAL_STATE_NOT_LEVEL (ambiguous); refused
#   (h) a second photon pair with the same final state: E_FINAL_STATE_SHARED; refused
#   (i) --describe --record rth: that record alone, its open energy window in words; an unknown name: exit 2
# "Refused" = an ERROR line of the run and no chiSquared.out.
#
#   ./tests/validate_run_refusals/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/validate_run_refusals.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

setup() {  # NAME MODEL: $WORK/NAME/t.azr and its data
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$HERE/models/data" "$WORK/$1/"
  sed 's/\r$//' "$HERE/models/$2.f3" > "$WORK/$1/t.azr"
}
validate() {  # NAME: --validate output in validate.log, its exit code in validate.exit
  (cd "$WORK/$1" && run_limited 120 "$AZURE2_BIN" --validate t.azr > validate.log 2>&1; echo $? > validate.exit)
}
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
codes() { grep -c "$2" "$WORK/$1/validate.log"; }
refused() { grep -q "$2" "$WORK/$1/run.log" && [ ! -s "$WORK/$1/output/chiSquared.out" ]; }
expect() {  # LABEL NAME CODE COUNT EXIT RUN-PATTERN REFUSED(1/0)
  local label=$1 name=$2 code=$3 count=$4 exit=$5 pattern=$6 refuse=$7
  validate "$name"; run "$name"
  local n; n=$(codes "$name" "$code")
  local e; e=$(cat "$WORK/$name/validate.exit")
  local r=0; refused "$name" "$pattern" && r=1
  if [ "$n" = "$count" ] && [ "$e" = "$exit" ] && grep -q "$pattern" "$WORK/$name/run.log" && [ "$r" = "$refuse" ]; then
    ok "$label"
  else
    bad "$label: $n x $code (want $count), --validate exit $e (want $exit), run: refused=$r (want $refuse)"
    grep -E "error|warning" "$WORK/$name/validate.log" | sed -n '1,3p'; grep -E "ERROR|WARNING" "$WORK/$name/run.log" | sed -n '1,3p'
  fi
}

echo "validate_run_refusals: $AZURE2_BIN"
setup a rth; validate a; run a
if [ "$(cat "$WORK/a/validate.exit")" = 0 ] && [ -s "$WORK/a/output/chiSquared.out" ]; then ok "(a) rth.f3: valid, runs"
else bad "(a) rth.f3: --validate exit $(cat "$WORK/a/validate.exit"), chiSquared.out $([ -s "$WORK/a/output/chiSquared.out" ] && echo written || echo absent)"; fi

setup b rth
awk 'NR == 3 { $2 = "2O" } NR == 6 { $0 = $1 " " $2 } { print }' "$HERE/models/data/rth.dat" > "$WORK/b/data/rth.dat"
expect "(b) two rows that do not parse: both reported (E_DATA_ROW); the run stops at the first and drops the data set" b \
  "E_DATA_ROW" 2 2 "Cannot parse line 3" 1
grep -q "line 3, column 2: '2O' is not a number" "$WORK/b/validate.log" && grep -q "line 6: 2 column(s)" "$WORK/b/validate.log" ||
  bad "(b) the rows' lines and columns: $(grep -c E_DATA_ROW "$WORK/b/validate.log") lines"

setup c rth
sed -i 's|^    norm           1 fixed|&\n    covariance     cov.txt|' "$WORK/c/t.azr"
expect "(c) a covariance file that does not exist: E_COVARIANCE_FILE, refused" c "E_COVARIANCE_FILE" 1 2 "ERROR: Covariance file" 1

setup d rth
sed -i 's|^    norm           1 fixed|&\n    covariance     cov.txt|' "$WORK/d/t.azr"
n=$(grep -cv '^#' "$HERE/models/data/rth.dat")
{ echo "file data/rth.dat $n"; printf 'diagonal'; for i in $(seq 1 "$n"); do [ "$i" = 3 ] && printf ' 0' || printf ' 0.0009'; done; echo; } > "$WORK/d/cov.txt"
expect "(d) a variance 0 on a row inside the window: E_COVARIANCE_FILE, refused" d "E_COVARIANCE_FILE" 1 2 "has a variance 0 <= 0" 1

setup e rth
awk 'NR == 5 { $4 = 0 } { print }' "$HERE/models/data/rth.dat" > "$WORK/e/data/rth.dat"
validate e; run e
if [ "$(codes e W_DATA_ERROR)" = 1 ] && [ "$(cat "$WORK/e/validate.exit")" = 0 ] && [ -s "$WORK/e/output/chiSquared.out" ]; then
  ok "(e) an uncertainty 0: W_DATA_ERROR, valid, the run computes"
else bad "(e) an uncertainty 0: $(codes e W_DATA_ERROR) W_DATA_ERROR, --validate exit $(cat "$WORK/e/validate.exit")"; fi

setup f capture
awk '!($1 == 0.5 && $2 == 1 && $3 == 0.4953)' "$HERE/models/capture.f3" > "$WORK/f/t.azr"
expect "(f) no final-state level: E_FINAL_STATE_NOT_LEVEL, refused" f "E_FINAL_STATE_NOT_LEVEL" 1 2 \
  "ERROR: external capture into pair 3 needs its final state" 1

setup g capture
extra="$(awk '$1 == 0.5 && $2 == 1 && $3 == 0.4953 { $3 = 0.4955; $9 = 9; print }' "$HERE/models/capture.f3")"
awk -v e="$extra" '/^<\/levels>/ { print e; print "" } { print }' "$HERE/models/capture.f3" > "$WORK/g/t.azr"
expect "(g) two levels within 1 keV of the final state: E_FINAL_STATE_NOT_LEVEL, refused" g "E_FINAL_STATE_NOT_LEVEL" 1 2 \
  "lie within 1 keV of its final-state energy" 1

setup h capture
awk '$1 == 0.5 && $2 == 1 && $3 == 0.4953 && $6 == 3 { print; $6 = 4; print; next }
     /^    norm / { print; print "    exit_pair      3"; next } { print }' "$HERE/models/capture.f3" > "$WORK/h/t.azr"
expect "(h) two photon pairs with one final state: E_FINAL_STATE_SHARED, refused" h "E_FINAL_STATE_SHARED" 1 2 \
  "ERROR: external capture into pairs 3 and 4" 1

# (i) --describe --record NAME: that record alone; an open energy window in words; an unknown name: exit 2
(cd "$WORK/a" && "$AZURE2_BIN" --describe t.azr --record rth > describe.txt 2>&1; echo $? > describe.exit;
 "$AZURE2_BIN" --describe t.azr --record nope > nope.txt 2>&1; echo $? > nope.exit)
if [ "$(cat "$WORK/a/describe.exit")" = 0 ] && [ "$(grep -c '^dataset ' "$WORK/a/describe.txt")" = 1 ] &&
   grep -q "window E_lab all energies (no energy_window)" "$WORK/a/describe.txt" && [ "$(cat "$WORK/a/nope.exit")" = 2 ]; then
  ok "(i) --describe --record: one record, the open window in words, an unknown name refused"
else bad "(i) --describe --record: exit $(cat "$WORK/a/describe.exit")/$(cat "$WORK/a/nope.exit"), $(head -c 300 "$WORK/a/describe.txt")"; fi

if [ "$failures" -eq 0 ]; then echo "validate_run_refusals: PASS"; else echo "validate_run_refusals: FAIL ($failures)"; exit 1; fi
