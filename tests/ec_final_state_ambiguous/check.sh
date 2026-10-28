#!/usr/bin/env bash
#
# An external-capture final state that matches more than one level is refused (audit 2026-10-25, 03_rmatrix_core F9 /
# checklist item 28; 2026-10-27).
#
# CNuc::ParseExternalCapture attaches a capture pair's final state to the level of the pair's J^pi group whose energy
# lies within 1 keV of the pair's E_x (JGroup::IsLevel, an energy-only match that returns the FIRST fit).  Since
# e239a12 two levels of one J^pi may sit within 1 keV of each other (they carry different level IDs), so the pair's
# E_x can match both: the first one in the file took the ANC silently and the other was ignored.  Now that is an ERROR
# naming the pair, the J^pi, E_x and the number of candidates, and the run writes nothing.
#
# Model: tests/13N with a second 1/2- level 0.8 keV above the ground state (the external-capture final state of its
# gamma pair), with its own ANC line.
#   (a) the plain 13N runs and reproduces expected/chiSquared.out;
#   (b) the duplicate-level model stops with the ERROR and writes no chiSquared.out;
#   (c) the duplicate moved to 2 keV (outside the window) runs: the final state is unambiguous again.
# dev b5faec4: (b) runs, attaching the capture to the first level.
#
#   ./tests/ec_final_state_ambiguous/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_final_state_ambiguous.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# the ground-state lines of 13N (J = 1/2-, E = 0): the proton ANC line and the gamma-pair line
GS="$(awk '/<levels>/ { f = 1; next } /<\/levels>/ { f = 0 } f && NF > 26 && $1 == 0.5 && $2 == -1 && $3 + 0 == 0' "$BASE/13N.azr")"
[ "$(printf '%s\n' "$GS" | wc -l)" -ge 2 ] || { echo "  FAIL  tests/13N/13N.azr no longer has the 1/2- ground-state lines"; exit 1; }
MAXID="$(awk '/<levels>/ { f = 1; next } /<\/levels>/ { f = 0 } f && NF > 26 { if ($9 + 0 > m) m = $9 + 0 } END { print m }' "$BASE/13N.azr")"

make_model() {  # NAME ENERGY-OF-THE-DUPLICATE ("" = none)
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$BASE/data" "$WORK/$1/"
  if [ -n "$2" ]; then
    printf '%s\n' "$GS" | awk -v e="$2" -v id="$((MAXID + 1))" '{ $3 = e; $9 = id; print }' > "$WORK/$1/extra"
    awk -v extra="$WORK/$1/extra" '{ sub(/\r$/, "") } /^<\/levels>$/ { while ((getline l < extra) > 0) print l; print "" } { print }' "$BASE/13N.azr" > "$WORK/$1/t.azr"
  else
    sed 's/\r$//' "$BASE/13N.azr" > "$WORK/$1/t.azr"
  fi
}
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
total() { awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2> /dev/null; }
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$BASE/expected/chiSquared.out")"
near() { awk -v a="$1" -v b="$2" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(a != "" && d < 1e-5) }'; }

make_model plain ""; make_model dup 0.0008; make_model far 0.002
run plain; run dup; run far
t="$(total plain)"
if near "$t" "$want"; then ok "(a) the plain 13N reproduces its expected total ($t)"; else bad "(a) the plain 13N gives '$t', expected $want"; fi
if [ ! -f "$WORK/dup/output/chiSquared.out" ] && grep -q "ERROR: external capture into pair 2: 2 levels of J^pi = 0.5- lie within 1 keV of its final-state energy" "$WORK/dup/run.log"; then ok "(b) a second 1/2- level 0.8 keV from the final state is refused: $(grep -o 'ERROR: external capture into pair 2: 2 levels[^,]*' "$WORK/dup/run.log" | head -1)"; else bad "(b) the duplicate final state was not refused (chiSquared.out present: $([ -f "$WORK/dup/output/chiSquared.out" ] && echo yes || echo no))"; grep -i "error" "$WORK/dup/run.log" | head -2; fi
t="$(total far)"
if [ -n "$t" ] && [ -f "$WORK/far/output/chiSquared.out" ]; then ok "(c) the duplicate moved to 2 keV runs ($t)"; else bad "(c) the duplicate at 2 keV did not run"; grep -i "error" "$WORK/far/run.log" | head -2; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
