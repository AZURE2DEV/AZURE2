#!/usr/bin/env bash
#
# An external-capture final state that matches more than one level is refused (audit 2026-10-25, 03_rmatrix_core F9 /
# checklist item 28; 2026-10-27).
#
# CompoundSystem::ParseExternalCapture attaches a capture pair's final state to the level of the pair's J^pi group whose energy
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
# The other way round, two photon pairs with one final state (2026-11-11): tests/13N with a copy of its gamma pair under
# a new key (its ground-state line) and a copy of its first capture data set into that pair.  The second pair took the
# level's capture parameters, and capture into the first ran without its external part and without a message
# (1.95e-9 against 4.32e-9 at the first point).  Now an ERROR names both pair keys and the level.
#   (d) the two-pair model (format 3, each capture data set naming its exit_pair) stops with that ERROR and writes no
#       chiSquared.out;
#   (e) so does the conversion of the same model written in the old grammar (AZURE2 --convert ... --to 3 of the
#       original 13N with the lines added; the converter writes it, the run refuses it).
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
# the two-pair model: the ground state's gamma-pair line again with a new pair key, and the first data set into the
# gamma pair again into the new one; the reaction label 12C(p,g)13N names both pairs now, so each capture data set says
# its exit pair (exit_pair <key>, project_format3.md section 3.5)
make_pairs() {  # NAME
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$BASE/data" "$WORK/$1/"
  awk '{ sub(/\r$/, "") }
       /<levels>/ { lv = 1 } /<\/levels>/ { lv = 0 }
       lv && NF > 26 { if ($6 + 0 > maxkey) maxkey = $6 + 0 }
       lv && NF > 26 && $1 == 0.5 && $2 == -1 && $3 + 0 == 0 && $27 == 10 { gamma = $6; line = $0 }
       { lines[++n] = $0; inlev[n] = lv }
       END {
         key = maxkey + 1; first = ""
         for (i = 1; i <= n; i++) {
           l = lines[i]
           if (l ~ /^dataset /) { rec = l; incap = 0 }
           if (l ~ /^[a-z_]+ [A-Za-z0-9_-]+$/ && l !~ /^dataset /) rec = ""
           if (rec != "" && l == "    reaction       12C(p,g)13N") {
             incap = 1; l = l "\n    exit_pair      " gamma
             if (first == "") { first = rec; copy = rec "_pair" key "\n" lines[i] "\n    exit_pair      " key; collecting = 1 }
           }
           if (collecting && i > 0 && lines[i] ~ /^    / && lines[i] !~ /reaction/ && rec == first) copy = copy "\n" lines[i]
           if (collecting && rec == first && (i == n || lines[i + 1] !~ /^    /)) { collecting = 0; l = l "\n" copy }
           print l
           if (lines[i] == line) { $0 = line; $6 = key; print }
         }
       }' "$BASE/13N.azr" > "$WORK/$1/t.azr"
}
# the same model in the old grammar (from the original of tests/13N, tests/convert/fixtures), for (e)
make_pairs_old() {  # NAME
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$BASE/data" "$WORK/$1/"
  awk '{ sub(/\r$/, "") }
       /<levels>/ { lv = 1 } /<\/levels>/ { lv = 0 } /<segmentsData>/ { sg = 1 } /<\/segmentsData>/ { sg = 0 }
       lv && NF > 26 { if ($6 + 0 > maxkey) maxkey = $6 + 0 }
       lv && NF > 26 && $1 == 0.5 && $2 == -1 && $3 + 0 == 0 && $27 == 10 { gamma = $6; line = $0 }
       { lines[++n] = $0; inlev[n] = lv; inseg[n] = sg }
       END {
         key = maxkey + 1; done = 0
         for (i = 1; i <= n; i++) {
           print lines[i]
           split(lines[i], f)
           if (lines[i] == line) { $0 = line; $6 = key; print }
           if (inseg[i] && !done && f[1] == 1 && f[3] == gamma) { $0 = lines[i]; $3 = key; print; done = 1 }
         }
       }' "$HERE/../convert/fixtures/tests/13N/13N.azr" > "$WORK/$1/t.azr"
}
run() { (cd "$WORK/$1" && run_limited 600 "$AZURE2_BIN" calculate --no-gui t.azr > run.log 2>&1); }
total() { awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2> /dev/null; }
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$BASE/expected/chiSquared.out")"
near() { awk -v a="$1" -v b="$2" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(a != "" && d < 1e-5) }'; }

make_model plain ""; make_model dup 0.0008; make_model far 0.002
run plain; run dup; run far
t="$(total plain)"
if near "$t" "$want"; then ok "(a) the plain 13N reproduces its expected total ($t)"; else bad "(a) the plain 13N gives '$t', expected $want"; fi
if [ ! -f "$WORK/dup/output/chiSquared.out" ] && grep -q "ERROR: external capture into pair 2: 2 levels of J^pi = 0.5- lie within 1 keV of its final-state energy" "$WORK/dup/run.log"; then ok "(b) a second 1/2- level 0.8 keV from the final state is refused: $(grep -o 'ERROR: external capture into pair 2: 2 levels[^,]*' "$WORK/dup/run.log" | sed -n 1p)"; else bad "(b) the duplicate final state was not refused (chiSquared.out present: $([ -f "$WORK/dup/output/chiSquared.out" ] && echo yes || echo no))"; grep -i "error" "$WORK/dup/run.log" | sed -n '1,2p'; fi
t="$(total far)"
if [ -n "$t" ] && [ -f "$WORK/far/output/chiSquared.out" ]; then ok "(c) the duplicate moved to 2 keV runs ($t)"; else bad "(c) the duplicate at 2 keV did not run"; grep -i "error" "$WORK/far/run.log" | sed -n '1,2p'; fi

make_pairs pairs; run pairs
PAIRS="ERROR: external capture into pairs 2 and 3: both have the final state J^pi = 0.5- at E_x = 0 MeV"
if [ ! -f "$WORK/pairs/output/chiSquared.out" ] && grep -q "$PAIRS" "$WORK/pairs/run.log"; then ok "(d) two photon pairs with one final state are refused: $(grep -o 'ERROR: external capture into pairs[^(]*' "$WORK/pairs/run.log" | sed -n 1p)"; else bad "(d) two photon pairs with one final state were not refused (chiSquared.out present: $([ -f "$WORK/pairs/output/chiSquared.out" ] && echo yes || echo no))"; grep -i "error" "$WORK/pairs/run.log" | sed -n '1,2p'; fi
make_pairs_old pairs_old
mkdir -p "$WORK/pairs3/output" "$WORK/pairs3/checks"; cp -r "$BASE/data" "$WORK/pairs3/"
(cd "$WORK/pairs_old" && run_limited 300 "$AZURE2_BIN" --convert t.azr "$WORK/pairs3/t.azr" --to 3 > convert.log 2>&1)
if [ -f "$WORK/pairs3/t.azr" ]; then
  run pairs3
  if [ ! -f "$WORK/pairs3/output/chiSquared.out" ] && grep -q "$PAIRS" "$WORK/pairs3/run.log"; then ok "(e) its format-3 conversion is refused the same way"; else bad "(e) the format-3 conversion of the two-pair model was not refused"; grep -i "error" "$WORK/pairs3/run.log" | sed -n '1,2p'; fi
else
  bad "(e) the two-pair model did not convert to format 3"; tail -3 "$WORK/pairs_old/convert.log"
fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
