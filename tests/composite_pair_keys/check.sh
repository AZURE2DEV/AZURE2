#!/usr/bin/env bash
#
# A sum of reactions in a project with two pairs of the same particles: examples/c12c12_tumino2018 has 12C + 12C as
# pair 1 (THM) and pair 6 (direct data), alike in particles and states, so the label 12C(12C,p)23Na names both and the
# pair keys of <levels> say which (docs/upgrade/project_format3.md section 3.5).  On a sum or a ratio entrance_pair and
# exit_pair apply to every term whose label names two pairs.  Checked, each extrapolation in a run of its own (one
# output file per entrance and exit pair):
#   1. format 3: sigma of 12C(12C,p)23Na[0.4402] + 0.5*12C(12C,a)20Ne[1.6337] with entrance_pair 6 equals, energy by
#      energy, sigma(p1) + 0.5 sigma(a1) of the two single reactions with entrance_pair 6 (to 1e-10 relative: the
#      output files print 11 digits), and
#      the same sum with entrance_pair 1 differs from it (the key reaches every term);
#   2. the old format: the same sum as an old composite test line (all its terms entering by pair 6) converts with no
#      question to the record of 1 with entrance_pair 6, and its run gives the same file as 1;
#   3. the old line with its terms entering by pairs 6 and 1 is refused by the converter (no format-3 form), naming
#      both pairs.
#
#   ./tests/composite_pair_keys/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/edit_guard.sh"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/composite_pair_keys.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

BASE="$ROOT/examples/c12c12_tumino2018/c12c12_tumino2018.azr"
OLD="$ROOT/tests/convert/fixtures/examples/c12c12_tumino2018/c12c12_tumino2018.azr"
SUM="12C(12C,p)23Na[0.4402] + 0.5*12C(12C,a)20Ne[1.6337]"

# only_extrapolation REACTION KEY FILE: the project with one extrapolation record (sigma, 2-5 MeV lab in 0.1 MeV)
only_extrapolation() {
  awk -v reaction="$1" -v key="$2" '
    /^<extrapolations>/ { print; print "extrapolation x"; print "    reaction       " reaction
                          print "    entrance_pair  " key; print "    quantity       sigma"; print "    unit           b"
                          print "    energies       range 2 5 MeV lab step 0.1"; skip = 1; next }
    /^<\/extrapolations>/ { skip = 0 }
    !skip { print }' "$3"
}
# old_test_line LINE FILE: the old project whose <segmentsTest> holds eight inactive lines and LINE as its ninth (the
# old <targetInt> lines 1 to 8 also fold the test lines of their indices; the ninth is not folded, as the records of 1)
old_test_line() {
  awk -v line="$1" '/^<segmentsTest>/ { print; for (i = 0; i < 8; i++) print "0  6  3  2  5  0.1  0  0  0  0  0"
                                        print line; skip = 1; next }
                    /^<\/segmentsTest>/ { skip = 0 } !skip { print }' "$2"
}
# run NAME PROJECT: extrapolate PROJECT in $WORK/NAME (the example's data beside it)
run() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  ln -s "$ROOT/examples/c12c12_tumino2018/data" "$WORK/$1/data"
  cp "$2" "$WORK/$1/m.azr"
  (cd "$WORK/$1" && $(guard_command 600) "$AZURE2_BIN" extrapolate m.azr > run.log 2>&1)
}

# ---- 1. format 3: the sum against its terms, and the key reaching every term
only_extrapolation "$SUM" 6 "$BASE" > "$WORK/sum6.azr" && must_differ "$BASE" "$WORK/sum6.azr" || bad "the sum's project"
only_extrapolation "12C(12C,p)23Na[0.4402]" 6 "$BASE" > "$WORK/p6.azr"
only_extrapolation "12C(12C,a)20Ne[1.6337]" 6 "$BASE" > "$WORK/a6.azr"
only_extrapolation "$SUM" 1 "$BASE" > "$WORK/sum1.azr"
for n in sum6 p6 a6 sum1; do run "$n" "$WORK/$n.azr"; done
s6="$WORK/sum6/output/AZUREOut_aa=6_R=3.extrap" p6="$WORK/p6/output/AZUREOut_aa=6_R=3.extrap"
a6="$WORK/a6/output/AZUREOut_aa=6_R=5.extrap" s1="$WORK/sum1/output/AZUREOut_aa=1_R=3.extrap"
if [ ! -s "$s6" ] || [ ! -s "$p6" ] || [ ! -s "$a6" ] || [ ! -s "$s1" ]; then
  bad "no extrapolation written ($(grep -h -m1 ERROR "$WORK"/*/run.log))"
else
  # column 4: sigma (b) at the energies of column 1; every row of the three files at the same energy
  if worst="$(paste "$s6" "$p6" "$a6" | awk -v re='^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$' '
        NF == 0 { next }
        { n++; if ($1 != $6 || $1 != $11 || $4 !~ re || $9 !~ re || $14 !~ re) { bad = 1; next }
          want = $9 + 0.5 * $14; d = $4 - want; if (d < 0) d = -d; s = want < 0 ? -want : want
          if (s > 0 && d / s > w) w = d / s; if ($4 + 0 > 0) positive++ }
        END { if (bad || n < 31 || positive < n) { print "layout or values wrong"; exit 1 } printf "%.1e", w
              exit !(w <= 1e-10) }')"; then
    ok "sum with entrance_pair 6 = sigma(p1) + 0.5 sigma(a1) of the single reactions, 31 energies, within $worst"
  else
    bad "sum with entrance_pair 6 against its terms: $worst"
  fi
  if differ="$(paste "$s6" "$s1" | awk 'NF == 0 { next } { d = $4 - $9; if (d < 0) d = -d; if (d > 1e-6 * ($4 < 0 ? -$4 : $4)) n++ }
                                         END { print n + 0; exit !(n > 20) }')"; then
    ok "the same sum with entrance_pair 1 differs at $differ of 31 energies (the key reaches every term)"
  else
    bad "the sum with entrance_pair 1 equals the one with 6 at most energies ($differ differ): the key does not bite"
  fi
fi

# ---- 2. the old format: an old composite test line whose terms all enter by pair 6
old_test_line "1  6  3  2  5  0.1  0  0  0  0  1  0  -1  1  6  5  -999  0.5" "$OLD" > "$WORK/old6.azr"
mkdir -p "$WORK/conv"
if ! (cd "$WORK/conv" && "$AZURE2_BIN" --convert "$WORK/old6.azr" "$WORK/conv/new6.azr" > convert.log 2>&1); then
  bad "the old composite line entering by pair 6 is not converted: $(grep -m1 -i "error\|refus" "$WORK/conv/convert.log")"
elif grep -q "question" <<< "$(grep -i "[1-9][0-9]* question" "$WORK/conv/convert.log")"; then
  bad "the conversion asks: $(cat "$WORK/conv/convert.log")"
else
  rec="$(awk '/^extrapolation x9$/ { on = 1 } on && /^$/ { exit } on' "$WORK/conv/new6.azr")"
  if grep -qF "    reaction       $SUM" <<< "$rec" && grep -q "^    entrance_pair  6$" <<< "$rec" &&
     ! grep -q "^    beam " <<< "$rec"; then
    ok "old composite line converted with no question: the sum with one entrance_pair 6"
  else
    bad "the converted record: $rec"
  fi
  run old6 "$WORK/conv/new6.azr"
  o6="$WORK/old6/output/AZUREOut_aa=6_R=3.extrap"
  if [ -s "$o6" ] && cmp -s "$o6" "$s6"; then
    ok "the converted old line gives the format-3 sum's file byte for byte"
  else
    bad "the converted old line's run differs from the format-3 sum ($(grep -m1 ERROR "$WORK/old6/run.log"))"
  fi
fi

# ---- 3. the old line whose terms enter by pairs 6 and 1: no format-3 form, refused
old_test_line "1  6  3  2  5  0.1  0  0  0  0  1  0  -1  1  1  5  -999  0.5" "$OLD" > "$WORK/old61.azr"
if (cd "$WORK/conv" && "$AZURE2_BIN" --convert "$WORK/old61.azr" "$WORK/conv/new61.azr" > convert61.log 2>&1); then
  bad "the old line entering by pairs 6 and 1 was converted"
elif grep -q "entrance pairs 1 and 6" "$WORK/conv/convert61.log" && [ ! -e "$WORK/conv/new61.azr" ]; then
  ok "the old line entering by pairs 6 and 1: refused, naming both pairs, nothing written"
else
  bad "the refusal of the old line entering by pairs 6 and 1: $(cat "$WORK/conv/convert61.log")"
fi

if [ "$failures" -eq 0 ]; then echo "composite_pair_keys: all checks passed"; else echo "composite_pair_keys: $failures FAILED"; fi
[ "$failures" -eq 0 ]
