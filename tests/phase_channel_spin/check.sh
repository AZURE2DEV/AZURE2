#!/usr/bin/env bash
#
# A phase-shift data set is one (J, l, s) partial wave: its field 's' (the old segment line's "spin <s>"), and an
# ambiguous (J, l) is refused (audit 2026-10-25, tier 3 / checklist t3-23; 2026-10-27).
#
# A phase-shift segment was specified by (J, l).  Up to dev 62c7ed2 the diagonal U-matrix elements of EVERY channel
# spin with that (J, l) were accumulated in one TempTMatrix and half the argument of the sum reported: for a target
# with spin (p + 3He 1P1 / 3P1, p + 7Be J = 1 l = 1 with s = 1 and 2) that is not the phase shift of either wave.  Now
# the phase branch keeps only the data set's channel spin; format 3 requires 's' of a phase shift, and the old
# grammar's segment whose (J, l) had several channel spins without "spin" (EData::Initialize refused it) is refused by
# the converter.
#
# Model: tests/scattering_length/be7p_fit.azr (p + 7Be, s = 1 or 2) plus two 1+ levels with l = 1 channels: level A
# (Ex 1.2) couples only to s = 1, level B (Ex 1.7) only to s = 2 (the other spin's width is 0), so the s = 1 wave is
# unaffected by B and the s = 2 one by A.  Truth for each spin: the same model with the other spin's channel lines
# removed, which is unambiguous and the same on every build.  Phase-shift data set J = 1, l = 1 at lab 0.3 - 2.0 MeV.
#   (a) s = 1 on the full model == the s = 1-only model at every energy (1e-9);  (b) the same for s = 2;
#   (c) the two phases differ from each other by more than 1 degree somewhere (teeth);
#   (d) the data set without 's' is refused (E_FIELD_REQUIRED) and writes no chiSquared.out; the old segment line
#       without "spin" (on the original fixture) is refused by AZURE2 --convert with the old engine's reason;
#   (e) s = 3/2 (no such channel) stops with an ERROR.
# dev 62c7ed2: (a) and (b) fail (the keyword is ignored, both spins summed), (d) and (e) fail (it runs).
#
#   ./tests/phase_channel_spin/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../scattering_length"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/phase_channel_spin.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# level lines (the file stores 2s and 2l): 1+ level A at Ex 1.2 (s = 1 l = 1 width 60 keV, s = 2 width 0),
# level B at Ex 1.7 (s = 1 width 0, s = 2 width 90 keV); pair columns copied from the fixture's p + 7Be lines
TAIL="0.5 1 1.5 -1 0 1.00783 7.01693 1 4 0.1375 0.1375 0 0 0.0 0 4.3 0 0 0"
A1="1 1 1.200000 0 1 1 2 2 4 1 0 60000 $TAIL"
A2="1 1 1.200000 0 1 1 4 2 4 1 0 0 $TAIL"
B1="1 1 1.700000 0 1 1 2 2 5 1 0 0 $TAIL"
B2="1 1 1.700000 0 1 1 4 2 5 1 0 90000 $TAIL"

# make_model NAME S LEVELS...: the fixture with its data sets replaced by one J = 1, l = 1 phase-shift data set (its
# field 's' S, none when S is empty), its extrapolations removed, and the given level lines added
make_model() {
  local name="$1" spin="$2"; shift 2
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks" "$WORK/$name/data"
  awk 'BEGIN { for (i = 0; i <= 34; i++) printf "%.3f 0.0 10.0 0.5\n", 0.3 + i * 0.05 }' > "$WORK/$name/data/p11.dat"
  printf '%s\n' "$@" > "$WORK/$name/levels"
  awk -v levelsfile="$WORK/$name/levels" -v spin="$spin" '
    { sub(/\r$/, "") }
    /^<\/levels>$/ { while ((getline l < levelsfile) > 0) print l }
    /^<data>$/ { print; print "dataset p11\n    reaction       7Be(p,p)7Be\n    quantity       phase_shift\n    unit           deg"
                 print "    J              1\n    l              1"; if (spin != "") print "    s              " spin
                 print "    file           data/p11.dat\n    columns        E_lab[MeV] _ value error"
                 print "    energy_window  0 10 MeV\n    norm           1 fixed"; skip = 1; next }
    /^<extrapolations>$/ { print; skip = 1; next }
    /^<\/data>$|^<\/extrapolations>$/ { skip = 0 }
    !skip { print }' "$BASE/be7p_fit.azr" > "$WORK/$name/t.azr"
}
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
phases() { awk '{ sub(/\r$/, "") } NF >= 5 { printf "%.10f %.10f\n", $1, $4 }' "$WORK/$1/output/AZUREOut_aa=1_R=1.out" 2> /dev/null; }
maxdiff() { paste "$1" "$2" | awk 'BEGIN { m = 0 } { d = $2 - $4; if (d < 0) d = -d; if (d > m) m = d } END { printf "%.3e", m }'; }

make_model full1 1 "$A1" "$A2" "$B1" "$B2"
make_model full2 2 "$A1" "$A2" "$B1" "$B2"
make_model only1 1 "$A1" "$B1"
make_model only2 2 "$A2" "$B2"
make_model nospin "" "$A1" "$A2" "$B1" "$B2"
make_model badspin 3/2 "$A1" "$A2" "$B1" "$B2"
for m in full1 full2 only1 only2 nospin badspin; do run "$m"; phases "$m" > "$WORK/$m.ph"; done
# (d) in the old grammar, through the converter: the original fixture with the old segment line without "spin"
mkdir -p "$WORK/oldnospin/data"
cp "$WORK/nospin/data/p11.dat" "$WORK/oldnospin/data/"
awk -v levelsfile="$WORK/nospin/levels" -v seg="1 1 1 0 10 0 180 2 1 1 1 0 0 0 0 0 data/p11.dat 0 0" '
  { sub(/\r$/, "") }
  /^<\/levels>$/ { while ((getline l < levelsfile) > 0) print l }
  /^<segmentsData>$/ { print; print seg; skip = 1; next }
  /^<segmentsTest>$/ { print; skip = 1; next }
  /^<\/segmentsData>$|^<\/segmentsTest>$/ { skip = 0 }
  !skip { print }' "$HERE/../convert/fixtures/tests/scattering_length/be7p_fit.azr" > "$WORK/oldnospin/t.azr"
(cd "$WORK/oldnospin" && "$AZURE2_BIN" --convert t.azr t3.azr --to 3 > convert.log 2>&1); oldexit=$?
n="$(wc -l < "$WORK/only1.ph" | tr -d ' ')"
if [ "$n" -eq 35 ]; then ok "the single-spin models wrote the 35 phase shifts"; else bad "only1 wrote $n phases, expected 35"; tail -3 "$WORK/only1/run.log"; fi

d1="$(maxdiff "$WORK/full1.ph" "$WORK/only1.ph")"; d2="$(maxdiff "$WORK/full2.ph" "$WORK/only2.ph")"
n1="$(wc -l < "$WORK/full1.ph" | tr -d ' ')"; n2="$(wc -l < "$WORK/full2.ph" | tr -d ' ')"
if [ "$n1" -eq 35 ] && awk -v d="$d1" 'BEGIN { exit !(d <= 1e-9) }'; then ok "(a) s = 1 on the full model == the s = 1 wave alone ($d1 deg)"; else bad "(a) s = 1 differs from the s = 1 wave alone by $d1 deg ($n1 phases; expected <= 1e-9)"; tail -2 "$WORK/full1/run.log"; fi
if [ "$n2" -eq 35 ] && awk -v d="$d2" 'BEGIN { exit !(d <= 1e-9) }'; then ok "(b) s = 2 on the full model == the s = 2 wave alone ($d2 deg)"; else bad "(b) s = 2 differs from the s = 2 wave alone by $d2 deg ($n2 phases; expected <= 1e-9)"; fi
d12="$(maxdiff "$WORK/only1.ph" "$WORK/only2.ph")"
if awk -v d="$d12" 'BEGIN { exit !(d > 1.0) }'; then ok "(c) the two waves differ by up to $d12 deg, so (a) and (b) have teeth"; else bad "(c) the two waves differ by only $d12 deg"; fi
if [ ! -f "$WORK/nospin/output/chiSquared.out" ] && grep -q "E_FIELD_REQUIRED: dataset 'p11': field 's' is missing" "$WORK/nospin/run.log"; then ok "(d) a phase-shift data set without its channel spin is refused (E_FIELD_REQUIRED)"; else bad "(d) the phase-shift data set without 's' was not refused"; grep -i "error" "$WORK/nospin/run.log" | head -2; fi
if [ "$oldexit" -ne 0 ] && [ ! -f "$WORK/oldnospin/t3.azr" ] && grep -q "a phase shift without 'spin s' whose (J, l) has several channel spins, which the old engine refused" "$WORK/oldnospin/convert.log"; then ok "(d) the old ambiguous segment without \"spin\" is refused by the converter with the old engine's reason"; else bad "(d) the old segment without \"spin\" was converted (exit $oldexit)"; head -3 "$WORK/oldnospin/convert.log"; fi
if [ ! -f "$WORK/badspin/output/chiSquared.out" ] && grep -q "ERROR: phase-shift segment #1 (J = 1, l = 1) asks for channel spin s = 1.5" "$WORK/badspin/run.log"; then ok "(e) \"s 3/2\" (no such channel) is refused"; else bad "(e) \"s 3/2\" was not refused"; grep -i "error" "$WORK/badspin/run.log" | head -2; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
