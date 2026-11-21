#!/usr/bin/env bash
#
# Capture from the interior of a hybrid potential through the whole engine (docs/upgrade/interior_potential.md
# section 6; the radial integrals themselves: ctest upgrade_interior_capture).  models/direct.azr: alpha + 12C with the
# Woods-Saxon of Koonin, Tombrello and Fox in every l and interior_boundary 0, E1 capture into the potential's own 0+
# bound state (2.427561 MeV below the threshold, given as a level with an ANC of 100 fm^-1/2), no level in 1-.  The
# cross section is then the potential model's direct capture, which knows no channel radius:
#
#   (a) with the channel radius 4.0 and 6.5 fm in place of 5.35 fm the cross section at 51 energies across the
#       potential's 1- resonance (E_cm 0.75 to 4.5 MeV) is the same to 5e-4 (the direct integrals' tolerance is 1e-4);
#       without the interior term it differs by more than a factor 1.1 somewhere (the check can fail);
#   (b) --no-ec-table gives the same numbers to 1e-9: these channels are integrated directly either way (ECTable.cpp);
#   (c) with the potential's record restricted to l = 1 the final channel (l = 0) has no interior: the run warns that
#       the final state gets no capture from the potential's interior, and the full model does not;
#       a final state on the potential's bound state, the intended use, draws a note instead of the warning that a
#       level describes a state of the potential twice;
#   (d) interior_charge_e1 changes the cross section, and is refused without interior_boundary
#       (E_POT_INTERIOR_CHARGE).
#
# An unbound final state (section 6.5, step 4).  models/unbound.azr: E1 capture from l = 2 and 4 into the potential's
# own l = 3 resonance (E_cm = 3.261 MeV, width 173 keV at 5.35 fm), given as a 3- level with that width, a narrow
# resonance whose outgoing wave continues inside the radius as the potential's regular solution:
#
#   (e) the cross section at 21 energies (E_cm 3.75 to 7.5 MeV) at 4.0 fm is that at 5.35 fm to 5 %, at 6.5 fm to
#       20 % (the level stays where it is while the potential's state moves by 15 keV with the radius, and the state
#       is 5 % of its energy wide: the narrow-resonance form is that good, no better); with the potential's record
#       restricted to l = 2, so that the final channel has no interior, the radii differ by more than a factor 1.5,
#       and that run is warned about while the full one is not;
#   (f) `line_window` on such a final state is refused: the line-shape integration is not built for a final channel
#       with a potential inside the radius.
#
# The rate mode:
#
#   (g) the rate of models/direct.azr at T9 = 0.5 to 10 (temperatures.dat) equals N_A <sigma v> integrated by
#       tests/reaction_rate/truth.py from a dense extrapolation of the same model (0.4 % steps, E_cm 0.05 to 30 MeV)
#       to 1e-4.  The rate's energy grid has a cluster on each resonance of the entrance pair's potential
#       (ReactionRate.cpp, AddPotentialResonances); without it the 412 keV wide 1- resonance of the potential sat
#       between backbone nodes and the rate was 0.2 % low from T9 = 3.
#
#   ./tests/interior_capture/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/interior_capture.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
OUT=output/AZUREOut_aa=1_R=2.extrap

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# model NAME SED-EXPRESSION [FLAGS]: the project with one substitution, run as an extrapolation
model() {
  local dir="$WORK/$1" flags="${3:-}"
  mkdir -p "$dir/output" "$dir/checks"
  sed "$2" "$HERE/models/direct.azr" > "$dir/t.azr"
  # shellcheck disable=SC2086
  (cd "$dir" && run_limited 600 "$AZURE2_BIN" extrapolate --use-brune $flags t.azr > run.log 2>&1)
}
# worst NAME1 NAME2: the largest relative difference of the two cross sections (column 4), "none" without rows
worst() {
  if [ ! -s "$WORK/$1/$OUT" ] || [ ! -s "$WORK/$2/$OUT" ]; then echo none; return; fi
  paste "$WORK/$1/$OUT" "$WORK/$2/$OUT" |
    awk 'NF >= 8 && $4 > 0 { n++; d = $(NF/2+4) / $4 - 1; if (d < 0) d = -d; if (d > w) w = d }
         END { if (n == 51) printf "%.2e", w; else print "none" }'
}
below() { awk -v d="$1" -v lim="$2" 'BEGIN { exit !(d ~ /^ *[-+]?\.?[0-9]/ && d < lim) }'; }
above() { awk -v d="$1" -v lim="$2" 'BEGIN { exit !(d ~ /^ *[-+]?\.?[0-9]/ && d > lim) }'; }

echo "interior_capture: $AZURE2_BIN"
model a535 's/^$//'
model a400 's/ 5\.350 0 0 0$/ 4.000 0 0 0/'
model a650 's/ 5\.350 0 0 0$/ 6.500 0 0 0/'
model direct535 's/^$//' --no-ec-table
# the final channel without the interior: the record for l = 1 alone
model l1_535 's/^    pair            a+12C$/    pair            a+12C\
    l               1/'
model l1_400 's/ 5\.350 0 0 0$/ 4.000 0 0 0/;s/^    pair            a+12C$/    pair            a+12C\
    l               1/'
model charge 's/^    interior_boundary 0$/    interior_boundary 0\
    interior_charge_e1 0.01/'

for r in a400 a650; do
  d=$(worst a535 $r)
  if below "$d" 5e-4; then ok "(a) channel radius ${r#a} (x 0.01 fm) against 5.35 fm: the cross section agrees to $d"
  else bad "(a) channel radius ${r#a} (x 0.01 fm) against 5.35 fm: the cross section differs by $d (limit 5e-4)"; fi
done
d=$(worst l1_535 l1_400)
if above "$d" 0.1; then ok "(a) without the final state's interior the two radii differ by $d (the check can fail)"
else bad "(a) without the final state's interior the two radii differ by only $d: the check proves nothing"; fi
d=$(worst a535 direct535)
if below "$d" 1e-9; then ok "(b) --no-ec-table: the same cross section ($d)"
else bad "(b) --no-ec-table differs by $d (limit 1e-9)"; fi

if grep -q "gets no capture from the potential's interior" "$WORK/l1_535/run.log"; then ok "(c) the final channel without an interior boundary constant is warned about"
else bad "(c) no warning for the final channel without an interior boundary constant"; fi
if grep -q "gets no capture from the potential's interior" "$WORK/a535/run.log"; then bad "(c) the full model is warned about"
else ok "(c) no such warning for the full model"; fi
if grep -q "NOTE: the external-capture final state .* lies on a bound state of the potential" "$WORK/a535/run.log" &&
   ! grep -q "J^pi=0+.* lies within one width of a bound state" "$WORK/a535/run.log"; then
  ok "(c) the final state on the potential's bound state draws a note, not the warning of a level that doubles it"
else bad "(c) the final state on the potential's bound state: no note, or the doubled-level warning"; fi

d=$(worst a535 charge)
if above "$d" 1.0; then ok "(d) interior_charge_e1 0.01 changes the cross section (relative change up to $d)"
else bad "(d) interior_charge_e1 0.01 leaves the cross section within $d"; fi
mkdir -p "$WORK/refused"
sed 's/^    interior_boundary 0$/    interior_charge_e1 0.01/' "$HERE/models/direct.azr" > "$WORK/refused/t.azr"
(cd "$WORK/refused" && run_limited 60 "$AZURE2_BIN" --validate t.azr > run.log 2>&1)
if grep -q "E_POT_INTERIOR_CHARGE" "$WORK/refused/run.log"; then ok "(d) interior_charge_e1 without interior_boundary is refused"
else bad "(d) interior_charge_e1 without interior_boundary is not refused"; fi

# (e), (f): the unbound final state
umodel() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  sed "$2" "$HERE/models/unbound.azr" > "$dir/t.azr"
  (cd "$dir" && run_limited 600 "$AZURE2_BIN" extrapolate --use-brune t.azr > run.log 2>&1)
}
# span NAME1 NAME2: the largest of sigma1/sigma2 and sigma2/sigma1 less 1 over the 21 energies, "none" without them
span() {
  if [ ! -s "$WORK/$1/$OUT" ] || [ ! -s "$WORK/$2/$OUT" ]; then echo none; return; fi
  paste "$WORK/$1/$OUT" "$WORK/$2/$OUT" |
    awk 'NF >= 8 && $4 > 0 && $(NF/2+4) > 0 { n++; d = $(NF/2+4) / $4; if (d < 1) d = 1 / d; d -= 1; if (d > w) w = d }
         END { if (n == 21) printf "%.2e", w; else print "none" }'
}
L2='s/^    pair            a+12C$/    pair            a+12C\
    l               2/'
umodel u535 's/^$//'
umodel u400 's/ 5\.350 0 0 0$/ 4.000 0 0 0/'
umodel u650 's/ 5\.350 0 0 0$/ 6.500 0 0 0/'
umodel n535 "$L2"
umodel n400 "$L2"';s/ 5\.350 0 0 0$/ 4.000 0 0 0/'
d=$(span u535 u400)
if below "$d" 0.05; then ok "(e) unbound final state, channel radius 4.0 against 5.35 fm: the cross section agrees to $d"
else bad "(e) unbound final state, channel radius 4.0 against 5.35 fm: the cross section differs by $d (limit 0.05)"; fi
d=$(span u535 u650)
if below "$d" 0.20; then ok "(e) unbound final state, channel radius 6.5 against 5.35 fm: the cross section agrees to $d"
else bad "(e) unbound final state, channel radius 6.5 against 5.35 fm: the cross section differs by $d (limit 0.20)"; fi
d=$(span n535 n400)
if above "$d" 0.5; then ok "(e) without the final channel's interior the radii 4.0 and 5.35 fm differ by $d (the check can fail)"
else bad "(e) without the final channel's interior the radii differ by only $d: the check proves nothing"; fi
if grep -q "is UNBOUND" "$WORK/u535/run.log" && ! grep -q "gets no capture from the potential's interior" "$WORK/u535/run.log" &&
   grep -q "gets no capture from the potential's interior" "$WORK/n535/run.log"; then
  ok "(e) the unbound final state is noted; only the run without its channel's interior is warned about"
else bad "(e) the unbound final state: no note, or the interior warning on the wrong run"; fi
mkdir -p "$WORK/window/output" "$WORK/window/checks"
sed 's/^    quantity       sigma$/    quantity       sigma\
    line_window    200 keV/' "$HERE/models/unbound.azr" > "$WORK/window/t.azr"
(cd "$WORK/window" && run_limited 600 "$AZURE2_BIN" extrapolate --use-brune t.azr > run.log 2>&1)
status=$?
if [ "$status" -ne 0 ] && grep -q "line window .* potential with interior_boundary" "$WORK/window/run.log" && [ ! -s "$WORK/window/$OUT" ]; then
  ok "(f) line_window on a final channel with an interior potential is refused"
else bad "(f) line_window on a final channel with an interior potential is not refused (exit status $status)"; fi

# (g): the rate against the integral of a dense extrapolation
if command -v python3 > /dev/null 2>&1; then
  mkdir -p "$WORK/rate/output" "$WORK/rate/checks" "$WORK/dense/output" "$WORK/dense/checks"
  cp "$HERE/models/direct.azr" "$WORK/rate/t.azr"
  cp "$HERE/temperatures.dat" "$WORK/rate/"
  energies="$(awk 'BEGIN { for (e = 0.05; e < 30; e *= 1.004) printf " %.9g", e * 16.0026 / 12 }')"
  sed "s/^    energies       range 1\.0 6\.0 MeV lab step 0\.1$/    energies       list$energies MeV lab/" \
    "$HERE/models/direct.azr" > "$WORK/dense/t.azr"
  (cd "$WORK/dense" && run_limited 900 "$AZURE2_BIN" extrapolate --use-brune t.azr > run.log 2>&1)
  (cd "$WORK/rate" && run_limited 600 "$AZURE2_BIN" rate --use-brune --rate-entrance=1 --rate-exit=2 \
     --temperature-file=temperatures.dat t.azr > run.log 2>&1)
  python3 "$HERE/../reaction_rate/truth.py" 4.0026 12.0 "$HERE/temperatures.dat" "$WORK/dense/$OUT" > "$WORK/truth.dat" 2> /dev/null
  d="$(paste <(awk 'NR > 1 && NF >= 2 { print $2 }' "$WORK/rate/output/reactionrates.out" 2> /dev/null) \
             <(awk 'NR > 1 && NF >= 2 { print $2 }' "$WORK/truth.dat") |
       awk 'NF == 2 && $2 > 0 { n++; d = $1 / $2 - 1; if (d < 0) d = -d; if (d > w) w = d }
            END { if (n == 6) printf "%.2e", w; else print "none" }')"
  if below "$d" 1e-4; then ok "(g) the rate at six temperatures equals the integral of the dense extrapolation to $d"
  else bad "(g) the rate differs from the integral of the dense extrapolation by $d (limit 1e-4)"; fi
else
  echo "  skip  (g) the rate against the dense extrapolation: no python3"
fi

if [ "$failures" -eq 0 ]; then echo "interior_capture: all checks passed"; else echo "interior_capture: $failures check(s) FAILED"; fi
exit "$failures"
