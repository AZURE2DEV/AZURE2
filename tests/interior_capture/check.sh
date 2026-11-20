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
#   (d) interior_charge_e1 changes the cross section, and is refused without interior_boundary
#       (E_POT_INTERIOR_CHARGE).
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

d=$(worst a535 charge)
if above "$d" 1.0; then ok "(d) interior_charge_e1 0.01 changes the cross section (relative change up to $d)"
else bad "(d) interior_charge_e1 0.01 leaves the cross section within $d"; fi
mkdir -p "$WORK/refused"
sed 's/^    interior_boundary 0$/    interior_charge_e1 0.01/' "$HERE/models/direct.azr" > "$WORK/refused/t.azr"
(cd "$WORK/refused" && run_limited 60 "$AZURE2_BIN" --validate t.azr > run.log 2>&1)
if grep -q "E_POT_INTERIOR_CHARGE" "$WORK/refused/run.log"; then ok "(d) interior_charge_e1 without interior_boundary is refused"
else bad "(d) interior_charge_e1 without interior_boundary is not refused"; fi

if [ "$failures" -eq 0 ]; then echo "interior_capture: all checks passed"; else echo "interior_capture: $failures check(s) FAILED"; fi
exit "$failures"
