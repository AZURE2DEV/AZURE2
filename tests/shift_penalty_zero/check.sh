#!/usr/bin/env bash
#
# Energy-shift penalties are measured from the prior's centre (ESegment::AddShiftPenalty, since 2026-10-25; in format 3
# the centre is written, "fit prior 0 <sigma>", as the converter writes the old lines): a free constant shift dE,
# sqrt(E) coefficient b and shift-function parameter a_j add (dE/sigma)^2, (b/sigma)^2, (a_j/sigma)^2; the value is
# only the starting value.  Model: tests/shift_fn/shift_fn.azr with two data sets whose free shift terms START away
# from zero:
#   xs_40: offset 0.004 fit prior 0 0.002, sqrt 0.003 fit prior 0 0.001                     -> 4 + 9
#   xs_90: offset 0.001 fit prior 0 0.002, energy_scale "a0*(x-1)", a0 0.002 fit prior 0 0.004 -> 0.25 + 0.25
#   (a) calculate with data: Total-Shift-Chi-Squared = 13.5 (the build before 2026-10-25 gives 0, the penalty
#       being measured from the start values);
#   (b) --use-gradient with AZURE_GRAD_CHECK=exit: every free shift component of the analytic gradient (penalty
#       rows included) agrees with a finite difference of the whole chi2;
#   (c) a fixed shift with an uncertainty adds nothing (only free terms are penalised).
#
#   ./tests/shift_penalty_zero/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
SRC="$HERE/../shift_fn"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/shift_penalty_zero.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
field() { awk -v f="$2" '/^Total-Chi-Squared:/ { for (i = 1; i < NF; i++) if ($i == f) print $(i+1) }' "$1/output/chiSquared.out" 2>/dev/null; }

# ds NAME FILE ANGLE [FIELD...]: a data set of the model (p + p elastic at ANGLE degrees) with these extra fields
ds() {
  local name="$1" file="$2" angle="$3"; shift 3
  printf 'dataset %s\n    reaction       1H(p,p)1H\n    quantity       dsigma_domega\n    unit           b/sr\n' "$name"
  printf '    file           data/%s\n    columns        E_lab[MeV] theta_cm[deg] value error\n' "$file"
  printf '    energy_window  0 100 MeV\n    angle_window   %s %s deg\n    norm           1 fixed\n' "$angle" "$angle"
  local f; for f in "$@"; do printf '    %s\n' "$f"; done
  echo
}
mk() {  # NAME RECORDS : shift_fn.azr with RECORDS as its <data> section
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -r "$SRC/data" "$dir/"
  # the records go through the environment: BSD awk (macOS) refuses a newline in -v ("newline in string")
  RECORDS="$2" awk '
    /^<data>/ { print; printf "%s\n", ENVIRON["RECORDS"]; skip = 1; next }
    /^<\/data>/ { skip = 0 }
    !skip { print }' "$SRC/shift_fn.azr" > "$dir/t.azr"
}

echo "shift_penalty_zero: $AZURE2_BIN"
mk calc "$(printf 'energy_scale xs_90_scale\n    function       "a0*(x-1)" MeV\n    param          a0 0.002 fit prior 0 0.004\n\n'
           ds xs_40 xs_40.dat 40 'shift          offset 0.004 MeV fit prior 0 0.002' 'shift          sqrt 0.003 MeV^1/2 fit prior 0 0.001'
           ds xs_90 xs_90.dat 90 'shift          offset 0.001 MeV fit prior 0 0.002' 'shift          xs_90_scale')"
(cd "$WORK/calc" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
s="$(field "$WORK/calc" Total-Shift-Chi-Squared:)"
if [ -n "$s" ] && awk -v s="$s" 'BEGIN { d = s - 13.5; if (d < 0) d = -d; exit !(d < 1e-6) }'; then
  ok "(a) penalties from zero: Total-Shift-Chi-Squared = $s (4 + 9 + 0.25 + 0.25)"
else
  bad "(a) Total-Shift-Chi-Squared '$s', expected 13.5 (0 = measured from the start values)"; tail -5 "$WORK/calc/run.log"
fi

cp -r "$WORK/calc" "$WORK/grad"
(cd "$WORK/grad" && printf '2\n\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > jacobian.log 2>&1)
line="$(grep -m1 'Jacobian path:' "$WORK/grad/jacobian.log")"
if grep -q "shift 0/2" <<< "$line" && grep -q "shift_fn 0/1" <<< "$line"; then
  ok "(b) analytic gradient = finite differences for dE (x2), b, a0 ($(grep -o 'worst relative difference [-0-9.e+]*' "$WORK/grad/jacobian.log" | sed -n 1p))"
else
  bad "(b) gradient check: '${line:-no Jacobian-path line}' $(grep 'components with relative' "$WORK/grad/jacobian.log" | sed -n '1,2p')"
fi

mk fixed "$(ds xs_40 xs_40.dat 40 'shift          offset 0.004 MeV fixed prior 0 0.002')"
(cd "$WORK/fixed" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
s="$(field "$WORK/fixed" Total-Shift-Chi-Squared:)"   # not written when no shift term is free
if [ -z "$(field "$WORK/fixed" Total-Chi-Squared:)" ]; then bad "(c) fixed shift: no chiSquared.out"
elif [ -z "$s" ] || awk -v s="$s" 'BEGIN { exit !(s + 0 == 0) }'; then ok "(c) a fixed shift adds no penalty (${s:-none reported})"
else bad "(c) fixed shift: Total-Shift-Chi-Squared '$s', expected none or 0"; fi

if [ "$failures" -eq 0 ]; then echo "shift_penalty_zero: PASS"; else echo "shift_penalty_zero: FAIL ($failures)"; exit 1; fi
