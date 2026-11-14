#!/usr/bin/env bash
#
# An energy shift is the true beam energy (include/ShiftedFold.h, the rule; README.md of this directory): the fold of
# a point measured at E_n with the shift s equals the fold of the same point without a shift at E_n + s(E_n).
#   1. models/identity.azr: 36 cases (thin target sigma, layer yield, thick-target yield) x (a formula stopping with
#      d ln eps/d ln E = -1.5, SRIM's evaluation, a table) x (no straggling, Bohr) x (no beam spread, 0.8 keV), each
#      with an offset (varied, with a prior), a sqrt or a function shift, against its twin without a shift at the
#      shifted energies in the same project (the same mesh): equal to 1e-9 (measured: equal in all 11 printed digits);
#   2. the same pairs in two projects (models/shifted.azr, models/twins.azr: two meshes, the nodes not shared): equal
#      to 3e-6, three times the 1.04e-6 measured (the beam spread's Gaussian on two different depth grids; 0 without a
#      spread on a layer, 6.5e-8 on the thick target);
#   3. models/fit.azr: the gradient of a fitted offset on a 72 keV layer with eps ~ E^-1.5 (the back edge's position
#      depends on the true energy) against a central difference of chi^2 (AZURE_GRAD_CHECK on the Jacobian path):
#      relative difference below 1e-5 (measured 3.9e-7; the existing gradient tests flag above 1e-2).
#
#   ./tests/shift_true_energy/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/shift_true_energy.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {  # NAME MODEL [ENV...]: menu 1 of models/MODEL in $WORK/NAME
  local name="$1" model="$2"; shift 2
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"
  ln -s "$HERE/models/data" "$WORK/$name/data"
  cp "$HERE/models/$model" "$WORK/$name/"
  (cd "$WORK/$name" && printf '1\n\n\n7\n' | env "$@" $(guard_command 900) "$AZURE2_BIN" --no-gui --no-readline "$model" > run.log 2>&1)
}

# ---- 1. the identity on one mesh
run identity identity.azr
if python3 -I "$HERE/compare.py" 1e-9 "$HERE/models/identity.azr" "$WORK/identity/output" > "$WORK/identity.txt" 2>&1; then
  ok "one project: every shifted data set = its twin at the shifted energies ($(tail -1 "$WORK/identity.txt"))"
else
  bad "one project: shifted data sets against their twins"; grep -v '^ok' "$WORK/identity.txt" | sed -n '1,12p'
fi

# ---- 2. the identity on two meshes
run shifted shifted.azr
run twins twins.azr
if python3 -I "$HERE/compare.py" 3e-6 "$HERE/models/shifted.azr" "$WORK/shifted/output" "$HERE/models/twins.azr" \
     "$WORK/twins/output" > "$WORK/apart.txt" 2>&1; then
  ok "two projects: every shifted data set = its twin at the shifted energies ($(tail -1 "$WORK/apart.txt"))"
else
  bad "two projects: shifted data sets against their twins"; grep -v '^ok' "$WORK/apart.txt" | sed -n '1,12p'
fi

# ---- 3. the gradient of a fitted shift against a central difference of chi^2
mkdir -p "$WORK/grad/output" "$WORK/grad/checks"
ln -s "$HERE/models/data" "$WORK/grad/data"
cp "$HERE/models/fit.azr" "$WORK/grad/"
(cd "$WORK/grad" && printf '2\n\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50 $(guard_command 900) "$AZURE2_BIN" \
   --use-gradient --no-gui --no-readline fit.azr > grad.log 2>&1)
line="$(grep -a '\[grad-check\] idx .* shift ' "$WORK/grad/grad.log" | sed -n 1p)"
relJ="$(echo "$line" | grep -o 'relJ [-0-9.e+]*' | cut -d' ' -f2)"
if [ -n "$relJ" ] && awk -v r="$relJ" 'BEGIN { exit !(r + 0 < 1e-5) }'; then
  ok "the shift's gradient (Jacobian path) = the central difference of chi^2 (relative $relJ, tolerance 1e-5)"
else
  bad "the shift's gradient against the central difference of chi^2: '${line:-no grad-check line}'"
fi

echo
if [ "$failures" -eq 0 ]; then echo "shift_true_energy: PASS"; else echo "shift_true_energy: FAIL ($failures)"; exit 1; fi
