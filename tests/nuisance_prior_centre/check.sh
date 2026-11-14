#!/usr/bin/env bash
#
# The centre of a nuisance parameter's prior is the number written (docs/upgrade/parameter_kinds.md section 4): with
# `<value> fit prior <centre> <width>` and value != centre, the penalty at the start is ((value - centre)/width)^2,
# for every kind that can carry a prior.  Calculations (mode 1) on copies of projects of the tree, each line moved so
# that value != centre; the terms of chiSquared.out against the numbers put together by hand:
#
#   (a) T2 (tests/upgrade/T2): norm 1 prior 1.05 5 % (0.907029); shift offset, sqrt and the shift function's a0 one
#       width off each (Total-Shift 3); isotope 15N 0.0037 prior 0.0036 0.0002 (0.25) and isotope_sum prior 0.9 1
#       on the sum 0.0037 (0.803354), Total-Abundance 1.05335 (six digits); areal density 2e17 prior 2.2e17 10 % (0.826446);
#   (b) screening 1000 eV prior 1200 200 (tests/electron_screening/grad.azr, 13N data): Total-Screening 1;
#   (c) stopping scale 1 prior 1.04 4 % (tests/stopping_yield/models/fit.azr, both targets): 0.924556;
#   (d) a kernel's width 14 keV prior 16 2 (tests/rectangular_resolution/models/fit.azr): Total-Kernel-Parameter 1;
#   (e) Doppler temperature 300 K prior 320 20 (tests/doppler_fgm models, dop): Total-Doppler 1;
#   (f) sample thickness 0.15 atoms/b prior 0.16 0.01 (tests/transmission_yield models, conv): Total-Thickness 1;
#   (g) where pyazr's module loads (PYAZR_ENGINE_DIR): the write-back of a fitted screening potential written in
#       keV, a Doppler temperature, an areal density and a sample thickness (_azure2.parameter_value_edits): the
#       value token only, in the unit written, the prior untouched.
#   (h) continued.sh: for an areal density, a Doppler temperature, a stopping scale and a kernel's parameter, the
#       project with the fitted value written as the start computes the fit's chi^2 and a second fit stays at the
#       minimum (the grids and limits are laid out from the prior's centre); a start outside the centre's range is
#       laid out from the start, with a note.
#
# Usage: tests/nuisance_prior_centre/check.sh path/to/AZURE2     (about a minute; (h) at two threads)

# CRLF-writing Windows build compares like the others.
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTS="$HERE/.."
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PY="$(command -v python3 || command -v python)"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/nuisance_prior_centre.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME SRC_DIR PROJECT SED: a copy of SRC_DIR (data linked), the project edited by SED, a calculation
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  for f in "$2"/*; do [ "$(basename "$f")" = output ] || ln -s "$f" "$dir/$(basename "$f")"; done
  rm -f "$dir/$3"; sed -e "$4" "$2/$3" > "$dir/$3"
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline "$3" > run.log 2>&1)
}
# term NAME LABEL: the value of " Total-LABEL-Chi-Squared:" in NAME's chiSquared.out ("" when absent)
term() {
  awk -v k="Total-$2-Chi-Squared:" '{ sub(/\r$/, "") } $1 == "Total-Chi-Squared:" { for (i = 1; i < NF; i++) if ($i == k) print $(i + 1) }' \
    "$WORK/$1/output/chiSquared.out" 2>/dev/null
}
# equal to the six digits chiSquared.out writes
near() { awk -v a="$1" -v b="$2" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(b ~ /^ *[-+]?\.?[0-9]/ && a != "" && d <= 1e-5 * (b < 0 ? -b : b) + 1e-12) }'; }
expect() {  # NAME LABEL WANT WHAT
  local got; got="$(term "$1" "$2")"
  if near "$got" "$3"; then ok "$4: Total-$2-Chi-Squared $got"; else bad "$4: Total-$2-Chi-Squared '$got', want $3"; fi
}

# (a) T2: every kind of its data sets and its natural target at once
T2SED='/^dataset s01_int_thin$/,/^$/s/norm           1 fit prior 1 5%/norm           1 fit prior 1.05 5%/
/^dataset s01_int_thin$/,/^$/s/shift          offset 0 MeV fit prior 0 0.003/shift          offset 0 MeV fit prior 0.003 0.003/
s/shift          sqrt 0 MeV^1\/2 fit prior 0 0.003/shift          sqrt 0 MeV^1\/2 fit prior 0.003 0.003/
s/param          a0 0 fit prior 0 0.01/param          a0 0 fit prior 0.01 0.01/
s/isotope        15N 0.0037 fit prior 0.0037 0.0002/isotope        15N 0.0037 fit prior 0.0036 0.0002/
s/isotope_sum    none/isotope_sum    prior 0.9 1/
s/areal_density  2e+17 atoms\/cm2 fit prior 2e+17 10%/areal_density  2e+17 atoms\/cm2 fit prior 2.2e+17 10%/'
run t2 "$TESTS/upgrade/T2" T2.azr "$T2SED"
expect t2 Norm 0.907029 "(a) a data set's norm"
expect t2 Shift 3 "(a) shift offset, sqrt and function parameter"
expect t2 Abundance 1.05335 "(a) isotope and isotope_sum"
expect t2 Density 0.826446 "(a) areal density (10 % of the centre)"

# (b) screening, on the 13N data
mkdir -p "$WORK/src_scr"; ln -s "$TESTS/13N/data" "$WORK/src_scr/data"; cp "$TESTS/electron_screening/grad.azr" "$WORK/src_scr/"
run scr "$WORK/src_scr" grad.azr 's/screening      1000 eV fit prior 1000 200/screening      1000 eV fit prior 1200 200/'
expect scr Screening 1 "(b) screening potential"

# (c) stopping scale (one parameter of two targets)
run stop "$TESTS/stopping_yield/models" fit.azr 's/scale 1 fit prior 1 4%/scale 1 fit prior 1.04 4%/'
expect stop Stopping-Scale 0.924556 "(c) stopping scale (4 % of the centre)"

# (d) a table kernel's width
run rect "$TESTS/rectangular_resolution/models" fit.azr 's/rectangular width 14 keV fit$/rectangular width 14 keV fit prior 16 2/'
expect rect Kernel-Parameter 1 "(d) a kernel's parameter"

# (e), (f) the generated models of tests/doppler_fgm and tests/transmission_yield
if "$PY" -I "$TESTS/doppler_fgm/make_models.py" "$TESTS/reaction_rate/models/n12C.azr" "$WORK/dop_models" > /dev/null 2>&1; then
  run dop "$WORK/dop_models" dop.azr 's/doppler        free_gas 300 K fixed/doppler        free_gas 300 K fit prior 320 20/'
  expect dop Doppler 1 "(e) Doppler temperature"
else bad "(e) the Doppler models were not built"; fi
if "$PY" -I "$TESTS/transmission_yield/make_models.py" "$TESTS/reaction_rate/models/n12C.azr" "$TESTS/natural_neon/natural_neon.azr" \
     "$WORK/tr_models" > /dev/null 2>&1; then
  run thk "$WORK/tr_models" conv.azr '0,/areal_density  0.15 atoms\/b fixed/s//areal_density  0.15 atoms\/b fit prior 0.16 0.01/'
  expect thk Thickness 1 "(f) sample thickness"
else bad "(f) the transmission models were not built"; fi

# (g) the write-back of the target's fitted values (needs the compiled module)
if [ -n "${PYAZR_ENGINE_DIR:-}" ]; then
  out="$("$PY" "$HERE/write_back.py" "$WORK" 2>&1)"; rc=$?
  echo "$out" | grep -E '^  (ok|FAIL) ' || true
  [ $rc -eq 0 ] || bad "(g) write-back: $(echo "$out" | tail -1)"
else
  echo "  skip  (g) write-back: PYAZR_ENGINE_DIR is not set"
fi

# (h) a fit continued from its written values computes the fit it continues (continued.sh: the four kinds whose value
#     lays out a grid; its own header lists the checks)
if OMP_NUM_THREADS=2 bash "$HERE/continued.sh" "$AZURE2_BIN" > "$WORK/continued.log" 2>&1; then
  grep -E '^  (ok|skip) ' "$WORK/continued.log"; ok "(h) continued fits"
else
  grep -E '^  (ok|FAIL|skip) |^[-<>] ' "$WORK/continued.log"; bad "(h) continued fits"
fi

if [ "$failures" -eq 0 ]; then echo "PASS: the prior centre of every nuisance kind"; else echo "FAIL: $failures check(s)"; exit 1; fi
