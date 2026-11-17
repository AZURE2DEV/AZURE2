#!/usr/bin/env bash
# Capture into a BROAD unbound final state: the line-shape integration (the field `line_window` of a capture data set
# or extrapolation; src/LineWindow.cpp).  The cross section is integrated over the final-state energy E_f within
# E_R +- W, with the photon energy E - E_x(E_f), the real scattering wave of the final channel at E_f as the exterior
# of external capture, and the one-level line shape rho(E_f)/N_f as the weight.  Without the field a capture keeps the
# narrow-resonance treatment (tests/ec_unbound_final) unchanged.
#
# Models (models/, dev's, written by AZURE2 --convert): a_* and b_*, 12C + p with capture pair 3 into the 13N 1/2+
# state at E_x = 2.3655 MeV (E_R = 422 keV, Gamma_p = 39 keV as in Rolfs & Azuma, NPA 227, 291 (1974)), a = 4.48 fm,
# every other width 0 (pure E1 direct capture from hard-sphere p waves; b_*: Gamma_p = 39 eV); d_*, g_*, h_*: the
# same with a 5/2+ level at E_R = 1.60 MeV (Gamma_p 60 keV) that decays by E2 (0.05 eV) into the 2.3655 state.
#   (a) direct capture, windowed over unwindowed, at four energies for W = 39 and 117 keV against an independent
#       reference (damped radial integrals + Richardson, mpmath/ODE Coulomb functions, quadrature over the window;
#       the numbers below): within 5e-5 / 3e-4;
#   (b) narrow limit: final width 39 eV, W = 20 Gamma: the reference (the window fraction of the line plus the
#       non-resonant continuum term) within 5e-5;
#   (c) internal capture alone (the E2 part, angle-integrated, so E1 and E2 add incoherently): windowed over
#       unwindowed = int rho/N_f ((E - E_f)/(E - E_R))^5 dE_f (scipy quad; the numbers below), within 2e-5;
#   (d) --no-brune gives the Brune result;
#   (e) set-up: a NOTE per windowed extrapolation; `line_window` on an elastic extrapolation is refused by the
#       validator (E_FIELD_NOT_APPLICABLE);
#   (f) internal capture alone, external capture into the pair OFF (the final level is then found by J^pi and energy
#       and its N_f formed from its channels): the ratio of (c), and the cross section itself equal to the internal
#       part of (c) within 2e-5; --no-brune the same; a pair whose excitation energy is 2 keV off the level stops with
#       an ERROR;
#   (g) the line shape follows the final level's CURRENT energy (Brune parametrizations): the level 20 keV higher
#       through a parameter file gives the cross sections of the project with the level there, within 1e-6.
#   (h) a capture cascade (`secondary_gamma` on the capture exit) with a line window is refused: the line-shape
#       integration has no secondary-gamma form (dev's rule; the two features meet on the combined branch only).
#   (i) a fit of a windowed data set with a fitter that reads the analytic gradient or the residual rows' Jacobian
#       (its row of Fitters(): lm, gsl-lm, migrad-gradient) is refused at set-up with no output: the adjoint does not
#       run over the line shape; plain MIGRAD (numerical gradient) fits it and writes its outputs.
#
#   ./tests/ec_unbound_broad/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/edit_guard.sh"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_unbound_broad.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MODEL [OPTIONS]: `extrapolate` of models/MODEL.azr in $WORK/NAME; returns the engine's exit status
run() {
  local name="$1" model="$2" options="${3:-}"
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"
  [ -f "$WORK/$name/t.azr" ] || cp "$HERE/models/$model.azr" "$WORK/$name/t.azr"
  (cd "$WORK/$name" && $(guard_command 900) "$AZURE2_BIN" extrapolate $options t.azr < /dev/null > run.log 2>&1)
}
CAP3="output/AZUREOut_aa=1_R=3.extrap"
# column 4 of NAME's capture output, one value per line (blank lines between extrapolations dropped)
values() { awk 'NF { print $4 }' "$WORK/$1/$CAP3" 2>/dev/null; }
# compare LABEL NUM DEN TOL REF...: (NUM/DEN)_i against REF_i, relative TOL
compare() {
  local label="$1" num="$2" den="$3" tol="$4"; shift 4
  local refs="$*"
  local out
  out=$(paste -d ' ' <(values "$num") <(values "$den") | awk -v refs="$refs" -v tol="$tol" '
    BEGIN { n = split(refs, r, " ") }
    NF == 2 && $2 + 0 != 0 { i++; v = $1 / $2; d = v / r[i] - 1; if (d < 0) d = -d; if (d > worst) worst = d
      printf "      %d: %.6f (ref %.6f)\n", i, v, r[i] }
    END { if (i != n) { print "      MISSING: " i + 0 " of " n " values"; exit 1 }
          printf "      worst relative deviation %.1e (tolerance %s)\n", worst, tol; exit !(worst <= tol + 0) }')
  local status=$?
  echo "$out"
  if [ $status -eq 0 ]; then ok "$label"; else bad "$label"; fi
}

echo "(a) direct capture through +-1 and +-3 Gamma windows vs the independent reference"
for m in a_none a_w1 a_w3; do run $m $m || bad "(a) $m: engine exit status $?"; done
compare "(a) W = 39 keV" a_w1 a_none 5e-5 0.70442 0.70120 0.69904 0.69831
compare "(a) W = 117 keV" a_w3 a_none 3e-4 1.84144 1.03167 0.91742 0.89875

echo "(b) narrow limit: final width 39 eV, W = 20 Gamma"
for m in b_none b_w20; do run $m $m || bad "(b) $m: engine exit status $?"; done
compare "(b) narrow limit" b_w20 b_none 5e-5 1.06495 0.99331

echo "(c) internal E2 capture alone, angle-integrated"
for m in d_none d_w1 d_ecnone d_ecw1; do run $m $m || bad "(c) $m: engine exit status $?"; done
# internal part = total - (external alone), windowed over unwindowed
mkdir -p "$WORK/d_int_w/output" "$WORK/d_int_n/output"
paste -d ' ' <(values d_w1) <(values d_ecw1) | awk '{ print 0, 0, 0, $1 - $2 }' > "$WORK/d_int_w/$CAP3"
paste -d ' ' <(values d_none) <(values d_ecnone) | awk '{ print 0, 0, 0, $1 - $2 }' > "$WORK/d_int_n/$CAP3"
compare "(c) internal capture over the line shape" d_int_w d_int_n 2e-5 0.698914 0.697808 0.697589

echo "(d) --no-brune"
run e_formal a_w1 "--no-brune" || bad "(d) engine exit status $?"
compare "(d) --no-brune = Brune" e_formal a_w1 1e-9 1 1 1 1

echo "(e) set-up"
if grep -q "NOTE: test segment #1: capture into a broad unbound final state" "$WORK/a_w1/run.log"; then
  ok "(e) NOTE for a windowed extrapolation"
else
  bad "(e) no NOTE for a windowed extrapolation"
fi
mkdir -p "$WORK/f_elastic/output" "$WORK/f_elastic/checks"
awk '!done && $1 == "reaction" { sub(/12C\(p,g\)13N\[[0-9.]*\]/, "12C(p,p)12C"); done = 1 } { print }' \
  "$HERE/models/a_w1.azr" > "$WORK/f_elastic/t.azr"
grep -q "reaction  *12C(p,p)12C" "$WORK/f_elastic/t.azr" || bad "(e) the elastic variant was not written"
run f_elastic a_w1
if grep -q "E_FIELD_NOT_APPLICABLE.*line_window" "$WORK/f_elastic/run.log" && [ ! -s "$WORK/f_elastic/output/AZUREOut_aa=1_R=1.extrap" ]; then
  ok "(e) line_window on an elastic extrapolation is refused"
else
  bad "(e) line_window on an elastic extrapolation was not refused"
  tail -5 "$WORK/f_elastic/run.log"
fi

echo "(f) internal capture alone, external capture into the pair off"
for m in g_none g_w1; do run $m $m || bad "(f) $m: engine exit status $?"; done
run g_formal g_w1 "--no-brune" || bad "(f) --no-brune: engine exit status $?"
compare "(f) internal capture alone over the line shape" g_w1 g_none 2e-5 0.698914 0.697808 0.697589
# (the internal part of (c) is a difference of two printed totals, kept to 6 digits)
compare "(f) = the internal part with external capture on" g_w1 d_int_w 2e-5 1 1 1
compare "(f) --no-brune = Brune" g_formal g_w1 1e-9 1 1 1
run g_off g_off
if grep -q "but its final state (J^pi = 0.5+, E_x = 2.3675 MeV) is not a level of the project" "$WORK/g_off/run.log" &&
  [ ! -s "$WORK/g_off/$CAP3" ]; then
  ok "(f) a pair energy off the level stops with an ERROR"
else
  bad "(f) a pair energy off the level was not refused"
  tail -5 "$WORK/g_off/run.log"
fi

echo "(g) the line follows the final level's current energy"
# the same parameters reach the engine twice: with the level (and its pair) 20 keV higher in the project file
# (h_file), and through that run's parameter file on the unchanged project (d_w1)
run h_file h_file || bad "(g) h_file: engine exit status $?"
mkdir -p "$WORK/h_par"
cp "$WORK/h_file/output/param.par" "$WORK/h_par/moved.par" 2>/dev/null
run h_par d_w1 "--parameter-file=moved.par" || bad "(g) h_par: engine exit status $?"
compare "(g) level energy from the parameter file = level energy in the project file" h_par h_file 1e-6 1 1 1
if paste -d ' ' <(values h_par) <(values d_w1) | awk 'NF == 2 && $2 + 0 != 0 { d = $1 / $2 - 1; if (d < 0) d = -d; n++; if (d > 0.02) moved++ } END { exit !(n > 0 && moved == n) }'; then
  ok "(g) the 20 keV move changes every cross section by more than 2 %"
else
  bad "(g) the moved level energy did not change the cross sections"
fi

echo "(h) not together with a capture cascade"
mkdir -p "$WORK/i_cascade/output" "$WORK/i_cascade/checks"
# the windowed project with one data set added: the secondary gamma ray of the same capture, one point
printf '0.704561  90  1e-6  1e-7\n' > "$WORK/i_cascade/cascade.dat"
awk '{ print } $1 == "<data>" {
  print "dataset cascade"
  print "    reaction        12C(p,g)13N[2.3655]"
  print "    quantity        dsigma_domega"
  print "    unit            b/sr"
  print "    line_window     39 keV"
  print "    secondary_gamma L 1 Ic 1/2 delta 0"
  print "    file            cascade.dat"
  print "    columns         E_lab[MeV] theta_cm[deg] value error"
  print "    energy_window   0.5 1.0 MeV"
  print "    angle_window    0 180 deg"
  print "    norm            1 fixed"
}' "$HERE/models/a_w1.azr" > "$WORK/i_cascade/t.azr"
# a data set is read by a calculation, not by an extrapolation
(cd "$WORK/i_cascade" && $(guard_command 900) "$AZURE2_BIN" calculate t.azr < /dev/null > run.log 2>&1)
if grep -q "ERROR: .*capture cascade.*not together with a line window" "$WORK/i_cascade/run.log" &&
  [ ! -s "$WORK/i_cascade/output/AZUREOut_aa=1_R=3.out" ]; then
  ok "(h) secondary_gamma with a line window is refused"
else
  bad "(h) secondary_gamma with a line window was not refused"
  tail -5 "$WORK/i_cascade/run.log"
fi

echo "(i) fitters that read the analytic derivatives are refused; plain MIGRAD fits"
# the windowed project with one data point; the final level's energy fixed (its widths stay free), so that MIGRAD's
# steps do not move the line centre and rebuild the line's capture integrals (seconds instead of minutes)
mkdir -p "$WORK/j_fit"
printf '0.704561  90  1e-6  1e-7\n' > "$WORK/j_fit/lw.dat"
edited "$HERE/models/a_w1.azr" "$WORK/j_fit/base.azr" awk '$1 == "<data>" {
  print
  print "dataset lw"
  print "    reaction       12C(p,g)13N[2.3655]"
  print "    quantity       dsigma_domega"
  print "    unit           b/sr"
  print "    line_window    39 keV"
  print "    file           lw.dat"
  print "    columns        E_lab[MeV] theta_lab[deg] value error"
  print "    norm           1 fixed"
  next
} { print }' || bad "(i) the data set was not added"
edited "$WORK/j_fit/base.azr" "$WORK/j_fit/t.azr" awk '$1 == "0.5" && $2 == "1" && $3 == "2.3655" && $4 == "0" { $4 = "1" }
  { print }' || bad "(i) the final level's energy was not fixed"
# fit_run NAME [SWITCH]: `fit` of $WORK/j_fit/t.azr in $WORK/NAME
fit_run() {
  local name="$1" switch="${2:-}"
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"
  cp "$WORK/j_fit/t.azr" "$WORK/j_fit/lw.dat" "$WORK/$name/"
  (cd "$WORK/$name" && $(guard_command 900) "$AZURE2_BIN" fit $switch --max-calls=40 t.azr < /dev/null > run.log 2>&1)
}
help_text="$("$AZURE2_BIN" --help 2>&1)"
for fitter in "lm --use-lm" "gsl-lm --use-gsl-lm" "migrad-gradient --use-gradient"; do
  name="${fitter%% *}" switch="${fitter#* }"
  if ! grep -q -- "$switch" <<< "$help_text"; then
    echo "  skip  (i) $switch: not in this build"
    continue
  fi
  fit_run "fit_$name" "$switch"
  status=$?
  if [ $status -ne 0 ] && grep -q "ERROR: segment #1 .*line window.*the fitter $name ($switch) is not supported" \
    "$WORK/fit_$name/run.log" && [ ! -s "$WORK/fit_$name/output/param.sav" ]; then
    ok "(i) $switch is refused, nothing written"
  else
    bad "(i) $switch was not refused (exit $status)"
    grep -a "ERROR" "$WORK/fit_$name/run.log" | sed -n '1,2p' | cut -c1-200
  fi
done
fit_run fit_migrad ""
status=$?
if [ $status -eq 0 ] && [ -s "$WORK/fit_migrad/output/param.sav" ] && [ -s "$WORK/fit_migrad/output/chiSquared.out" ] &&
  ! grep -q "ERROR" "$WORK/fit_migrad/run.log"; then
  ok "(i) plain MIGRAD fits the windowed data set"
else
  bad "(i) plain MIGRAD did not fit (exit $status)"
  tail -5 "$WORK/fit_migrad/run.log"
fi

if [ $failures -eq 0 ]; then
  echo "PASS: capture into a broad unbound final state"
  exit 0
fi
echo "FAIL ($failures)"
exit 1
