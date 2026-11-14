#!/usr/bin/env bash
#
# Shift functions (an energy_scale record with a function and its parameters, section 3.6 of project_format3.md):
# the shift of every data set that names the record becomes dE(x) = a + f(x; a0, a1, ...), x the point's original
# lab energy in MeV, a the record's offset, every aj a parameter (value, prior, fixed or fit).  The data sets that name
# one record share it (a shift group).
#
#   1. function "a0*sqrt(x)" is the sqrt shift term: identical model values (1e-10), point by point, also through
#      target integration (sub-points); the centred form a0*(sqrt(x)-sqrt(1.7)) with the offset moved by
#      a0*sqrt(1.7) too;
#   2. a0*x^1.5 equals the same data written at the shifted energies;
#   3. a second data set of a scale uses its function: equal to the same function in a record of its own;
#   4. a fit recovers a and a0 of a linear calibration error shared by two data sets (chi2 = 0 at the truth, only
#      reported with Minuit's numerical gradient), with --use-gradient and --use-lm, and writes the parameters
#      G.offset and G.a0, Total-Shift-Chi-Squared and the shifts.out function line;
#   5. --use-gradient's group-local shift columns agree with full-dataset finite differences
#      (AZURE_SHIFT_GRAD_CHECK);
#   6. invalid input stops with a message: a function plus a sqrt term on one data set, an undeclared or an unused
#      parameter, a term with no effect; two indistinguishable terms give a warning.  The refusals of the old block
#      grammar (a follower with a function of its own, sqrtshift and shiftfn on one line, a short block) are the
#      converter's: the old lines through AZURE2 --convert, refused with the old engine's reason;
#   7. a .sav value does not free a fixed function parameter (design note D4).
#
#   ./tests/shift_fn/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/shift_fn.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
seg() { awk -F, -v k="$2" '$1 == k { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
field() { awk -v f="$2" '/^Total-Chi-Squared:/ { for (i = 1; i < NF; i++) if ($i == f) print $(i+1) }' "$1/output/chiSquared.out" 2>/dev/null; }
par() { awk -v n="$2" '$1 == n { print $2 }' "$1/output/param.sav" 2>/dev/null; }
par_file() { awk -v n="$2" '$1 == n { print $2 }' "$1" 2>/dev/null; }
relclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m == 0) m = 1e-300; exit !(b ~ /^ *[-+]?\.?[0-9]/ && d <= t * m) }'; }
absclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d < t) }'; }
# Largest relative difference between the model columns (4) of two data sets' blocks in AZUREOut_aa=1_R=1.out
# (blocks are separated by blank lines, in the order of the data sets).
blockdiff() {  # DIR SET_A SET_B
  awk -v A="$2" -v B="$3" '
    BEGIN { blk = 1 }
    NF == 0 { if (prev) blk++; prev = 0; next }
    { prev = 1; if (blk == A) a[na++] = $4; if (blk == B) b[nb++] = $4 }
    END {
      if (na == 0 || na != nb) { print "mismatch"; exit }
      m = 0; for (i = 0; i < na; i++) { d = a[i] - b[i]; if (d < 0) d = -d; s = b[i] < 0 ? -b[i] : b[i]; if (s > 0 && d / s > m) m = d / s }
      printf "%.3g\n", m
    }' "$1/output/AZUREOut_aa=1_R=1.out"
}

# ds NAME FILE ANGLE [FIELD...]: a data set of the model (p + p elastic at ANGLE degrees, data/FILE) with these extra
# fields ("shift offset 0.004 MeV fixed", "shift G", ...); TARGET=1 adds the thin target of section 1
ds() {
  local name="$1" file="$2" angle="$3"; shift 3
  printf 'dataset %s\n    reaction       1H(p,p)1H\n    quantity       dsigma_domega\n    unit           b/sr\n' "$name"
  printf '    file           data/%s\n    columns        E_lab[MeV] theta_cm[deg] value error\n' "$file"
  printf '    energy_window  0 100 MeV\n    angle_window   %s %s deg\n    norm           1 fixed\n' "$angle" "$angle"
  local f; for f in "$@"; do printf '    %s\n' "$f"; done
  [ "${TARGET:-0}" = 1 ] && printf '    target         thin\n    quadrature     points 20\n'
  echo
}
# scale NAME FUNCTION OFFSET PARAM...: an energy_scale record; OFFSET "" for none, each PARAM "a0 0.003 fixed"
scale() {
  local name="$1" fn="$2" off="$3"; shift 3
  printf 'energy_scale %s\n' "$name"
  [ -n "$off" ] && printf '    offset         %s\n' "$off"
  printf '    function       "%s" MeV\n' "$fn"
  local p; for p in "$@"; do printf '    param          %s\n' "$p"; done
  echo
}
THIN='target thin
    kind           layer
    active         1H
    areal_density  1e18 atoms/cm2 fixed
    stopping       formula "a0" a0=1e-20 MeV*cm2/atom
    straggling     none'

# mk NAME LEVELS_FIXED RECORDS -> $WORK/NAME/t.azr: tests/shift_fn/shift_fn.azr with RECORDS as its <data> section
# (and the thin target in <targets> when a record names it); LEVELS_FIXED 1 fixes every level energy and width
mk() {
  local name="$1" fix="$2" dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/data" "$dir/data"
  # the records go through the environment: BSD awk (macOS) refuses a newline in -v ("newline in string")
  RECORDS="$3" TARGETS="$(case "$3" in *"target         thin"*) printf '%s' "$THIN" ;; esac)" awk -v fix="$fix" '
    /<levels>/ { inlev = 1; print; next }
    /<\/levels>/ { inlev = 0 }
    inlev && fix && NF > 12 { $4 = 1; $11 = 1; print; next }   # fix every energy and width
    /^<data>/ { print; printf "%s\n", ENVIRON["RECORDS"]; skip = 1; next }
    /^<\/data>/ { skip = 0 }
    /^<targets>/ { print; if (ENVIRON["TARGETS"] != "") print ENVIRON["TARGETS"]; skip = 1; next }
    /^<\/targets>/ { skip = 0 }
    skip { next }
    { print }' "$HERE/shift_fn.azr" > "$dir/t.azr"
}
run() {  # NAME MODE EXTRA_FLAGS [SAV]
  local input
  if [ "$2" = 1 ]; then input='1\n\n\n7\n'; else input="2\n\n${4:-}\n\n7\n"; fi
  (cd "$WORK/$1" && printf "$input" | "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -s "$WORK/$1/output/chiSquared.out" ]; then
    bad "$1: no chiSquared.out"; tail -15 "$WORK/$1/run.log"; return 1
  fi
}
fails_with() {  # NAME MESSAGE_FRAGMENT RECORDS
  mk "$1" 0 "$3"
  (cd "$WORK/$1" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  if [ -s "$WORK/$1/output/chiSquared.out" ]; then bad "$1: ran, should have stopped"
  elif grep -qF "$2" "$WORK/$1/run.log"; then ok "$1: stops with \"$2\""
  else bad "$1: no \"$2\" in the log"; tail -5 "$WORK/$1/run.log"; fi
}
# The old block grammar, through the converter: the original project (tests/convert/fixtures) with these
# <segmentsData> lines, converted; the old engine refused each, so the converter must refuse it or the conversion must
# be refused by the run, with MESSAGE either way.
OLD="$HERE/../convert/fixtures/tests/shift_fn/shift_fn.azr"
convert_refuses() {  # NAME MESSAGE_FRAGMENT SEGMENT_LINE...
  local name="$1" msg="$2" dir="$WORK/$1"; shift 2
  mkdir -p "$dir/output" "$dir/checks"; cp -R "$HERE/data" "$dir/data"
  SEGS="$(printf '%s\n' "$@")" awk '
    /<segmentsData>/ { print; print ENVIRON["SEGS"]; skip = 1; next }
    /<\/segmentsData>/ { skip = 0 }
    skip { next }
    { print }' "$OLD" > "$dir/old.azr"
  (cd "$dir" && "$AZURE2_BIN" --convert old.azr t.azr --to 3 > convert.log 2>&1)
  local where="the converter"
  if [ -f "$dir/t.azr" ]; then
    where="the run of the conversion"
    (cd "$dir" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr >> convert.log 2>&1)
  fi
  if [ -s "$dir/output/chiSquared.out" ]; then bad "$name: the conversion ran, the old engine refused it"
  elif grep -qF "$msg" "$dir/convert.log"; then ok "$name: refused by $where (\"$msg\")"
  else bad "$name: no \"$msg\" from $where"; tail -3 "$dir/convert.log"; fi
}

# ---- 1. sqrt(E): function "a0*sqrt(x)" == the sqrt term ---------------------------
for t in plain target; do
  tg=0; [ "$t" = target ] && tg=1
  mk "sqrt_$t" 0 "$(scale s2 'a0*sqrt(x)' '' 'a0 0.003 fixed'
                    scale s3 'a0*(sqrt(x)-sqrt(1.7))' '' 'a0 0.003 fixed'
                    TARGET=$tg ds xs_40 xs_40.dat 40 'shift          offset 0.004 MeV fixed' 'shift          sqrt 0.003 MeV^1/2 fixed'
                    TARGET=$tg ds xs_40_2 xs_40.dat 40 'shift          offset 0.004 MeV fixed' 'shift          s2'
                    TARGET=$tg ds xs_40_3 xs_40.dat 40 'shift          offset 0.0079115 MeV fixed' 'shift          s3')"
  # 0.0079115 = 0.004 + 0.003*sqrt(1.7) to 7 digits; data set 3 is checked to 1e-6 only.
  if run "sqrt_$t" 1 ""; then
    d12="$(blockdiff "$WORK/sqrt_$t" 1 2)"; d13="$(blockdiff "$WORK/sqrt_$t" 1 3)"
    if awk -v d="$d12" 'BEGIN { exit !(d ~ /^ *[-+]?\.?[0-9]/ && d != "mismatch" && d + 0 <= 1e-10) }'; then ok "$t: a0*sqrt(x) = the sqrt term, model values to $d12"; else bad "$t: a0*sqrt(x) vs the sqrt term differ by $d12"; fi
    if awk -v d="$d13" 'BEGIN { exit !(d ~ /^ *[-+]?\.?[0-9]/ && d != "mismatch" && d + 0 <= 1e-6) }'; then ok "$t: centred form agrees ($d13)"; else bad "$t: centred form differs by $d13"; fi
    s1="$(seg "$WORK/sqrt_$t" 1)"
    if ! relclose "$s1" 24 1e-3; then ok "$t: the shift is really applied (chi2 $s1, unshifted 24)"; else bad "$t: chi2 $s1 = the unshifted 24"; fi
  fi
done
if [ -n "$(seg "$WORK/sqrt_plain" 1)" ] && ! relclose "$(seg "$WORK/sqrt_target" 1)" "$(seg "$WORK/sqrt_plain" 1)" 1e-3; then
  ok "target integration really on (chi2 $(seg "$WORK/sqrt_target" 1) vs $(seg "$WORK/sqrt_plain" 1))"
else bad "target integration made no difference"; fi

# ---- 2. a0*x^1.5 == the same data written at the shifted energies -------------
mk x15 0 "$(scale s1 'a0*x^1.5' '' 'a0 0.003 fixed'
            ds xs_40 xs_40.dat 40 'shift          offset 0.002 MeV fixed' 'shift          s1'
            ds xs_40_x15_preshifted xs_40_x15_preshifted.dat 40)"
if run x15 1 ""; then
  d="$(blockdiff "$WORK/x15" 1 2)"
  if awk -v d="$d" 'BEGIN { exit !(d ~ /^ *[-+]?\.?[0-9]/ && d != "mismatch" && d + 0 <= 1e-8) }'; then ok "a0*x^1.5 = pre-shifted data ($d)"; else bad "a0*x^1.5 vs pre-shifted data: $d"; fi
fi

# ---- 3. a second data set of a scale uses its function --------------------------
mk follower 0 "$(scale G 'a0*x' '0.001 MeV fixed' 'a0 0.003 fixed'
                 scale H 'a0*x' '0.001 MeV fixed' 'a0 0.003 fixed'
                 ds xs_40 xs_40.dat 40 'shift          G'
                 ds xs_90 xs_90.dat 90 'shift          G'
                 ds xs_90_2 xs_90.dat 90 'shift          H')"
if run follower 1 ""; then
  d="$(blockdiff "$WORK/follower" 2 3)"
  if awk -v d="$d" 'BEGIN { exit !(d ~ /^ *[-+]?\.?[0-9]/ && d != "mismatch" && d + 0 <= 1e-10) }'; then ok "second data set = the same function given in its own record ($d)"; else bad "second data set vs own record: $d"; fi
  if grep -q "shift function on segment #1 (group 'G', segments #1 2)" "$WORK/follower/run.log"; then ok "function and group noted in the log"; else bad "no shift-function note in the log"; fi
fi

# ---- 4. fit a shared linear calibration error -----------------------------------
# Both files carry labels L with L + a + a0*(L - 1) = E, a = -0.004 (the shift at
# 1 MeV), a0 = 0.006, and the model's exact values as yields: chi2 = 0 at the truth.
# With Minuit's numerical gradient the chi2 is only reported: on exact data it is
# very sensitive to where MIGRAD stops (macOS CI 2026-10-21: chi2 808068 with a and
# a0 within 5e-8 of the truth, then 0 on a rerun), so a and a0 decide there.
GROUP="$(scale G 'a0*(x-1)' '0 MeV fit prior 0 1' 'a0 0 fit prior 0 1'
         ds xs_40_lin_mislabelled xs_40_lin_mislabelled.dat 40 'shift          G'
         ds xs_90_lin_mislabelled xs_90_lin_mislabelled.dat 90 'shift          G')"
fit_linear() {  # NAME FLAGS [CHI2: strict (chi2 < 1e-3 required, default) | info]
  mk "$1" 1 "$GROUP"
  run "$1" 2 "$2" || return
  local dir="$WORK/$1" a a0 t f
  a="$(par "$dir" G.offset)"; a0="$(par "$dir" G.a0)"; t="$(total "$dir")"
  if absclose "$a" -0.004 1e-5 && absclose "$a0" 0.006 1e-5; then ok "$1: a = $a, a0 = $a0 (truth -0.004, 0.006)"; else bad "$1: a = $a, a0 = $a0, expected -0.004, 0.006"; fi
  if [ "${3:-strict}" = info ]; then
    echo "  info  $1: chi2 $t (not checked: a and a0 above decide)"
  elif awk -v t="$t" 'BEGIN { exit !(t ~ /^ *[-+]?\.?[0-9]/ && t + 0 < 1e-3) }'; then ok "$1: chi2 = 0 at the truth (chi2 $t)"; else bad "$1: chi2 $t, expected ~0"; fi
  f="$(field "$dir" Total-Shift-Chi-Squared:)"
  if [ -n "$f" ]; then ok "$1: chiSquared.out reports Total-Shift-Chi-Squared ($f)"; else bad "$1: no Total-Shift-Chi-Squared"; fi
  if grep -q "^# 1,1 2,a0\*(x-1)," "$dir/output/shifts.out"; then ok "$1: shifts.out has the function line"; else bad "$1: no function line in shifts.out"; fi
}
fit_linear fit_numerical "" info
fit_linear fit_analytic "--use-gradient"
fit_linear fit_lm "--use-lm"

# ---- 5. group-local shift gradient == full-dataset finite difference ------------
mk gradcheck 1 "$(scale G 'a0*(x-1)+a1*(x-1)^2' '0.001 MeV fit prior 0 1' 'a0 0.001 fit prior 0 1' 'a1 0.0002 fit prior 0 1'
                  ds xs_40_lin_mislabelled xs_40_lin_mislabelled.dat 40 'shift          G'
                  ds xs_90_lin_mislabelled xs_90_lin_mislabelled.dat 90 'shift          G')"
# Adjoint path: each shift column against a full-dataset difference (AZURE_SHIFT_GRAD_CHECK).
(cd "$WORK/gradcheck" && printf '2\n\n\n\n7\n' | AZURE_GRAD_ADJOINT=1 AZURE_SHIFT_GRAD_CHECK=1 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > adjoint.log 2>&1)
nchk="$(grep -c 'shift-grad-check' "$WORK/gradcheck/adjoint.log")"
worst="$(grep -o 'reldiff=[-0-9.e+]*' "$WORK/gradcheck/adjoint.log" | cut -d= -f2 | sort -g | tail -1)"
if [ "$nchk" -ge 3 ] && awk -v w="$worst" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 < 1e-6) }'; then ok "adjoint path: group-local shift columns = full differences ($nchk checks, worst $worst)"
else bad "adjoint-path shift gradient check: $nchk checks, worst reldiff '$worst'"; fi
# Jacobian path (the --use-gradient default): AZURE_GRAD_CHECK compares every free
# component (here a, a0, a1) with a finite difference of the whole chi2.
(cd "$WORK/gradcheck" && printf '2\n\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > jacobian.log 2>&1)
if grep -q "Jacobian path: shift 0/1 shift_fn 0/2" "$WORK/gradcheck/jacobian.log"; then ok "Jacobian path: a, a0, a1 agree with finite differences ($(grep -o 'worst relative difference [-0-9.e+]*' "$WORK/gradcheck/jacobian.log"))"
else bad "Jacobian-path gradient check: $(grep 'components with relative' "$WORK/gradcheck/jacobian.log")"; fi

# ---- 6. invalid input ---------------------------------------------------------
fails_with with_sqrt "two descriptions of one energy calibration" \
  "$(scale G 'a0*x' '' 'a0 0 fixed'; ds xs_40 xs_40.dat 40 'shift          sqrt 0.001 MeV^1/2 fixed' 'shift          G')"
fails_with undeclared "uses a1 but declares 1 parameter" \
  "$(scale G 'a0*x+a1' '' 'a0 0 fixed'; ds xs_40 xs_40.dat 40 'shift          G')"
fails_with unused "declares a1 but never uses it" \
  "$(scale G 'a0*x' '' 'a0 0 fixed' 'a1 0 fixed'; ds xs_40 xs_40.dat 40 'shift          G')"
fails_with no_effect "has no effect at this group's energies" \
  "$(scale G 'a0*0*x' '' 'a0 0 fit'; ds xs_40 xs_40.dat 40 'shift          G')"
mk degenerate 0 "$(scale G 'a0+0*x' '0 MeV fit prior 0 0.01' 'a0 0 fit prior 0 0.01'; ds xs_40 xs_40.dat 40 'shift          G')"
(cd "$WORK/degenerate" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
if grep -q "WARNING: shift function on segment #1: its free terms cannot be told apart" "$WORK/degenerate/run.log"; then ok "a constant function next to a free a: degeneracy warning"
else bad "no degeneracy warning for a0 + 0*x with a free"; fi
convert_refuses follower_fn "but it is not the first segment of shift group 'G'" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftgroup G" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0 shiftgroup G shiftfn \"a0*x\" 1 0 0 0"
convert_refuses old_with_sqrt "two descriptions of one energy calibration" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 sqrtshift 0.001 0 0 shiftfn \"a0*x\" 1 0 0 0"
convert_refuses short_block "declares 2 parameters but gives" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftfn \"a0*x+a1\" 2 0 0 0"

# ---- 7. a .sav value does not free a fixed function parameter -------------------
mk sav_fixed 1 "$(scale G 'a0*(x-1)' '0 MeV fit prior 0 1' 'a0 0 fixed'; ds xs_40_lin_mislabelled xs_40_lin_mislabelled.dat 40 'shift          G')"
if run sav_fixed 2 ""; then
  # POSIX BRE only: BSD sed (macOS) reads "\+" as a literal "+", which left the .sav unedited.
  sed 's/^\( *G[.]a0  *\)[-0-9.e+]*/\10.002/' "$WORK/sav_fixed/output/param.sav" > "$WORK/sav_fixed/start.sav"
  if ! absclose "$(par_file "$WORK/sav_fixed/start.sav" G.a0)" 0.002 1e-12; then
    bad "test setup: start.sav does not carry a0 = 0.002"
  elif run sav_fixed 2 "" start.sav; then
    v="$(par "$WORK/sav_fixed" G.a0)"
    if absclose "$v" 0.002 1e-12; then ok "fixed a0 stays at the .sav value ($v)"; else bad "fixed a0 moved to $v (.sav 0.002)"; fi
  fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: shift functions"; else echo "FAIL: $failures check(s)"; exit 1; fi
