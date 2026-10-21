#!/usr/bin/env bash
#
# Shift functions ("shiftfn" block at the end of a segment line): the shift becomes
# dE(x) = a + f(x; a0, a1, ...), x the point's original lab energy in MeV, with
# every aj a fit parameter (value, penalty width, vary flag).  On a shift group the
# function goes on the group's first segment and applies to every member.
#
#   1. shiftfn "a0*sqrt(x)" is the existing sqrtshift term: identical model values
#      (1e-10), point by point, also through target integration (sub-points);
#      the centred form a0*(sqrt(x)-sqrt(1.7)) with a moved by a0*sqrt(1.7) too;
#   2. a0*x^1.5 equals the same data written at the shifted energies;
#   3. a follower uses its leader's function;
#   4. a fit recovers a and a0 of a linear calibration error shared by two
#      segments (chi2 = 0 at the truth; only reported for the numerical gradient),
#      with Minuit's numerical gradient,
#      --use-gradient and --use-lm, and writes segment_<key>_energy_shift_fn0,
#      Total-Shift-Chi-Squared and the shifts.out function line;
#   5. --use-gradient's group-local shift columns agree with full-dataset finite
#      differences (AZURE_SHIFT_GRAD_CHECK);
#   6. invalid input stops with a message: a follower declaring a function, a
#      function plus sqrtshift, an undeclared or unused parameter, a short block;
#      a term with no effect is an error, two indistinguishable terms a warning;
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
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
seg() { awk -F, -v k="$2" '$1 == k { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
field() { awk -v f="$2" '/^Total-Chi-Squared:/ { for (i = 1; i < NF; i++) if ($i == f) print $(i+1) }' "$1/output/chiSquared.out" 2>/dev/null; }
par() { awk -v n="$2" '$1 == n { print $2 }' "$1/output/param.sav" 2>/dev/null; }
par_file() { awk -v n="$2" '$1 == n { print $2 }' "$1" 2>/dev/null; }
relclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m == 0) m = 1e-300; exit !(d <= t * m) }'; }
absclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d < t) }'; }
# Largest relative difference between the model columns (4) of two segments' blocks in
# AZUREOut_aa=1_R=1.out (blocks are separated by blank lines, in segment order).
blockdiff() {  # DIR SEG_A SEG_B
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

# mk NAME LEVELS_FIXED TARGETINT_LINE SEGMENT_LINE... -> $WORK/NAME/t.azr
mk() {
  local name="$1" fix="$2" tint="$3"; shift 3
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/data" "$dir/data"
  printf '%s\n' "$@" > "$dir/segs.txt"
  printf '%s\n' "$tint" > "$dir/tint.txt"
  awk -v fix="$fix" -v segs="$dir/segs.txt" -v tint="$dir/tint.txt" '
    /<levels>/ { inlev = 1; print; next }
    /<\/levels>/ { inlev = 0 }
    inlev && fix && NF > 12 { $4 = 1; $11 = 1; print; next }   # fix every energy and width
    /<targetInt>/ { print; while ((getline l < tint) > 0) if (l != "") print l; next }
    /<segmentsData>/ { print; while ((getline l < segs) > 0) print l; skip = 1; next }
    /<\/segmentsData>/ { skip = 0 }
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
fails_with() {  # NAME MESSAGE_FRAGMENT SEGMENT_LINE...
  local name="$1" msg="$2"; shift 2
  mk "$name" 0 "" "$@"
  (cd "$WORK/$name" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  if [ -s "$WORK/$name/output/chiSquared.out" ]; then bad "$name: ran, should have stopped"
  elif grep -qF "$msg" "$WORK/$name/run.log"; then ok "$name: stops with \"$msg\""
  else bad "$name: no \"$msg\" in the log"; tail -5 "$WORK/$name/run.log"; fi
}

# ---- 1. sqrt(E): shiftfn "a0*sqrt(x)" == sqrtshift ---------------------------
TINT='1 "1,2,3" 20 0 0 1 1e18 "a0" 1 1e-20 0 0 0 "" 0'
for t in plain target; do
  tint=""; [ "$t" = target ] && tint="$TINT"
  mk "sqrt_$t" 0 "$tint" \
    "1 1 1 0 100 40 40 4 1 0 0 0.004 0 0 data/xs_40.dat 0 0 sqrtshift 0.003 0 0" \
    "1 1 1 0 100 40 40 4 1 0 0 0.004 0 0 data/xs_40.dat 0 0 shiftfn \"a0*sqrt(x)\" 1 0.003 0 0" \
    "1 1 1 0 100 40 40 4 1 0 0 0.0079115 0 0 data/xs_40.dat 0 0 shiftfn \"a0*(sqrt(x)-sqrt(1.7))\" 1 0.003 0 0"
  # 0.0079115 = 0.004 + 0.003*sqrt(1.7) to 7 digits; segment 3 is checked to 1e-6 only.
  if run "sqrt_$t" 1 ""; then
    d12="$(blockdiff "$WORK/sqrt_$t" 1 2)"; d13="$(blockdiff "$WORK/sqrt_$t" 1 3)"
    if awk -v d="$d12" 'BEGIN { exit !(d != "mismatch" && d + 0 <= 1e-10) }'; then ok "$t: a0*sqrt(x) = sqrtshift, model values to $d12"; else bad "$t: a0*sqrt(x) vs sqrtshift differ by $d12"; fi
    if awk -v d="$d13" 'BEGIN { exit !(d != "mismatch" && d + 0 <= 1e-6) }'; then ok "$t: centred form agrees ($d13)"; else bad "$t: centred form differs by $d13"; fi
    s1="$(seg "$WORK/sqrt_$t" 1)"
    if ! relclose "$s1" 24 1e-3; then ok "$t: the shift is really applied (chi2 $s1, unshifted 24)"; else bad "$t: chi2 $s1 = the unshifted 24"; fi
  fi
done
if [ -n "$(seg "$WORK/sqrt_plain" 1)" ] && ! relclose "$(seg "$WORK/sqrt_target" 1)" "$(seg "$WORK/sqrt_plain" 1)" 1e-3; then
  ok "target integration really on (chi2 $(seg "$WORK/sqrt_target" 1) vs $(seg "$WORK/sqrt_plain" 1))"
else bad "target integration made no difference"; fi

# ---- 2. a0*x^1.5 == the same data written at the shifted energies -------------
mk x15 0 "" \
  "1 1 1 0 100 40 40 4 1 0 0 0.002 0 0 data/xs_40.dat 0 0 shiftfn \"a0*x^1.5\" 1 0.003 0 0" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40_x15_preshifted.dat 0 0"
if run x15 1 ""; then
  d="$(blockdiff "$WORK/x15" 1 2)"
  if awk -v d="$d" 'BEGIN { exit !(d != "mismatch" && d + 0 <= 1e-8) }'; then ok "a0*x^1.5 = pre-shifted data ($d)"; else bad "a0*x^1.5 vs pre-shifted data: $d"; fi
fi

# ---- 3. a follower uses its leader's function ---------------------------------
mk follower 0 "" \
  "1 1 1 0 100 40 40 4 1 0 0 0.001 0 0 data/xs_40.dat 0 0 shiftgroup G shiftfn \"a0*x\" 1 0.003 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0 shiftgroup G" \
  "1 1 1 0 100 90 90 4 1 0 0 0.001 0 0 data/xs_90.dat 0 0 shiftfn \"a0*x\" 1 0.003 0 0"
if run follower 1 ""; then
  d="$(blockdiff "$WORK/follower" 2 3)"
  if awk -v d="$d" 'BEGIN { exit !(d != "mismatch" && d + 0 <= 1e-10) }'; then ok "follower = the same function given directly ($d)"; else bad "follower vs direct: $d"; fi
  if grep -q "shift function on segment #1 (group 'G', segments #1 2)" "$WORK/follower/run.log"; then ok "function and group noted in the log"; else bad "no shift-function note in the log"; fi
fi

# ---- 4. fit a shared linear calibration error -----------------------------------
# Both files carry labels L with L + a + a0*(L - 1) = E, a = -0.004 (the shift at
# 1 MeV), a0 = 0.006, and the model's exact values as yields: chi2 = 0 at the truth.
# With Minuit's numerical gradient the chi2 is only reported: on exact data it is
# very sensitive to where MIGRAD stops (macOS CI 2026-10-21: chi2 808068 with a and
# a0 within 5e-8 of the truth, then 0 on a rerun), so a and a0 decide there.
GROUP=("1 1 1 0 100 40 40 4 1 0 0 0 1 1 data/xs_40_lin_mislabelled.dat 0 0 shiftgroup G shiftfn \"a0*(x-1)\" 1 0 1 1"
       "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90_lin_mislabelled.dat 0 0 shiftgroup G")
fit_linear() {  # NAME FLAGS [CHI2: strict (chi2 < 1e-3 required, default) | info]
  mk "$1" 1 "" "${GROUP[@]}"
  run "$1" 2 "$2" || return
  local dir="$WORK/$1" a a0 t f
  a="$(par "$dir" segment_1_energy_shift)"; a0="$(par "$dir" segment_1_energy_shift_fn0)"; t="$(total "$dir")"
  if absclose "$a" -0.004 1e-5 && absclose "$a0" 0.006 1e-5; then ok "$1: a = $a, a0 = $a0 (truth -0.004, 0.006)"; else bad "$1: a = $a, a0 = $a0, expected -0.004, 0.006"; fi
  if [ "${3:-strict}" = info ]; then
    echo "  info  $1: chi2 $t (not checked: a and a0 above decide)"
  elif awk -v t="$t" 'BEGIN { exit !(t + 0 < 1e-3) }'; then ok "$1: chi2 = 0 at the truth (chi2 $t)"; else bad "$1: chi2 $t, expected ~0"; fi
  f="$(field "$dir" Total-Shift-Chi-Squared:)"
  if [ -n "$f" ]; then ok "$1: chiSquared.out reports Total-Shift-Chi-Squared ($f)"; else bad "$1: no Total-Shift-Chi-Squared"; fi
  if grep -q "^# 1,1 2,a0\*(x-1)," "$dir/output/shifts.out"; then ok "$1: shifts.out has the function line"; else bad "$1: no function line in shifts.out"; fi
}
fit_linear fit_numerical "" info
fit_linear fit_analytic "--use-gradient"
fit_linear fit_lm "--use-lm"

# ---- 5. group-local shift gradient == full-dataset finite difference ------------
mk gradcheck 1 "" \
  "1 1 1 0 100 40 40 4 1 0 0 0.001 1 1 data/xs_40_lin_mislabelled.dat 0 0 shiftgroup G shiftfn \"a0*(x-1)+a1*(x-1)^2\" 2 0.001 1 1 0.0002 1 1" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90_lin_mislabelled.dat 0 0 shiftgroup G"
# Adjoint path: each shift column against a full-dataset difference (AZURE_SHIFT_GRAD_CHECK).
(cd "$WORK/gradcheck" && printf '2\n\n\n\n7\n' | AZURE_GRAD_ADJOINT=1 AZURE_SHIFT_GRAD_CHECK=1 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > adjoint.log 2>&1)
nchk="$(grep -c 'shift-grad-check' "$WORK/gradcheck/adjoint.log")"
worst="$(grep -o 'reldiff=[-0-9.e+]*' "$WORK/gradcheck/adjoint.log" | cut -d= -f2 | sort -g | tail -1)"
if [ "$nchk" -ge 3 ] && awk -v w="$worst" 'BEGIN { exit !(w + 0 < 1e-6) }'; then ok "adjoint path: group-local shift columns = full differences ($nchk checks, worst $worst)"
else bad "adjoint-path shift gradient check: $nchk checks, worst reldiff '$worst'"; fi
# Jacobian path (the --use-gradient default): AZURE_GRAD_CHECK compares every free
# component (here a, a0, a1) with a finite difference of the whole chi2.
(cd "$WORK/gradcheck" && printf '2\n\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > jacobian.log 2>&1)
if grep -q "Jacobian path: shift 0/1 shift_fn 0/2" "$WORK/gradcheck/jacobian.log"; then ok "Jacobian path: a, a0, a1 agree with finite differences ($(grep -o 'worst relative difference [-0-9.e+]*' "$WORK/gradcheck/jacobian.log"))"
else bad "Jacobian-path gradient check: $(grep 'components with relative' "$WORK/gradcheck/jacobian.log")"; fi

# ---- 6. invalid input ---------------------------------------------------------
fails_with follower_fn "is not the first segment of shift group 'G'" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftgroup G" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0 shiftgroup G shiftfn \"a0*x\" 1 0 0 0"
fails_with with_sqrt "two descriptions of one energy calibration" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 sqrtshift 0.001 0 0 shiftfn \"a0*x\" 1 0 0 0"
fails_with undeclared "uses a1 but declares 1 parameter" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftfn \"a0*x+a1\" 1 0 0 0"
fails_with unused "declares a1 but never uses it" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftfn \"a0*x\" 2 0 0 0 0 0 0"
fails_with short_block "declares 2 parameters but gives" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftfn \"a0*x+a1\" 2 0 0 0"
fails_with no_effect "has no effect at this group's energies" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftfn \"a0*0*x\" 1 0 0 1"
mk degenerate 0 "" "1 1 1 0 100 40 40 4 1 0 0 0 0.01 1 data/xs_40.dat 0 0 shiftfn \"a0+0*x\" 1 0 0.01 1"
(cd "$WORK/degenerate" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
if grep -q "WARNING: shift function on segment #1: its free terms cannot be told apart" "$WORK/degenerate/run.log"; then ok "a constant function next to a free a: degeneracy warning"
else bad "no degeneracy warning for a0 + 0*x with a free"; fi

# ---- 7. a .sav value does not free a fixed function parameter -------------------
mk sav_fixed 1 "" "1 1 1 0 100 40 40 4 1 0 0 0 1 1 data/xs_40_lin_mislabelled.dat 0 0 shiftfn \"a0*(x-1)\" 1 0 0 0"
if run sav_fixed 2 ""; then
  # POSIX BRE only: BSD sed (macOS) reads "\+" as a literal "+", which left the .sav unedited.
  sed 's/^\(segment_1_energy_shift_fn0  *\)[-0-9.e+]*/\10.002/' "$WORK/sav_fixed/output/param.sav" > "$WORK/sav_fixed/start.sav"
  if ! absclose "$(par_file "$WORK/sav_fixed/start.sav" segment_1_energy_shift_fn0)" 0.002 1e-12; then
    bad "test setup: start.sav does not carry a0 = 0.002"
  elif run sav_fixed 2 "" start.sav; then
    v="$(par "$WORK/sav_fixed" segment_1_energy_shift_fn0)"
    if absclose "$v" 0.002 1e-12; then ok "fixed a0 stays at the .sav value ($v)"; else bad "fixed a0 moved to $v (.sav 0.002)"; fi
  fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: shift functions"; else echo "FAIL: $failures check(s)"; exit 1; fi
