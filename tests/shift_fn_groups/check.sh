#!/usr/bin/env bash
#
# Shift-function groups (an energy_scale record with a function and no offset, named by several data sets): every
# data set naming it applies the record's function with the record's parameters, while each keeps its OWN constant
# energy shift ("shift offset ..." of the data set; a record with an offset shares the constant as well).  Model and
# data: tests/shift_fn.
#
#   1. a member applies the leader's function: the same as declaring that function
#      directly, with the member's own constant (not the leader's);
#   2. a fit recovers a shared a0 and two different constants, a1 = -0.004 and
#      a2 = -0.007018 (the 90 deg labels moved by +3 keV), with Minuit's numerical
#      gradient, --use-gradient and --use-lm, with one a0 parameter (the leader's) in
#      param.sav and the shifts.out function line listing both members;
#   3. the analytic gradient of the member's constant and of the shared a0 agrees with
#      finite differences (AZURE_GRAD_CHECK, Jacobian path; AZURE_SHIFT_GRAD_CHECK,
#      adjoint path);
#   4. invalid input of the old grammar is refused by the converter (the old lines through AZURE2 --convert, the
#      original tests/shift_fn from tests/convert/fixtures): a segment in a shiftgroup and a shiftfngroup, a first
#      member without a shiftfn, a later member declaring one;
#   5. a composite (SUM) member: the same as the composite declaring the function
#      directly, and the component is included (components of a shifted segment share its shift).
#
#   ./tests/shift_fn_groups/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL="$HERE/../shift_fn"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/shift_fn_groups.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
par() { awk -v n="$2" '$1 == n { print $2 }' "$1/output/param.sav" 2>/dev/null; }
absclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d < t) }'; }
blockdiff() {  # DIR SEG_A SEG_B: largest relative difference of the model columns
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
# fields; COMP=1 makes it the sum of itself and 0.5 x the 40 deg elastic cross section
ds() {
  local name="$1" file="$2" angle="$3" reaction="1H(p,p)1H"; shift 3
  [ "${COMP:-0}" = 1 ] && reaction="1H(p,p)1H + 0.5*1H(p,p)1H @ 40 deg cm"
  printf 'dataset %s\n    reaction       %s\n    quantity       dsigma_domega\n    unit           b/sr\n' "$name" "$reaction"
  printf '    file           data/%s\n    columns        E_lab[MeV] theta_cm[deg] value error\n' "$file"
  printf '    energy_window  0 100 MeV\n    angle_window   %s %s deg\n    norm           1 fixed\n' "$angle" "$angle"
  local f; for f in "$@"; do printf '    %s\n' "$f"; done
  echo
}
# scale NAME FUNCTION PARAM...: an energy_scale record with a function and no offset; each PARAM "a0 0.003 fixed"
scale() {
  local name="$1" fn="$2"; shift 2
  printf 'energy_scale %s\n    function       "%s" MeV\n' "$name" "$fn"
  local p; for p in "$@"; do printf '    param          %s\n' "$p"; done
  echo
}
# mk NAME LEVELS_FIXED RECORDS -> $WORK/NAME/t.azr (tests/shift_fn's model with RECORDS as its <data> section)
mk() {
  local name="$1" fix="$2" dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$MODEL/data" "$dir/data"
  # The 90 deg mislabelled data with every label moved by +3 keV: its constant becomes
  # a - 0.003*(1 + a0) = -0.007018 with the shared a0 = 0.006 unchanged.
  awk 'NF >= 4 { $1 = sprintf("%.12f", $1 + 0.003) } { print }' OFS='\t' \
    "$MODEL/data/xs_90_lin_mislabelled.dat" > "$dir/data/xs_90_lin_mislabelled_p3.dat"
  # the records go through the environment: BSD awk (macOS) refuses a newline in -v ("newline in string")
  RECORDS="$3" awk -v fix="$fix" '
    /<levels>/ { inlev = 1; print; next }
    /<\/levels>/ { inlev = 0 }
    inlev && fix && NF > 12 { $4 = 1; $11 = 1; print; next }   # fix every energy and width
    /^<data>/ { print; printf "%s\n", ENVIRON["RECORDS"]; skip = 1; next }
    /^<\/data>/ { skip = 0 }
    skip { next }
    { print }' "$MODEL/shift_fn.azr" > "$dir/t.azr"
}
run() {  # NAME MODE EXTRA_FLAGS
  local input
  if [ "$2" = 1 ]; then input='1\n\n\n7\n'; else input='2\n\n\n\n7\n'; fi
  (cd "$WORK/$1" && printf "$input" | "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -s "$WORK/$1/output/chiSquared.out" ]; then
    bad "$1: no chiSquared.out"; tail -15 "$WORK/$1/run.log"; return 1
  fi
}
OLD="$HERE/../convert/fixtures/tests/shift_fn/shift_fn.azr"
convert_refuses() {  # NAME MESSAGE_FRAGMENT SEGMENT_LINE...: the original with these lines is refused by the converter
  local name="$1" msg="$2" dir="$WORK/$1"; shift 2
  mkdir -p "$dir"; cp -R "$MODEL/data" "$dir/data"
  SEGS="$(printf '%s\n' "$@")" awk '
    /<segmentsData>/ { print; print ENVIRON["SEGS"]; skip = 1; next }
    /<\/segmentsData>/ { skip = 0 }
    skip { next }
    { print }' "$OLD" > "$dir/old.azr"
  (cd "$dir" && "$AZURE2_BIN" --convert old.azr t.azr --to 3 > convert.log 2>&1)
  if [ -f "$dir/t.azr" ]; then bad "$name: converted, the old engine refused it"
  elif grep -qF "$msg" "$dir/convert.log"; then ok "$name: the converter refuses it (\"$msg\")"
  else bad "$name: no \"$msg\" in the converter's log"; tail -3 "$dir/convert.log"; fi
}

# ---- 1. a member applies the leader's function with its own constant ------------
mk member 0 "$(scale F 'a0*x' 'a0 0.003 fixed'
                scale H 'a0*x' 'a0 0.003 fixed'
                ds xs_40 xs_40.dat 40 'shift          offset 0.001 MeV fixed' 'shift          F'
                ds xs_90 xs_90.dat 90 'shift          offset 0.002 MeV fixed' 'shift          F'
                ds xs_90_2 xs_90.dat 90 'shift          offset 0.002 MeV fixed' 'shift          H'
                ds xs_90_3 xs_90.dat 90 'shift          offset 0.001 MeV fixed' 'shift          H')"
if run member 1 ""; then
  d23="$(blockdiff "$WORK/member" 2 3)"; d24="$(blockdiff "$WORK/member" 2 4)"
  if awk -v d="$d23" 'BEGIN { exit !(d ~ /^ *[-+]?\.?[0-9]/ && d != "mismatch" && d + 0 <= 1e-10) }'; then ok "member = the function given directly with its own constant ($d23)"; else bad "member vs direct: $d23"; fi
  if awk -v d="$d24" 'BEGIN { exit !(d != "mismatch" && d + 0 > 1e-6) }'; then ok "member does not take the leader's constant ($d24 from it)"; else bad "member uses the leader's constant ($d24)"; fi
  if grep -q "shift function on segment #1 (shift-function group 'F', segments #1 2; each keeps its own constant shift)" "$WORK/member/run.log"; then ok "group noted in the log"; else bad "no shift-function-group note in the log"; grep -i "shift function" "$WORK/member/run.log" | sed -n '1,3p'; fi
fi

# ---- 2. fit a shared a0 with two different constants --------------------------------
GROUP="$(scale F 'a0*(x-1)' 'a0 0 fit prior 0 1'
         ds xs_40_lin_mislabelled xs_40_lin_mislabelled.dat 40 'shift          offset 0 MeV fit prior 0 1' 'shift          F'
         ds xs_90_lin_mislabelled_p3 xs_90_lin_mislabelled_p3.dat 90 'shift          offset 0 MeV fit prior 0 1' 'shift          F')"
fit() {  # NAME FLAGS [strict|info]
  mk "$1" 1 "$GROUP"
  run "$1" 2 "$2" || return
  local dir="$WORK/$1" a1 a2 a0 t
  a1="$(par "$dir" xs_40_lin_mislabelled.shift.offset)"; a2="$(par "$dir" xs_90_lin_mislabelled_p3.shift.offset)"
  a0="$(par "$dir" F.a0)"; t="$(total "$dir")"
  if absclose "$a1" -0.004 1e-5 && absclose "$a2" -0.007018 1e-5 && absclose "$a0" 0.006 1e-5; then
    ok "$1: a1 = $a1, a2 = $a2, a0 = $a0 (truth -0.004, -0.007018, 0.006)"
  else bad "$1: a1 = $a1, a2 = $a2, a0 = $a0, expected -0.004, -0.007018, 0.006"; fi
  if [ "${3:-strict}" = info ]; then echo "  info  $1: chi2 $t (not checked: the parameters above decide)"
  elif awk -v t="$t" 'BEGIN { exit !(t ~ /^ *[-+]?\.?[0-9]/ && t + 0 < 1e-3) }'; then ok "$1: chi2 = 0 at the truth ($t)"; else bad "$1: chi2 $t, expected ~0"; fi
  if [ "$(awk '$1 ~ /a0/' "$dir/output/param.sav" | wc -l | tr -d ' ')" = 1 ]; then ok "$1: one a0 parameter, the record's"; else bad "$1: a0 rows: $(awk '$1 ~ /a0/ { printf "%s ", $1 }' "$dir/output/param.sav")"; fi
  if grep -q "^# 1,1 2,a0\*(x-1)," "$dir/output/shifts.out"; then ok "$1: shifts.out function line lists both members"; else bad "$1: shifts.out function line: $(grep '^# 1' "$dir/output/shifts.out")"; fi
}
fit fit_numerical "" info
fit fit_analytic "--use-gradient"
fit fit_lm "--use-lm"

# ---- 3. analytic gradient of the member's constant and the shared a0 ----------------
mk gradcheck 1 "$(scale F 'a0*(x-1)+a1*(x-1)^2' 'a0 0.001 fit prior 0 1' 'a1 0.0002 fit prior 0 1'
                  ds xs_40_lin_mislabelled xs_40_lin_mislabelled.dat 40 'shift          offset 0.001 MeV fit prior 0 1' 'shift          F'
                  ds xs_90_lin_mislabelled_p3 xs_90_lin_mislabelled_p3.dat 90 'shift          offset -0.002 MeV fit prior 0 1' 'shift          F')"
(cd "$WORK/gradcheck" && printf '2\n\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > jacobian.log 2>&1)
if grep -q "Jacobian path: shift 0/2 shift_fn 0/2" "$WORK/gradcheck/jacobian.log"; then ok "Jacobian path: a1, a2, a0, a1' agree with finite differences ($(grep -o 'worst relative difference [-0-9.e+]*' "$WORK/gradcheck/jacobian.log"))"
else bad "Jacobian-path gradient check: $(grep -E 'Jacobian path|components with relative' "$WORK/gradcheck/jacobian.log" | sed -n '1,3p')"; fi
(cd "$WORK/gradcheck" && printf '2\n\n\n\n7\n' | AZURE_GRAD_ADJOINT=1 AZURE_SHIFT_GRAD_CHECK=1 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > adjoint.log 2>&1)
nchk="$(grep -c 'shift-grad-check' "$WORK/gradcheck/adjoint.log")"
worst="$(grep -o 'reldiff=[-0-9.e+]*' "$WORK/gradcheck/adjoint.log" | cut -d= -f2 | sort -g | tail -1)"
if [ "$nchk" -ge 4 ] && awk -v w="$worst" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 < 1e-6) }'; then ok "adjoint path: shift columns = full differences ($nchk checks, worst $worst)"
else bad "adjoint-path shift gradient check: $nchk checks, worst reldiff '$worst'"; fi

# ---- 4. invalid input ---------------------------------------------------------------
convert_refuses both_groups "is in shift group 'G' and in shift-function group 'F'" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftgroup G shiftfngroup F shiftfn \"a0*x\" 1 0 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0 shiftgroup G"
convert_refuses no_fn "the first line of shift-function group F declares no shiftfn" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftfngroup F" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0 shiftfngroup F shiftfn \"a0*x\" 1 0 0 0"
convert_refuses two_fns "but it is not the first segment of shift-function group 'F'" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftfngroup F shiftfn \"a0*x\" 1 0 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0 shiftfngroup F shiftfn \"a0*x\" 1 0 0 0"

# ---- 5. a composite member (SUM: 90 deg elastic + 0.5 x 40 deg elastic) --------------
mk composite 0 "$(scale F 'a0*x' 'a0 0.003 fixed'
                   scale H 'a0*x' 'a0 0.003 fixed'
                   ds xs_40 xs_40.dat 40 'shift          offset 0.001 MeV fixed' 'shift          F'
                   COMP=1 ds xs_90 xs_90.dat 90 'shift          offset 0.002 MeV fixed' 'shift          F'
                   COMP=1 ds xs_90_2 xs_90.dat 90 'shift          offset 0.002 MeV fixed' 'shift          H'
                   ds xs_90_3 xs_90.dat 90 'shift          offset 0.002 MeV fixed' 'shift          H')"
if run composite 1 ""; then
  d23="$(blockdiff "$WORK/composite" 2 3)"; d24="$(blockdiff "$WORK/composite" 2 4)"
  if awk -v d="$d23" 'BEGIN { exit !(d ~ /^ *[-+]?\.?[0-9]/ && d != "mismatch" && d + 0 <= 1e-10) }'; then ok "composite member = the composite with the function given directly ($d23)"; else bad "composite member vs direct: $d23"; fi
  if awk -v d="$d24" 'BEGIN { exit !(d != "mismatch" && d + 0 > 1e-6) }'; then ok "the component is included ($d24 from the plain segment)"; else bad "composite equals the plain segment ($d24)"; fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: shift-function groups"; else echo "FAIL: $failures check(s)"; exit 1; fi
