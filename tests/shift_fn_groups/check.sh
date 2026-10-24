#!/usr/bin/env bash
#
# Shift-function groups ("shiftfngroup <label>"): every segment with the label applies
# the shift function (shiftfn) declared on the first of them, with that segment's
# parameters, while each keeps its OWN constant energy shift (a shiftgroup shares the
# constant as well).  Model and data: tests/shift_fn.
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
#   4. invalid input stops with a message: a segment in a shiftgroup and a
#      shiftfngroup, a first member without a shiftfn, a later member declaring one;
#   5. a composite (SUM) member: the same as the composite declaring the function
#      directly, and every output identical with and without the component mapping
#      (AZURE_NO_COMPONENT_MAPPING=1; components of a shifted segment share its shift).
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
trap 'rm -rf "$WORK"' EXIT

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

# mk NAME LEVELS_FIXED SEGMENT_LINE... -> $WORK/NAME/t.azr (tests/shift_fn's model)
mk() {
  local name="$1" fix="$2"; shift 2
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$MODEL/data" "$dir/data"
  # The 90 deg mislabelled data with every label moved by +3 keV: its constant becomes
  # a - 0.003*(1 + a0) = -0.007018 with the shared a0 = 0.006 unchanged.
  awk 'NF >= 4 { $1 = sprintf("%.12f", $1 + 0.003) } { print }' OFS='\t' \
    "$MODEL/data/xs_90_lin_mislabelled.dat" > "$dir/data/xs_90_lin_mislabelled_p3.dat"
  printf '%s\n' "$@" > "$dir/segs.txt"
  awk -v fix="$fix" -v segs="$dir/segs.txt" '
    /<levels>/ { inlev = 1; print; next }
    /<\/levels>/ { inlev = 0 }
    inlev && fix && NF > 12 { $4 = 1; $11 = 1; print; next }   # fix every energy and width
    /<segmentsData>/ { print; while ((getline l < segs) > 0) print l; skip = 1; next }
    /<\/segmentsData>/ { skip = 0 }
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
fails_with() {  # NAME MESSAGE_FRAGMENT SEGMENT_LINE...
  local name="$1" msg="$2"; shift 2
  mk "$name" 0 "$@"
  (cd "$WORK/$name" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  if [ -s "$WORK/$name/output/chiSquared.out" ]; then bad "$name: ran, should have stopped"
  elif grep -qF "$msg" "$WORK/$name/run.log"; then ok "$name: stops with \"$msg\""
  else bad "$name: no \"$msg\" in the log"; tail -5 "$WORK/$name/run.log"; fi
}

# ---- 1. a member applies the leader's function with its own constant ------------
mk member 0 \
  "1 1 1 0 100 40 40 4 1 0 0 0.001 0 0 data/xs_40.dat 0 0 shiftfngroup F shiftfn \"a0*x\" 1 0.003 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0.002 0 0 data/xs_90.dat 0 0 shiftfngroup F" \
  "1 1 1 0 100 90 90 4 1 0 0 0.002 0 0 data/xs_90.dat 0 0 shiftfn \"a0*x\" 1 0.003 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0.001 0 0 data/xs_90.dat 0 0 shiftfn \"a0*x\" 1 0.003 0 0"
if run member 1 ""; then
  d23="$(blockdiff "$WORK/member" 2 3)"; d24="$(blockdiff "$WORK/member" 2 4)"
  if awk -v d="$d23" 'BEGIN { exit !(d != "mismatch" && d + 0 <= 1e-10) }'; then ok "member = the function given directly with its own constant ($d23)"; else bad "member vs direct: $d23"; fi
  if awk -v d="$d24" 'BEGIN { exit !(d != "mismatch" && d + 0 > 1e-6) }'; then ok "member does not take the leader's constant ($d24 from it)"; else bad "member uses the leader's constant ($d24)"; fi
  if grep -q "shift function on segment #1 (shift-function group 'F', segments #1 2; each keeps its own constant shift)" "$WORK/member/run.log"; then ok "group noted in the log"; else bad "no shift-function-group note in the log"; grep -i "shift function" "$WORK/member/run.log" | head -3; fi
fi

# ---- 2. fit a shared a0 with two different constants --------------------------------
GROUP=("1 1 1 0 100 40 40 4 1 0 0 0 1 1 data/xs_40_lin_mislabelled.dat 0 0 shiftfngroup F shiftfn \"a0*(x-1)\" 1 0 1 1"
       "1 1 1 0 100 90 90 4 1 0 0 0 1 1 data/xs_90_lin_mislabelled_p3.dat 0 0 shiftfngroup F")
fit() {  # NAME FLAGS [strict|info]
  mk "$1" 1 "${GROUP[@]}"
  run "$1" 2 "$2" || return
  local dir="$WORK/$1" a1 a2 a0 t
  a1="$(par "$dir" segment_1_energy_shift)"; a2="$(par "$dir" segment_2_energy_shift)"
  a0="$(par "$dir" segment_1_energy_shift_fn0)"; t="$(total "$dir")"
  if absclose "$a1" -0.004 1e-5 && absclose "$a2" -0.007018 1e-5 && absclose "$a0" 0.006 1e-5; then
    ok "$1: a1 = $a1, a2 = $a2, a0 = $a0 (truth -0.004, -0.007018, 0.006)"
  else bad "$1: a1 = $a1, a2 = $a2, a0 = $a0, expected -0.004, -0.007018, 0.006"; fi
  if [ "${3:-strict}" = info ]; then echo "  info  $1: chi2 $t (not checked: the parameters above decide)"
  elif awk -v t="$t" 'BEGIN { exit !(t + 0 < 1e-3) }'; then ok "$1: chi2 = 0 at the truth ($t)"; else bad "$1: chi2 $t, expected ~0"; fi
  if [ -z "$(par "$dir" segment_2_energy_shift_fn0)" ]; then ok "$1: one a0 parameter, on the leader"; else bad "$1: the member has its own a0"; fi
  if grep -q "^# 1,1 2,a0\*(x-1)," "$dir/output/shifts.out"; then ok "$1: shifts.out function line lists both members"; else bad "$1: shifts.out function line: $(grep '^# 1' "$dir/output/shifts.out")"; fi
}
fit fit_numerical "" info
fit fit_analytic "--use-gradient"
fit fit_lm "--use-lm"

# ---- 3. analytic gradient of the member's constant and the shared a0 ----------------
mk gradcheck 1 \
  "1 1 1 0 100 40 40 4 1 0 0 0.001 1 1 data/xs_40_lin_mislabelled.dat 0 0 shiftfngroup F shiftfn \"a0*(x-1)+a1*(x-1)^2\" 2 0.001 1 1 0.0002 1 1" \
  "1 1 1 0 100 90 90 4 1 0 0 -0.002 1 1 data/xs_90_lin_mislabelled_p3.dat 0 0 shiftfngroup F"
(cd "$WORK/gradcheck" && printf '2\n\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > jacobian.log 2>&1)
if grep -q "Jacobian path: shift 0/2 shift_fn 0/2" "$WORK/gradcheck/jacobian.log"; then ok "Jacobian path: a1, a2, a0, a1' agree with finite differences ($(grep -o 'worst relative difference [-0-9.e+]*' "$WORK/gradcheck/jacobian.log"))"
else bad "Jacobian-path gradient check: $(grep -E 'Jacobian path|components with relative' "$WORK/gradcheck/jacobian.log" | head -3)"; fi
(cd "$WORK/gradcheck" && printf '2\n\n\n\n7\n' | AZURE_GRAD_ADJOINT=1 AZURE_SHIFT_GRAD_CHECK=1 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > adjoint.log 2>&1)
nchk="$(grep -c 'shift-grad-check' "$WORK/gradcheck/adjoint.log")"
worst="$(grep -o 'reldiff=[-0-9.e+]*' "$WORK/gradcheck/adjoint.log" | cut -d= -f2 | sort -g | tail -1)"
if [ "$nchk" -ge 4 ] && awk -v w="$worst" 'BEGIN { exit !(w + 0 < 1e-6) }'; then ok "adjoint path: shift columns = full differences ($nchk checks, worst $worst)"
else bad "adjoint-path shift gradient check: $nchk checks, worst reldiff '$worst'"; fi

# ---- 4. invalid input ---------------------------------------------------------------
fails_with both_groups "is in shift group 'G' and in shift-function group 'F'" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftgroup G shiftfngroup F shiftfn \"a0*x\" 1 0 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0 shiftgroup G"
fails_with no_fn "is the first segment of shift-function group 'F' but declares no shiftfn" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftfngroup F" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0 shiftfngroup F shiftfn \"a0*x\" 1 0 0 0"
fails_with two_fns "is not the first segment of shift-function group 'F'" \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0 shiftfngroup F shiftfn \"a0*x\" 1 0 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0 shiftfngroup F shiftfn \"a0*x\" 1 0 0 0"

# ---- 5. a composite member (SUM: 90 deg elastic + 0.5 x 40 deg elastic) --------------
COMP="1 0 -1 1 1 1 40 0.5 0"
mk composite 0 \
  "1 1 1 0 100 40 40 4 1 0 0 0.001 0 0 data/xs_40.dat 0 0 shiftfngroup F shiftfn \"a0*x\" 1 0.003 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0.002 0 0 data/xs_90.dat $COMP shiftfngroup F" \
  "1 1 1 0 100 90 90 4 1 0 0 0.002 0 0 data/xs_90.dat $COMP shiftfn \"a0*x\" 1 0.003 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0.002 0 0 data/xs_90.dat 0 0 shiftfn \"a0*x\" 1 0.003 0 0"
if run composite 1 ""; then
  d23="$(blockdiff "$WORK/composite" 2 3)"; d24="$(blockdiff "$WORK/composite" 2 4)"
  if awk -v d="$d23" 'BEGIN { exit !(d != "mismatch" && d + 0 <= 1e-10) }'; then ok "composite member = the composite with the function given directly ($d23)"; else bad "composite member vs direct: $d23"; fi
  if awk -v d="$d24" 'BEGIN { exit !(d != "mismatch" && d + 0 > 1e-6) }'; then ok "the component is included ($d24 from the plain segment)"; else bad "composite equals the plain segment ($d24)"; fi
  mkdir -p "$WORK/composite_nomap"; cp -R "$WORK/composite/data" "$WORK/composite/t.azr" "$WORK/composite_nomap/"
  mkdir -p "$WORK/composite_nomap/output" "$WORK/composite_nomap/checks"
  (cd "$WORK/composite_nomap" && printf '1\n\n\n7\n' | AZURE_NO_COMPONENT_MAPPING=1 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  nf=0; nd=0
  for f in "$WORK/composite/output/"*.out; do
    nf=$((nf + 1))
    awk 'NR == FNR { x[FNR] = $0; m = FNR; next } { k = FNR; if (x[FNR] != $0) e++ } END { exit !(e == 0 && k == m) }' \
      "$f" "$WORK/composite_nomap/output/$(basename "$f")" 2>/dev/null || nd=$((nd + 1))
  done
  if [ "$nf" -ge 2 ] && [ "$nd" = 0 ]; then ok "composite: all $nf output files identical with and without the component mapping"
  else bad "composite: $nd of $nf output files differ with AZURE_NO_COMPONENT_MAPPING=1"; fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: shift-function groups"; else echo "FAIL: $failures check(s)"; exit 1; fi
