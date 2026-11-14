#!/usr/bin/env bash
#
# Normalization groups ("normgroup <label>" at the end of a segment line): the
# segments carrying one label share ONE normalization, that of the group's first
# segment (its value, vary flag and penalty).
#
#   1. a follower takes its leader's norm, whatever its own line says: chi2
#      identical to the same segment given that norm directly;
#   2. a fit recovers a common norm from two segments scaled by the same factor,
#      and writes the follower's (fixed) parameter as the shared value -- with
#      Minuit's numerical gradient, --use-gradient and --use-lm;
#   3. two segments scaled by DIFFERENT factors get one compromise norm,
#      strictly between the two (independent norms would reach both);
#   4. the group's penalty is counted once, on the leader, although the
#      follower's own line carries the same uncertainty;
#   5. the analytic norm gradient (adjoint and Jacobian paths) agrees with a
#      finite difference of the whole chi2.
#
#   ./tests/norm_groups/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/norm_groups.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
seg() { awk -F, -v k="$2" '$1 == k { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
segnormchi() { awk -F, -v k="$2" '$1 == k { print $5 }' "$1/output/chiSquared.out" 2>/dev/null; }
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
normtotal() { awk '/^Total-Chi-Squared:/ { print $4 }' "$1/output/chiSquared.out" 2>/dev/null; }
par() { awk -v n="$2" '$1 == n { print $2 }' "$1/output/param.sav" 2>/dev/null; }
# relclose A B TOL -> 0 when |A-B| <= TOL*max(|B|,1e-300)
relclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m == 0) m = 1e-300; exit !(d <= t * m) }'; }
absclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d < t) }'; }

# mk NAME LEVELS_FIXED RECORD... -> $WORK/NAME/t.azr with those records as its <data> section (format 3)
mk() {
  local name="$1" fix="$2"; shift 2
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/data" "$dir/data"
  printf '%s\n' "$@" > "$dir/segs.txt"
  awk -v fix="$fix" -v segs="$dir/segs.txt" '
    /<levels>/ { inlev = 1; print; next }
    /<\/levels>/ { inlev = 0 }
    inlev && fix && NF > 12 { $4 = 1; $11 = 1; print; next }   # fix every energy and width
    /^<data>$/ { print; while ((getline l < segs) > 0) print l; skip = 1; next }
    /^<\/data>$/ { skip = 0 }
    skip { next }
    { print }' "$HERE/norm_groups.azr" > "$dir/t.azr"
}
# ds FILE-STEM ANGLE NORM [SHIFT]: a data set of 1H(p,p) dsigma/dOmega at one c.m. angle (deg), named after its file;
# NORM a nuisance ("1.2 fixed") or the name of a normalization record; SHIFT the value of a shift field
ds() {
  printf 'dataset %s\n    reaction       1H(p,p)1H\n    quantity       dsigma_domega\n    unit           b/sr\n' "$1"
  printf '    file           data/%s.dat\n    columns        E_lab[MeV] theta_cm[deg] value error\n' "$1"
  printf '    energy_window  0 100 MeV\n    angle_window   %s %s deg\n    norm           %s' "$2" "$2" "$3"
  [ -n "${4:-}" ] && printf '\n    shift          %s' "$4"
  printf '\n'
}
# nr LABEL VALUE: a normalization record shared by the data sets that name it (format 3; old "normgroup LABEL", the
# value of the group's first line)
nr() { printf 'normalization %s\n    value          %s\n' "$1" "$2"; }
run() {  # NAME MODE EXTRA_FLAGS
  local input
  if [ "$2" = 1 ]; then input='1\n\n\n7\n'; else input='2\n\n\n\n7\n'; fi
  (cd "$WORK/$1" && printf "$input" | "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -s "$WORK/$1/output/chiSquared.out" ]; then
    bad "$1: no chiSquared.out"; tail -15 "$WORK/$1/run.log"; return 1
  fi
}

# Every data file is the model's exact cross section times a factor (0.8 or
# 0.9), with 1% errors, so a norm of 1/factor is the chi2 = 0 point.

# ---- 1. a follower takes the leader's norm -----------------------------------
# Segment 2's own norm (1, and "free with 10%") is ignored: it uses segment 1's 1.2.
mk group 0 "$(nr D "1.2 fixed")" "$(ds xs_40_s08 40 D)" "$(ds xs_90_s09 90 D)"
mk direct 0 "$(ds xs_40_s08 40 "1.2 fixed")" "$(ds xs_90_s09 90 "1.2 fixed")"
mk plain 0 "$(ds xs_40_s08 40 "1.2 fixed")" "$(ds xs_90_s09 90 "1 fixed")"
if run group 1 "" && run direct 1 "" && run plain 1 ""; then
  g1="$(seg "$WORK/group" 1)"; g2="$(seg "$WORK/group" 2)"; d1="$(seg "$WORK/direct" 1)"; d2="$(seg "$WORK/direct" 2)"; p2="$(seg "$WORK/plain" 2)"
  if relclose "$g1" "$d1" 1e-10; then ok "leader keeps its own norm (chi2 $g1 = $d1)"; else bad "leader: grouped $g1 != direct $d1"; fi
  if relclose "$g2" "$d2" 1e-10; then ok "follower uses the leader's norm (chi2 $g2 = $d2)"; else bad "follower: grouped $g2 != direct $d2"; fi
  if relclose "$g2" "$p2" 1e-3; then bad "follower at its own norm (chi2 $g2 = $p2)"; else ok "follower really rescaled (chi2 $g2 vs own norm $p2)"; fi
  if [ "$(normtotal "$WORK/group")" = 0 ]; then ok "no penalty from the follower's own 10%"; else bad "norm penalty $(normtotal "$WORK/group"), expected 0"; fi
  if grep -q "normalization group 'D': segments #1 + 2 share one normalization (segment #1's, fixed)" "$WORK/group/run.log"; then ok "group noted in the log"; else bad "no normalization-group note in the log"; fi
fi

# ---- 1b. a fixed leader keeps the group fixed ----------------------------------
# The group's value fixed at 1.2: no free norm may appear.
# Only the leader's shift is free (the data want 0), so the fit has one parameter.
mk fixedleader 1 "$(nr D "1.2 fixed")" "$(ds xs_40_s08 40 D "offset 0.001 MeV fit prior 0 0.01")" "$(ds xs_90_s09 90 D)"
if run fixedleader 2 ""; then
  # format 3 states the group's value once: a fixed group has no norm parameter at all (the old file kept the
  # follower's own "free" flag as a fixed follower row at the leader's 1.2; that flag cannot be written now)
  nrows="$(awk '$1 ~ /^D[.]value/' "$WORK/fixedleader/output/param.sav" 2>/dev/null | wc -l | tr -d ' ')"
  sh="$(par "$WORK/fixedleader" xs_40_s08.shift.offset)"
  if [ -s "$WORK/fixedleader/output/param.sav" ] && [ "$nrows" = 0 ] && [ -n "$sh" ]; then ok "fixed leader: no norm parameter of the group, the leader's shift fitted ($sh)"
  else bad "fixed leader: $nrows rows of D.value in param.sav, shift '$sh'"; fi
  if grep -q "share one normalization (segment #1's, fixed)" "$WORK/fixedleader/run.log"; then ok "fixed leader: the log calls the group fixed"; else bad "fixed leader: no 'fixed' group note"; fi
fi

# ---- 2. fit one common norm --------------------------------------------------
# Both files are the model times 0.8: n = 1.25 is the chi2 = 0 minimum.  The
# follower's own line says "free, 10%": ignored (its parameter is fixed).
fit_common() {  # NAME EXTRA_FLAGS
  mk "$1" 1 "$(nr F "1 fit")" "$(ds xs_40_s08 40 F)" "$(ds xs_90_s08 90 F)"
  run "$1" 2 "$2" || return
  local n1 n2 t dir="$WORK/$1"
  n1="$(par "$dir" F.value)"; n2="$(par "$dir" "F.value[xs_90_s08]")"; t="$(total "$dir")"
  if absclose "$n1" 1.25 1e-4; then ok "$1: shared norm n = $n1 (truth 1.25)"; else bad "$1: shared norm n = $n1, expected 1.25"; fi
  if relclose "$n2" "$n1" 1e-12; then ok "$1: follower written with the shared value ($n2)"; else bad "$1: follower written as $n2, leader $n1"; fi
  if awk -v t="$t" 'BEGIN { exit !(t + 0 < 1e-3) }'; then ok "$1: chi2 = 0 at the truth (chi2 $t)"; else bad "$1: chi2 $t, expected ~0"; fi
}
fit_common fit_numerical ""
fit_common fit_analytic "--use-gradient"
fit_common fit_lm "--use-lm"

# ---- 3. different scales -> one compromise -----------------------------------
mk compromise 1 "$(nr F "1 fit")" "$(ds xs_40_s08 40 F)" "$(ds xs_90_s09 90 F)"
if run compromise 2 "--use-gradient"; then
  n1="$(par "$WORK/compromise" F.value)"; n2="$(par "$WORK/compromise" "F.value[xs_90_s09]")"
  if awk -v n="$n1" 'BEGIN { exit !(n > 1.1112 && n < 1.2499) }' && relclose "$n2" "$n1" 1e-12; then
    ok "one norm for both, between the two truths (n = $n1; truths 1.111 / 1.25)"
  else
    bad "expected one norm strictly between 1.111 and 1.25, got $n1 / $n2"
  fi
fi

# ---- 4. the penalty is counted once ------------------------------------------
# Leader and follower both say "1 +- 10%"; data want 1.25, so the penalty pulls.
# Grouped: chi2 + ((n-1)/0.1)^2 once.  Total-Norm-Chi-Squared must be that one
# term, and the follower's own row must show no penalty.
mk penalty 1 "$(nr F "1 fit prior 1 10%")" "$(ds xs_40_s08 40 F)" "$(ds xs_90_s08 90 F)"
if run penalty 2 "--use-gradient"; then
  n="$(par "$WORK/penalty" F.value)"; pt="$(normtotal "$WORK/penalty")"; p2="$(segnormchi "$WORK/penalty" 2)"
  want="$(awk -v n="$n" 'BEGIN { printf "%.10g", ((n - 1) / 0.1) ^ 2 }')"
  if relclose "$pt" "$want" 1e-4 && [ "$p2" = 0 ]; then ok "one penalty ((n-1)/0.1)^2 = $pt at n = $n; follower row 0"
  else bad "norm penalty $pt (follower row $p2), expected one term $want at n = $n"; fi
  if awk -v n="$n" 'BEGIN { exit !(n > 1.0 && n < 1.25) }'; then ok "penalty pulls the shared norm below the data's 1.25 (n = $n)"; else bad "n = $n not between 1 and 1.25"; fi
fi

# ---- 5. analytic norm gradient == finite difference ----------------------------
# Leader's norm free, follower's fixed: the follower's points must enter the
# leader's column.  The follower's shift is free too, so the check is not
# about a single component.
mk gradcheck 1 "$(nr G "1.1 fit prior 1.1 10%")" "$(ds xs_40_s08 40 G)" "$(ds xs_90_s09 90 G "offset 0.001 MeV fit prior 0 0.01")"
(cd "$WORK/gradcheck" && printf '2\n\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > check.log 2>&1)
line="$(grep 'components with relative' "$WORK/gradcheck/check.log")"
if grep -q "norm 0/1 .*Jacobian path: norm 0/1" <<< "$line"; then ok "leader's norm: adjoint and Jacobian gradients = finite difference"
else bad "norm gradient check: $line"; fi
if grep -q "Jacobian-path vs adjoint gradient over 2 free parameters" "$WORK/gradcheck/check.log"; then ok "two free parameters (the follower's norm is fixed)"
else bad "free-parameter count: $(grep 'Jacobian-path vs adjoint' "$WORK/gradcheck/check.log")"; fi

if [ "$failures" -eq 0 ]; then echo "PASS: normalization groups"; else echo "FAIL: $failures check(s)"; exit 1; fi
