#!/usr/bin/env bash
#
# Explicit prior centres.  Old files had rows "segment_N_norm prior_centre c" and "segment_N_energy_shift prior_centre
# c" of <parameterSettings> (EData::ApplyPriorCentres) that moved the centre of that segment's prior away from the
# norm / shift field of its <segmentsData> line, which stayed the start value.  Format 3 writes the centre in the field
# ("norm <start> fit prior <centre> <width>"); the rows are the converter's (it folds them into the fields).
#
# Model: tests/15N_p_a, data set 3 (schardt) has norm 0.9950877877494786 with a 15 % prior; calculate with data, so
# every norm and shift stays at its field.
#   none  the file as it is                     -> chi2_0, norm term 0
#   same  "prior 0.9950877877494786 15%" written again (centre = the start value)
#                                               -> stdout and chiSquared.out byte-identical to none
#   one   centre 1                              -> chi2_0 + ((n - 1)/0.15)^2
#   shift data set 3's shift free, error 0.01 MeV, centre 0.005 MeV
#                                               -> + ((0 - 0.005)/0.01)^2 = 0.25
#   stale the original with data set 2's line inactive and a row for it, converted: the row is dropped, the chi2 of
#         the same file without the row
#   zero  centre 0 ("norm <n> fit prior 0 15%"): refused by --validate and by the run with E_PRIOR_CENTRE
#   unknown, typo, bad, nonpos: rows of the original naming no line (segment_9_norm), misspelt (segment_3_nrom),
#         malformed (abc) or a norm centre 0: the conversion is refused (a prior that is silently not applied is a
#         wrong fit that looks right)
# The "Total Chi-Squared:" line on stdout includes the priors; chiSquared.out gives the norm term separately
# (Total-Norm-Chi-Squared).
#
#   ./tests/prior_centre/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"
SRC="$HERE/../15N_p_a"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/prior_centre.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

# run NAME "ROWS (\n-separated) or empty" [free-shift] [inactive2] -> the run directory
# run NAME CENTRE [SHIFT-CENTRE] [inactive2]: data set 3's norm prior centred at CENTRE (its 15 % width), its shift
# free (prior width 0.01 MeV) centred at SHIFT-CENTRE when given, data set 2 inactive when the fourth argument is set
run() {
  local d="$WORK/$1" centre="$2" shift3="${3:-}" off2="${4:-}"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  awk -v c="$centre" -v shift3="$shift3" -v off2="$off2" '
    { sub(/\r$/, "") }
    /^dataset / { n++; print; if (n == 2 && off2 != "") print "    active         no"; if (n == 3 && shift3 != "") print "    shift          offset 0 MeV fit prior " shift3 " 0.01"; next }
    n == 3 && $1 == "norm" { $0 = "    norm           " $2 " fit prior " c " 15%" }
    {print}
  ' "$SRC/15N_p_a.azr" > "$d/run.azr"
  (cd "$d" && $RUN "$AZURE2_BIN" calculate --no-gui run.azr 2>&1 | head -c 2000000 > run.log)
  echo "$d"
}
# old NAME ROW [inactive2]: the original (tests/convert/fixtures) with a prior_centre row, segment 2 inactive when the
# third argument is set, converted to $WORK/NAME/run.azr; the converter's log in convert.log
OLD="$HERE/../convert/fixtures/tests/15N_p_a/15N_p_a.azr"
old() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  awk -v row="$2" -v off2="${3:-}" '{ sub(/\r$/, "") }
    /^<segmentsData>/ { seg = 1; n = 0; print; next } /^<\/segmentsData>/ { seg = 0 }
    seg && NF { n++; if (n == 2 && off2 != "") $1 = "0" }
    /^<\/parameterSettings>/ { print row } { print }' "$OLD" > "$d/old.azr"
  (cd "$d" && $RUN "$AZURE2_BIN" --convert old.azr run.azr --to 3 > convert.log 2>&1)
  echo "$d"
}
total() { grep -oE 'Total Chi-Squared: [0-9.eE+-]+' "$1/run.log" | awk '{print $3}'; }
normchi() { awk '/^Total-Chi-Squared:/ {print $4}' "$1/output/chiSquared.out" 2>/dev/null; }

bad=0
pass() { echo "  ok    $1"; }
fail() { echo "  FAIL  $1"; bad=1; }
close() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN{d=a-b; if(d<0)d=-d; exit !(d<=t)}'; }

grep -q "<parameterSettings>" "$SRC/15N_p_a.azr" || { echo "  FAIL  tests/15N_p_a has no <parameterSettings>"; exit 1; }

d0="$(run none 0.9950877877494786)";   c0="$(total "$d0")"
ds="$(run same 0.9950877877494786)"
d1="$(run one 1)"; c1="$(total "$d1")"
dz="$(run shift0 0.9950877877494786 0)"; cz="$(total "$dz")"
dh="$(run shift 0.9950877877494786 0.005)"; ch="$(total "$dh")"
dy="$(run off2 0.9950877877494786 "" inactive2)"; cy="$(total "$dy")"
dx="$(old stale "segment_2_norm prior_centre 1" inactive2)"
(cd "$dx" && $RUN "$AZURE2_BIN" calculate --no-gui run.azr 2>&1 | head -c 2000000 > run.log); cx="$(total "$dx")"
du="$(old unknown "segment_9_norm prior_centre 1")"
dt="$(old typo "segment_3_nrom prior_centre 1")"
db="$(old bad "segment_3_norm prior_centre abc")"
dn="$(old nonpos "segment_3_norm prior_centre 0")"

echo "  chi2: classic $c0   centre 1 $c1   free shift $cz   its centre 0.005 $ch   line 2 inactive $cy, with its row $cx"
for v in "$c0" "$c1" "$cz" "$ch" "$cy" "$cx"; do
  [ -n "$v" ] || { echo "  FAIL  a run produced no chi-squared"; exit 1; }
done

# Byte-identical files without cmp (not in every MSYS2 image); \r stripped for CRLF builds.
same() { [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
if same "$d0/run.log" "$ds/run.log" && same "$d0/output/chiSquared.out" "$ds/output/chiSquared.out"; then
  pass "a centre equal to the field changes nothing (stdout and chiSquared.out byte-identical)"
else
  fail "a centre equal to the field changed the output"
fi
[ "$(normchi "$d0")" = 0 ] && pass "classic: norm term 0" || fail "classic: norm term $(normchi "$d0")"
want="$(awk 'BEGIN{n=0.9950877877494786; printf "%.12g", ((n-1)/0.15)^2}')"
if close "$(awk -v a="$c1" -v z="$c0" 'BEGIN{print a-z}')" "$want" 1e-6; then
  pass "centre 1: the prior term is ((n - 1)/0.15)^2 = $want"
else
  fail "centre 1 added $(awk -v a="$c1" -v z="$c0" 'BEGIN{print a-z}'), expected $want"
fi
close "$(normchi "$d1")" "$want" 1e-5 && pass "and chiSquared.out reports it" || fail "chiSquared.out norm term $(normchi "$d1")"
if close "$(awk -v a="$ch" -v z="$cz" 'BEGIN{print a-z}')" 0.25 1e-6; then
  pass "shift centre 0.005 MeV: + 0.25"
else
  fail "shift centre added $(awk -v a="$ch" -v z="$cz" 'BEGIN{print a-z}'), expected 0.25"
fi
if [ "$cx" = "$cy" ] && grep -q "prior_centre row(s) folded" "$dx/convert.log"; then
  pass "a row for an inactive segment line: dropped by the converter, the chi2 without it"
else
  fail "a row for an inactive segment line: chi2 $cx (without the row $cy), converter: $(tail -1 "$dx/convert.log")"
fi
for c in "unknown:$du:names no parameter a prior_centre row can move" "typo:$dt:cannot be read" "bad:$db:cannot be read" \
         "nonpos:$dn:must be positive"; do
  IFS=: read -r what dir msg <<< "$c"
  if [ ! -f "$dir/run.azr" ] && grep -q "not converted: <parameterSettings> row .*$msg" "$dir/convert.log"; then
    pass "the row of the original ($what) is refused by the converter"
  else
    fail "the row of the original ($what) was not refused ($(grep -m1 ERROR "$dir/convert.log"))"
  fi
done
# zero: the format-3 counterpart of nonpos, "norm <n> fit prior 0 15%": refused by --validate and by the run
# (E_PRIOR_CENTRE; a width relative to a centre 0 means nothing, and a prior silently not applied is a wrong fit)
dzero="$(run zero 0)"
(cd "$dzero" && $RUN "$AZURE2_BIN" --validate run.azr > validate.log 2>&1; echo "exit $?" >> validate.log)
if grep -q "E_PRIOR_CENTRE" "$dzero/validate.log" && ! grep -q "^exit 0$" "$dzero/validate.log"; then
  pass "a fitted norm with a prior centred at 0: refused by --validate (E_PRIOR_CENTRE)"
else
  fail "a prior centred at 0: --validate gave $(tail -1 "$dzero/validate.log"), $(grep -c E_PRIOR_CENTRE "$dzero/validate.log") E_PRIOR_CENTRE"
fi
if grep -q "E_PRIOR_CENTRE" "$dzero/run.log" && [ ! -s "$dzero/output/chiSquared.out" ]; then
  pass "and by the run (E_PRIOR_CENTRE, no chiSquared.out)"
else
  fail "a prior centred at 0 ran: $(grep -c E_PRIOR_CENTRE "$dzero/run.log") E_PRIOR_CENTRE in the log"
fi
[ "$bad" -eq 0 ] && echo "  ok    explicit prior centres"
exit "$bad"
