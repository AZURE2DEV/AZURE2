#!/usr/bin/env bash
#
# Explicit prior centres: rows "segment_N_norm prior_centre c" and
# "segment_N_energy_shift prior_centre c" of <parameterSettings>
# (EData::ReadPriorCentres) move the centre of that segment's prior away from
# the norm / shift field of its <segmentsData> line, which stays the start
# value.  Without a row the field is both, as it always was.
#
# Model: tests/15N_p_a, segment 3 (Schardt) has norm 0.9950877877494786 with a
# 15 % prior; calculate with data, so every norm and shift stays at its field.
#   none  classic file                          -> chi2_0, norm term 0
#   same  centre = the field                    -> stdout and chiSquared.out
#                                                  byte-identical to none
#   one   centre 1                              -> chi2_0 + ((n - 1)/0.15)^2
#   shift segment 3 shift free, error 0.01 MeV, centre 0.005 MeV
#                                               -> + ((0 - 0.005)/0.01)^2 = 0.25
#   stale a row for segment 9 (no such segment) -> warning, chi2_0
#   bad   "segment_3_norm prior_centre abc", and a centre 0 for a norm
#                                               -> refused, no chi-squared
# The "Total Chi-Squared:" line on stdout includes the priors; chiSquared.out
# gives the norm term separately (Total-Norm-Chi-Squared).
#
#   ./tests/prior_centre/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../15N_p_a"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/prior_centre.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

# run NAME "ROWS (\n-separated) or empty" [free-shift] -> the run directory
run() {
  local d="$WORK/$1" rows="$2" shift3="${3:-}"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  awk -v rows="$rows" -v shift3="$shift3" '
    { sub(/\r$/, "") }
    /^<segmentsData>/ {seg = 1; n = 0; print; next}
    /^<\/segmentsData>/ {seg = 0}
    seg && NF {
      n++
      if (n == 3 && shift3 != "") { $12 = "0"; $13 = "0.01"; $14 = "1" }
    }
    /^<\/parameterSettings>/ && rows != "" {print rows}
    {print}
  ' "$SRC/15N_p_a.azr" > "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline run.azr > run.log 2>&1)
  echo "$d"
}
total() { grep -oE 'Total Chi-Squared: [0-9.eE+-]+' "$1/run.log" | awk '{print $3}'; }
normchi() { awk '/^Total-Chi-Squared:/ {print $4}' "$1/output/chiSquared.out" 2>/dev/null; }

bad=0
pass() { echo "  ok    $1"; }
fail() { echo "  FAIL  $1"; bad=1; }
close() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN{d=a-b; if(d<0)d=-d; exit !(d<=t)}'; }

grep -q "<parameterSettings>" "$SRC/15N_p_a.azr" || { echo "  FAIL  tests/15N_p_a has no <parameterSettings>"; exit 1; }

d0="$(run none "")";             c0="$(total "$d0")"
ds="$(run same "segment_3_norm prior_centre 0.9950877877494786")"
d1="$(run one "segment_3_norm prior_centre 1")"; c1="$(total "$d1")"
dz="$(run shift0 "" free)";      cz="$(total "$dz")"
dh="$(run shift "segment_3_energy_shift prior_centre 0.005" free)"; ch="$(total "$dh")"
dx="$(run stale "segment_9_norm prior_centre 1")"; cx="$(total "$dx")"
db="$(run bad "segment_3_norm prior_centre abc")"; cb="$(total "$db")"
dn="$(run nonpos "segment_3_norm prior_centre 0")"; cn="$(total "$dn")"

echo "  chi2: classic $c0   centre 1 $c1   free shift $cz   its centre 0.005 $ch   stale $cx"
for v in "$c0" "$c1" "$cz" "$ch" "$cx"; do
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
if [ "$cx" = "$c0" ] && grep -q "names no active data segment" "$dx/run.log"; then
  pass "a row for a missing segment: warned, ignored"
else
  fail "a row for a missing segment: chi2 $cx, warning $(grep -c 'names no active' "$dx/run.log")"
fi
if [ -z "$cb" ] && grep -q "is not \"segment_N_norm prior_centre value\"" "$db/run.log"; then
  pass "a malformed row is refused"
else
  fail "a malformed row was not refused (chi2 '$cb')"
fi
if [ -z "$cn" ] && grep -q "must be positive" "$dn/run.log"; then
  pass "a norm centre 0 is refused"
else
  fail "a norm centre 0 was not refused (chi2 '$cn')"
fi
[ "$bad" -eq 0 ] && echo "  ok    explicit prior centres"
exit "$bad"
