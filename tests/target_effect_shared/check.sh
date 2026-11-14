#!/usr/bin/env bash
#
# One <targetInt> line shared by several segments must act exactly like one
# line per segment.
#
# Every segment listed on a <targetInt> line (and every component of an
# advanced SUM/RATIO segment) points to the same TargetEffect object.  Until
# 2026-10-01 EData::ReadTargetEffectsFile converted that object's Gaussian
# sigma from lab to c.m. once per segment, so a line naming N segments applied
# sigma * (m_t/(m_p+m_t))^N to all of them, and a SUM segment's component
# narrowed the kernel of its parent (22Ne+a, 10-1-26: Harms "4-6" chi2 2037
# shared vs 1379 split).  TargetEffect::ConvertSigmaToCM now converts once.
#
# Model: tests/beam_spread_ay (two Deineko data sets naming one beam record, beam_targetInt1, as the old
# "1,2" line); the split variant gives the second its own copy of the record.
#
#   ./tests/target_effect_shared/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../beam_spread_ay"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
# A hang fails the check instead of stalling the suite (tests/lib/guard.sh).
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/edit_guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-600}")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/te_shared.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# make_variant NAME SPLIT: the project as it is; SPLIT=1 gives deineko_160 a beam record of its own, beam_split, a copy
# of the shared beam_targetInt1 (the same spread), so that each data set has its own record.  Every edit must change
# the file (tests/lib/edit_guard.sh).
SPLIT='$1 == "dataset" { ds = $2 } /^<\/data>/ { ds = "" }
  ds == "deineko_160" && $1 == "beam" { $0 = "    beam           beam_split" }
  /^<\/targets>/ { print ""; print "beam beam_split"; print "    spread         0.030 MeV lab" }
  { print }'
make_variant() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$SRC/data" "$dir/data"
  cp "$SRC/beam_spread_ay.azr" "$dir/t.azr"
  [ "$2" = 1 ] || return 0
  edited "$dir/t.azr" "$dir/t.azr" awk "$SPLIT" || bad "$1: the beam record was not split"
  [ "$(grep -c '^beam beam_' "$dir/t.azr")" -eq 2 ] && [ "$(grep -c '^    beam           beam_split$' "$dir/t.azr")" -eq 1 ] \
    || bad "$1: not two beam records, one data set on each"
}

run() { (cd "$WORK/$1" && $RUN "$AZURE2_BIN" calculate --no-gui t.azr > run.log 2>&1); }
seg() { awk -F, -v k="$2" '$1 == k { print $2 }' "$WORK/$1/output/chiSquared.out" 2>/dev/null; }

make_variant shared 0
make_variant split 1
run shared
run split

for k in 1 2; do
  a="$(seg shared $k)"; b="$(seg split $k)"
  if [ -z "$a" ] || [ -z "$b" ]; then bad "segment $k: no chi-squared (shared '$a', split '$b')"
  elif [ "$a" = "$b" ]; then ok "segment $k: shared line $a = separate lines $b"
  else bad "segment $k: shared line $a, separate lines $b"; fi
done

[ "$failures" -eq 0 ] && echo "  PASS" || echo "  FAIL ($failures)"
[ "$failures" -eq 0 ]
