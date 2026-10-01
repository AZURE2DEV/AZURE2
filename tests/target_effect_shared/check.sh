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
# Model: tests/target_effect_ranges (two Deineko segments sharing "1,2").
#
#   ./tests/target_effect_shared/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../target_effect_ranges"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/te_shared.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# make_variant NAME SPLIT: SPLIT=1 replaces the "1,2" line by identical "1" and
# "2" lines.  The energy-range tokens are dropped so the whole segment is
# convolved and the comparison does not depend on the automatic decision.
make_variant() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$SRC/data" "$dir/data"
  awk -v sep="$2" '
    /<targetInt>/ { inte = 1; print; next }
    /<\/targetInt>/ { inte = 0 }
    inte && $2 == "\"1,2\"" {
      NF = NF - 3
      if (sep) { $2 = "\"1\""; print; $2 = "\"2\""; print; next }
    }
    { print }' "$SRC/target_effect_ranges.azr" > "$dir/t.azr"
}

run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
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
