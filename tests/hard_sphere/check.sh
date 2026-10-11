#!/usr/bin/env bash
#
# Hard-sphere scattering in every J^pi (includeAllJPi): a model whose only way to
# get every partial wave's potential scattering was a fixed zero-width "dummy"
# level per J^pi must give the same elastic cross sections as the same real
# levels with includeAllJPi=1 and no dummies -- on the A-matrix and the R-matrix
# routes, with and without --use-brune.  And with the option off, dropping the
# dummies must change the result (the partial waves really were missing).
#
#   ./tests/hard_sphere/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/hard_sphere.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MODEL ROUTE(true|false) FLAGS  -> extrapolation in $WORK/NAME/output
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  sed "1,2s/^\(true\|false\) #Perform A-Matrix/$3 #Perform A-Matrix/" "$2" > "$dir/t.azr"
  (cd "$dir" && printf '3\n\n\n7\n' | timeout 300 "$AZURE2_BIN" $4 --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -s "$dir/output/AZUREOut_aa=1_R=1.extrap" ]; then
    bad "$1: no extrapolation written"; tail -15 "$dir/run.log"; return 1
  fi
}
# maxrel A B -> largest relative difference between the cross-section columns of two .extrap files
maxrel() {
  paste "$1/output/AZUREOut_aa=1_R=1.extrap" "$2/output/AZUREOut_aa=1_R=1.extrap" | awk '
    NF >= 8 { n = NF / 2; a = $4; b = $(n + 4); d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m < 1e-300) m = 1e-300;
              r = d / m; if (r > worst) worst = r; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
same() { awk -v r="$1" 'BEGIN { exit !(r + 0 <= 1e-10) }'; }

# Both models write the same grid: elastic p+12C, 1-3 MeV, differential 60-160 deg and angle-integrated.
sed -n '/<segmentsTest>/,/<\/segmentsTest>/p' "$HERE/models/no_dummies.azr" | grep -q "^1 1 1" || { bad "test segments missing"; exit 1; }
# option off, no dummies: strip the two config keys
sed '/^includeAllJPi=/d; /^maxOrbitalMomentum=/d' "$HERE/models/no_dummies.azr" > "$WORK/no_dummies_off.azr"

for route in true false; do
  for flags in "" "--use-brune"; do
    tag="A$( [ "$route" = true ] && echo mat || echo rmat )${flags:+_brune}"
    run "dummies_$tag" "$HERE/models/with_dummies.azr" "$route" "$flags" || continue
    run "option_$tag"  "$HERE/models/no_dummies.azr"   "$route" "$flags" || continue
    run "neither_$tag" "$WORK/no_dummies_off.azr"      "$route" "$flags" || continue
    r="$(maxrel "$WORK/dummies_$tag" "$WORK/option_$tag")"
    if same "$r"; then ok "$tag: dummy levels == includeAllJPi (max rel diff $r)"; else bad "$tag: dummy levels vs includeAllJPi differ (max rel diff $r)"; fi
    r2="$(maxrel "$WORK/dummies_$tag" "$WORK/neither_$tag")"
    if same "$r2"; then bad "$tag: option off without dummies gives the same result ($r2): the missing partial waves had no effect"; else ok "$tag: option off without dummies differs (max rel diff $r2), so the partial waves matter"; fi
  done
done

# The hard-sphere groups are reported when the compound-nucleus check is on.
mkdir -p "$WORK/report/output" "$WORK/report/checks"
sed '5s/^none #Compond/screen #Compond/' "$HERE/models/no_dummies.azr" > "$WORK/report/t.azr"
(cd "$WORK/report" && printf '3\n\n\n7\n' | timeout 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
n="$(grep -c 'Hard-sphere group added' "$WORK/report/run.log")"
if [ "$n" = 6 ]; then ok "six hard-sphere groups reported (1/2- 3/2+ 5/2- 7/2- 7/2+ 9/2+)"; else bad "expected 6 hard-sphere groups reported, got $n"; grep 'Hard-sphere' "$WORK/report/run.log"; fi

if [ "$failures" -eq 0 ]; then echo "PASS: hard-sphere groups"; else echo "FAIL: $failures check(s)"; exit 1; fi
