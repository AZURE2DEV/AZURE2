#!/usr/bin/env bash
#
# External capture across a charged-particle threshold (models from make_models.py): 16O(n,g)17O(g.s.) with channel
# capture through alpha + 13C, which opens at E_cm(n) = 2.21556 MeV.  Within ~1 keV above that threshold the alpha +
# 13C Coulomb functions over- and underflow (eta > 100); before 2026-10-24 the channel-capture integral there came out
# as 0, as an overflow, as NaN or 30 % off, and the cross section of this model reached 1e177 b (direct integrals) and
# 1e305 b (energy tables).  Now the integral is bridged linearly from the threshold value (closed-channel form) to
# the open form at eta = 100 (ECIntegral::CalculateDirect).  With the energy tables (default) and with --no-ec-table:
#
#   (a) every cross section is finite and positive;
#   (b) the excitation function is smooth: each point lies within 1e-5 (relative) of the straight line through its
#       neighbours (grid 0.01-0.25 keV across the threshold; the physical curvature is far below that);
#   (c) the two routes agree to 1e-6.
#
#   ./tests/ec_threshold/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_threshold.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
OUT=output/AZUREOut_aa=1_R=3.out

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {   # run DIR [FLAGS]
  local dir="$WORK/$1" flags="${2:-}"
  mkdir -p "$dir/output" "$dir/checks"
  cp -r "$HERE/models/"* "$dir/"
  # shellcheck disable=SC2086
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-long-wavelength --use-brune $flags \
    --no-gui --no-readline ng.azr > run.log 2>&1)
}

echo "ec_threshold: $AZURE2_BIN"
run tables
run direct --no-ec-table
for r in tables direct; do
  f="$WORK/$r/$OUT"
  if [ ! -s "$f" ]; then bad "($r) no output"; continue; fi
  res=$(awk 'NF >= 4 { n++; e[n] = $1; s[n] = $4 }
    END {
      bad = 0; worst = 0
      for (i = 1; i <= n; i++) if (!(s[i] > 0 && s[i] < 1e300)) bad++
      for (i = 2; i < n; i++) {
        lin = s[i-1] + (s[i+1] - s[i-1]) * (e[i] - e[i-1]) / (e[i+1] - e[i-1])
        d = (s[i] - lin) / s[i]; if (d < 0) d = -d; if (d > worst) worst = d
      }
      printf "%d %d %.2e", n, bad, worst
    }' "$f")
  set -- $res
  if [ "$2" -eq 0 ]; then ok "(a, $r) all $1 cross sections finite and positive"; else bad "(a, $r) $2 of $1 cross sections not finite/positive (or > 1e300)"; fi
  if awk -v w="$3" 'BEGIN { exit !(w < 1e-5) }'; then ok "(b, $r) smooth: largest deviation from the neighbours' line $3"
  else bad "(b, $r) not smooth: largest deviation from the neighbours' line $3"; fi
done
if [ -s "$WORK/tables/$OUT" ] && [ -s "$WORK/direct/$OUT" ]; then
  dev=$(paste "$WORK/tables/$OUT" "$WORK/direct/$OUT" |
    awk 'NF >= 8 && $4 != 0 { n++; d = ($4 - $(NF/2+4)) / $4; if (d < 0) d = -d; if (d > w) w = d }
         END { if (n) printf "%.2e", w; else print "none" }')
  if [ "$dev" != none ] && awk -v d="$dev" 'BEGIN { exit !(d < 1e-6) }'; then ok "(c) energy tables = direct integrals ($dev)"
  else bad "(c) energy tables differ from the direct integrals ($dev)"; fi
fi

if [ "$failures" -eq 0 ]; then echo "ec_threshold: PASS"; else echo "ec_threshold: FAIL ($failures)"; exit 1; fi
