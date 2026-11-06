#!/usr/bin/env bash
#
# Regenerate expected/rates_*.dat: a dense mode-3 extrapolation of sigma(E) for each truth
# model, integrated by truth.py.  Mode 3 does not go through ReactionRate.cpp, so the
# reference is independent of the code under test.  Also writes expected/moved.par, the
# parameter file of moved_azr (the narrow level at E_x + 5 keV), which check.sh feeds to
# the moved model.
#
#   ./tests/reaction_rate/make_truth.sh path/to/AZURE2

set -euo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: make_truth.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/rate_truth.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

python3 "$HERE/make_models.py"
mkdir -p "$HERE/expected"
for m in p12C:1.00728 n12C:1.008665 moved_azr:1.00728; do
  name="${m%%:*}"; m1="${m##*:}"
  d="$WORK/$name"; mkdir -p "$d/output" "$d/checks"
  cp "$HERE/models/${name}_truth.azr" "$d/t.azr"
  (cd "$d" && printf '3\n\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  python3 "$HERE/truth.py" "$m1" 12.0 "$HERE/temperatures.dat" "$d"/output/AZUREOut_aa=1_R=2.extrap \
    > "$HERE/expected/rates_${name}.dat"
  echo "$name: $(wc -l < "$d"/output/AZUREOut_aa=1_R=2.extrap) extrapolation points"
  [ "$name" = moved_azr ] && cp "$d/output/param.par" "$HERE/expected/moved.par"
done
cat "$HERE"/expected/rates_*.dat
