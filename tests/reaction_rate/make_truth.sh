#!/usr/bin/env bash
#
# Regenerate expected/rates_*.dat: a dense mode-3 extrapolation of sigma(E) for each truth
# model (models/*_truth.azr, committed; README.md: how the models were made), integrated by truth.py.  Mode 3 does not go through ReactionRate.cpp, so the
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

mkdir -p "$HERE/expected"
for m in p12C:1.00728 n12C:1.008665 moved_azr:1.00728; do
  name="${m%%:*}"; m1="${m##*:}"
  d="$WORK/$name"; mkdir -p "$d/output" "$d/checks"
  cp "$HERE/models/${name}_truth.azr" "$d/t.azr"
  (cd "$d" && "$AZURE2_BIN" extrapolate --no-gui t.azr > run.log 2>&1)
  python3 "$HERE/truth.py" "$m1" 12.0 "$HERE/temperatures.dat" "$d"/output/AZUREOut_aa=1_R=2.extrap \
    > "$HERE/expected/rates_${name}.dat"
  echo "$name: $(wc -l < "$d"/output/AZUREOut_aa=1_R=2.extrap) extrapolation points"
  [ "$name" = moved_azr ] && cp "$d/output/param.par" "$HERE/expected/moved.par"
done
# p12C_closed: the exit pair 3 opens at E_cm = 48 MeV, so the rate's grid has one point of sigma > 0 and the rate
# takes its per-temperature quadrature; the truth is its dense extrapolation (closed_model.py --truth), at temperatures_closed.dat
d="$WORK/p12C_closed"; mkdir -p "$d/output" "$d/checks"
python3 "$HERE/closed_model.py" "$HERE/models/p12C.azr" --truth > "$d/t.azr"
(cd "$d" && "$AZURE2_BIN" extrapolate --no-gui t.azr > run.log 2>&1)
python3 "$HERE/truth.py" 1.00728 12.0 "$HERE/temperatures_closed.dat" "$d"/output/AZUREOut_aa=1_R=3.extrap \
  > "$HERE/expected/rates_p12C_closed.dat"
echo "p12C_closed: $(cat "$d"/output/AZUREOut_aa=1_R=3.extrap | wc -l) extrapolation points"
cat "$HERE"/expected/rates_*.dat
