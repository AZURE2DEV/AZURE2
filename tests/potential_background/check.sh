#!/usr/bin/env bash
#
# The interior background of a hybrid pair (docs/upgrade/potential_background.md; README.md here):
#   1. the identity: models/identity_a35.azr (a_c = 3.5 fm, B_bg = 0) and identity_a50.azr (a_c = 5.0 fm, B_bg = -1),
#      3He + 4He with the Gaussian potential VD in both halves and levels without widths, give the potential's phase
#      shifts in seven partial waves (data/phases.csv, make_phases.py: scipy and mpmath, not the engine) to 1e-4 deg,
#      and each other's to 1e-4 deg;
#   2. the background acts: exterior_only.azr (the same potential, no background) misses the s-wave phase by >= 1 deg at some energy;
#   3. a level below threshold in a partial wave with a background is refused with the reason (below_threshold.azr);
#   4. the analytic gradient of fit_level.azr (a free 1/2+ level on top of the s-wave background) equals finite
#      differences (AZURE_GRAD_CHECK);
#   5. identity_a35.azr at 1 and 3 threads: every output file identical;
#   6. (1b) scatteringLengths.out does not compute a0 for an s wave with a background (its threshold formula reads
#      the point-Coulomb functions).
#
#   ./tests/potential_background/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/potential_background.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MODEL COMMAND THREADS [env ...] [-- flags]: in $WORK/NAME; returns the program's exit status
run() {
  local name="$1" model="$2" command="$3" threads="$4"; shift 4
  local envs=() flags=()
  while [ $# -gt 0 ] && [ "$1" != "--" ]; do envs+=("$1"); shift; done
  [ $# -gt 0 ] && shift && flags=("$@")
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"
  [ -e "$WORK/$name/data" ] || ln -s "$HERE/models/data" "$WORK/$name/data"
  cp "$model" "$WORK/$name/m.azr"
  (cd "$WORK/$name" && env OMP_NUM_THREADS="$threads" ${envs[@]+"${envs[@]}"} $(guard_command 600) "$AZURE2_BIN" \
     "$command" ${flags[@]+"${flags[@]}"} m.azr > run.log 2>&1)
}
extrap() { echo "$WORK/$1/output/AZUREOut_aa=1_R=1.extrap"; }

# ---- 1. the identity
run a35 "$HERE/models/identity_a35.azr" extrapolate 1
run a50 "$HERE/models/identity_a50.azr" extrapolate 1
for n in a35 a50; do
  if res="$(python3 -I "$HERE/phase_check.py" "$(extrap $n)" "$HERE/data/phases.csv" 1e-4 2>&1)"; then
    ok "$n: phase shifts = the potential's (make_phases.py): $res"
  else
    bad "$n: phase shifts against make_phases.py: $res"; tail -3 "$WORK/$n/run.log"
  fi
done
if res="$(python3 -I "$HERE/phase_check.py" "$(extrap a35)" "$(extrap a50)" 1e-4 --pair 2>&1)"; then
  ok "a_c = 3.5 fm (B_bg 0) and 5.0 fm (B_bg -1) agree: $res"
else
  bad "a_c = 3.5 and 5.0 fm: $res"
fi

# ---- 1b. the threshold formula of the scattering length does not know the background: refused in the report
if grep -q "not computed: interior background" "$WORK/a35/output/scatteringLengths.out"; then
  ok "scatteringLengths.out: the s wave with a background is not computed (the formula reads point-Coulomb functions)"
else
  bad "scatteringLengths.out should refuse the s wave with a background"
fi

# ---- 2. without the background the exterior alone misses the potential's phase
run ext "$HERE/models/exterior_only.azr" extrapolate 1
if res="$(python3 -I "$HERE/phase_check.py" "$(extrap ext)" "$HERE/data/phases.csv" 1.0 --differs 2>&1)"; then
  ok "the background acts: without it $res"
else
  bad "without the background the s wave should miss the potential's phase: $res"
fi

# ---- 3. refused below threshold
if run below "$HERE/models/below_threshold.azr" extrapolate 1; then
  bad "a bound 3/2- level with a background in l = 1, J = 3/2 ran"
elif grep -q "has an interior background and is evaluated at or below threshold" "$WORK/below/run.log"; then
  ok "a level below threshold in a background partial wave is refused with the reason"
else
  bad "below threshold: refused without the reason"; tail -3 "$WORK/below/run.log"
fi

# ---- 4. gradient
run grad "$HERE/models/fit_level.azr" fit 2 AZURE_GRAD_CHECK=exit -- --use-gradient
worst="$(awk '/^\[grad-check\] idx / { for (i = 1; i <= NF; i++) if ($i == "rel") print $(i + 1) }' "$WORK/grad/run.log" |
         sort -g | tail -1)"
n="$(grep -c '^\[grad-check\] idx ' "$WORK/grad/run.log")"
if [ "$n" = 2 ] && awk -v w="$worst" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 < 1e-6) }'; then
  ok "analytic gradient = finite differences on 2 parameters (worst relative $worst)"
else
  bad "gradient check: $n parameters, worst relative ${worst:-none}"; tail -3 "$WORK/grad/run.log"
fi

# ---- 5. threads
run t1 "$HERE/models/identity_a35.azr" extrapolate 1
run t3 "$HERE/models/identity_a35.azr" extrapolate 3
if [ -n "$(ls "$WORK/t1/output")" ] && diff -r "$WORK/t1/output" "$WORK/t3/output" > /dev/null; then
  ok "identity_a35: every output file identical at 1 and 3 threads"
else
  bad "identity_a35: outputs differ between 1 and 3 threads"
fi

[ "$failures" -eq 0 ] && echo "potential_background: PASS" || echo "potential_background: $failures FAILED"
exit $((failures > 0))
