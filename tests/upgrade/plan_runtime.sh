#!/bin/bash
# The calculation plan at run time (include/CalculationPlan.h, docs/upgrade/data_model.md S5): with AZURE_PLAN_CHECK=1
# the point pass runs its energy tasks and folds one after the other and checks that each writes only its own slots
# and that no slot is written by two of them (a breach stops the run); then every unit is solved alone and compared
# bit for bit with the solve it shares at its compound energy (azure2::CheckEnergyTasks; a difference stops the run).
# Mode-1 calculations of T1, T2 (with a SUM and a RATIO), T3, 13N, natural_neon (sample sums), 18O_p_a_thm (THM),
# thick_target_resonance (seven capture exits under a thick target), target_effect_ranges (automatic tolerance,
# windows, probes), electron_screening's grad.azr (see below) and tests/transmission_yield's conv (a transmission and
# a yield under a Gaussian spread), with and
# without the check: each run exits 0, the checked run reports "0 differ", and every output file is byte for byte the
# same.  Shown failing with a task that also calculated the next task's host (injected once: "writes slot ...,
# outside its records"), with a T-matrix kept across exit pairs (injected once: "T-matrix differs"), and on the
# screening model by T-matrix groups without the amplitude source (about 1,600 units differ).
#   plan_runtime.sh path/to/AZURE2
set -u
export LC_ALL=C OMP_NUM_THREADS=2
[ $# -ge 1 ] || { echo "usage: plan_runtime.sh path/to/AZURE2"; exit 2; }
BIN="$1"
TESTS="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail=0
ok()  { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; fail=1; }

mkdir -p "$WORK/ty"
python3 -I "$TESTS/transmission_yield/make_models.py" "$TESTS/reaction_rate/models/n12C.azr" \
  "$TESTS/natural_neon/natural_neon.azr" "$WORK/ty" > "$WORK/ty.log" 2>&1 || { bad "transmission_yield models"; tail -3 "$WORK/ty.log"; }
cp -r "$TESTS/natural_neon/data" "$WORK/ty/" 2>/dev/null
# 13N's capture data under Gaussian spreads with electron screening (tests/electron_screening/grad.azr): points and
# sub-points of two data sets at one energy hold set-up external-capture amplitudes that differ in the last bits.
mkdir -p "$WORK/scr" && cp "$TESTS/electron_screening/grad.azr" "$WORK/scr/" && cp -r "$TESTS/13N/data" "$WORK/scr/"
# run NAME SRC_DIR AZR CHECK -> $WORK/NAME.CHECK (exit code in .rc)
run() {
  local d="$WORK/$1.$4"
  mkdir -p "$d" && cp -r "$2/." "$d/" && rm -rf "$d/output" "$d/checks" && mkdir -p "$d/output" "$d/checks"
  (cd "$d" && printf '1\n\n\n7\n' | AZURE_PLAN_CHECK=$4 "$BIN" --no-gui --no-readline "$3" > run.log 2>&1; echo $? > rc)
}
for p in "T1:$TESTS/upgrade/T1:T1.azr" "T2:$TESTS/upgrade/T2:T2.azr" "T3:$TESTS/upgrade/T3:T3.azr" "13N:$TESTS/13N:13N.azr" \
         "natural_neon:$TESTS/natural_neon:natural_neon.azr" "thm:$TESTS/18O_p_a_thm:18O_p_a_thm.azr" "conv:$WORK/ty:conv.azr" \
         "thick:$TESTS/thick_target_resonance:thick_target_resonance.azr" \
         "ranges:$TESTS/target_effect_ranges:target_effect_ranges.azr" "screening:$WORK/scr:grad.azr"; do
  IFS=: read -r name src azr <<< "$p"
  run "$name" "$src" "$azr" 0
  run "$name" "$src" "$azr" 1
  r0="$(cat "$WORK/$name.0/rc")"; r1="$(cat "$WORK/$name.1/rc")"
  n="$(ls "$WORK/$name.0/output" | wc -l | tr -d ' ')"
  solve="$(grep -a "AZURE_PLAN_CHECK: [0-9]" "$WORK/$name.1/run.log" | tail -1)"
  if [ "$r0" = 0 ] && [ "$r1" = 0 ] && [ "$n" -gt 0 ] && [[ "$solve" == *" 0 differ" ]] &&
     diff -rq "$WORK/$name.0/output" "$WORK/$name.1/output" > /dev/null; then
    ok "$name: the checked run passes (${solve#AZURE_PLAN_CHECK: }), its $n output files identical"
  else
    bad "$name: exit $r0 / $r1 (unchecked / checked), $n files; $(grep -a AZURE_PLAN_CHECK "$WORK/$name.1/run.log" | head -1)"
  fi
done
[ $fail -eq 0 ] && echo "plan_runtime: PASS" || echo "plan_runtime: FAIL"
exit $fail
