#!/usr/bin/env bash
#
# A capture final state above a particle threshold is refused (E_FINAL_STATE_UNBOUND, experiment/src/FinalStates.cpp
# UnboundFinalStates): the engine computes external capture into it as bound (|E_f - S| as a binding energy) and a
# point below it at a negative photon energy (NaN).  Model: tests/upgrade/T1, whose 3/2- final state of 12C(p,g)13N[1.0]
# (photon pair 3) is moved, with its level, to E_f; S_p = 1.94349 MeV.
#   (a) E_f = 2.5 MeV, external capture on (E1, M1, E2): --validate names the code on the capture data set cap_g1 with
#       its exit code 2; the run prints the same line and is refused (no chiSquared.out, exit code not 0)
#   (b) E_f = 2.5 MeV, external capture off (column 31 = 0): the same from the entrance threshold, on cap_g1 and on the
#       extrapolation x2 (both capture from p+12C into pair 3)
#   (c) E_f = 1.94 MeV, 3.5 keV below the threshold: valid, runs, and gives the chi^2 the start commit gave
#       (ENGINE_CHI2 below, "Total Chi-Squared" of the build of 1bd1a0d5, 12 digits; 1e-9 relative)
#
#   ./tests/ec_final_state_unbound/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_final_state_unbound.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"
ENGINE_CHI2=219429.963689  # "Total Chi-Squared" of the run of (c) with the start commit
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

variant() {  # NAME E_F MASK: T1 with the 3/2- final state (level and photon pair 3) at E_F, pair 3's column 31 = MASK
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp -R "$HERE/../upgrade/T1/data" "$d/"
  awk -v e="$2" -v mask="$3" '
    /^<levels>/ { inside = 1 } /^<\/levels>/ { inside = 0 }
    inside && NF >= 31 {
      if ($1 == 1.5 && $2 == -1 && $3 == 1.0) $3 = e
      if ($6 == 3) { $17 = e; $31 = mask }
    }
    { gsub(/13N\[1\.0\]/, "13N[" e "]"); print }' "$HERE/../upgrade/T1/T1.azr" > "$d/t.azr"
}
validate() { (cd "$WORK/$1" && $RUN "$AZURE2_BIN" --validate t.azr > validate.log 2>&1; echo $? > validate.exit); }
run() { (cd "$WORK/$1" && $RUN "$AZURE2_BIN" calculate --no-gui t.azr > run.log 2>&1; echo $? > run.exit); }

refused() {  # LABEL NAME COUNT: COUNT x E_FINAL_STATE_UNBOUND, one on cap_g1, by --validate (exit 2); the run refused
  local n e r
  validate "$2"; run "$2"
  n=$(grep -c "E_FINAL_STATE_UNBOUND" "$WORK/$2/validate.log")
  e=$(cat "$WORK/$2/validate.exit"); r=$(cat "$WORK/$2/run.exit")
  local said="E_FINAL_STATE_UNBOUND.*Capture into an unbound final state is not computed by this version"
  if [ "$n" = "$3" ] && [ "$e" = 2 ] && grep -q "cap_g1.reaction: error E_FINAL_STATE_UNBOUND" "$WORK/$2/validate.log" &&
     grep -q "^ERROR: .*$said" "$WORK/$2/run.log" &&
     [ "$r" != 0 ] && [ ! -s "$WORK/$2/output/chiSquared.out" ]; then
    ok "$1"
  else
    bad "$1: $n x E_FINAL_STATE_UNBOUND, --validate exit $e, run exit $r"
    sed -n '1,3p' "$WORK/$2/validate.log"; grep -E "ERROR" "$WORK/$2/run.log" | sed -n '1,3p'
  fi
}

echo "ec_final_state_unbound: $AZURE2_BIN"
variant a 2.5 7
refused "(a) final state at 2.5 MeV above S_p, external capture: refused by --validate and the run" a 1
grep -q "level ID 2), lies at or above the threshold of pair 1 (p+12C), E_x = 1.94349 MeV" "$WORK/a/validate.log" ||
  bad "(a) the message names the level and both energies: $(grep E_FINAL_STATE_UNBOUND "$WORK/a/validate.log")"

variant b 2.5 0
refused "(b) final state at 2.5 MeV, internal capture only: refused from the entrance threshold" b 2
grep -q "photon energy is negative" "$WORK/b/validate.log" ||
  bad "(b) the message: $(grep E_FINAL_STATE_UNBOUND "$WORK/b/validate.log")"

variant c 1.94 7
validate c; run c
chi2=$(awk '/^Total Chi-Squared:/ { print $3 }' "$WORK/c/run.log" 2>/dev/null)
if [ "$(cat "$WORK/c/validate.exit")" = 0 ] && [ "$(cat "$WORK/c/run.exit")" = 0 ] &&
   awk -v a="$chi2" -v b="$ENGINE_CHI2" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(a == a + 0 && d <= 1e-9 * b) }'
then
  ok "(c) final state 3.5 keV below S_p: valid, runs, chi^2 $chi2 as the start commit's"
else
  bad "(c) final state 3.5 keV below S_p: --validate exit $(cat "$WORK/c/validate.exit"),\
 run exit $(cat "$WORK/c/run.exit"), chi^2 '$chi2' (want $ENGINE_CHI2)"
fi

if [ "$failures" -eq 0 ]; then echo "ec_final_state_unbound: PASS"
else echo "ec_final_state_unbound: FAIL ($failures)"; exit 1; fi
