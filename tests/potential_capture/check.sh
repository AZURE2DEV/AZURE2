#!/usr/bin/env bash
#
# The interior direct capture of a potential background (docs/upgrade/potential_background.md section 5; README.md
# here).  models/capture_a40.azr.in: 3He(alpha,gamma)7Be with the Gaussian potential VD of Tursunov et al. in both halves,
# the entrance s1/2, d3/2, d5/2 waves with `background interior ... capture interior`, the p3/2 and p1/2 bound states
# as the final states at the potential's own energies and ANCs, levels without widths, E1:
#   1. the capture cross sections into both states at E_cm = 0.02-9 MeV equal the direct potential-model calculation
#      of make_reference.py (scipy and mpmath, 0 to infinity, not the engine; data/reference.csv) at a_c = 3.5, 4.0,
#      4.87 and 5.5 fm;
#   2. and each other's (a_c independence);
#   3. without `capture interior` the same model depends on a_c and misses the reference (the test can fail);
#  3b. through the poles of R_pot (to 9 MeV) the capture does not depend on B_bg or on the energy table;
#   4. a final level away from the potential's eigenvalue (the p3/2 depth changed) is refused, with the depth that
#      would bind it; the set-up note prints the potential's own ANC;
#   5. the analytic gradient of models/fit_level.azr.in (a 1/2+ level on top of the s-wave background, its particle and
#      gamma widths and the ground state's ANC free, capture data) equals finite differences (AZURE_GRAD_CHECK);
#   6. 1 and 3 threads: identical outputs.
#
#   ./tests/potential_capture/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/edit_guard.sh"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/potential_capture.XXXXXX")"
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
  (cd "$WORK/$name" && env OMP_NUM_THREADS="$threads" ${envs[@]+"${envs[@]}"} $(guard_command 900) "$AZURE2_BIN" \
     "$command" ${flags[@]+"${flags[@]}"} m.azr > run.log 2>&1)
}

# the channel radius of every particle row (column 28 of <levels>, pair type 0 in column 27)
radius() {
  edited "$HERE/models/capture_a40.azr.in" "$WORK/models/$2" awk -v a="$1" '
    /^<levels>/ { inl = 1 } /^<\/levels>/ { inl = 0 }
    inl && NF >= 31 && $27 == 0 { $28 = a } { print }'
}
mkdir -p "$WORK/models"
radius 3.5 a35.azr
cp "$HERE/models/capture_a40.azr.in" "$WORK/models/a40.azr"
radius 4.87 a487.azr
radius 5.5 a55.azr
for m in a35 a487 a55; do
  edited "$WORK/models/$m.azr" "$WORK/models/${m}_nocap.azr" sed 's/ capture interior$//'
done
edited "$WORK/models/a40.azr" "$WORK/models/a40_nocap.azr" sed 's/ capture interior$//'

# ---- 1, 2. the identity and the a_c independence
for m in a35 a40 a487 a55; do
  if run "$m" "$WORK/models/$m.azr" extrapolate 2; then :; else bad "$m: the run failed"; tail -3 "$WORK/$m/run.log"; fi
done
if res="$(python3 -I "$HERE/compare.py" "$HERE/data/reference.csv" "$WORK/a35" "$WORK/a40" "$WORK/a487" "$WORK/a55" --tol 1e-5 2>&1)"; then
  ok "sigma(3/2-), sigma(1/2-) at a_c = 3.5, 4.0, 4.87, 5.5 fm = the potential model's (make_reference.py): $res"
else
  bad "against make_reference.py: $res"
fi
if res="$(python3 -I "$HERE/compare.py" --spread "$WORK/a35" "$WORK/a40" "$WORK/a487" "$WORK/a55" --tol 1e-5 2>&1)"; then
  ok "a_c independent: $res"
else
  bad "a_c dependence: $res"
fi

# ---- 3. without the interior capture: a_c dependent, and off the reference
for m in a35 a55; do run "${m}_nocap" "$WORK/models/${m}_nocap.azr" extrapolate 2; done
if res="$(python3 -I "$HERE/compare.py" --spread "$WORK/a35_nocap" "$WORK/a55_nocap" --at-least 0.1 2>&1)"; then
  ok "without capture interior the result depends on a_c (3.5 against 5.5 fm): $res"
else
  bad "without capture interior the result should depend on a_c: $res"
fi
if res="$(python3 -I "$HERE/compare.py" "$HERE/data/reference.csv" "$WORK/a35_nocap" --at-least 0.1 2>&1)"; then
  ok "without capture interior (a_c = 3.5 fm) the reference is missed: $res"
else
  bad "without capture interior the reference should be missed: $res"
fi

# ---- 3b. the energy table through the poles of R_pot (0.02-9 MeV; R_pot has an s-wave pole at E_cm = 7.59 MeV for
#          a_c = 4.87 fm, B_bg = 0, and s, d5/2, d3/2 poles at 5.14, 7.02, 7.53 MeV for 5.5 fm, B_bg = 0): B_bg = 0 and
#          B_bg = -0.05 (s), -0.27 (d3/2), -0.25 (d5/2) give the same capture, with the energy table and without it
#          (--no-ec-table), with and without the interior capture.  Before 2026-11-20 the table divided by the dressed
#          sqrt(P~) and exp(-i phi~), which vanish or jump at the poles: B_bg = 0 against B_bg != 0 differed by up to a
#          factor 50 above 6 MeV.
bnonzero() {
  edited "$1" "$2" awk '
    /wave            l 0 J 1\/2 / { sub(/boundary 0.0/, "boundary -0.05") }
    /wave            l 2 J 3\/2 / { sub(/boundary 0.0/, "boundary -0.27") }
    /wave            l 2 J 5\/2 / { sub(/boundary 0.0/, "boundary -0.25") } { print }'
}
for m in a487 a55 a487_nocap a55_nocap; do bnonzero "$WORK/models/$m.azr" "$WORK/models/${m}_b.azr"; done
run a487_b "$WORK/models/a487_b.azr" extrapolate 2
run a55_b "$WORK/models/a55_b.azr" extrapolate 2
run a487_b_notab "$WORK/models/a487_b.azr" extrapolate 2 -- --no-ec-table
run a487_nocap "$WORK/models/a487_nocap.azr" extrapolate 2
run a487_nocap_b "$WORK/models/a487_nocap_b.azr" extrapolate 2
run a55_nocap_b "$WORK/models/a55_nocap_b.azr" extrapolate 2
run a55_nocap_notab "$WORK/models/a55_nocap.azr" extrapolate 2 -- --no-ec-table
if res="$(python3 -I "$HERE/compare.py" "$HERE/data/reference.csv" "$WORK/a487_b" "$WORK/a55_b" "$WORK/a487_b_notab" --tol 1e-5 2>&1)"; then
  ok "B_bg != 0 (4.87 and 5.5 fm; 4.87 fm without the energy table) = the potential model's: $res"
else
  bad "B_bg != 0 against make_reference.py: $res"
fi
if res="$(python3 -I "$HERE/compare.py" --spread "$WORK/a487_nocap" "$WORK/a487_nocap_b" --tol 1e-5 2>&1)" &&
   res2="$(python3 -I "$HERE/compare.py" --spread "$WORK/a55_nocap" "$WORK/a55_nocap_b" "$WORK/a55_nocap_notab" --tol 1e-5 2>&1)"; then
  ok "without capture interior, B_bg = 0 and != 0 (and table off at 5.5 fm) agree: 4.87 fm $res; 5.5 fm $res2"
else
  bad "without capture interior, B_bg = 0 against != 0 or table on against off: $res ${res2:-}"
fi

# ---- 4. the bound-state rule: with the p3/2 depth 83.65 MeV (written 83.589554) the potential binds 3/2- about 15 keV
#         away from the level: refused, with the depth that binds at the level's energy; the note of the base run
edited "$WORK/models/a40.azr" "$WORK/models/off.azr" sed 's/l 1 J 3\/2 depth 83.589554 MeV/l 1 J 3\/2 depth 83.65 MeV/'
if run off "$WORK/models/off.azr" extrapolate 1; then
  bad "a final state away from the potential's eigenvalue should be refused"
elif depth="$(grep -o 'The depth that binds it at E_f is [0-9.]* MeV' "$WORK/off/run.log")" &&
     awk -v d="$(awk '{ print $(NF - 1) }' <<< "$depth")" 'BEGIN { exit !(d ~ /^[0-9.]+$/ && (d - 83.589554) ^ 2 < 1e-8) }'; then
  ok "depth 83.65 MeV: refused, ${depth%%$'\n'*} (the written 83.589554 to 1e-4 MeV)"
else
  bad "depth 83.65 MeV: refused without the right depth"; tail -3 "$WORK/off/run.log"
fi
note="$(grep 'NOTE: interior capture into the final state of photon pair 2' "$WORK/a40/run.log" || true)"
if [ -n "$note" ]; then
  ok "set-up note: ${note#NOTE: }"
else
  bad "no set-up note of the potential's ANC"
fi

# ---- 5. gradient
run grad "$HERE/models/fit_level.azr.in" fit 2 AZURE_GRAD_CHECK=exit -- --use-gradient
worst="$(awk '/^\[grad-check\] idx / { for (i = 1; i <= NF; i++) if ($i == "rel") print $(i + 1) }' "$WORK/grad/run.log" |
         sort -g | tail -1)"
n="$(grep -c '^\[grad-check\] idx ' "$WORK/grad/run.log")"
if [ "$n" = 3 ] && awk -v w="$worst" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 < 1e-5) }'; then
  ok "analytic gradient = finite differences on 3 parameters (worst relative $worst)"
else
  bad "gradient check: $n parameters, worst relative ${worst:-none}"; tail -3 "$WORK/grad/run.log"
fi

# ---- 6. threads
run t3 "$WORK/models/a40.azr" extrapolate 3
run t1 "$WORK/models/a40.azr" extrapolate 1
if [ -n "$(ls "$WORK/t1/output")" ] && diff -r "$WORK/t1/output" "$WORK/t3/output" > /dev/null; then
  ok "capture_a40: every output file identical at 1 and 3 threads"
else
  bad "capture_a40: outputs differ between 1 and 3 threads"
fi

[ "$failures" -eq 0 ] && echo "potential_capture: PASS" || echo "potential_capture: $failures FAILED"
exit $((failures > 0))
