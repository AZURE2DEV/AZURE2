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
#   6. 1 and 3 threads: identical outputs;
#   7. `only bound_state` on the p waves: their elastic channel functions are those without a potential (identical),
#      the capture equals the one with the potential in the p waves where it is negligible beyond a_c (8 fm), the
#      bound-state rule and note, and the record through --validate, --describe and pyazr.
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

# ---- 7. `only bound_state` (potential_background.md section 5.9): the p waves' potential gives the bound-state shape
#         of the interior capture only, their channel functions are the point-Coulomb ones
#   7a. elastic: models/bound_only.azr.in (p waves `only bound_state`, a 3/2- and a 1/2- level with widths on top of the
#       bound ones, p-wave phase shifts) against the same levels with no potential at all, both without the Coulomb
#       table (one Coulomb routine): the p-wave phase shifts and the p-wave levels of parameters.out are identical
mkdir -p "$WORK/models"
cp "$HERE/models/bound_only.azr.in" "$WORK/models/bso.azr"
edited "$WORK/models/bso.azr" "$WORK/models/bso_nopot.azr" awk '/^<potentials>/ { skip = 1 } !skip { print } /^<\/potentials>/ { skip = 0 }'
run bso "$WORK/models/bso.azr" extrapolate 2 -- --no-coulomb-table
run bso_nopot "$WORK/models/bso_nopot.azr" extrapolate 2 -- --no-coulomb-table
pw() { awk 'NF >= 4 { print $1, $4 }' "$WORK/$1/output/AZUREOut_aa=1_R=1.extrap"; }
plev() { awk '/^J = / { p = ($3 ~ /-$/) } p' "$WORK/$1/output/parameters.out"; }
if [ -s "$WORK/bso/output/AZUREOut_aa=1_R=1.extrap" ] && [ "$(pw bso)" = "$(pw bso_nopot)" ] &&
   [ -n "$(plev bso)" ] && [ "$(plev bso)" = "$(plev bso_nopot)" ]; then
  ok "only bound_state: p-wave phase shifts ($(pw bso | wc -l) values, 0.1-9 MeV) and p-wave levels identical to the model without a potential"
else
  bad "only bound_state: the p-wave channel functions differ from those without a potential"; tail -3 "$WORK/bso/run.log"
fi
# the check bites: with the potential acting in the p waves (no `only`) their phase shifts differ
edited "$WORK/models/bso.azr" "$WORK/models/bso_full.azr" sed 's/ only bound_state$//'
run bso_full "$WORK/models/bso_full.azr" extrapolate 2 -- --no-coulomb-table
if [ -s "$WORK/bso_full/output/AZUREOut_aa=1_R=1.extrap" ] && [ "$(pw bso_full)" != "$(pw bso_nopot)" ]; then
  ok "without only bound_state the p-wave phase shifts differ from those without a potential (the check bites)"
else
  bad "with the potential in the p waves their phase shifts should differ"
fi
#   7b. capture: p waves `only bound_state` against the potential in them (section 1's models): equal where the p-wave
#       potential beyond a_c is negligible (a_c = 8 fm, V(8 fm) = 3 keV), different where it is not (4.87 fm)
for a in 4.87 8.0; do
  edited "$HERE/models/capture_a40.azr.in" "$WORK/models/f$a.azr" awk -v a="$a" '
    /^<levels>/ { inl = 1 } /^<\/levels>/ { inl = 0 } inl && NF >= 31 && $27 == 0 { $28 = a } { print }'
  edited "$WORK/models/f$a.azr" "$WORK/models/b$a.azr" sed 's/^\(    wave            l 1 J .* MeV\)$/\1 only bound_state/'
  run "f$a" "$WORK/models/f$a.azr" extrapolate 2
  run "b$a" "$WORK/models/b$a.azr" extrapolate 2
done
if res="$(python3 -I "$HERE/compare.py" --spread "$WORK/f8.0" "$WORK/b8.0" --tol 2e-4 2>&1)" &&
   res2="$(python3 -I "$HERE/compare.py" --spread "$WORK/f4.87" "$WORK/b4.87" --at-least 1e-2 2>&1)"; then
  ok "capture, only bound_state against the potential in the p waves: a_c = 8 fm $res; 4.87 fm $res2"
else
  bad "capture, only bound_state against the potential in the p waves: 8 fm ${res:-}; 4.87 fm ${res2:-}"
fi
#   7c. the bound-state rule with only bound_state: a depth off the level is refused with the binding depth; the note
#       gives the potential's own ANC (its own exterior: 4.0931604450 fm^-1/2 by make_reference.py) and the jump of
#       the log-derivative at a_c
edited "$WORK/models/bso.azr" "$WORK/models/bso_off.azr" sed 's/l 1 J 3\/2 depth 83.589554 MeV/l 1 J 3\/2 depth 83.65 MeV/'
if run bso_off "$WORK/models/bso_off.azr" extrapolate 1; then
  bad "only bound_state: a depth off the level should be refused"
elif depth="$(grep -o 'The depth that binds it at E_f is [0-9.]* MeV' "$WORK/bso_off/run.log")" &&
     awk -v d="$(awk '{ print $(NF - 1) }' <<< "$depth")" 'BEGIN { exit !(d ~ /^[0-9.]+$/ && (d - 83.589554) ^ 2 < 1e-8) }'; then
  ok "only bound_state, depth 83.65 MeV: refused, ${depth%%$'\n'*}"
else
  bad "only bound_state, depth 83.65 MeV: refused without the right depth"; tail -3 "$WORK/bso_off/run.log"
fi
note="$(grep 'only bound_state, E_f' "$WORK/bso/run.log" || true)"
anc="$(sed -n 's/.*its own ANC \([0-9.]*\) fm.*/\1/p' <<< "$note")"
anc="${anc%%$'\n'*}"
if awk -v c="$anc" 'BEGIN { exit !(c ~ /^[0-9.]+$/ && (c - 4.0931604450) ^ 2 < 1e-8) }' && grep -q "a W'/W = " <<< "$note"; then
  ok "only bound_state note: own ANC $anc fm^-1/2 (make_reference.py 4.0931604450), with the log-derivative jump at a_c"
else
  bad "only bound_state note: '${note:-none}'"
fi
#   7d. the record: validate, describe, and pyazr's unedited round trip keep the qualifier
if "$AZURE2_BIN" --validate "$WORK/models/bso.azr" > "$WORK/v.log" 2>&1 &&
   "$AZURE2_BIN" --describe "$WORK/models/bso.azr" > "$WORK/d.log" 2>&1 &&
   [ "$(grep -c 'only bound_state (point-Coulomb' "$WORK/d.log")" = 2 ]; then
  ok "only bound_state: --validate accepts it, --describe shows both p waves"
else
  bad "only bound_state: --validate or --describe"; tail -3 "$WORK/v.log" "$WORK/d.log"
fi
if python3 -c "import pyazr" 2> /dev/null; then
  cp "$WORK/models/bso.azr" "$WORK/rt.azr"
  if (cd "$WORK" && python3 -c "import pyazr; pyazr.Project.open('rt.azr').save('rt2.azr')") 2> "$WORK/rt.log" &&
     cmp -s "$WORK/models/bso.azr" "$WORK/rt2.azr"; then
    ok "only bound_state: pyazr opens and saves the project byte for byte"
  else
    bad "only bound_state: pyazr round trip"; tail -3 "$WORK/rt.log"
  fi
else
  echo "  skip  pyazr round trip (pyazr not importable)"
fi

[ "$failures" -eq 0 ] && echo "potential_capture: PASS" || echo "potential_capture: $failures FAILED"
exit $((failures > 0))
