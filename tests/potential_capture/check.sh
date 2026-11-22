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
#      bound-state rule and note, and the record through --validate, --describe and pyazr;
#   8. a narrow resonance of the background's potential: table on = table off = make_reference.py;
#   9. capture through a closed entrance wave with a background, continuous through its threshold;
#  10. a subtracted box pole put back as a level with the box state's particle and photon amplitudes reproduces the
#      unsubtracted potential (capture and elastic);
#  11. E2 and M1 through the interior capture against an independent calculation.
#  12. the interior charge (`charge_e1 q` on a wave with capture interior, section 14): with q minus the pair's own
#      E1 charge number the interior operator vanishes and the capture is that without `capture interior`; the cross
#      section is a quadratic in q; `charge_e1` without `capture interior` is refused.
#  13. an unbound final state (section 14.3), models/unbound.azr.in: alpha + 12C with the Woods-Saxon of Koonin,
#      Tombrello and Fox, E1 from the d and g waves (background and capture interior) into the potential's own l = 3
#      resonance (E_cm 3.261 MeV, 173 keV wide at 5.35 fm), given as a 3- level of that width: the cross section at
#      a_c = 4.0 and 6.5 fm is that at 5.35 fm to 1 %; without capture interior the radii differ by more than a
#      factor 1.5; the set-up note; `line_window` on it is refused.
#  14. a level on a state of the potential (section 14.5): a 1- level of 400 keV at E_x = 9.60 MeV in a p wave with
#      an interior background lies on the potential's own 1- resonance (9.561 MeV, 412 keV): the set-up warns that it
#      describes the state twice; with `subtract 1` on the wave it is a note; models/unbound.azr.in itself draws
#      neither.
#  15. tools/convert_interior_boundary.py (section 14.6): models/old_layout.azr.in, a project of the branch
#      feature/upgrade-6-interior-potential (`interior_boundary 0` on every l, E1 capture into the potential's own 0+
#      bound state), becomes a `wave l 1 background interior boundary 0 capture interior` (the final state's s wave
#      gets no background), and its capture is the same at a_c = 4.0 and 5.35 fm to 5e-4 (2.3e-4); a second record
#      of the pair with another Woods-Saxon radius becomes a wave with its own radius, one with another Coulomb
#      radius is refused.
#  16. a wave's own Woods-Saxon radius and diffuseness (section 14.8): the p wave with `radius 3.2 fm diffuseness
#      0.55 fm` on a central term of 3.55 and 0.6 fm gives the cross sections of the central term 3.2 and 0.55 fm with
#      the s wave given 3.55 and 0.6 fm, to the digits written; and they differ from the run without them.
#  17. a bound state of the background put back as a level by its observed values (section 14.9):
#      models/bound_level.azr.in (make_bound_level.py; s1/2 depth 104 MeV, the 2-node s state bound by 0.126 MeV,
#      subtracted, and a 1/2+ level with the bound state's energy, ANC and photon widths from scipy alone) gives the
#      capture of the unsubtracted model into both final states to 1e-4 (2e-5); without the level's photon widths it
#      does not.
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

# ---- 8. a narrow resonance of the background's potential (review 2026-11-20, finding 1): the s1/2 depth 100 MeV puts
#         an s-wave resonance about 20 keV wide at E_cm = 0.30 MeV; the energy table, which tabulated a quotient of the
#         dressed phase, was up to 21 % off across it.  A background wave is now integrated directly: table on and off
#         are the same numbers, and both equal make_reference.py --s-depth=100 (data/reference_narrow.csv)
edited "$HERE/models/capture_a40.azr.in" "$WORK/models/narrow.azr" awk '
  /wave            l 0 J 1\/2 depth 67.5 MeV/ { sub(/depth 67.5 MeV/, "depth 100.0 MeV") }
  /^    energies       list / { $0 = "    energies       list 0.27 0.28 0.29 0.295 0.3 0.305 0.31 0.32 MeV cm" } { print }'
run narrow "$WORK/models/narrow.azr" extrapolate 2
run narrow_notab "$WORK/models/narrow.azr" extrapolate 2 -- --no-ec-table
if res="$(python3 -I "$HERE/compare.py" "$HERE/data/reference_narrow.csv" "$WORK/narrow" "$WORK/narrow_notab" --tol 2e-4 2>&1)" &&
   res2="$(python3 -I "$HERE/compare.py" --spread "$WORK/narrow" "$WORK/narrow_notab" --tol 1e-12 2>&1)"; then
  ok "narrow s-wave resonance at 0.30 MeV: table on = table off ($res2), = make_reference.py ($res)"
else
  bad "narrow s-wave resonance: ${res:-} ${res2:-}"
fi

# ---- 9. a closed entrance wave with a background (review 2026-11-20, finding 2): models/closed_channel.azr.in, a
#         second pair 3He + 4He(2.0 MeV) with an s-wave background and capture interior, the final state's ANC in pair 1
#         set to 0 so that all capture runs through pair 2's channel (a 1/2+ level at E_x = 3.6 MeV in both pairs,
#         --no-ignore-externals): the capture is continuous through pair 2's threshold at E_cm = 2.0 MeV, closed below
#         (the decaying wave over the dressed W~_d(a)) and open above (the dressed outgoing wave).  With the undressed
#         W~(a) it jumped by a factor 2.44 there.
cp "$HERE/models/closed_channel.azr.in" "$WORK/models/closed.azr"
run closed "$WORK/models/closed.azr" extrapolate 2 -- --no-ignore-externals
sig() { awk -v E="$1" '$1 + 0 == E + 0 { print $4 }' "$WORK/closed/output/AZUREOut_aa=1_R=3.extrap"; }
lo="$(sig 1.9999)"; hi="$(sig 2.0001)"; lo2="$(sig 1.999)"; hi2="$(sig 2.001)"
if awk -v a="$lo" -v b="$hi" -v c="$lo2" -v d="$hi2" 'BEGIN { exit !(a ~ /^[0-9]/ && b ~ /^[0-9]/ && c ~ /^[0-9]/ && d ~ /^[0-9]/ &&
       a > 1e-14 && (b - a) / a < 1e-3 && (b - a) / a > -1e-3 && ((b - a) - (d - c) / 10) ^ 2 < (0.2 * (d - c) / 10) ^ 2) }'; then
  ok "capture through a closed background wave continuous through its threshold: sigma(1.9999) = $lo, sigma(2.0001) = $hi b (the slope of 1.999 - 2.001 MeV)"
else
  bad "capture through a closed background wave not continuous at its threshold: $lo2 $lo | $hi $hi2"
fi
grep -q "Could not\|ERROR" "$WORK/closed/run.log" && bad "closed_channel run: $(grep -m1 'ERROR' "$WORK/closed/run.log")"

# ---- 10. subtract + capture end to end (potential_background.md section 5.11): the s-wave box state with 2 nodes of
#          section 9's model at 4.87 fm (E_p = 7.59 MeV) subtracted and put back as a 1/2+ level whose particle and
#          photon amplitudes are the box state's (models/subtract_level.azr.in, make_subtract_level.py: Brune's E_b and
#          observed width from the formal (E_p, gamma_p, B_bg), the photon amplitude from the box state's overlap with
#          the bound state by scipy) reproduces the unsubtracted model: capture into both states and the s-wave phase
#          shift.  Without the level's photon widths, or without the level, it does not.
ph='
extrapolation ph_s
    reaction       4He(h,h)4He
    quantity       phase_shift
    unit           deg
    J              1/2
    l              0
    s              1/2
    energies       list 0.1 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 7.5 7.59 8.0 9.0 MeV cm
'
addph() { awk -v ph="$ph" '/^<\/extrapolations>/ { printf "%s", ph } { print }' "$1"; }
edited "$WORK/models/a487.azr" "$WORK/models/sub_U.azr" addph
edited "$HERE/models/subtract_level.azr.in" "$WORK/models/sub_S.azr" addph
edited "$WORK/models/sub_S.azr" "$WORK/models/sub_S0g.azr" awk '$1 == "0.5" && $2 == "1" && ($6 == "2" || $6 == "3") && NF >= 31 { $12 = "0" } { print }'
edited "$WORK/models/sub_S.azr" "$WORK/models/sub_Snolev.azr" awk '$1 == "0.5" && $2 == "1" && NF >= 31 && $3 != "40" { $3 = "40"; $12 = "0" } { print }'
for m in sub_U sub_S sub_S0g sub_Snolev; do run "$m" "$WORK/models/$m.azr" extrapolate 2 -- --no-ignore-externals; done
phs() { awk 'NF >= 4 { print $1, $4 }' "$WORK/$1/output/AZUREOut_aa=1_R=1.extrap"; }
dph="$(paste <(phs sub_U) <(phs sub_S) | awk 'NF == 4 { d = $2 - $4; if (d < 0) d = -d; if (d > m) m = d } END { if (NR != 13) print "x"; else printf "%.3g", m + 0 }')"
if res="$(python3 -I "$HERE/compare.py" --spread "$WORK/sub_U" "$WORK/sub_S" --tol 1e-5 2>&1)" &&
   awk -v d="$dph" 'BEGIN { exit !(d ~ /^[0-9.e+-]+$/ && d + 0 < 1e-4) }'; then
  ok "subtract + its level = the unsubtracted potential: capture $res; s-wave phase to $dph deg"
else
  bad "subtract + its level against the unsubtracted potential: capture ${res:-}; s-wave phase ${dph:-}"
fi
if res="$(python3 -I "$HERE/compare.py" --spread "$WORK/sub_U" "$WORK/sub_S0g" --at-least 0.05 2>&1)" &&
   res2="$(python3 -I "$HERE/compare.py" --spread "$WORK/sub_U" "$WORK/sub_Snolev" --at-least 0.05 2>&1)"; then
  ok "the level carries the subtracted pole's capture: without its photon widths $res; without the level $res2"
else
  bad "without the level's photon widths or the level the capture should differ: ${res:-} ${res2:-}"
fi

# ---- 11. E2 and M1 through the interior capture (potential_background.md section 5.12): models/e2m1.azr.in, the
#          final state VD's bound 0d5/2 state (d5/2 `only bound_state`), the s1/2 and d3/2 waves with capture interior;
#          the photon pair's multipole E2 (s1/2 -> d5/2, d3/2 -> d5/2 through the interior, d5/2 -> d5/2 exterior
#          only), then M1 (d3/2 -> d5/2 interior, d5/2 -> d5/2), orbital (the pair's g factors are 0).  The final
#          level's ANC is 1e-3 fm^-1/2, so that its own pole in the d5/2 channel is negligible: the capture is then
#          (1e-3/C)^2 times the independent calculation of make_reference_e2m1.py (data/reference_e2m1.csv, at the
#          state's own ANC C).  Without capture interior on d3/2 (M1) or s1/2 (E2) it is not.
cp "$HERE/models/e2m1.azr.in" "$WORK/models/e2.azr"
edited "$WORK/models/e2.azr" "$WORK/models/m1.azr" awk '$27 == "10" && NF >= 31 { $8 = "2"; $31 = "2" } { print }'
edited "$WORK/models/e2.azr" "$WORK/models/e2_nos.azr" sed 's/\(l 0 J 1\/2 depth 67.5 MeV background interior boundary 0.0\) capture interior$/\1/'
edited "$WORK/models/m1.azr" "$WORK/models/m1_nod.azr" sed 's/\(l 2 J 3\/2 depth 66.0 MeV background interior boundary 0.0\) capture interior$/\1/'
for m in e2 m1 e2_nos m1_nod; do run "$m" "$WORK/models/$m.azr" extrapolate 2 -- --no-ignore-externals; done
e2m1() {  # e2m1 RUN COLUMN: the largest |sigma/sigma_ref - 1| over the energies, the reference scaled to the ANC 1e-3
  python3 -I - "$HERE/data/reference_e2m1.csv" "$WORK/$1/output/AZUREOut_aa=1_R=2.extrap" "$2" <<'PY'
import csv, re, sys
lines = open(sys.argv[1]).read().split("\n")
C = float(re.search(r"ANC = ([0-9.eE+-]+)", lines[0]).group(1))
ref = {round(float(r["E_cm_MeV"]), 4): float(r[sys.argv[3]]) for r in csv.DictReader(l for l in lines if l and not l.startswith("#"))}
got = {round(float(l.split()[0]), 4): float(l.split()[3]) for l in open(sys.argv[2]) if len(l.split()) >= 5}
d = [abs(got[E] / (ref[E] * (1e-3 / C) ** 2) - 1.0) for E in ref if E in got]
print(f"{max(d):.3e}" if len(d) == len(ref) else "missing")
PY
}
d2="$(e2m1 e2 sigma_E2_b)"; d1="$(e2m1 m1 sigma_M1_b)"; n2="$(e2m1 e2_nos sigma_E2_b)"; n1="$(e2m1 m1_nod sigma_M1_b)"
if awk -v a="$d2" -v b="$d1" 'BEGIN { exit !(a ~ /^[0-9.e+-]+$/ && b ~ /^[0-9.e+-]+$/ && a + 0 < 1e-4 && b + 0 < 1e-4) }'; then
  ok "E2 and M1 with interior capture = make_reference_e2m1.py: largest relative difference E2 $d2, M1 $d1"
else
  bad "E2 ($d2) or M1 ($d1) with interior capture against make_reference_e2m1.py"
fi
if awk -v a="$n2" -v b="$n1" 'BEGIN { exit !(a ~ /^[0-9.e+-]+$/ && b ~ /^[0-9.e+-]+$/ && a + 0 > 1e-2 && b + 0 > 1e-2) }'; then
  ok "without capture interior on s1/2 (E2) or d3/2 (M1) the reference is missed: E2 $n2, M1 $n1"
else
  bad "without capture interior the E2 ($n2) and M1 ($n1) references should be missed"
fi

# ---- 12. the interior charge (potential_background.md section 14).  In the long-wavelength form the E1 operator of
#          the relative motion is e (Z1 m2/M - Z2 m1/M) r; `charge_e1 q` adds e q r inside a_c.  q = -(Z1 m2 - Z2 m1)/M
#          (0.2811 for 3He + 4He with the model's masses) removes the interior operator: the run equals the one
#          without `capture interior`.  The amplitude is linear in q, so sigma(-c) + sigma(c) - 2 sigma(0) =
#          sigma(2c) + sigma(0) - 2 sigma(c).
charge() {  # charge NAME Q: capture_a40 at five energies, charge_e1 Q on every wave with capture interior ("": none)
  edited "$HERE/models/capture_a40.azr.in" "$WORK/models/$1.azr" awk -v q="$2" '
    / capture interior$/ && q != "" { $0 = $0 " charge_e1 " q }
    /^    energies       list / { $0 = "    energies       list 0.1 0.5 1.0 2.0 4.0 MeV cm" } { print }' &&
    run "$1" "$WORK/models/$1.azr" extrapolate 2
}
OWN="$(awk 'BEGIN { printf "%.17g", 2 * (4.0026032 - 3.0160293) / (4.0026032 + 3.0160293) }')"
charge q_zero ""; charge q_minus "-$OWN"; charge q_plus "$OWN"; charge q_two "$(awk -v o="$OWN" 'BEGIN { printf "%.17g", 2 * o }')"
edited "$WORK/models/q_zero.azr" "$WORK/models/q_nocap.azr" sed 's/ capture interior$//' && run q_nocap "$WORK/models/q_nocap.azr" extrapolate 2
if res="$(python3 -I - "$WORK" <<'PY'
import sys
def sigma(name):
    out = []
    for key in (2, 3):
        out += [float(l.split()[3]) for l in open(f"{sys.argv[1]}/{name}/output/AZUREOut_aa=1_R={key}.extrap") if len(l.split()) >= 5]
    return out
z, m, p, t, n = (sigma(x) for x in ("q_zero", "q_minus", "q_plus", "q_two", "q_nocap"))
assert len(z) == 10 and all(len(x) == 10 for x in (m, p, t, n)), "ten cross sections each"
same = max(abs(a / b - 1) for a, b in zip(m, n))
moved = max(abs(a / b - 1) for a, b in zip(z, n))
quad = max(abs((a + b - 2 * c) / (d + c - 2 * b) - 1) for a, b, c, d in zip(m, p, z, t))
print(f"q = -own against no interior capture {same:.1e}, q = 0 against it {moved:.1e}, quadratic to {quad:.1e}")
sys.exit(0 if same < 1e-9 and moved > 1e-2 and quad < 1e-5 else 1)
PY
)"; then ok "interior charge: $res"
else bad "interior charge: ${res:-no result}"; fi
edited "$WORK/models/q_nocap.azr" "$WORK/models/q_refused.azr" sed 's/\(background interior boundary 0.0\)$/\1 charge_e1 0.1/'
mkdir -p "$WORK/q_refused"
if "$AZURE2_BIN" --validate "$WORK/models/q_refused.azr" > "$WORK/q_refused/run.log" 2>&1; then
  bad "charge_e1 without capture interior is accepted"
elif grep -q "E_POTENTIAL_WAVE.*charge_e1 is the charge of the interior capture" "$WORK/q_refused/run.log"; then
  ok "charge_e1 without capture interior is refused (E_POTENTIAL_WAVE)"
else bad "charge_e1 without capture interior: refused, but not in these words"; fi

# ---- 13. an unbound final state (potential_background.md section 14.3).  The level stays where it is while the
#          potential's state moves by 15 keV with the radius, and the state is 5 % of its energy wide: 1 %, not 1e-4.
ub() {  # ub NAME RADIUS [nocap]
  edited "$HERE/models/unbound.azr.in" "$WORK/models/$1.azr" awk -v a="$2" -v nocap="${3:-}" '
    $27 == "0" && NF >= 31 { $28 = a } nocap != "" { sub(/ capture interior$/, "") } { print }' &&
    run "$1" "$WORK/models/$1.azr" extrapolate 2
}
cp "$HERE/models/unbound.azr.in" "$WORK/models/ub535.azr"; run ub535 "$WORK/models/ub535.azr" extrapolate 2
ub ub400 4.000; ub ub650 6.500; ub ubn400 4.000 nocap
edited "$HERE/models/unbound.azr.in" "$WORK/models/ubn535.azr" sed 's/ capture interior$//' && run ubn535 "$WORK/models/ubn535.azr" extrapolate 2
if res="$(python3 -I - "$WORK" <<'PY'
import sys
def sigma(name):
    return [float(l.split()[3]) for l in open(f"{sys.argv[1]}/{name}/output/AZUREOut_aa=1_R=2.extrap") if len(l.split()) >= 5]
a, b, c, n, m = (sigma(x) for x in ("ub535", "ub400", "ub650", "ubn535", "ubn400"))
assert len(a) == 21 and all(len(x) == 21 and min(x) > 0 for x in (b, c, n, m)), "21 cross sections each"
near = max(max(abs(x / y - 1) for x, y in zip(b, a)), max(abs(x / y - 1) for x, y in zip(c, a)))
far = max(max(x / y, y / x) for x, y in zip(m, n))
print(f"4.0 and 6.5 fm against 5.35 fm to {near:.1e}; without capture interior the radii differ by a factor {far:.2f}")
sys.exit(0 if near < 1e-2 and far > 1.5 else 1)
PY
)"; then ok "unbound final state: $res"
else bad "unbound final state: ${res:-no result}"; fi
if grep -q "NOTE: interior capture into the UNBOUND final state of photon pair 2" "$WORK/ub535/run.log" &&
   ! grep -q "interior capture into the UNBOUND" "$WORK/ubn535/run.log"; then
  ok "unbound final state: the set-up note, only with capture interior"
else bad "unbound final state: no set-up note, or one without capture interior"; fi
edited "$HERE/models/unbound.azr.in" "$WORK/models/ubwin.azr" sed 's/^    quantity       sigma$/    quantity       sigma\
    line_window    200 keV/'
if run ubwin "$WORK/models/ubwin.azr" extrapolate 2; then bad "line_window with interior capture ran"
elif grep -q "potential with interior capture, which the line-shape integration is not built" "$WORK/ubwin/run.log"; then
  ok "line_window on a final state whose pair has interior capture is refused"
else bad "line_window with interior capture: stopped, but not in these words"; fi

# ---- 14. a level on a state of the potential (potential_background.md section 14.5)
edited "$HERE/models/unbound.azr.in" "$WORK/models/double.azr" awk '
  /wave            l 4/ { print; print "    wave            l 1 background interior boundary 0"; next }
  /^3 -1 10.4229 / { print; print ""; l = $0; sub(/^3 -1 10.4229 1 1 1 0 6 1 1 1 173042.0/, "1 -1 9.60 1 1 1 0 2 1 1 1 400000.0", l); print l; next }
  { print }' && run double "$WORK/models/double.azr" extrapolate 2
edited "$WORK/models/double.azr" "$WORK/models/double_sub.azr" sed 's/l 1 background interior boundary 0$/l 1 background interior boundary 0 subtract 1/' &&
  run double_sub "$WORK/models/double_sub.azr" extrapolate 2
W14="lies within one width of a resonance of the potential of pair 1 (l = 1)"
if grep -q "WARNING: J^pi=1-.*$W14" "$WORK/double/run.log" && grep -q "a level on top of it describes it twice" "$WORK/double/run.log"; then
  ok "a level on the potential's 1- resonance is warned about"
else bad "a level on the potential's 1- resonance: no warning"; fi
if grep -q "NOTE: J^pi=1-.*$W14" "$WORK/double_sub/run.log" && ! grep -q "WARNING: J^pi=1-.*$W14" "$WORK/double_sub/run.log"; then
  ok "with subtract on the wave it is a note"
else bad "with subtract on the wave: no note, or still a warning"; fi
if [ -s "$WORK/ub535/run.log" ] && ! grep -q "lies within one width of a" "$WORK/ub535/run.log"; then
  ok "no such message for a model without a level on a state of a background wave"
else bad "the unbound model draws the doubled-level message"; fi

# ---- 15. the converter from the layout of feature/upgrade-6-interior-potential (section 14.6)
CONVERT="$HERE/../../tools/convert_interior_boundary.py"
if python3 -I "$CONVERT" "$HERE/models/old_layout.azr.in" "$WORK/models/conv535.azr" 2> "$WORK/convert.log" &&
   grep -q "^    wave            l 1 background interior boundary 0 capture interior$" "$WORK/models/conv535.azr" &&
   [ "$(grep -c "^    wave " "$WORK/models/conv535.azr")" = 1 ] && ! grep -q "interior_boundary" "$WORK/models/conv535.azr" &&
   grep -q "l = 0: holds an external-capture final state" "$WORK/convert.log"; then
  ok "the old layout converts: one wave (l = 1), the final state's s wave without a background, said"
else bad "the converter: $(tail -1 "$WORK/convert.log" 2> /dev/null)"; fi
edited "$WORK/models/conv535.azr" "$WORK/models/conv400.azr" awk '$27 == "0" && NF >= 31 { $28 = "4.000" } { print }'
run conv535 "$WORK/models/conv535.azr" extrapolate 2; run conv400 "$WORK/models/conv400.azr" extrapolate 2
if res="$(python3 -I - "$WORK" <<'PY'
import sys
def sigma(name):
    return [float(l.split()[3]) for l in open(f"{sys.argv[1]}/{name}/output/AZUREOut_aa=1_R=2.extrap") if len(l.split()) >= 5]
a, b = sigma("conv535"), sigma("conv400")
assert len(a) == 51 and len(b) == 51 and min(a) > 0, "51 cross sections each"
d = max(abs(x / y - 1) for x, y in zip(a, b))
print(f"the same to {d:.1e}")
sys.exit(0 if d < 5e-4 else 1)
PY
)"; then
  ok "the converted project's capture at 4.0 and 5.35 fm: $res"
else bad "the converted project's capture depends on the radius: ${res:-}"; fi
second() {  # second NAME COULOMB_RADIUS: old_layout with a second record for l = 1 (radius 3.70 fm) of the pair
  awk -v rc="$2" '/^<\/potentials>$/ { print "potential a12C_p"; print "    pair            a+12C"; print "    coulomb_radius  " rc " fm"
       print "    central         woods_saxon depth 40 MeV radius 3.70 fm diffuseness 0.6 fm"; print "    matching_radius 15 fm"
       print "    step            0.0275 fm"; print "    l               1"; print "    interior_boundary 0" } { print }' \
    "$HERE/models/old_layout.azr.in" > "$WORK/models/$1.azr"
}
second old_two 3.55; second old_rc 3.40
if python3 -I "$CONVERT" "$WORK/models/old_two.azr" "$WORK/models/conv_two.azr" 2> "$WORK/convert_two.log" &&
   grep -q "^    wave            l 1 depth 40 MeV radius 3.70 fm background interior boundary 0 capture interior$" "$WORK/models/conv_two.azr" &&
   "$AZURE2_BIN" --validate "$WORK/models/conv_two.azr" > "$WORK/convert_two.validate" 2>&1; then
  ok "a second record with another Woods-Saxon radius becomes a wave with its own depth and radius, and validates"
else bad "a second record with another radius: $(tail -1 "$WORK/convert_two.log") $(grep -m1 error "$WORK/convert_two.validate" 2> /dev/null)"; fi
if python3 -I "$CONVERT" "$WORK/models/old_rc.azr" "$WORK/models/conv_rc.azr" 2> "$WORK/convert_rc.log"; then
  bad "two records with different Coulomb radii were converted"
elif grep -q "differ in coulomb_radius" "$WORK/convert_rc.log" && [ ! -e "$WORK/models/conv_rc.azr" ]; then
  ok "two records of one pair with different Coulomb radii are refused, nothing written"
else bad "two records with different Coulomb radii: $(tail -1 "$WORK/convert_rc.log")"; fi

# ---- 16. a wave's own Woods-Saxon geometry (section 14.8): one potential written two ways
edited "$WORK/models/conv535.azr" "$WORK/models/geo_wave.azr" sed 's/^\(    wave            l 1\) background/\1 radius 3.2 fm diffuseness 0.55 fm background/' &&
  run geo_wave "$WORK/models/geo_wave.azr" extrapolate 2
edited "$WORK/models/conv535.azr" "$WORK/models/geo_central.azr" awk '
  /^    central / { print "    central         woods_saxon depth 45.891 MeV radius 3.2 fm diffuseness 0.55 fm"
                    print "    wave            l 0 radius 3.55 fm diffuseness 0.6 fm"; next } { print }' &&
  run geo_central "$WORK/models/geo_central.azr" extrapolate 2
if res="$(python3 -I - "$WORK" <<'PY'
import sys
def sigma(name):
    return [float(l.split()[3]) for l in open(f"{sys.argv[1]}/{name}/output/AZUREOut_aa=1_R=2.extrap") if len(l.split()) >= 5]
a, b, c = sigma("geo_wave"), sigma("geo_central"), sigma("conv535")
assert len(a) == 51 and len(b) == 51 and len(c) == 51 and min(a) > 0, "51 cross sections each"
same = max(abs(x / y - 1) for x, y in zip(a, b))
moved = max(abs(x / y - 1) for x, y in zip(a, c))
print(f"written two ways: the same to {same:.1e}; against the central geometry in the p wave: {moved:.1e}")
sys.exit(0 if same < 1e-9 and moved > 1e-2 else 1)
PY
)"; then ok "a wave's own radius and diffuseness: $res"
else bad "a wave's own radius and diffuseness: ${res:-no result}"; fi
mkdir -p "$WORK/geo_refused"
sed 's/depth 67.5 MeV background interior boundary 0.0 capture interior$/depth 67.5 MeV radius 3.0 fm background interior boundary 0.0 capture interior/' \
  "$HERE/models/capture_a40.azr.in" > "$WORK/models/geo_refused.azr"
if "$AZURE2_BIN" --validate "$WORK/models/geo_refused.azr" > "$WORK/geo_refused/run.log" 2>&1; then bad "a wave radius on a Gaussian term is accepted"
elif grep -q "E_POTENTIAL_WAVE.*radius belongs to a woods_saxon central term" "$WORK/geo_refused/run.log"; then
  ok "a wave radius on a Gaussian central term is refused (E_POTENTIAL_WAVE)"
else bad "a wave radius on a Gaussian term: refused, but not in these words"; fi

# ---- 17. a bound state of the background as a level given by its observed values (section 14.9).  The level's ANC
#          and photon widths are the normalised bound state's, from make_bound_level.py (scipy); the engine has to
#          turn the ANC into the reduced width of the dressed channel and split the observed photon width into the
#          internal amplitude and the channel's (interior and exterior) part.
edited "$WORK/models/a487.azr" "$WORK/models/bound_U.azr" sed 's/l 0 J 1\/2 depth 67.5 MeV/l 0 J 1\/2 depth 104.0 MeV/'
edited "$HERE/models/bound_level.azr.in" "$WORK/models/bound_S0g.azr" awk '$1 == "0.5" && $2 == "1" && ($6 == "2" || $6 == "3") && NF >= 31 { $12 = "0" } { print }'
run bound_U "$WORK/models/bound_U.azr" extrapolate 2 -- --no-ignore-externals
run bound_S "$HERE/models/bound_level.azr.in" extrapolate 2 -- --no-ignore-externals
run bound_S0g "$WORK/models/bound_S0g.azr" extrapolate 2 -- --no-ignore-externals
if res="$(python3 -I - "$WORK" <<'PY'
import sys
def sigma(name, pair):
    return [float(l.split()[3]) for l in open(f"{sys.argv[1]}/{name}/output/AZUREOut_aa=1_R={pair}.extrap") if len(l.split()) >= 5]
same = moved = 0.0
for pair in (2, 3):
    u, s, n = sigma("bound_U", pair), sigma("bound_S", pair), sigma("bound_S0g", pair)
    assert len(u) == 24 and len(s) == 24 and len(n) == 24 and min(u) > 0, "24 cross sections each"
    same = max(same, max(abs(x / y - 1) for x, y in zip(s, u)))
    moved = max(moved, max(abs(x / y - 1) for x, y in zip(n, u)))
print(f"the level = the potential's bound state to {same:.1e}; without its photon widths {moved:.1e}")
sys.exit(0 if same < 1e-4 and moved > 0.2 else 1)
PY
)"; then ok "a subtracted bound state put back by its observed energy, ANC and photon widths: $res"
else bad "a subtracted bound state as a level: ${res:-no result}"; fi
grep -q "Could not\|ERROR" "$WORK/bound_S/run.log" && bad "bound_level run: $(grep -m1 'ERROR' "$WORK/bound_S/run.log")"

[ "$failures" -eq 0 ] && echo "potential_capture: PASS" || echo "potential_capture: $failures FAILED"
exit $((failures > 0))
