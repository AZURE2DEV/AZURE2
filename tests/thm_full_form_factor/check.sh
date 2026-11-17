#!/usr/bin/env bash
#
# Full form factor of a THM experiment: `vertexModel=ff` on an `experiment[<name>]` line of the <thm> section
# (ThmFormFactor.h; docs/source/theory/thm_implementation.rst, "Full form factor"; docs/upgrade/thm_form_factor.md).
# The vertex of every level becomes M_l + Delta_lambda, Delta the volume term of a Woods-Saxon level function.
#
# Model: tests/18O_p_a_thm (p + 18O, two 1/2+ levels at 609 and 812 keV, entrance l = 0, energies fixed), with
# vertex=perlevel.
#
#   (a) ffVolume=off is byte-identical to vertexModel=pw (the surface term alone), with and without a ps window;
#   (b) ff changes the model; the set-up report lists both levels with the Wildermuth node count of a proton on 18O
#       (l = 0: 1s, n_r = 1) and R = 1 + Delta/M_l at E_lambda as the thm-review implementation (4882ad04) gave them on
#       this project (2.77, 4.63; that code is the source of this one, so the check guards the port, and the physics
#       itself is checked against independent values by the ctests physics_transfer_form_factor and
#       pyazr_physics_table); ffNodes=1/2+:0 changes them;
#   (c) refusals: vertex=constant, a free level energy, an ff key without vertexModel=ff, an ffNodes term that matches
#       no level, an unknown vertexModel;
#   (d) the output is byte-identical on 1 and 4 threads;
#   (e) the analytic gradient (d chi2 / d gamma of the four widths, entrance and exit) agrees with finite differences
#       (AZURE_GRAD_CHECK), with and without a ps window: the entrance widths' derivative carries the volume term.
#       Shown to fail on a build whose gradient took the surface term alone (relative error 0.25).
#
#   ./tests/thm_full_form_factor/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
O18="$HERE/../18O_p_a_thm"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_full_form_factor.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
. "$HERE/../lib/thm_records.sh"
. "$HERE/../lib/edit_guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"
OUT="AZUREOut_aa=1_R=2.out"
BOTH="lc723_thm_points lc723_band_mid"
KIN="beam=18O target=d spectator=n Ebeam=54"
PER="vertex=perlevel"
fail=0

# Every width of <levels> made free (column 11 of the level lines).
FREE_W='/^<levels>/ { L = 1 } /^<\/levels>/ { L = 0 } L && NF > 20 { $11 = 0 } { print }'
# The first level's energy made free (column 4 of its two <levels> lines).
FREE_E='/^<levels>/ { L = 1 } /^<\/levels>/ { L = 0 } L && $3 == "8.602600" { $4 = 0 } { print }'

# run NAME BLOCK [EDIT] -- a "calculate with data" run of a copy of tests/18O_p_a_thm with BLOCK as its <thm> section,
# the .azr first passed through the awk program named EDIT (checked to change the file).  Status in NAME/status.
run() {
  local d="$WORK/$1"
  rm -rf "$d"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$O18/data" "$d/"
  cp "$O18/18O_p_a_thm.azr" "$d/base.azr"
  if [ -n "${3:-}" ]; then
    edited "$d/base.azr" "$d/run.azr" awk "${!3}" || { echo 1 > "$d/status"; return; }
  else
    cp "$d/base.azr" "$d/run.azr"
  fi
  thm_set "$d/run.azr" "$2" $BOTH || { echo 1 > "$d/status"; return; }
  (cd "$d" && $RUN "$AZURE2_BIN" calculate --no-gui run.azr 2>&1 | head -c 1000000 > log;
   echo "${PIPESTATUS[0]}" > status)
}

# (a) ----------------------------------------------------------------------
echo "(a) ffVolume=off == vertexModel=pw"
run pw "$PER
experiment[A] $KIN"
run off "$PER
experiment[A] $KIN vertexModel=ff ffVolume=off"
run pw_ps "$PER
experiment[A] $KIN ps=hulthen:0-30 psNodes=6"
run off_ps "$PER
experiment[A] $KIN ps=hulthen:0-30 psNodes=6 vertexModel=ff ffVolume=off"
if ran pw && ran off && ran pw_ps && ran off_ps; then
  for f in chiSquared.out normalizations.out "$OUT"; do
    same "$WORK/pw/output/$f" "$WORK/off/output/$f" && ok "$f byte-identical" || bad "$f differs"
    same "$WORK/pw_ps/output/$f" "$WORK/off_ps/output/$f" && ok "$f byte-identical with a ps window" \
      || bad "$f differs with a ps window"
  done
fi

# (b) ----------------------------------------------------------------------
echo "(b) ff: the volume term"
run ff "$PER
experiment[A] $KIN vertexModel=ff"
run ff_ps "$PER
experiment[A] $KIN ps=hulthen:0-30 psNodes=6 vertexModel=ff"
run ff0 "$PER
experiment[A] $KIN vertexModel=ff ffNodes=1/2+:0"
if ran ff && ran ff_ps && ran ff0; then
  same "$WORK/pw/output/$OUT" "$WORK/ff/output/$OUT" && bad "ff gives the pw model" || ok "ff changes the model"
  same "$WORK/pw_ps/output/$OUT" "$WORK/ff_ps/output/$OUT" && bad "ff gives the pw model with a ps window" \
    || ok "ff changes the model with a ps window"
  rows="$(tr -d '\r' < "$WORK/ff/log" | awk '$1 == "1/2+" && $3 == "l" && $6 == "n_r" { print $2, $8, $NF }')"
  n="$(printf '%s\n' "$rows" | grep -c .)"
  [ "$n" = 2 ] && ok "set-up report: two level rows" || bad "set-up report rows: $n"
  printf '%s\n' "$rows" | awk '{ print }' | while read -r lev nodes R; do
    want=$([ "$lev" = 1 ] && echo 2.77 || echo 4.63)
    awk -v r="$R" -v w="$want" -v n="$nodes" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && w ~ /^ *[-+]?\.?[0-9]/ && n == 1 && r - w < 0.01 && w - r < 0.01) }' \
      && ok "level $lev: n_r = $nodes (Wildermuth), R = $R (thm-review: $want)" \
      || { echo "  FAIL  level $lev: n_r = $nodes, R = $R (want 1, $want)"; touch "$WORK/badrow"; }
  done
  [ -f "$WORK/badrow" ] && fail=1
  grep -q 'Wildermuth, proton' "$WORK/ff/log" && ok "node counts from the Wildermuth rule" || bad "no Wildermuth rows"
  grep -q '(ffNodes)' "$WORK/ff0/log" && ok "ffNodes= takes over" || bad "ffNodes= not used"
  same "$WORK/ff/output/$OUT" "$WORK/ff0/output/$OUT" && bad "the node count changes nothing" \
    || ok "the node count changes the model"
fi

# (c) ----------------------------------------------------------------------
echo "(c) refused"
refuse() {  # refuse NAME MESSAGE-FRAGMENT BLOCK [EDIT]
  run "$1" "$3" "${4:-}"
  if [ "$(cat "$WORK/$1/status")" != 0 ] && grep -q "ERROR" "$WORK/$1/log" && grep -q -- "$2" "$WORK/$1/log"; then
    ok "$1: $(grep -m1 'ERROR' "$WORK/$1/log" | tr -d '\r' | cut -c1-150)"
  else
    bad "$1: not refused with '$2' (status $(cat "$WORK/$1/status"))"; tail -3 "$WORK/$1/log" | sed 's/^/        /'
  fi
}
refuse constant "vertexModel=ff needs vertex=perlevel" "experiment[A] $KIN vertexModel=ff"
refuse free_energy "has a free energy" "$PER
experiment[A] $KIN vertexModel=ff" FREE_E
refuse key_alone "ffWell= belongs to the full form factor" "$PER
experiment[A] $KIN ffWell=1.2,0.6"
refuse nodes_unmatched "matches no level" "$PER
experiment[A] $KIN vertexModel=ff ffNodes=3/2-:1"
refuse value "expected pw, dw or ff" "$PER
experiment[A] $KIN vertexModel=full"

# (d) ----------------------------------------------------------------------
echo "(d) the same numbers on 1 and 4 threads"
threads_saved="$OMP_NUM_THREADS"
export OMP_NUM_THREADS=1
run t1 "$PER
experiment[A] $KIN vertexModel=ff"
export OMP_NUM_THREADS=4
run t4 "$PER
experiment[A] $KIN vertexModel=ff"
export OMP_NUM_THREADS="$threads_saved"
if ran t1 && ran t4; then
  for f in chiSquared.out "$OUT"; do
    same "$WORK/t1/output/$f" "$WORK/t4/output/$f" && ok "$f byte-identical on 1 and 4 threads" \
      || bad "$f differs between 1 and 4 threads"
  done
fi

# (e) ----------------------------------------------------------------------
echo "(e) the analytic gradient against finite differences"
gradcheck() {  # gradcheck NAME BLOCK: AZURE_GRAD_CHECK of a fit of the widths
  local d="$WORK/$1"
  rm -rf "$d"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$O18/data" "$d/"
  edited "$O18/18O_p_a_thm.azr" "$d/run.azr" awk "$FREE_W" || { bad "$1: no width made free"; return; }
  thm_set "$d/run.azr" "$2" $BOTH || { bad "$1: thm_set"; return; }
  (cd "$d" && AZURE_GRAD_CHECK=exit $RUN "$AZURE2_BIN" fit --use-gradient run.azr < /dev/null > log 2>&1)
  local line
  line="$(grep "components with relative disagreement" "$d/log" | sed -n 1p)"
  if grep -q "gamma 0/4" <<< "$line"; then ok "$1: ${line#*analytic block used; }"
  else bad "$1: ${line:-no [grad-check] line}"; grep "grad-check\] idx" "$d/log" | sed 's/^/        /'; fi
}
gradcheck grad "$PER
experiment[A] $KIN vertexModel=ff"
gradcheck grad_ps "$PER
experiment[A] $KIN ps=hulthen:0-30 psNodes=6 vertexModel=ff"

[ "$fail" = 0 ] && echo "thm_full_form_factor: all checks passed" || echo "thm_full_form_factor: FAILED"
exit "$fail"
