#!/bin/bash
# The node tables at run time (include/NodeTable.h, docs/upgrade/energy_grid.md section 2.6): with AZURE_NODE_CHECK=1
# every hand-over and every point pass checks that each unit of each plan task is at the node its slot names -- on the
# stage table or on the working copy's overflow table -- with the node's channel rows bit for bit (a breach stops the
# run).  Runs: a short MIGRAD fit of T2 (re-placed grids and shifted folded segments: the overflow table is used), the
# calculations of shift_sqrt_subpoints/models/sqrt.azr (folded points under a sqrt(E) shift: on the fixed nodes of the
# shared mesh, so none on the overflow table), 18O_p_a_thm (THM nodes)
# and natural_neon (folded points of several compound nuclei), each with and without the check: each run exits 0,
# every check line reports 0 units off the tables, 0 on a wrong node and 0 rows that differ, the fit uses the overflow
# table and the shifted calculation does not, and the calculations' output files are byte for byte the same.  Shown failing with a slot of a
# re-placed grid pointing at the next node (injected once: "on a wrong node").  Each check also checks the sharing
# groups of the slots (include/SlotGroup.h: every range in one block, the writer an earlier reader, one group per
# block, the slots of a point on a shared mesh at the nodes its plan gives), reported on its own line; shown failing
# with one point's slot range shifted by one site (injected once: "reads of another group", "off their mesh plan").
# Each check also checks that every reader of slots without an angle holds its angular row and Coulomb amplitudes at
# each slot's energy (WorkingCopy::CheckReaderRows; hybrid_potential: folded elastic data); shown failing with
# the rows of a Q-coefficient particle exit taken at the point's energy instead of the slot's (injected once: the T2
# fit reports rows that differ).  THM data sets of one experiment with weight tables of equal rows share their
# blocks, and with other rows do not (shown failing with the tables compared by pointer: "blocks 9 (T2), 10 (equal
# tables), 10 (other rows)").
#   node_check.sh path/to/AZURE2
set -u
export LC_ALL=C OMP_NUM_THREADS=2
[ $# -ge 1 ] || { echo "usage: node_check.sh path/to/AZURE2"; exit 2; }
BIN="$1"
TESTS="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail=0
ok()  { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; fail=1; }
# run NAME SRC_DIR AZR CHECK MENU [FLAGS] -> $WORK/NAME.CHECK (exit code in rc)
run() {
  local d="$WORK/$1.$4"
  mkdir -p "$d" && cp -r "$2/." "$d/" && rm -rf "$d/output" "$d/checks" && mkdir -p "$d/output" "$d/checks"
  (cd "$d" && printf "$5" | AZURE_NODE_CHECK=$4 "$BIN" --no-gui --no-readline ${6:-} "$3" > run.log 2>&1; echo $? > rc)
}
# clean NAME: the checked run's check lines all clean, and how many of them used the overflow table
clean() {
  local log="$WORK/$1.1/run.log" dirty
  lines="$(grep -ac 'AZURE_NODE_CHECK: [0-9]' "$log")"
  dirty="$(grep -a 'AZURE_NODE_CHECK: [0-9]' "$log" | grep -vc ' 0 off the tables, 0 on a wrong node; channel rows 0 differ')"
  overflow="$(grep -a 'AZURE_NODE_CHECK: [0-9]' "$log" | grep -vc ' 0 on the overflow table')"
  # the sharing groups of the slots (WorkingCopy::CheckSlotGroups), one line per check
  groups="$(grep -ac 'AZURE_NODE_CHECK: slot groups:' "$log")"
  grouped="$(grep -a 'AZURE_NODE_CHECK: slot groups:' "$log" | grep -c ' 0 ranges outside a block, 0 slots with a wrong writer, 0 reads of another group, 0 ranges off their mesh plan')"
  # the readers' rows and Coulomb amplitudes at their slots (WorkingCopy::CheckReaderRows), one line per check
  readers="$(grep -a 'AZURE_NODE_CHECK: reader rows:' "$log" | awk '{ n += $4 } END { print n + 0 }')"
  rowlines="$(grep -ac 'AZURE_NODE_CHECK: reader rows:' "$log")"
  rowsclean="$(grep -a 'AZURE_NODE_CHECK: reader rows:' "$log" | grep -c ' slots, 0 differ$')"
  [ "$lines" -gt 0 ] && [ "$dirty" -eq 0 ] && [ "$groups" -eq "$lines" ] && [ "$grouped" -eq "$groups" ] &&
    [ "$rowlines" -eq "$lines" ] && [ "$rowsclean" -eq "$rowlines" ]
}
CALC='1\n\n\n7\n'
FIT='2\n\n\n\n7\n'
for p in "sqrt:$TESTS/shift_sqrt_subpoints/models:sqrt.azr:none" "thm:$TESTS/18O_p_a_thm:18O_p_a_thm.azr:0" \
         "natural_neon:$TESTS/natural_neon:natural_neon.azr:0" "elastic:$TESTS/hybrid_potential:hybrid_potential.azr:0" \
         "ay:$TESTS/13N:13N.azr:0"; do
  IFS=: read -r name src azr wantOverflow <<< "$p"
  run "$name" "$src" "$azr" 0 "$CALC"
  run "$name" "$src" "$azr" 1 "$CALC"
  r0="$(cat "$WORK/$name.0/rc")"; r1="$(cat "$WORK/$name.1/rc")"
  n="$(ls "$WORK/$name.0/output" | wc -l | tr -d ' ')"
  if [ "$r0" = 0 ] && [ "$r1" = 0 ] && [ "$n" -gt 0 ] && clean "$name" && { [ "$wantOverflow" = 0 ] || { [ "$wantOverflow" = none ] && [ "$overflow" -eq 0 ]; } ||
                                           { [ "$wantOverflow" = 1 ] && [ "$overflow" -gt 0 ]; }; } &&
     diff -rq "$WORK/$name.0/output" "$WORK/$name.1/output" > /dev/null; then
    ok "$name: $lines checks clean ($overflow with the overflow table, $readers reader checks of slots), its $n output files identical"
  else
    bad "$name: exit $r0 / $r1, $n files; $(grep -a AZURE_NODE_CHECK "$WORK/$name.1/run.log" | grep -v ' 0 off the tables, 0 on a wrong node; channel rows 0 differ' | head -2)"
  fi
done
# Two data sets of one reaction on one mesh share their sites (SlotGroup.h: the key holds the observable and the
# target's line, not the segment): hybrid_potential's meyer_84 and meyer_114 read one block of 2413 slots.
if grep -aq 'AZURE_NODE_CHECK: slot groups: 1 blocks, 2413 slots' "$WORK/elastic.1/run.log"; then
  ok "elastic: two data sets of one reaction share one block"
else
  bad "elastic: two data sets of one reaction do not share their block: $(grep -a -m1 'slot groups' "$WORK/elastic.1/run.log")"
fi
# The THM data sets of one experiment share their sites (TermSettings.h, SameObservable: the experiment's attachments
# by identity, a data set's weight table by value): T2's s11 and s12 (experiment LC09, one beam line) weighted by two
# tables with equal rows keep T2's blocks; with other rows the second reads a block of its own.
. "$TESTS/lib/thm_records.sh"
for v in same other; do
  mkdir -p "$WORK/w$v" && cp -r "$TESTS/upgrade/T2/." "$WORK/w$v/"
  printf '0 2\n1 3\n' > "$WORK/w$v/wa.dat"
  if [ $v = same ]; then printf '0 2\n1 3\n' > "$WORK/w$v/wb.dat"; else printf '0 2\n1 4\n' > "$WORK/w$v/wb.dat"; fi
  thm_weight "$WORK/w$v/T2.azr" s11_thm_hoes_a wa.dat && thm_weight "$WORK/w$v/T2.azr" s12_thm_hoes_b wb.dat ||
    bad "weights: the data sets are not in T2"
  run "w$v" "$WORK/w$v" T2.azr 1 "$CALC"
done
b0="$(grep -a -m1 'AZURE_NODE_CHECK: slot groups:' "$WORK/T2.1/run.log" 2>/dev/null | awk '{print $4}')"
[ -n "$b0" ] || { run T2 "$TESTS/upgrade/T2" T2.azr 1 "$CALC"; b0="$(grep -a -m1 'AZURE_NODE_CHECK: slot groups:' "$WORK/T2.1/run.log" | awk '{print $4}')"; }
bs="$(grep -a -m1 'AZURE_NODE_CHECK: slot groups:' "$WORK/wsame.1/run.log" | awk '{print $4}')"
bo="$(grep -a -m1 'AZURE_NODE_CHECK: slot groups:' "$WORK/wother.1/run.log" | awk '{print $4}')"
if [ -n "$b0" ] && [ "$bs" = "$b0" ] && [ "$bo" = "$((b0 + 1))" ] && clean wsame && clean wother; then
  ok "THM weights: equal tables share the blocks ($bs = $b0), other rows take one more ($bo)"
else
  bad "THM weights: blocks $b0 (T2), $bs (equal tables), $bo (other rows)"
fi
run T2fit "$TESTS/upgrade/T2" T2.azr 1 "$FIT" --max-calls=30
r1="$(cat "$WORK/T2fit.1/rc")"
if [ "$r1" = 0 ] && clean T2fit && [ "$overflow" -gt 0 ]; then
  ok "T2 fit (30 calls): $lines checks clean, $overflow with the overflow table, $readers reader checks of slots"
else
  bad "T2 fit: exit $r1; $(grep -a AZURE_NODE_CHECK "$WORK/T2fit.1/run.log" | grep -v ' 0 off the tables, 0 on a wrong node; channel rows 0 differ' | head -2)"
fi
[ $fail -eq 0 ] && echo "node_check: PASS" || echo "node_check: FAIL"
exit $fail
