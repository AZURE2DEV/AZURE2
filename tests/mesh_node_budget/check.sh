#!/usr/bin/env bash
#
# A mesh too large for memory is refused before it is made (kMeshNodeBudget, include/MeshPlanner.h): the message names
# the data set and point that ask for the most nodes, their kernel's width, the planned width and the step, and what
# to do.  Model: tests/thin_layer_limit's (a 15N layer of 1e14 atoms/cm^2, a kernel of 6e-8 MeV) with a varied energy
# shift on its two folded data sets:
#   1. a 3 keV prior on the shift (menu 1): the planned support, 3 prior widths each side, at the layer grid's step
#      asks for about 3e7 nodes; refused at the set-up in about a second (exit 1, "Calculation was aborted.").  Before,
#      it took 28 GB in 3 minutes and did not finish;
#   2. no prior, the data a flat chi^2 (placeholders 1 +- 1, the model 1e-9), 1e16 atoms/cm^2 (menu 2): MIGRAD steps the
#      shift by MeV and the re-plan asks for 1.7e7 nodes; the fit ends with the message (the hand-over runs inside
#      Minuit's parallel gradient, so the process ends by std::terminate, exit 134) instead of std::bad_alloc after
#      tens of GB.
#
#   ./tests/mesh_node_budget/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
ulimit -c 0  # case 2 ends by std::terminate: no core file
WORK="$(mktemp -d "${TMPDIR:-/tmp}/mesh_node_budget.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
MESSAGE="the mesh of a compound nucleus would have more than 1e+06 nodes"

# case NAME DENSITY SHIFT MENU: the model in $WORK/NAME, run with the menu answers MENU (at most 30 s); sets $rc, $secs
case_run() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$HERE/../thin_layer_limit/models/data" "$d/data"
  python3 -I - "$HERE/../thin_layer_limit/models/thin_limit.azr" "$d/m.azr" "$2" "$3" <<'PY'
import sys
text, out, done = open(sys.argv[1]).read().replace("areal_density  1e+14 ", "areal_density  %s " % sys.argv[3]), [], 0
for line in text.split("\n"):
    out.append(line)
    if line.startswith("    norm ") and done < 2:
        out.append("    shift          " + sys.argv[4])
        done += 1
open(sys.argv[2], "w").write("\n".join(out))
PY
  local t0=$SECONDS
  (cd "$d" && printf "$4" | $(guard_command 30) "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
  rc=$?
  secs=$((SECONDS - t0))
}

case_run prior 1e+14 "offset 0.02 MeV fit prior 0 0.003" '1\n\n\n7\n'
if [ "$rc" = 1 ] && grep -q "^ERROR: $MESSAGE" "$WORK/prior/run.log" &&
   grep -q "data set 'thin', the point at E_lab 0.600000 MeV: its kernel is .* wide, its planned support .* at a step of" "$WORK/prior/run.log" &&
   grep -q "Give the energy shift a narrower prior, or fix it" "$WORK/prior/run.log" && [ "$secs" -le 10 ]; then
  ok "a 3 keV prior on a 1e14 layer: refused at the set-up in ${secs} s"
else
  bad "a 3 keV prior: exit $rc after ${secs} s: $(grep -E 'ERROR|what' "$WORK/prior/run.log" | sed -n '1,2p' | cut -c1-200)"
fi
case_run flat 1e+16 "offset 0.02 MeV fit" '2\nn\n\n\n7\n'
if grep -q "$MESSAGE" "$WORK/flat/run.log" && ! grep -q "bad_alloc" "$WORK/flat/run.log" && [ "$secs" -le 20 ]; then
  ok "a free shift on a flat chi^2: the fit ends with the message in ${secs} s (exit $rc)"
else
  bad "a flat fit: exit $rc after ${secs} s: $(grep -E 'ERROR|what|bad_alloc' "$WORK/flat/run.log" | sed -n '1,2p' | cut -c1-200)"
fi
if [ "$failures" -ne 0 ]; then echo "  FAIL ($failures)"; exit 1; fi
echo "  PASS"
