#!/usr/bin/env bash
#
# A parameter file's data-set entries are not applied in a run without data (audit item 12, 2026-10-26).
#
# An extrapolation numbered its <segmentsTest> grids 1..N and emitted segment_<k>_energy_shift / _norm parameters for
# them, so a param.sav written by a fit handed DATA segment k's fitted energy shift to TEST segment k by name
# (AZUREParams::ReadUserParameters) and moved the requested grid -- the energy-shift twin of the <targetInt> trap.
# Now those entries are skipped with a note.  In format 3 the extrapolations' parameters are named by their records
# (x1.shift.offset, ...).  Model: tests/sum_lab_frame/models/extrap.azr (13C + alpha, no EC).
#
#   (a) mode 3 with a param.sav that sets x1.shift.offset = 0.01 gives the same grid as mode 3 without it (first
#       energies and values identical); the run log carries the note;
#   (b) a data set's offset (one.shift.offset = 0.01) in a run WITH data (mode 1 on a two-point data set) does apply
#       the shift, as before;
#   (c) pyazr (skipped without the extension): load_sav in extrapolation mode leaves the energies alone.
#
#   ./tests/extrap_sav_segments/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/extrap_sav_segments.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
MODEL="$HERE/../sum_lab_frame/models/extrap.azr"

prep() { mkdir -p "$WORK/$1/output" "$WORK/$1/checks"; cp "$MODEL" "$WORK/$1/t.azr"; }
first() { awk 'NF >= 4 { printf "%.10e %.10e\n", $1, $4; exit }' "$1"; }

# ---- (a) extrapolation with and without the .sav
prep plain
(cd "$WORK/plain" && printf '3\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --use-brune --no-gui --no-readline t.azr > run.log 2>&1)
f="$WORK/plain/output/AZUREOut_aa=2_R=3.extrap"
[ -s "$f" ] && [ -s "$WORK/plain/output/param.par" ] || { bad "plain extrapolation wrote no output"; tail -5 "$WORK/plain/run.log"; exit 1; }
grep -q "^ *x1.shift.offset " "$WORK/plain/output/param.par" || { bad "param.par has no x1.shift.offset line"; head -3 "$WORK/plain/output/param.par"; exit 1; }
awk '$1 == "x1.shift.offset" { $2 = "0.01"; print; next } { print }' "$WORK/plain/output/param.par" > "$WORK/shift.sav"
prep sav; cp "$WORK/shift.sav" "$WORK/sav/"
# mode 3 asks for the uncertainty band (n) before the parameter file
(cd "$WORK/sav" && printf '3\nn\nshift.sav\n\n7\n' | run_limited 300 "$AZURE2_BIN" --use-brune --no-gui --no-readline t.azr > run.log 2>&1)
g="$WORK/sav/output/AZUREOut_aa=2_R=3.extrap"
[ -s "$g" ] || { bad "extrapolation with the .sav wrote no output"; tail -5 "$WORK/sav/run.log"; exit 1; }
grep -q "Reading User Parameter File" "$WORK/sav/run.log" || bad "the .sav was not read"
# energies identical; the values to 1e-6 (the file carries the level parameters with 7 digits)
cmp_out="$(paste <(awk 'NF >= 4 { print $1, $4 }' "$f") <(awk 'NF >= 4 { print $1, $4 }' "$g") | awk '
  { if ($1 != $3) { print "GRID"; exit } d = $2 - $4; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300; r = d / m; if (r > w) w = r; n++ }
  END { printf "%.2e %d", w, n }')"
if [ "$cmp_out" = "GRID" ]; then bad "(a) the .sav moved the extrapolation grid: $(first "$f") vs $(first "$g")"
elif awk -v w="${cmp_out%% *}" 'BEGIN { exit !(w + 0 <= 1e-6) }'; then ok "(a) the grid is unchanged with x1.shift.offset = 0.01 in the file (first energy $(first "$g" | cut -d' ' -f1) MeV; values to ${cmp_out%% *}, ${cmp_out##* } points)"
else bad "(a) values differ by ${cmp_out%% *} with the .sav"; fi
if grep -q "segment_\* entries of shift.sav (data-segment normalizations and energy shifts) are not applied" "$WORK/sav/run.log"; then ok "(a) the run log says the segment_* entries were skipped"; else bad "(a) no note about the skipped entries"; fi

# ---- (b) with data the shift does apply
prep data; mkdir -p "$WORK/data/data"
printf '6.0  30  1.0  0.1\n6.2  30  1.0  0.1\n' > "$WORK/data/data/one.dat"
python3 - "$WORK/data/t.azr" <<'EOF'
import sys
p = sys.argv[1]; s = open(p).read()
one = ("dataset one\n    reaction       13C(a,n)16O[6.0494]\n    quantity       dsigma_domega\n"
       "    unit           b/sr\n    file           data/one.dat\n    columns        E_lab[MeV] theta_lab[deg] value error\n"
       "    energy_window  5.9 6.4 MeV\n    angle_window   0 180 deg\n    norm           1.0 fixed\n")
assert s.count("<data>\n</data>") == 1
s = s.replace("<data>\n</data>", "<data>\n" + one + "</data>")
open(p, "w").write(s)
EOF
(cd "$WORK/data" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --use-brune --no-gui --no-readline t.azr > run.log 2>&1)
awk '$1 == "one.shift.offset" { $2 = "0.01"; print; next } { print }' "$WORK/data/output/param.par" > "$WORK/data/shift.sav"
e0="$(awk 'NF >= 5 { print $1; exit }' "$WORK/data/output/AZUREOut_aa=2_R=3.out")"
(cd "$WORK/data" && printf '1\nshift.sav\n\n7\n' | run_limited 300 "$AZURE2_BIN" --use-brune --no-gui --no-readline t.azr > run2.log 2>&1)
e1="$(awk 'NF >= 5 { print $1; exit }' "$WORK/data/output/AZUREOut_aa=2_R=3.out")"
if awk -v a="$e0" -v b="$e1" 'BEGIN { d = b - a; exit !(d > 0.005 && d < 0.009) }'; then ok "(b) with data the file's shift moves the point as before (E_cm $e0 -> $e1)"; else bad "(b) with data the shift did not apply as expected ($e0 -> $e1)"; fi

# ---- (c) pyazr load_sav in extrapolation mode
PY="$(command -v python3 || true)"
if [ -n "$PY" ] && (cd /tmp && PYTHONPATH="$ROOT" "$PY" -c "import pyazr._azure2" > /dev/null 2>&1); then
  prep py; cp "$WORK/shift.sav" "$WORK/py/"
  out="$(cd "$WORK/py" && PYTHONPATH="$ROOT" "$PY" - <<'EOF' 2>&1
import numpy as np
from pyazr import Session
with Session("t.azr", cwd=".", data_mode=False) as m:
    e0 = np.asarray(m.calculate_energies(m.params_rwa)[0], float)
    x = m.load_sav("shift.sav")
    e1 = np.asarray(m.calculate_energies(x)[0], float)
print("MAXDIFF", np.max(np.abs(e1 - e0)))
EOF
)"
  d="$(echo "$out" | awk '/^MAXDIFF/ { print $2 }')"
  if [ -n "$d" ] && awk -v d="$d" 'BEGIN { exit !(d + 0 < 1e-12) }'; then ok "(c) pyazr load_sav in extrapolation mode leaves the grid alone (max |dE| $d)"; else bad "(c) pyazr load_sav in extrapolation mode: $(echo "$out" | tail -3)"; fi
else
  echo "  skip  (c) pyazr extension not importable"
fi

echo
if [ "$failures" -eq 0 ]; then echo "extrap_sav_segments: PASS"; else echo "extrap_sav_segments: FAIL ($failures)"; exit 1; fi
