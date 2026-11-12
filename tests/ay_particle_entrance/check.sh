#!/usr/bin/env bash
#
# An analyzing power needs a particle entrance: an analyzing_power data set on a beta-delayed reaction
# (16N(b-,a)12C, tests/beta_delayed_alpha) is refused (experiment/src/Build.cpp, E_REACTION_NOT_ALLOWED, the message
# names the data set and says a particle entrance is needed), by --validate (exit 2) and by a run (exit 1, no
# chiSquared.out).  Before, the label rule let it through (a beta is no d, a or g), the engine computed the analyzing
# power as 0, and the analytic gradient treated its points as the beta-delayed spectrum (wrong by about 1e-3).
#
#   ./tests/ay_particle_entrance/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ay_particle_entrance.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

mkdir -p "$WORK/output" "$WORK/checks"
cp -r "$HERE/../beta_delayed_alpha/data" "$WORK/data"
printf '1.0 90 0.1 0.1\n1.2 90 -0.1 0.1\n' > "$WORK/data/ay_beta.dat"
python3 -I - "$HERE/../beta_delayed_alpha/16N.azr" "$WORK/t.azr" <<'PY'
import sys
s = open(sys.argv[1]).read()
add = """dataset ay_beta
    reaction       16N(b-,a)12C
    quantity       analyzing_power
    unit           1
    file           data/ay_beta.dat
    columns        E_lab[MeV] theta_cm[deg] value error
    norm           1 fixed
"""
open(sys.argv[2], "w").write(s.replace("</data>", add + "</data>", 1))
PY
RUN="$(guard_command 120)"
out="$(cd "$WORK" && $RUN "$AZURE2_BIN" --validate t.azr 2>&1)"; rc=$?
if [ "$rc" = 2 ] && printf '%s\n' "$out" | grep -q "ay_beta.*E_REACTION_NOT_ALLOWED.*an analyzing power needs a particle entrance"; then
  ok "--validate refuses it (exit 2): $(printf '%s\n' "$out" | grep -m1 'particle entrance' | cut -c1-160)"
else
  bad "--validate: exit $rc, $(printf '%s\n' "$out" | grep -E 'error|valid' | head -3 | tr '\n' ' ')"
fi
(cd "$WORK" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); rc=$?
if [ "$rc" = 1 ] && grep -q "an analyzing power needs a particle entrance" "$WORK/run.log" && [ ! -f "$WORK/output/chiSquared.out" ]; then
  ok "the run refuses it (exit 1, no chiSquared.out)"
else
  bad "run: exit $rc, $(grep -iE 'error|chi' "$WORK/run.log" | head -3 | tr '\n' ' ')"
fi
if [ "$failures" -ne 0 ]; then echo "  FAIL ($failures)"; exit 1; fi
echo "  PASS"
