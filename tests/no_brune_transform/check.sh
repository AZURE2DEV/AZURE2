#!/usr/bin/env bash
#
# --no-brune on a project with no standard-R-matrix equivalent: a clean refusal, not a crash.
#
# The standard formalism needs the input (observed) parameters transformed to formal ones (CNuc::TransformIn), a
# generalized eigenproblem N a = E M a that has a real solution only for a positive-definite M = 1 - sum gamma
# gamma' dS/dE (Brune, PRC 66, 044611).  tests/7Li_p_ay (a Brune fit, 0+ background poles with partial widths up to
# 135 MeV) makes M indefinite.  GSL reported the failure, AZURE2's error handler ignores that code, and before
# 2026-11-01 the unset eigenvalues were used (formal energies ~1e-320); TransformOut then threw "z is not finite in
# log_Gamma" and the CLI aborted with a core dump (pyazr: RuntimeError at load).  Checked:
#   1. 7Li_p_ay --no-brune, mode 1: the run ends without a crash, says that no standard R-matrix corresponds to the
#      0+ levels and aborts the calculation (no chiSquared.out);
#   2. the same project with Brune (default) still runs;
#   3. 13N --no-brune (a valid transform) gives its usual chi2, 336480.
#
#   ./tests/no_brune_transform/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/no_brune.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME PROJECT [FLAGS...]: mode 1 of tests/PROJECT under $WORK/NAME; sets $rc.
run() {
  local dir="$WORK/$1" proj="$2"; shift 2
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/../$proj/data" "$dir/data"
  cp "$HERE/../$proj/$proj.azr" "$dir/"
  (cd "$dir" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline "$@" "$proj.azr" > run.log 2>&1)
  rc=$?
}
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2>/dev/null; }

run li7 7Li_p_ay --no-brune
[ "$rc" -lt 128 ] && ok "7Li_p_ay --no-brune: no crash (exit $rc)" || bad "7Li_p_ay --no-brune: killed by a signal (exit $rc)"
grep -q "No standard R-matrix (constant boundary conditions) corresponds to the levels of J^pi = 0+" "$WORK/li7/run.log" \
  && ok "7Li_p_ay --no-brune: says no standard R-matrix corresponds to the 0+ levels" \
  || bad "7Li_p_ay --no-brune: no explanation in the log"
grep -q "Calculation was aborted" "$WORK/li7/run.log" && [ ! -s "$WORK/li7/output/chiSquared.out" ] \
  && ok "7Li_p_ay --no-brune: calculation aborted, no chiSquared.out" \
  || bad "7Li_p_ay --no-brune: not aborted cleanly"

run li7b 7Li_p_ay
t="$(total li7b)"
[ -n "$t" ] && [ "$t" != nan ] && ok "7Li_p_ay (Brune): runs, chi2 $t" || bad "7Li_p_ay (Brune): chi2 '$t'"

run n13 13N --no-brune
t="$(total n13)"
[ "$t" = 336480 ] && ok "13N --no-brune: chi2 336480 as before" || bad "13N --no-brune: chi2 '$t' (expected 336480)"

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
