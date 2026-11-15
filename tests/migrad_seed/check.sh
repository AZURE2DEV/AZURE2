#!/usr/bin/env bash
#
# The start of a --use-gradient MIGRAD fit.
#
# MIGRAD starts from the inverse Hessian diag(1 / G2) and takes the step -g / G2.  With the analytic Gauss-Newton G2
# that first step was unbounded for a parameter the data hardly constrain (G2 near zero) and for a limited parameter
# near its limit (Minuit2 leaves the term g d2x/du2 out of the internal G2).  On a 547-parameter 12C+a fit
# (2026-11-15) the unit step went from chi^2 1.2e5 to 4.4e6; when the 12-evaluation line search then finds no
# improvement, Minuit2 starts a numerical HESSE of n(n-1)/2 function evaluations with no gradient call, which looks
# like a stalled fit.  Now (AZURECalc::BoundSeedCurvature) the seed curvature is raised so that no parameter moves
# more than ten of its initial step sizes in the first step, and if the first line search still fails the fit is
# restarted once with Minuit2's numerical seed.  Checked on tests/13N (14 free parameters):
#
#   1. the fit converges, says which share of the first steps was bounded, and reaches the chi^2 of the unbounded
#      seed (AZURE_SEED_BOUND=0);
#   2. with the first line search made to fail (AZURE_SEED_STALL_TEST=1) the run restarts with the numerical seed
#      and converges to the same chi^2.
#
#   ./tests/migrad_seed/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/migradseed.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# fit NAME: a mode-2 gradient fit of a private copy of tests/13N under $WORK/NAME (environment from the caller).
fit() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/../13N/data" "$dir/data"
  cp "$HERE/../13N/13N.azr" "$dir/t.azr"
  (cd "$dir" && printf '2\n\n\n7\n' | "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
}
chi() { awk '/Total-Chi-Squared:/ { print $2; exit }' "$WORK/$1/output/chiSquared.out" 2>/dev/null; }
same() { awk -v a="$1" -v b="$2" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(a > 0 && b > 0 && d < 0.05) }'; }
converged() { grep -q "^MIGRAD converged: valid minimum" "$WORK/$1/run.log"; }

# ---- 1. bounded seed against the unbounded one ------------------------------
fit bound
AZURE_SEED_BOUND=0 fit free
converged bound && ok "bounded seed: MIGRAD converged" || bad "bounded seed: no convergence line"
converged free && ok "AZURE_SEED_BOUND=0: MIGRAD converged" || bad "AZURE_SEED_BOUND=0: no convergence line"
grep -q "^MIGRAD seed: the first step of [0-9]* of [0-9]* free parameters was bounded" "$WORK/bound/run.log" \
  && ok "bounded seed: the log says how many first steps were bounded" \
  || bad "bounded seed: no 'MIGRAD seed:' line"
grep -q "^MIGRAD seed:" "$WORK/free/run.log" && bad "AZURE_SEED_BOUND=0 still bounds the seed" \
  || ok "AZURE_SEED_BOUND=0 switches the bound off"
cb="$(chi bound)"; cf="$(chi free)"
same "$cb" "$cf" && ok "same minimum: chi^2 $cb (bounded) and $cf (unbounded)" \
  || bad "different minima: chi^2 '$cb' (bounded) and '$cf' (unbounded)"

# ---- 2. a failed first line search restarts with the numerical seed ---------
AZURE_SEED_STALL_TEST=1 fit stall
grep -q "restarting with Minuit2's numerical seed" "$WORK/stall/run.log" \
  && ok "failed first step: the run restarts with the numerical seed" \
  || bad "failed first step: no restart line in the log"
converged stall && ok "failed first step: MIGRAD converged after the restart" \
  || bad "failed first step: no convergence line"
cs="$(chi stall)"
same "$cs" "$cb" && ok "failed first step: same minimum, chi^2 $cs" \
  || bad "failed first step: chi^2 '$cs' against '$cb'"

echo
if [ "$failures" -eq 0 ]; then echo "PASS: MIGRAD seed"; else echo "FAIL: $failures check(s)"; fi
exit "$failures"
