#!/usr/bin/env bash
#
# MCMC from a starting point with Park's J <= 0.
#
# Under --use-park the posterior is zero where a level has J = 1 - sum gamma^2
# dS/dE < 0 (AZURECalcMCMC rejects such points).  A run started there used to
# say nothing and sample with acceptance 0: every walker stayed at its start.
# Now (AZURECalcMCMC::RunMCMCSampling):
#   - a start with J <= 0 for some level is reported, the level named with
#     its J ("WARNING: the MCMC starting point is outside Park's parameter
#     space");
#   - if no walker of the initial ensemble has a finite posterior, the run is
#     refused ("MCMC not started"), writes no samples.mcmc and the CLI exits
#     non-zero;
#   - if some walkers are inside, it runs (with the warning).
#
# Model: tests/identical_pp_res with the 0+ level's width raised past the
# bound (as tests/park_formalism section 4): 51 MeV (J = -0.54, the whole
# 5 % ball outside: refused, in physical parameters and in RWA) and 33.5 MeV
# (J = -0.009, the ball straddles J = 0: runs); the project as it is (J > 0:
# no warning) as the control.  24 walkers, 100 steps, a few seconds a run.
#
#   ./tests/mcmc_park_start/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-600}")"
export OMP_NUM_THREADS=1

WORK="$(mktemp -d "${TMPDIR:-/tmp}/mcmc_park.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
fail=0
SRC="$HERE/../identical_pp_res"

# stage NAME WIDTH: the project with the 0+ level's width WIDTH (eV), or as
# it is for an empty WIDTH.
stage() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -R "$SRC/data" "$WORK/$1/data"
  # First " 300000.0 " only (the 0+ level); awk, not sed -i (BSD/macOS).
  tr -d '\r' < "$SRC/identical_pp_res.azr" |
    awk -v w="$2" 'w != "" && !done && index($0, " 300000.0 ") { sub(/ 300000\.0 /, " " w " "); done = 1 } { print }' \
    > "$WORK/$1/run.azr"
}
# mcmc NAME RWA(yes|no): the sampler under --use-park; NAME/status, NAME/log.
mcmc() {
  local d="$WORK/$1"
  (cd "$d" && printf '6\n\n24\n100\n\n0.5\n1\n%s\nyes\n' "$2" |
     $RUN "$AZURE2_BIN" --no-gui --no-readline --use-park run.azr 2>&1 |
     head -c 20000000 | tr '\r' '\n' > log
   echo "${PIPESTATUS[1]}" > status)
}
WARN="WARNING: the MCMC starting point is outside Park's parameter space"
LEVEL="J^pi=0+, E=0.0000 MeV (J-group 1, level 1): J = -"

refused() {  # NAME
  local d="$WORK/$1"
  grep -qF "$WARN" "$d/log" && grep -qF "$LEVEL" "$d/log" &&
    ok "$1: start reported, the level named ($(grep -F "$LEVEL" "$d/log" | head -1 | sed 's/^ *//'))" ||
    bad "$1: no warning naming the level with J <= 0"
  grep -q "MCMC not started" "$d/log" && ok "$1: refused (no walker inside the support)" ||
    bad "$1: not refused"
  [ "$(cat "$d/status")" != 0 ] && ok "$1: exit status $(cat "$d/status")" || bad "$1: exit status 0"
  [ ! -f "$d/output/samples.mcmc" ] && ok "$1: no samples.mcmc written" || bad "$1: samples.mcmc written"
  ! grep -q "Acceptance fraction" "$d/log" && ok "$1: no step taken" || bad "$1: the sampler ran"
}
ran() {  # NAME WARNED(1|0)
  local d="$WORK/$1"
  if [ "$(cat "$d/status")" = 0 ] && grep -q "MCMC sampling completed successfully" "$d/log" &&
     [ -f "$d/output/samples.mcmc" ]; then ok "$1: runs (status 0, samples.mcmc)"
  else bad "$1: did not run (status $(cat "$d/status"))"; tail -5 "$d/log" | sed 's/^/        /'; fi
  local acc
  acc="$(awk '/Acceptance fraction:/ { print $NF }' "$d/log" | tail -1)"
  awk -v a="$acc" 'BEGIN { exit !(a + 0 > 0) }' && ok "$1: acceptance fraction $acc" ||
    bad "$1: acceptance fraction '$acc'"
  if [ "$2" = 1 ]; then
    grep -qF "$WARN" "$d/log" && grep -qF "$LEVEL" "$d/log" && ok "$1: start with J <= 0 reported" ||
      bad "$1: no start warning"
  else
    ! grep -qF "$WARN" "$d/log" && ok "$1: no start warning (J > 0)" || bad "$1: spurious start warning"
  fi
}

stage wall 51000000;   mcmc wall no
stage wall_rwa 51000000; mcmc wall_rwa yes
stage edge 33500000;   mcmc edge no
stage inside "";       mcmc inside no

echo "== J = -0.54 at the start, physical parameters"
refused wall
echo "== J = -0.54 at the start, reduced width amplitudes"
refused wall_rwa
echo "== J = -0.009 at the start (part of the ensemble inside)"
ran edge 1
echo "== J > 0 (the project as it is)"
ran inside 0

[ "$fail" -eq 0 ] && echo "  PASS" || echo "  FAIL"
[ "$fail" -eq 0 ]
