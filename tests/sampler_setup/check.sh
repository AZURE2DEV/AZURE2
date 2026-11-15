#!/usr/bin/env bash
#
# The sampler (the mcmc command) sets up the compound system and the data as the program's run does
# (SetUpProgramCompound, src/RunSetUp.cpp; docs/upgrade/cnuc_split.md section 12.3):
#   - the check files: tests/upgrade/T1 with the compound-nucleus check (line 4 of <config>) and the data check
#     (line 6) set to "file"; the sampler writes checks/compoundnucleus.chk and checks/data.chk, equal to those the
#     calculate command writes;
#   - a refused last part (tests/7Li_p_ay with --no-brune: no standard R-matrix): the sampler's screen text from
#     "Filling Compound Nucleus..." to "Calculation was aborted." equals the calculate command's, the blank line before
#     "Calculation was aborted." included.
#
#   ./tests/sampler_setup/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/sampler_setup.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
MCMC="mcmc --walkers=8 --steps=100 --parameter-spread=0.5 --sampler-threads=1 --overwrite-samples --basis=rwa"

# case NAME SOURCE AZR: $WORK/NAME/run and $WORK/NAME/mcmc, copies of SOURCE with empty output/ and checks/
case_dirs() {
  local d
  for d in run mcmc; do
    mkdir -p "$WORK/$1/$d"
    cp -R "$2/." "$WORK/$1/$d/"
    rm -rf "$WORK/$1/$d/output" "$WORK/$1/$d/checks"
    mkdir -p "$WORK/$1/$d/output" "$WORK/$1/$d/checks"
  done
}

echo "-- the check files (T1, compound nucleus and data checks to file)"
case_dirs checks "$HERE/../upgrade/T1"
for d in run mcmc; do
  awk '/^<config>/ { inside = 1; n = 0; print; next }
       /^<\/config>/ { inside = 0 }
       inside { n++; if (n == 4 || n == 6) sub(/^none/, "file") }
       { print }' "$HERE/../upgrade/T1/T1.azr" > "$WORK/checks/$d/T1.azr"
done
(cd "$WORK/checks/run" && $RUN "$AZURE2_BIN" calculate --no-gui T1.azr < /dev/null > screen.txt 2>&1)
(cd "$WORK/checks/mcmc" && $RUN "$AZURE2_BIN" $MCMC T1.azr < /dev/null > screen.txt 2>&1)
[ -s "$WORK/checks/mcmc/output/samples.mcmc" ] || { bad "the sampler wrote no samples.mcmc"; tail -20 "$WORK/checks/mcmc/screen.txt"; }
for f in compoundnucleus.chk data.chk; do
  if [ ! -s "$WORK/checks/run/checks/$f" ]; then bad "the calculation wrote no checks/$f"
  elif [ ! -f "$WORK/checks/mcmc/checks/$f" ]; then bad "the sampler wrote no checks/$f"
  elif cmp -s "$WORK/checks/run/checks/$f" "$WORK/checks/mcmc/checks/$f"; then
    ok "the sampler writes checks/$f, equal to the calculation's ($(wc -l < "$WORK/checks/mcmc/checks/$f" | tr -d ' ') lines)"
  else bad "the sampler's checks/$f differs from the calculation's"; fi
done

echo "-- a refused last part (7Li_p_ay, --no-brune)"
case_dirs finish "$HERE/../7Li_p_ay"
(cd "$WORK/finish/run" && $RUN "$AZURE2_BIN" calculate --no-gui --no-brune 7Li_p_ay.azr < /dev/null > screen.txt 2>&1)
(cd "$WORK/finish/mcmc" && $RUN "$AZURE2_BIN" $MCMC --no-brune 7Li_p_ay.azr < /dev/null > screen.txt 2>&1; echo $? > exit.txt)
# the set-up's lines: from "Filling Compound Nucleus..." to "Calculation was aborted."
setup_text() { awk '/^Filling Compound Nucleus\.\.\./ { on = 1 } on { print } /^Calculation was aborted\./ { exit }' "$1"; }
setup_text "$WORK/finish/run/screen.txt" > "$WORK/finish/run.txt"
setup_text "$WORK/finish/mcmc/screen.txt" > "$WORK/finish/mcmc.txt"
if ! grep -q '^Calculation was aborted\.' "$WORK/finish/run.txt"; then
  bad "the calculation did not abort"; tail -20 "$WORK/finish/run/screen.txt"
elif cmp -s "$WORK/finish/run.txt" "$WORK/finish/mcmc.txt"; then
  ok "the sampler's text equals the calculation's ($(wc -l < "$WORK/finish/mcmc.txt" | tr -d ' ') lines, a blank line before \"Calculation was aborted.\")"
else
  bad "the sampler's text differs from the calculation's:"; diff "$WORK/finish/run.txt" "$WORK/finish/mcmc.txt"
fi
[ "$(cat "$WORK/finish/mcmc/exit.txt")" = 1 ] && ok "the sampler's exit status is 1" \
  || bad "the sampler's exit status is $(cat "$WORK/finish/mcmc/exit.txt"), not 1"

if [ "$failures" -eq 0 ]; then echo "PASS: the sampler's set-up is the run's"; exit 0; fi
echo "FAIL: $failures check(s)"; exit 1
