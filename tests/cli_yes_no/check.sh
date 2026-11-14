#!/usr/bin/env bash
# An answer cannot reach the wrong question: every answer of a run is a named option of its command.
#
# The text menu asked "Calculate cross-section uncertainty band? (y/n)" BEFORE the parameter file in modes 2 and 3, and
# a script that left that line out ('3\nparam.sav\n\n7') gave the parameter-file name to the y/n question; the run then
# used the <levels> of the .azr instead of the file, without an error (found on a 7Be A_y test-segment calculation and
# in the 12C+a "run_again" fits of July-August 2026).  The menu is gone: the band is the flag --covariance-band, the
# parameter file --parameter-file=FILE, and a command line that leaves out a required answer, gives a value to a flag
# or names a file that does not exist is refused with exit 2 and the three parts (what, what the command takes, how to
# fix) before anything is calculated.
#
# Checked on tests/7Li_p_ay (one A_y test segment) with a parameter file written by a calculation:
#   (a) extrapolate --parameter-file=fit.par reads the file and writes the .extrap, exit 0;
#   (b) the file name where the band's answer was (--covariance-band=fit.par) is refused: the band takes no value;
#   (c) the file name without its option (a second project file) is refused;
#   (d) a missing required answer is refused, naming it: minos without --delta-chi2, rate without its pairs and
#       without its temperatures, mcmc without the sampler's numbers;
#   (e) a run without a command is refused with the list of commands, whatever its standard input holds (the old
#       answers, nothing);
#   (f) none of the refused runs calculated anything.
#
# Usage:
#   ./tests/cli_yes_no/check.sh path/to/AZURE2
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cli_yes_no.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

PROJ=7Li_p_ay
mkdir -p "$WORK/base/output" "$WORK/base/checks"
cp -R "$HERE/../$PROJ/data" "$WORK/base/data"
cp "$HERE/../$PROJ/$PROJ.azr" "$WORK/base/"
# A parameter file to point at: the param.par a calculation writes.
(cd "$WORK/base" && "$AZURE2_BIN" calculate "$PROJ.azr" < /dev/null > calculate.log 2>&1)
[ -s "$WORK/base/output/param.par" ] && cp "$WORK/base/output/param.par" "$WORK/fit.par" \
  || bad "the calculation wrote no param.par (see the log)"

# run NAME INPUT ARGS...: AZURE2 ARGS on a fresh copy, INPUT (printf format) on its standard input; sets $rc
run() {
  local name="$1" input="$2" dir="$WORK/$1"; shift 2
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$WORK/base/data" "$WORK/base/$PROJ.azr" "$dir/"
  cp "$WORK/fit.par" "$dir/fit.par" 2>/dev/null
  (cd "$dir" && printf "$input" | "$AZURE2_BIN" "$@" > run.log 2>&1)
  rc=$?
}
said() { grep -qF -- "$2" "$WORK/$1/run.log"; }
extrap() { ls "$WORK/$1"/output/*.extrap 2>/dev/null | head -1; }
# refused NAME WANT: exit 2, the ERROR containing WANT, its fix lines, and nothing calculated
refused() {
  if [ "$rc" = 2 ] && said "$1" "ERROR: " && said "$1" "$2" && [ "$(grep -c '^       ' "$WORK/$1/run.log")" -ge 2 ] &&
     [ -z "$(find "$WORK/$1/output" -type f)" ] && ! said "$1" "Calling AZURE2"; then
    ok "$3: refused (exit 2, three parts), nothing calculated"
  else
    bad "$3: exit $rc, $(grep -m1 ERROR "$WORK/$1/run.log" || echo 'no ERROR'); outputs: $(ls "$WORK/$1/output" | wc -l)"
  fi
}

# (a)
run right "" extrapolate --parameter-file=fit.par "$PROJ.azr"
said right "Reading User Parameter File" && [ "$rc" = 0 ] && [ -n "$(extrap right)" ] \
  && ok "(a) extrapolate --parameter-file=fit.par: the parameter file is read, .extrap written, exit 0" \
  || bad "(a) extrapolate --parameter-file=fit.par: exit $rc, $(said right 'Reading User Parameter File' && echo read || echo 'not read')"

# (b), (c)
run band_value "" extrapolate --covariance-band=fit.par "$PROJ.azr"
refused band_value "--covariance-band=fit.par: --covariance-band takes no value." \
  "(b) the file name given to the band (--covariance-band=fit.par)"
run bare_file "" extrapolate fit.par "$PROJ.azr"
refused bare_file "Two project files: fit.par and $PROJ.azr." "(c) the file name without its option"

# (d)
run minos "" minos "$PROJ.azr"
refused minos "The command minos needs --delta-chi2=X." "(d) minos without --delta-chi2"
run rate_pairs "" rate --temperature-file=fit.par "$PROJ.azr"
refused rate_pairs "The command rate needs --rate-entrance=N, --rate-exit=N." "(d) rate without its pairs"
run rate_temperatures "" rate --rate-entrance=1 --rate-exit=-1 "$PROJ.azr"
refused rate_temperatures "needs --temperature-min=T, --temperature-max=T, --temperature-step=T" \
  "(d) rate without its temperatures"
if grep -q '^	mcmc:' <<< "$("$AZURE2_BIN" --help 2>&1)"; then
  run mcmc "" mcmc --walkers=4 "$PROJ.azr"
  refused mcmc "needs --steps=N, --sampler-threads=N, --basis=rwa|physical" "(d) mcmc without the sampler's numbers"
fi

# (e)
run no_command '3\nn\nfit.par\n\n7\n' "$PROJ.azr"
refused no_command "No command:" "(e) the old menu answers on the input of a run without a command"
said no_command "calculate, fit, extrapolate, minos, rate" && ok "  ... the refusal lists the commands" \
  || bad "  ... the refusal does not list the commands"
run no_command_empty "" "$PROJ.azr"
refused no_command_empty "No command:" "(e) a run without a command and an empty input"

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
