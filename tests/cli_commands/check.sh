#!/usr/bin/env bash
#
# AZURE2 <command> [options] FILE runs its calculation, and a wrong command line is refused.
#
# The commands and their options are one table (src/RunRequest.cpp); a run goes through one preparation
# (src/RunPreparation.cpp).  Checked on the binary as built:
#
#   1. calculate and extrapolate (T1), fit with the band and its scaling, minos (identical_pp_res), calculate from
#      the fit's parameter file (relative to the project's directory: the same chi-squared as the fit), rate (the
#      closed 12C(p,g) model, both temperature forms) and, in a build with the sampler, mcmc (13N_fast, both bases):
#      exit 0 and the files of the calculation, with the old menu's answers on standard input (a command reads none);
#   2. refused with exit 2 and the three parts (what; "<command> takes:"; how to fix), before anything is
#      calculated: an unknown command, an option of another command, a missing required answer, a value out of
#      range (--delta-chi2=0, a negative --rate-entrance), a file that does not exist, a command option without a
#      command, a run without a command;
#   3. --help, --describe --json and docs/source/reference/command_line.rst name the same commands and options;
#      --help lists the commands first, prints each range once (from the table, not the help text) and marks the
#      rate's temperatures as one of two ways;
#   4. tools/menu_to_command.py turns an old script's menu answers into the command (and refuses the answers that
#      made no run or that depend on the project).
#
# A run by command wrote the same files as the run by menu with the same answers while the menu existed (the branch
# that removed the menu compared every project in the three bases and the eight fits with the menu's reference
# outputs).  Shown failing: the range of --walkers written into its help text again (printed twice); --delta-chi2
# with its range back to "at least 0"; the parameter file dropped from the preparation (ApplyAnswers).
#
#   ./tests/cli_commands/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C OMP_NUM_THREADS=1

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"
TESTS="$HERE/.."
DOC="$HERE/../../docs/source/reference/command_line.rst"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cli_commands.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
fail=0

# fixture NAME SRC_DIR FILE...: a fresh project folder $WORK/NAME with output/ and checks/
fixture() {
  local name="$1" src="$2"; shift 2
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"
  for f in "$@"; do cp -R "$src/$f" "$WORK/$name/"; done
}

# command NAME AZR WANT COMMAND...: the run by command, the old menu's answers on its input (it must not read them);
# exit 0 and every file of WANT (globs relative to the folder) written
command() {
  local name="$1" azr="$2" want="$3"; shift 3
  (cd "$WORK/$name" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" "$@" "$azr" > run.log 2>&1; echo $? > status)
  local st missing="" f
  st="$(cat "$WORK/$name/status")"
  for f in $want; do compgen -G "$WORK/$name/$f" > /dev/null || missing="$missing $f"; done
  if [ "$st" = 0 ] && [ -z "$missing" ] && ! grep -q "Please select\|azure2: " "$WORK/$name/run.log"; then
    ok "$name: $* ran, exit 0, wrote $want"
  else
    bad "$name: $*: exit $st, missing:$missing"
    tail -5 "$WORK/$name/run.log" | sed 's/^/        /'
  fi
}
total() { tr -d '\r' < "$WORK/$1/output/chiSquared.out" 2> /dev/null | awk '/^Total-Chi-Squared:/ { print $2 }'; }

# ---- 1. the commands run ---------------------------------------------------------------------------------------
fixture calc "$TESTS/upgrade/T1" T1.azr data
command calc T1.azr "output/chiSquared.out output/AZUREOut_*.out output/param.par" calculate
fixture extrapolate "$TESTS/upgrade/T1" T1.azr data
command extrapolate T1.azr "output/AZUREOut_*.extrap" extrapolate --covariance-band
fixture fit "$TESTS/identical_pp_res" identical_pp_res.azr data
command fit identical_pp_res.azr "output/param.sav output/covariance.dat" fit --covariance-band --scale-covariance
fixture minos "$TESTS/identical_pp_res" identical_pp_res.azr data
command minos identical_pp_res.azr "output/param.sav" minos --delta-chi2=1
grep -q "Performing parameter error analysis" "$WORK/minos/run.log" && ok "minos: the intervals were computed" ||
  bad "minos: no MINOS analysis in the log"
# the fit's parameters as the start of a calculation, relative to the project's directory: the fit's chi-squared
fixture start "$TESTS/identical_pp_res" identical_pp_res.azr data
cp "$WORK/fit/output/param.sav" "$WORK/start/fitted.sav" 2> /dev/null
command start identical_pp_res.azr "output/chiSquared.out" calculate --parameter-file=fitted.sav
if grep -q "Reading User Parameter File" "$WORK/start/run.log" &&
   awk -v a="$(total start)" -v b="$(total fit)" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(a != "" && d <= 1e-5 * b) }'; then  # the file has six digits
  ok "start: the parameter file was read: chi2 $(total start), the fit's $(total fit)"
else
  bad "start: chi2 '$(total start)' from --parameter-file, the fit's '$(total fit)'"
fi
fixture rate "$TESTS/reaction_rate/models" p12C.azr
cp "$TESTS/reaction_rate/temperatures.dat" "$WORK/rate/"
command rate p12C.azr "output/reactionrates.out" rate --rate-entrance=1 --rate-exit=2 \
  --temperature-file=temperatures.dat
fixture rate_grid "$TESTS/reaction_rate/models" p12C.azr
command rate_grid p12C.azr "output/reactionrates.out" rate --rate-entrance=1 --rate-exit=-1 \
  --temperature-min=0.01 --temperature-max=2 --temperature-step=0.1
help="$("$AZURE2_BIN" --help 2>&1 | tr -d '\r')"
if grep -q '^	mcmc:' <<< "$help"; then
  for b in physical rwa; do
    fixture "mcmc_$b" "$TESTS/mcmc_physical_norm" 13N_fast.azr
    cp -R "$TESTS/13N/data" "$WORK/mcmc_$b/"
    command "mcmc_$b" 13N_fast.azr "output/samples.mcmc" mcmc --walkers=4 --steps=100 --parameter-spread=0.5 \
      --sampler-threads=1 --basis=$b --overwrite-samples
  done
else
  echo "  (no sampler in this build: mcmc not checked)"
fi

# ---- 2. refusals -------------------------------------------------------------------------------------------------
mkdir -p "$WORK/refuse/output" "$WORK/refuse/checks"
cp -R "$TESTS/upgrade/T1/T1.azr" "$TESTS/upgrade/T1/data" "$WORK/refuse/"
# refused LABEL WANT ARGS...: exit 2, "ERROR: ..." containing WANT, the line of what the command takes (or the
# commands), a fix line, and no calculation
refused() {
  local label="$1" want="$2" out st; shift 2
  out="$(cd "$WORK/refuse" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" "$@" 2>&1 | tr -d '\r')"
  st=$(cd "$WORK/refuse" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" "$@" > /dev/null 2>&1; echo $?)
  if [ "$st" = 2 ] && grep -qF -- "ERROR: " <<< "$out" && grep -qF -- "$want" <<< "$out" &&
     grep -qE -- "( takes: |The commands)" <<< "$out" && [ "$(grep -c '^       ' <<< "$out")" -ge 2 ] &&
     ! grep -q "Calling AZURE2" <<< "$out"; then
    ok "$label: refused, exit 2, in three parts"
  else
    bad "$label: exit $st, not the refusal: $out"
  fi
}
refused "unknown command" "fitt is not a command" fitt T1.azr
refused "option of another command" "--walkers is not an option of the command fit" fit --walkers=4 T1.azr
refused "missing required answer" "needs --delta-chi2=X" minos T1.azr
refused "missing temperatures" "needs --temperature-min=T" rate --rate-entrance=1 --rate-exit=2 T1.azr
refused "value out of range" "--walkers takes a whole number from 2 to 1000" mcmc --walkers=1 --steps=100 \
  --sampler-threads=1 --basis=rwa T1.azr
refused "a file that does not exist" "--parameter-file=no_such.sav: cannot read" fit --parameter-file=no_such.sav \
  T1.azr
refused "elastic rate" "Cannot calculate rate for elastic scattering" rate --rate-entrance=1 --rate-exit=1 \
  --temperature-min=1 --temperature-max=2 --temperature-step=1 T1.azr
refused "a command option without a command" "--delta-chi2 is an option of a command" --delta-chi2=1 T1.azr
refused "no project file" "needs the project file" fit
refused "no command (the old menu's answers on the input)" "No command:" T1.azr
refused "--delta-chi2=0" "--delta-chi2 takes a number above 0" minos --delta-chi2=0 T1.azr
refused "a negative entrance pair" "--rate-entrance takes a whole number of at least 1" rate --rate-entrance=-1 \
  --rate-exit=2 --temperature-file=T1.azr T1.azr
[ -z "$(find "$WORK/refuse/output" -type f)" ] && ok "the refused command lines calculated nothing" ||
  bad "a refused command line wrote outputs"

# ---- 3. --help, --describe --json and the manual name the same commands and options ------------------------------
json="$(cd "$WORK/refuse" && "$AZURE2_BIN" --describe T1.azr --json 2> /dev/null | tr -d '\r')"
described="$(python3 -c '
import json, sys
d = json.loads(sys.stdin.read())
for c in d["commands"]:
    print("command", c["name"])
    for o in c["options"]:
        print("option", c["name"], o["option"])
' <<< "$json" 2>&1)"
helped="$(awk '/^Commands/ { on = 1 } on && /^\t[a-z]+:/ { sub(/^\t/, ""); sub(/:.*/, ""); c = $0; print "command", c }
  on && /^ +Options: / { sub(/^ +Options: /, ""); sub(/\.$/, ""); n = split($0, o, ", ");
  for (i = 1; i <= n; i++) { sub(/=.*/, "", o[i]); print "option", c, o[i] } }' <<< "$help")"
if [ -n "$described" ] && [ "$described" = "$helped" ]; then
  ok "--help and --describe --json list the same $(grep -c '^command' <<< "$helped") commands and their options"
else
  bad "--help and --describe --json differ: $(diff <(echo "$helped") <(echo "$described"))"
fi
doc="$(tr -d '\r' < "$DOC")"
missing=""
for c in $(awk '$1 == "command" { print $2 }' <<< "$helped"); do
  grep -qF -- "\`\`$c\`\`" <<< "$doc" || missing="$missing $c"
done
for o in $(awk '$1 == "option" { print $3 }' <<< "$helped" | sort -u); do
  grep -qF -- "\`\`$o" <<< "$doc" || missing="$missing $o"
done
[ -z "$missing" ] && ok "command_line.rst names every command and option" ||
  bad "command_line.rst does not name:$missing"

# the help: the commands first, each range once, the temperatures as one of two ways
[ "$(grep -m1 -n '^Commands' <<< "$help" | cut -d: -f1)" -lt "$(grep -m1 -n -- '--use-brune' <<< "$help" | cut -d: -f1)" ] &&
  ok "--help lists the commands before the program's options" || bad "--help does not list the commands first"
twice=""
while IFS= read -r line; do
  # the row's help text (before the table's "; <range>; required|default"): no range of its own
  text="${line#*:}"; text="${text%%; *}"
  grep -qE '(from|to|above|least|between) -?[0-9]' <<< "$text" && twice="$twice ${line%%:*}"
done < <(awk '/^The options of the commands/ { on = 1; next } on && /^\t--/ { print } on && !/^\t/ { on = 0 }' <<< "$help")
[ -z "$twice" ] && ok "--help prints each option's range once (from the table)" || bad "--help: a range in the help text of:$twice"
tf="$(grep -m1 '^	--temperature-file:' <<< "$help")"
grep -q "required, or instead --temperature-min, --temperature-max and --temperature-step" <<< "$tf" &&
  ok "--help: --temperature-file is one of two ways" || bad "--help: --temperature-file: $tf"

# ---- 4. tools/menu_to_command.py: an old script's menu answers -> the command ------------------------------------
TOOL="$HERE/../../tools/menu_to_command.py"
# converts ANSWERS WANT [TOOL OPTIONS]: the tool prints WANT, exit 0
converts() {
  local answers="$1" want="$2" got; shift 2
  got="$(python3 -I "$TOOL" "$@" "$answers" 2> /dev/null)"
  [ "$got" = "$want" ] && ok "menu_to_command '$answers' -> $want" || bad "menu_to_command '$answers': '$got', not '$want'"
}
converts '1\n\n\n7\n' calculate
converts '3\nn\nparam.sav\n\n7\n' "extrapolate --parameter-file=param.sav"
converts '2\ny\ny\nfit.sav\n\n7\n' "fit --covariance-band --scale-covariance --parameter-file=fit.sav"
converts '2\nfit.sav\n\n7\n' "fit --parameter-file=fit.sav" --flags=--covariance-band
converts '4\n1\nn\n\n\n7\n' "minos --delta-chi2=1"
converts '5\n\n1\n-1\nno\n0.01\n2\n0.01\n' \
  "rate --rate-entrance=1 --rate-exit=-1 --temperature-min=0.01 --temperature-max=2 --temperature-step=0.01"
converts '6\n\n24\n100\n\n0.5\n1\nno\nyes\n' \
  "mcmc --walkers=24 --steps=100 --energy-spread=0.5 --sampler-threads=1 --basis=physical --overwrite-samples"
converts '1\n\nintEC.dat\n7\n' "calculate --capture-amplitudes=intEC.dat" --capture
# refuses ANSWERS WHAT: exit 1 and the reason on standard error
refuses() {
  local st err; err="$(python3 -I "$TOOL" "$1" 2>&1 > /dev/null)"; st=$?
  [ "$st" = 1 ] && grep -q -- "$2" <<< "$err" && ok "menu_to_command '$1': refused ($2)" ||
    bad "menu_to_command '$1': exit $st, $err"
}
refuses '3\nparam.sav\n\n7\n' "scripted-answer trap"
refuses '1\n\nintEC.dat\n7\n' "with and without external capture"
refuses '7\n' "Exit"

[ "$fail" -eq 0 ] && echo "  PASS" || echo "  FAIL"
[ "$fail" -eq 0 ]
