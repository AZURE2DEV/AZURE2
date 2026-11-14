#!/usr/bin/env bash
#
# AZURE2 <command> [options] FILE runs what the menu runs with the same answers, and refuses a wrong command line.
#
# The command line and the menu fill one run request from one table (src/RunRequest.cpp); both runs go through the
# same preparation (src/RunPreparation.cpp).  Checked on the binary as built:
#
#   1. calculate (T1), extrapolate (T1), fit and minos (identical_pp_res), rate (the closed 12C(p,g) model, both
#      temperature forms) and, in a build with the sampler, mcmc (13N_fast, both bases): every file the run by
#      command writes is identical to the run by menu with the same answers, and the exit statuses are equal;
#   2. a command reads nothing from standard input (a run with a menu script on its input still runs the command);
#   3. refused with exit 2 and the three parts (what; "<command> takes:"; how to fix), before anything is
#      calculated: an unknown command, an option of another command, a missing required answer, a value out of
#      range, a file that does not exist, a command option without a command;
#   4. --help, --describe --json and docs/source/reference/command_line.rst name the same commands and options.
#
# First passed on the build that added the commands.  Shown failing: with the command line's defaults left unset
# (ReadRunCommand), the mcmc runs differ from the menu's; with the parameter file dropped from the shared preparation
# (ApplyAnswers), the check that the log read it fails (the two runs stay equal there, as both take that path).
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

# fixture NAME SRC_DIR FILE...: a fresh project folder $WORK/NAME/{menu,command} with output/ and checks/
fixture() {
  local name="$1" src="$2"; shift 2
  for r in menu command; do
    mkdir -p "$WORK/$name/$r/output" "$WORK/$name/$r/checks"
    for f in "$@"; do cp -R "$src/$f" "$WORK/$name/$r/"; done
  done
}

# pair NAME AZR MENU COMMAND...: the run by menu (the answers MENU on standard input) and by command (the menu's
# answers on its input too: a command must not read them); the same files and exit status
pair() {
  local name="$1" azr="$2" menu="$3"; shift 3
  (cd "$WORK/$name/menu" && printf "$menu" | $RUN "$AZURE2_BIN" --no-readline "$azr" > ../menu.log 2>&1
   echo $? > ../menu.status)
  (cd "$WORK/$name/command" && printf "$menu" | $RUN "$AZURE2_BIN" "$@" "$azr" > ../command.log 2>&1
   echo $? > ../command.status)
  local files n=0 differ="" sm sc
  files="$(cd "$WORK/$name/menu" && find . -type f | sort)"
  for f in $files; do
    n=$((n + 1))
    cmp -s "$WORK/$name/menu/$f" "$WORK/$name/command/$f" || differ="$differ $f"
  done
  [ "$files" = "$(cd "$WORK/$name/command" && find . -type f | sort)" ] || differ="$differ (the lists of files)"
  sm="$(cat "$WORK/$name/menu.status")"; sc="$(cat "$WORK/$name/command.status")"
  if [ -z "$differ" ] && [ "$sm" = "$sc" ] && [ "$(find "$WORK/$name/menu/output" -type f | wc -l)" -gt 0 ]; then
    ok "$name: $* gives the menu's $n files, exit $sc"
  else
    bad "$name: $*: differs from the menu (exit $sm / $sc):$differ"
    tail -5 "$WORK/$name/command.log" | sed 's/^/        /'
  fi
}

# ---- 1. the run by command is the run by menu ----------------------------------------------------------------------
fixture calc "$TESTS/upgrade/T1" T1.azr data
pair calc T1.azr '1\n\n\n7\n' calculate
fixture extrapolate "$TESTS/upgrade/T1" T1.azr data
pair extrapolate T1.azr '3\ny\n\n\n7\n' extrapolate --covariance-band
fixture fit "$TESTS/identical_pp_res" identical_pp_res.azr data
pair fit identical_pp_res.azr '2\ny\ny\n\n\n7\n' fit --covariance-band --scale-covariance
fixture minos "$TESTS/identical_pp_res" identical_pp_res.azr data
pair minos identical_pp_res.azr '4\n1\nn\n\n\n7\n' minos --delta-chi2=1
# the fit's parameters as the start of a calculation, relative to the project's directory
fixture start "$TESTS/identical_pp_res" identical_pp_res.azr data
for r in menu command; do cp "$WORK/fit/menu/output/param.sav" "$WORK/start/$r/fitted.sav" 2> /dev/null; done
pair start identical_pp_res.azr '1\nfitted.sav\n\n7\n' calculate --parameter-file=fitted.sav
# both runs share the preparation, so a parameter file dropped there would leave them equal: the log says it was read
grep -q "Reading User Parameter File" "$WORK/start/command.log" && ok "start: the parameter file was read" ||
  bad "start: the run did not read --parameter-file"
fixture rate "$TESTS/reaction_rate/models" p12C.azr
for r in menu command; do cp "$TESTS/reaction_rate/temperatures.dat" "$WORK/rate/$r/"; done
pair rate p12C.azr '5\n\n1\n2\nyes\ntemperatures.dat\n' rate --rate-entrance=1 --rate-exit=2 \
  --temperature-file=temperatures.dat
fixture rate_grid "$TESTS/reaction_rate/models" p12C.azr
pair rate_grid p12C.azr '5\n\n1\n-1\nno\n0.01\n2\n0.1\n' rate --rate-entrance=1 --rate-exit=-1 \
  --temperature-min=0.01 --temperature-max=2 --temperature-step=0.1
help="$("$AZURE2_BIN" --help 2>&1 | tr -d '\r')"
if grep -q '^	mcmc:' <<< "$help"; then
  for b in physical rwa; do
    fixture "mcmc_$b" "$TESTS/mcmc_physical_norm" 13N_fast.azr
    for r in menu command; do cp -R "$TESTS/13N/data" "$WORK/mcmc_$b/$r/"; done
  done
  pair mcmc_physical 13N_fast.azr '6\n\n\n4\n100\n0.5\n\n1\nno\nyes\n' mcmc --walkers=4 --steps=100 \
    --parameter-spread=0.5 --sampler-threads=1 --basis=physical --overwrite-samples
  pair mcmc_rwa 13N_fast.azr '6\n\n\n4\n100\n\n0.5\n1\nyes\nyes\n' mcmc --walkers=4 --steps=100 \
    --energy-spread=0.5 --sampler-threads=1 --basis=rwa --overwrite-samples
else
  echo "  (no sampler in this build: mcmc not checked)"
fi

# ---- 2. a command reads nothing from standard input ----------------------------------------------------------------
if grep -q "Please select" "$WORK/calc/command.log" || grep -q "azure2: " "$WORK/calc/command.log"; then
  bad "a command printed the menu"
else
  ok "a command prints no menu and leaves standard input alone (the menu's script on its input changed nothing)"
fi

# ---- 3. refusals -------------------------------------------------------------------------------------------------
mkdir -p "$WORK/refuse/output" "$WORK/refuse/checks"
cp -R "$TESTS/upgrade/T1/T1.azr" "$TESTS/upgrade/T1/data" "$WORK/refuse/"
# refused LABEL WANT ARGS...: exit 2, "ERROR: ..." containing WANT, the line of what the command takes (or the
# commands), a fix line, and no calculation
refused() {
  local label="$1" want="$2" out st; shift 2
  out="$(cd "$WORK/refuse" && $RUN "$AZURE2_BIN" calculate "$@" 2>&1 | tr -d '\r')"
  st=$(cd "$WORK/refuse" && $RUN "$AZURE2_BIN" calculate "$@" > /dev/null 2>&1; echo $?)
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
[ -z "$(find "$WORK/refuse/output" -type f)" ] && ok "the refused command lines calculated nothing" ||
  bad "a refused command line wrote outputs"

# ---- 4. --help, --describe --json and the manual name the same commands and options ------------------------------
json="$(cd "$WORK/refuse" && "$AZURE2_BIN" --describe T1.azr --json 2> /dev/null | tr -d '\r')"
described="$(python3 -c '
import json, sys
d = json.loads(sys.stdin.read())
for c in d["commands"]:
    print("command", c["name"])
    for o in c["options"]:
        print("option", c["name"], o["option"])
' <<< "$json" 2>&1)"
helped="$(awk '/^Commands:/ { on = 1 } on && /^\t[a-z]+:/ { sub(/^\t/, ""); sub(/:.*/, ""); c = $0; print "command", c }
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

[ "$fail" -eq 0 ] && echo "  PASS" || echo "  FAIL"
[ "$fail" -eq 0 ]
