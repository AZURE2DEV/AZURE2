#!/usr/bin/env bash
#
# AZURE2 --help lists exactly the flags the command line accepts.
#
# The help text (src/AZURE2.cpp) and the parser (parseOptions, src/ProgramMain.cpp) are
# two lists kept by hand; they had drifted: --help printed
# "--no-long-wavelenth" (the parser takes --no-long-wavelength) and left out
# --covariance-band, --scale-covariance and --use-api.  Checked here, on the
# binary as built (a flag behind a build option, --use-nlopt, counts only
# where it is compiled in):
#
#   1. every flag the help prints is accepted (no "Unknown option" error);
#   2. every flag literal of src/AZURE2.cpp and src/ProjectCommands.cpp (the --convert, --validate and --schema
#      commands), of the option table and of the switch table (support/src/Environment.cpp) is either accepted and
#      printed,
#      or neither (compiled out);
#   3. the flag table of docs/source/reference/command_line.rst names every
#      printed flag, and no flag the source does not have.
#
# Each probe runs the binary with the flag and a project file that does not
# exist, so nothing is calculated; the whole check takes about a second.
#
#   ./tests/cli_help/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-60}")"
export OMP_NUM_THREADS=1
SRC="$HERE/../../src/AZURE2.cpp"
SRC_COMMANDS="$HERE/../../src/ProjectCommands.cpp"
SRC_PROGRAM="$HERE/../../src/ProgramMain.cpp"  # parseOptions: the run of a command (main's and the front ends' facade)
SRC_OPTIONS="$HERE/../../io/src/ConfigSection.cpp"  # the flags of the options of the calculation (their table)
SRC_SWITCHES="$HERE/../../support/src/Environment.cpp"  # the options of the switches that move numbers (their table)
SRC_FITTERS="$HERE/../../src/Fitters.cpp"  # the fitters' switches (the fitter table, with the fitters not built)
SRC_RUN="$HERE/../../src/RunRequest.cpp"  # the options of the commands (the table of a run's questions)
DOC="$HERE/../../docs/source/reference/command_line.rst"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cli_help.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
fail=0

if [ ! -f "$SRC" ] || [ ! -f "$SRC_COMMANDS" ] || [ ! -f "$DOC" ]; then
  bad "src/AZURE2.cpp or docs/source/reference/command_line.rst not found next to tests/"
  echo "  FAIL"; exit 1
fi

# The flags as printed: "<tab>--flag:" at the start of a line.
(cd "$WORK" && $RUN "$AZURE2_BIN" --help < /dev/null 2>&1 | head -c 100000 | tr -d '\r' > help.txt)
# A flag that takes a value prints as "--flag[=FILE]:".
printed="$(sed -n 's/^[[:space:]]*\(--[a-z][a-z0-9-]*\)\(\[=[A-Z]*\]\)\{0,1\}:.*/\1/p' "$WORK/help.txt" | sort -u)"
# The flags of the source: every "--flag" string literal of src/AZURE2.cpp, src/ProjectCommands.cpp and the option table.
literals="$(cat "$SRC" "$SRC_COMMANDS" "$SRC_PROGRAM" "$SRC_OPTIONS" "$SRC_SWITCHES" "$SRC_FITTERS" "$SRC_RUN" | tr -d '\r' |
  grep -o '"--[a-z][a-z0-9-]*"' | tr -d '"' | sort -u)"
# The flags of the documentation's table: "* - ``--flag``" rows.
# and the switch options without a value: the "Option" column of the list "Environment switches".
documented="$( (tr -d '\r' < "$DOC" | sed -n 's/^ *\* - ``\(--[a-z][a-z0-9-]*\)\(\[=[A-Z]*\]\)\{0,1\}``.*/\1/p';
  tr -d '\r' < "$DOC" | sed -n '/^Environment switches$/,$s/^     - ``\(--[a-z][a-z0-9-]*\)``$/\1/p') | sort -u)"

if [ -z "$printed" ] || [ -z "$literals" ] || [ -z "$documented" ]; then
  bad "could not read the flags (help: $(echo $printed | wc -w), source: $(echo $literals | wc -w), docs: $(echo $documented | wc -w))"
  echo "  FAIL"; exit 1
fi

# The lists are searched from here-strings, never through a pipe: under pipefail a quiet grep that leaves at its
# first match can end the writer with SIGPIPE, and the search then counts as failed (seen as random failures of
# this test on a loaded machine).
# accepted FLAG: the binary takes FLAG without the parser's error.  A
# missing project file ends the run right after the options are parsed.
accepted() {
  local out
  [ "$1" = --help ] && { grep -q "^Syntax: AZURE2" "$WORK/help.txt"; return; }
  out="$(cd "$WORK" && $RUN "$AZURE2_BIN" --no-gui "$1" no_such_project.azr < /dev/null 2>&1 | head -c 100000 | tr -d '\r')"
  ! grep -qF "Unknown option $1." <<< "$out" && ! grep -qF "which this build of AZURE2 does not have" <<< "$out"
}

# A fitter's switch whose fitter this build leaves out (the absent rows of src/Fitters.cpp) is refused in three
# parts: what (the fitter is not in this build), the switches allowed, how to build it; exit 2.
absent="$(sed -n '/AbsentFitters()/,/^}/s/.*{"[a-z-]*", "\(--[a-z-]*\)".*/\1/p' "$SRC_FITTERS")"
if grep -qF "#ifndef USE_NLOPT" "$SRC_FITTERS" && [ -z "$absent" ]; then
  bad "could not read the absent fitters' switches from src/Fitters.cpp (check this test)"
fi
for f in $absent; do
  out="$(cd "$WORK" && $RUN "$AZURE2_BIN" --no-gui "$f" no_such_project.azr < /dev/null 2>&1 | head -c 100000 |
    tr -d '\r')"
  status=$(cd "$WORK" && $RUN "$AZURE2_BIN" --no-gui "$f" no_such_project.azr < /dev/null > /dev/null 2>&1; echo $?)
  grep -qF -- "$f" "$WORK/help.txt" && continue  # built in: accepted and printed (checked below)
  if [ "$status" = 2 ] && grep -qF "this build of AZURE2 does not have" <<< "$out" &&
     grep -qF "The fitter switches of this build: --use-" <<< "$out" && grep -qF "build AZURE2 with" <<< "$out"; then
    ok "$f: a fitter not built, refused (exit 2) with what, the switches allowed and how to build it"
  else
    bad "$f: a fitter not built: exit $status, not the three-part refusal: $out"
  fi
done

for f in $(printf '%s\n%s\n' "$printed" "$literals" | sort -u); do
  in_help=no; in_src=no
  grep -qxF -- "$f" <<< "$printed" && in_help=yes
  grep -qxF -- "$f" <<< "$literals" && in_src=yes
  if accepted "$f"; then acc=yes; else acc=no; fi
  if [ "$in_help" = yes ] && [ "$acc" = no ]; then
    bad "$f is printed by --help but the parser refuses it"
  elif [ "$acc" = yes ] && [ "$in_help" = no ]; then
    bad "$f is accepted but --help does not print it"
  elif [ "$in_src" = no ]; then
    bad "$f is printed and accepted but is no literal of src/AZURE2.cpp, src/ProgramMain.cpp or src/ProjectCommands.cpp (check this test)"
  elif [ "$acc" = yes ]; then
    ok "$f: accepted and printed"
  else
    ok "$f: compiled out of this build, and not printed"
  fi
done

# The documentation's table: every printed flag, and nothing the source lacks.
for f in $printed; do
  grep -qxF -- "$f" <<< "$documented" || bad "$f is printed but not in the table of command_line.rst"
done
for f in $documented; do
  grep -qxF -- "$f" <<< "$literals" || bad "command_line.rst documents $f, which src/AZURE2.cpp and src/ProgramMain.cpp do not have"
done
[ "$fail" -eq 0 ] && ok "command_line.rst documents the $(echo $printed | wc -w) printed flags and no other"

[ "$fail" -eq 0 ] && echo "  PASS" || echo "  FAIL"
[ "$fail" -eq 0 ]
