#!/usr/bin/env bash
#
# AZURE2 --help lists exactly the flags the command line accepts.
#
# The help text and the parser (parseOptions and main in src/AZURE2.cpp) are
# two lists kept by hand; they had drifted: --help printed
# "--no-long-wavelenth" (the parser takes --no-long-wavelength) and left out
# --covariance-band, --scale-covariance and --use-api.  Checked here, on the
# binary as built (flags behind a build option, --no-readline and
# --use-nlopt, count only where they are compiled in):
#
#   1. every flag the help prints is accepted (no "Unknown option" warning);
#   2. every flag literal of src/AZURE2.cpp is either accepted and printed,
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
DOC="$HERE/../../docs/source/reference/command_line.rst"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cli_help.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
fail=0

if [ ! -f "$SRC" ] || [ ! -f "$DOC" ]; then
  bad "src/AZURE2.cpp or docs/source/reference/command_line.rst not found next to tests/"
  echo "  FAIL"; exit 1
fi

# The flags as printed: "<tab>--flag:" at the start of a line.
(cd "$WORK" && $RUN "$AZURE2_BIN" --help < /dev/null 2>&1 | head -c 100000 | tr -d '\r' > help.txt)
printed="$(sed -n 's/^[[:space:]]*\(--[a-z][a-z-]*\):.*/\1/p' "$WORK/help.txt" | sort -u)"
# The flags of the source: every "--flag" string literal of src/AZURE2.cpp.
literals="$(tr -d '\r' < "$SRC" | grep -o '"--[a-z][a-z-]*"' | tr -d '"' | sort -u)"
# The flags of the documentation's table: "* - ``--flag``" rows.
documented="$(tr -d '\r' < "$DOC" | sed -n 's/^ *\* - ``\(--[a-z][a-z-]*\)``.*/\1/p' | sort -u)"

if [ -z "$printed" ] || [ -z "$literals" ] || [ -z "$documented" ]; then
  bad "could not read the flags (help: $(echo $printed | wc -w), source: $(echo $literals | wc -w), docs: $(echo $documented | wc -w))"
  echo "  FAIL"; exit 1
fi

# accepted FLAG: the binary takes FLAG without the parser's warning.  A
# missing project file ends the run right after the options are parsed.
accepted() {
  local out
  [ "$1" = --help ] && { grep -q "^Syntax: AZURE2" "$WORK/help.txt"; return; }
  out="$(cd "$WORK" && $RUN "$AZURE2_BIN" --no-gui "$1" no_such_project.azr < /dev/null 2>&1 | head -c 100000 | tr -d '\r')"
  ! printf '%s\n' "$out" | grep -qF "Unknown option $1."
}

for f in $(printf '%s\n%s\n' "$printed" "$literals" | sort -u); do
  in_help=no; in_src=no
  printf '%s\n' "$printed" | grep -qxF -- "$f" && in_help=yes
  printf '%s\n' "$literals" | grep -qxF -- "$f" && in_src=yes
  if accepted "$f"; then acc=yes; else acc=no; fi
  if [ "$in_help" = yes ] && [ "$acc" = no ]; then
    bad "$f is printed by --help but the parser refuses it"
  elif [ "$acc" = yes ] && [ "$in_help" = no ]; then
    bad "$f is accepted but --help does not print it"
  elif [ "$in_src" = no ]; then
    bad "$f is printed and accepted but is no literal of src/AZURE2.cpp (check this test)"
  elif [ "$acc" = yes ]; then
    ok "$f: accepted and printed"
  else
    ok "$f: compiled out of this build, and not printed"
  fi
done

# The documentation's table: every printed flag, and nothing the source lacks.
for f in $printed; do
  printf '%s\n' "$documented" | grep -qxF -- "$f" || bad "$f is printed but not in the table of command_line.rst"
done
for f in $documented; do
  printf '%s\n' "$literals" | grep -qxF -- "$f" || bad "command_line.rst documents $f, which src/AZURE2.cpp does not have"
done
[ "$fail" -eq 0 ] && ok "command_line.rst documents the $(echo $printed | wc -w) printed flags and no other"

[ "$fail" -eq 0 ] && echo "  PASS" || echo "  FAIL"
[ "$fail" -eq 0 ]
