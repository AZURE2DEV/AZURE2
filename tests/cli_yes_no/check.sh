#!/usr/bin/env bash
# The y/n questions of the text menu (2026-11-01).
#
# Modes 2 and 3 ask "Calculate cross-section uncertainty band? (y/n)" BEFORE the parameter file.  A script that left
# that line out ('3\nparam.sav\n\n7') gave the parameter-file name to the y/n question, and before 2026-11-01 any answer
# not starting with y, 1 or t was "no": the blank line after it then reached the parameter-file prompt ("leave blank
# for new file"), and the run used the <levels> of the .azr instead of the file, without an error.  Found on a 7Be A_y
# test-segment calculation whose curves disagreed with the data-segment ones, and in the 12C+a "run_again" fits of
# July-August 2026 (meant to continue from output/param.sav, they started from <levels>).  Now the answers are
# y, yes, 1, t, true and n, no, 0, f, false (any case) and a blank line ("no"); job scripts answer 1.  Anything else
# stops a piped run with an ERROR naming the answer and exit status 1, before anything is calculated.
#
# Checked in mode 3 on tests/7Li_p_ay (one A_y test segment) with a parameter file written by a mode-1 run.
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
# A parameter file to point at: the param.par a mode-1 run writes.
(cd "$WORK/base" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline "$PROJ.azr" > mode1.log 2>&1)
[ -s "$WORK/base/output/param.par" ] && cp "$WORK/base/output/param.par" "$WORK/fit.par" \
  || bad "mode 1 wrote no param.par (see the log)"

# run NAME INPUT: mode-3 input INPUT (printf format) in a fresh copy; sets $rc.
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$WORK/base/data" "$WORK/base/$PROJ.azr" "$dir/"
  cp "$WORK/fit.par" "$dir/fit.par" 2>/dev/null
  (cd "$dir" && printf "$2" | "$AZURE2_BIN" --no-gui --no-readline "$PROJ.azr" > run.log 2>&1)
  rc=$?
}
said() { grep -q "$2" "$WORK/$1/run.log"; }
extrap() { ls "$WORK/$1"/output/*.extrap 2>/dev/null | head -1; }

run missed '3\nfit.par\n\n\n7\n'
said missed 'ERROR: "Calculate cross-section uncertainty band?" takes y or n, not "fit.par"' && [ "$rc" = 1 ] \
  && ok "parameter-file name given to the y/n question: refused, answer named, exit 1" \
  || bad "parameter-file name given to the y/n question: exit $rc, ERROR $(said missed 'takes y or n' && echo printed || echo missing)"
! said missed "Creating New param.par File" && [ -z "$(extrap missed)" ] \
  && ok "  ... nothing calculated (no new param.par, no .extrap)" \
  || bad "  ... the run went on: $(said missed 'Creating New param.par File' && echo 'new param.par from <levels>,') $(extrap missed)"
said missed "BEFORE the parameter file" && ok "  ... the message says the band question comes first" \
  || bad "  ... no hint about the prompt order"

run right '3\nn\nfit.par\n\n7\n'
said right "Reading User Parameter File" && [ "$rc" = 0 ] && [ -n "$(extrap right)" ] \
  && ok "'3 n fit.par': the parameter file is read, .extrap written, exit 0" \
  || bad "'3 n fit.par': exit $rc, parameter file $(said right 'Reading User Parameter File' && echo read || echo 'not read')"

run upper '3\n No \nfit.par\n\n7\n'
said upper "Reading User Parameter File" && [ "$rc" = 0 ] \
  && ok "' No ' (case, spaces) is an answer" || bad "' No ': exit $rc"

run blank '3\n\nfit.par\n\n7\n'
said blank "Reading User Parameter File" && [ "$rc" = 0 ] \
  && ok "a blank line is 'no'" || bad "blank answer: exit $rc"

run zero '3\n0\nfit.par\n\n7\n'
said zero "Reading User Parameter File" && [ "$rc" = 0 ] \
  && ok "'0' is 'no' (job scripts answer 1 / 0)" || bad "'0': exit $rc"

run one_name '3\n1x\nfit.par\n\n7\n'
said one_name 'takes y or n, not "1x"' && [ "$rc" = 1 ] \
  && ok "'1x' is refused (only whole answers count)" || bad "'1x': exit $rc"

# Byte-identical files without cmp (not in every MSYS2 image: the Windows runner had none, 2026-11-05);
# \r stripped so a CRLF-writing Windows build compares like the others.
same_file() { [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
[ -n "$(extrap right)" ] && same_file "$(extrap right)" "$(extrap upper)" && same_file "$(extrap right)" "$(extrap blank)" \
  && ok "  ... 'n', ' No ' and blank give the same .extrap" || bad "  ... the .extrap files differ"

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
