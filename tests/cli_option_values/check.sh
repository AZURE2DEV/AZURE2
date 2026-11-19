#!/usr/bin/env bash
#
# A value-taking flag given without its value, or with a value that does not parse, is an error: "ERROR: ..." naming
# the flag, the usage line, exit status 2, and nothing is run (include/CommandLineValues.h, parseOptions in
# src/ProgramMain.cpp).  Before, --max-calls N (without "=") was an unknown option with a warning and the fit ran without a
# limit, --max-calls=x a warning, --refreeze-covariance=x read as 0, --coulomb-table= an empty file name, and an empty
# AZURE_REFINE_QUADRATURE no refinement.  Each case below, and that the good values are still accepted (the run then
# stops only at the missing project file, exit 255).  The same AZURE_REFINE_QUADRATURE in pyazr raises ValueError
# (when the compiled module can be imported, PYAZR_ENGINE_DIR; else that case is skipped).  An unknown option is the same
# kind of error (before, "WARNING: Unknown option" and the run went on: a mistyped --use-gradiant fitted without the
# gradient), on the command line and in AZURE_OPTIONS_FILE; --no-readline, the option of the old menu's scripts, is
# refused in three parts (the program has no menu: drop it and give the command).
#
#   ./tests/cli_option_values/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS=1
. "$HERE/../lib/guard.sh"
RUN="$(guard_command 60)"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cli_option_values.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# probe ENV ARGS...: run the binary (calculate) on a project that does not exist; sets $out and $rc
probe() {
  local env="$1"; shift
  (cd "$WORK" && env $env $RUN "$AZURE2_BIN" calculate --no-gui "$@" no_such_project.azr < /dev/null > out.txt 2>&1)
  rc=$?
  out="$(head -c 20000 "$WORK/out.txt")"
}
# refused NAME TEXT ENV ARGS...: exit 2, an ERROR line with TEXT, the usage line, and the project never opened
refused() {
  local name="$1" text="$2"; shift 2
  probe "$@"
  if [ "$rc" = 2 ] && grep -q "^ERROR: .*$text" <<< "$out" && grep -q "Syntax: AZURE2" <<< "$out" &&
     ! grep -q "Could not open" <<< "$out"; then
    ok "$name: refused, exit 2 ($(grep -m1 '^ERROR' <<< "$out"))"
  else
    bad "$name: exit $rc, $(printf '%s\n' "$out" | sed -n '1,3p' | tr '\n' ' ')"
  fi
}
# accepted NAME ENV ARGS...: no ERROR, the run goes on to the missing project (exit 255; 127 in the MSYS2 shell of the
# Windows CI, where the same exit -1 of a Windows program reads so: seen there, 2026-11-15)
accepted() {
  local name="$1"; shift
  probe "$@"
  if { [ "$rc" = 255 ] || [ "$rc" = 127 ]; } && ! grep -q "^ERROR" <<< "$out" && grep -q "Could not open" <<< "$out"; then
    ok "$name: accepted"
  else
    bad "$name: exit $rc, $(printf '%s\n' "$out" | sed -n '1,3p' | tr '\n' ' ')"
  fi
}

refused "--max-calls without a value" "--max-calls needs a value: --max-calls=N" "" --max-calls
refused "--max-calls 50 (a space, not =)" "--max-calls needs a value: --max-calls=50" "" --max-calls 50
refused "--max-calls=x" "--max-calls needs a whole number from 0 to 1000000000, not 'x'" "" --max-calls=x
refused "--max-calls=" "--max-calls needs a whole number" "" --max-calls=
refused "--max-calls=-1" "--max-calls needs a whole number" "" --max-calls=-1
refused "--max-calls=10x" "--max-calls needs a whole number" "" --max-calls=10x
accepted "--max-calls=50" "" --max-calls=50
refused "--refreeze-covariance=x" "--refreeze-covariance needs a whole number from 0 to 1000" "" --refreeze-covariance=x
refused "--refreeze-covariance=-2" "--refreeze-covariance needs a whole number" "" --refreeze-covariance=-2
accepted "--refreeze-covariance (default 5)" "" --refreeze-covariance
accepted "--refreeze-covariance=3" "" --refreeze-covariance=3
refused "--coulomb-table= (no file)" "--coulomb-table= needs a value" "" --coulomb-table=
accepted "--coulomb-table=FILE" "" --coulomb-table=some_table.bin
accepted "--coulomb-table (the default file)" "" --coulomb-table
# (the default table file in the scratch folder: an engine that took the empty FILE as the default would build there)
refused "--build-coulomb-table= (no file)" "--build-coulomb-table= needs a FILE" \
  "HOME=$WORK AZURE2_COULOMB_TABLE=$WORK/table.bin" --build-coulomb-table=
refused "--use-gradiant (unknown)" "Unknown option --use-gradiant\." "" --use-gradiant
refused "--use-gradient=1 (unknown)" "Unknown option --use-gradient=1\." "" --use-gradient=1
refused "an unknown option after good ones" "Unknown option --no-such-flag\." "" --use-brune --no-such-flag --use-lm
accepted "--use-gradient --use-brune" "" --use-gradient --use-brune
refused "--no-readline (the menu is gone)" "Unknown option --no-readline\.  AZURE2 has no menu" "" --no-readline
grep -q "Drop --no-readline and give the command" <<< "$out" && ok "  ... the refusal says how to fix it" ||
  bad "  ... no fix in the refusal"
printf -- '--no-readline\n' > "$WORK/readline.txt"
refused "--no-readline in AZURE_OPTIONS_FILE" "Unknown option --no-readline" "AZURE_OPTIONS_FILE=$WORK/readline.txt"
printf -- '--use-gradiant\n' > "$WORK/unknown.txt"
refused "an unknown option in AZURE_OPTIONS_FILE" "Unknown option --use-gradiant\." "AZURE_OPTIONS_FILE=$WORK/unknown.txt"
printf -- '--max-calls\n' > "$WORK/options.txt"
refused "--max-calls without a value in AZURE_OPTIONS_FILE" "--max-calls needs a value" "AZURE_OPTIONS_FILE=$WORK/options.txt"
refused "AZURE_REFINE_QUADRATURE=x" "AZURE_REFINE_QUADRATURE needs a whole number from 1 to 64, not 'x'" "AZURE_REFINE_QUADRATURE=x"
refused "AZURE_REFINE_QUADRATURE= (empty)" "AZURE_REFINE_QUADRATURE needs a whole number" "AZURE_REFINE_QUADRATURE="
refused "AZURE_REFINE_QUADRATURE=0" "AZURE_REFINE_QUADRATURE needs a whole number" "AZURE_REFINE_QUADRATURE=0"
refused "AZURE_REFINE_QUADRATURE=65" "AZURE_REFINE_QUADRATURE needs a whole number" "AZURE_REFINE_QUADRATURE=65"
accepted "AZURE_REFINE_QUADRATURE=2" "AZURE_REFINE_QUADRATURE=2"

# pyazr: a session refuses a bad AZURE_REFINE_QUADRATURE with ValueError (the process goes on); a good one opens
py="$(cd "$WORK" && AZURE_REFINE_QUADRATURE=x python3 -I - "$HERE/../.." "$HERE/../thin_layer_limit/models/thin_limit.azr" <<'PY' 2>&1
import os, sys
sys.path.insert(0, sys.argv[1])
try:
    from pyazr import Session
except Exception as err:
    print("skip", type(err).__name__, err)
    sys.exit(0)
try:
    Session(sys.argv[2], cwd=os.getcwd())
    print("opened")
except ValueError as err:
    print("ValueError", err)
PY
)"
case "$py" in
  skip*) echo "  skip  pyazr not importable: ${py#skip }" ;;
  *"ValueError AZURE_REFINE_QUADRATURE needs a whole number"*) ok "pyazr: AZURE_REFINE_QUADRATURE=x raises ValueError" ;;
  *) bad "pyazr: AZURE_REFINE_QUADRATURE=x: $(printf '%s\n' "$py" | tail -2 | tr '\n' ' ')" ;;
esac

if [ "$failures" -ne 0 ]; then echo "  FAIL ($failures)"; exit 1; fi
echo "  PASS"
