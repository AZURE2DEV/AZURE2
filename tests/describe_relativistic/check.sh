#!/usr/bin/env bash
#
# --validate and --describe read the project's <config> relativistic setting, as a run does.
#
# The setting reaches the data-set build through the project reader (experiment/ProjectReader.cpp, RelativisticFlag);
# a THM data set is refused with it (E_THM_RELATIVISTIC, the THM vertex is non-relativistic).  Model: tests/18O_p_a_thm
# with relativistic=1 added to <config>.
#   (a) without the line both commands accept the project (exit 0);
#   (b) with it, --validate reports E_THM_RELATIVISTIC and exits 2;
#   (c) with it, --describe reports E_THM_RELATIVISTIC and exits 2.
#
#   ./tests/describe_relativistic/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../18O_p_a_thm"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/describe_relativistic.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

for v in plain rel; do
  mkdir -p "$WORK/$v"
  cp -R "$SRC/data" "$WORK/$v/data"
  awk -v rel="$([ "$v" = rel ] && echo 1 || echo 0)" '{ sub(/\r$/, "") } /^<\/config>$/ && rel == 1 { print "relativistic=1" } { print }' \
    "$SRC/18O_p_a_thm.azr" > "$WORK/$v/t.azr"
done
cmd() { (cd "$WORK/$1" && "$AZURE2_BIN" "$2" t.azr > "$2.log" 2>&1); echo $?; }

a1="$(cmd plain --validate)"; a2="$(cmd plain --describe)"
if [ "$a1" = 0 ] && [ "$a2" = 0 ]; then ok "(a) the project without the line validates and is described"; else bad "(a) exits $a1 / $a2 without the line"; fi
b="$(cmd rel --validate)"
if [ "$b" = 2 ] && grep -q E_THM_RELATIVISTIC "$WORK/rel/--validate.log"; then ok "(b) --validate refuses the THM data set under relativistic=1"; else bad "(b) --validate exit $b, E_THM_RELATIVISTIC printed: $(grep -c E_THM_RELATIVISTIC "$WORK/rel/--validate.log")"; fi
c="$(cmd rel --describe)"
if [ "$c" = 2 ] && grep -q E_THM_RELATIVISTIC "$WORK/rel/--describe.log"; then ok "(c) --describe refuses the THM data set under relativistic=1"; else bad "(c) --describe exit $c, E_THM_RELATIVISTIC printed: $(grep -c E_THM_RELATIVISTIC "$WORK/rel/--describe.log")"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
