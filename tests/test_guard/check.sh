#!/usr/bin/env bash
#
# The time limit the checks put on each AZURE2 run (tests/lib/guard.sh) must
# hold where GNU timeout is missing (macOS without coreutils): the shell
# fallback run_guard is exercised here directly, on every platform.
#
#   ./tests/test_guard/check.sh [path/to/AZURE2]   (the binary is not used)

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
fail=0

run_guard 5 true && ok "a quick command: status 0" || bad "a quick command failed"
run_guard 5 false; [ $? -ne 0 ] && ok "a failing command keeps its status" || bad "a failure was passed"
t0=$(date +%s)
run_guard 1 sleep 30; status=$?
t1=$(date +%s)
[ "$status" -ne 0 ] && [ $((t1 - t0)) -lt 10 ] && ok "a hang is stopped after the limit ($((t1 - t0)) s, status $status)" \
  || bad "a hang was not stopped (status $status, $((t1 - t0)) s)"
out="$(printf 'one\ntwo\n' | run_guard 5 cat | head -c 100)"
[ "$out" = "$(printf 'one\ntwo')" ] && ok "stdin reaches the command, stdout the pipe" || bad "pipe: '$out'"
case "$(guard_command 7)" in
  "timeout 7" | "gtimeout 7" | "run_guard 7") ok "guard_command: $(guard_command 7)" ;;
  *) bad "guard_command: $(guard_command 7)" ;;
esac

[ "$fail" -eq 0 ] && echo "PASSED" || { echo "FAILED"; exit 1; }
