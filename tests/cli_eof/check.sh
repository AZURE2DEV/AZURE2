#!/usr/bin/env bash
#
# The program reads nothing from standard input.
#
# The text menu read its answers from standard input and re-asked until it got a valid one; at end of file such a loop
# printed its prompt forever (one MCMC run wrote 10 GB of "Number of Threads (1 or more):" and filled the disk), and a
# job started without `< /dev/null` waited on an open input.  The menu is gone: a run is a command with named options
# (AZURE2 <command> [options] FILE) and reads no input at all.
#
# Checked on tests/identical_pp_res (p+p, one pair, no external capture), copied:
#   (a) calculate and fit run to the end, and rate (an elastic pair) is refused at once, with an input that never ends
#       (`yes`) and with a closed input, exit 0, in a few seconds and with a bounded log;
#   (b) a run without a command ends at once with exit 2 and the refusal, whatever its input: none, the old menu
#       answers, an input that never ends;
#   (c) every required option of every command left out in turn: exit 2, the refusal names it, nothing calculated
#       (what the menu's end-of-input stops at a required prompt became).
#
#   ./tests/cli_eof/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
# A run that waits on its input is stopped here; a correct one ends in seconds.
LIMIT="${EOF_TIME_LIMIT:-60}"
RUN="$(guard_command "$LIMIT")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
CAP=2000000
MAXLOG=100000

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cli_eof.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
fail=0

SRC="$HERE/../identical_pp_res"
n=0
# go NAME INPUT ARGS...: AZURE2 ARGS on a fresh copy with INPUT on standard input ("never": `yes`, an input that
# never ends; "closed": < /dev/null; else printf %b of it); NAME/status, NAME/log, NAME/seconds
go() {
  local name="$1" input="$2" d="$WORK/$1"; shift 2
  mkdir -p "$d/output" "$d/checks"
  cp -R "$SRC/data" "$d/data"
  cp "$SRC/identical_pp_res.azr" "$d/run.azr"
  local t0=$SECONDS
  case "$input" in
    never) (cd "$d" && { yes 2> /dev/null || true; } | $RUN "$AZURE2_BIN" "$@" 2>&1 | head -c "$CAP" > log
            echo "${PIPESTATUS[1]}" > status) ;;
    closed) (cd "$d" && $RUN "$AZURE2_BIN" "$@" < /dev/null 2>&1 | head -c "$CAP" > log
             echo "${PIPESTATUS[0]}" > status) ;;
    *) (cd "$d" && printf '%b' "$input" | $RUN "$AZURE2_BIN" "$@" 2>&1 | head -c "$CAP" > log
        echo "${PIPESTATUS[1]}" > status) ;;
  esac
  echo $((SECONDS - t0)) > "$d/seconds"
}
bounded() { [ "$(wc -c < "$WORK/$1/log")" -le "$MAXLOG" ] && [ "$(cat "$WORK/$1/seconds")" -lt "$LIMIT" ]; }

# (a) the commands do not read their input
for input in never closed; do
  for args in "calculate" "fit" \
              "rate --rate-entrance=1 --rate-exit=1 --temperature-min=0.1 --temperature-max=1 --temperature-step=0.1"; do
    n=$((n + 1)); go "a$n" "$input" $args run.azr
    set -- $args
    if [ "$1" = rate ]; then  # a single pair: the rate of elastic scattering is refused before anything is read
      [ "$(cat "$WORK/a$n/status")" = 2 ] && bounded "a$n" \
        && ok "(a) rate (elastic) with input '$input': refused at once, exit 2" \
        || bad "(a) rate (elastic) with input '$input': status $(cat "$WORK/a$n/status"), $(cat "$WORK/a$n/seconds") s"
      continue
    fi
    [ "$(cat "$WORK/a$n/status")" = 0 ] && bounded "a$n" && [ -s "$WORK/a$n/output/param.par" ] \
      && ok "(a) $1 with input '$input': ran to the end, exit 0, $(wc -c < "$WORK/a$n/log") bytes of log" \
      || bad "(a) $1 with input '$input': status $(cat "$WORK/a$n/status"), $(cat "$WORK/a$n/seconds") s, $(wc -c < "$WORK/a$n/log") bytes"
  done
done

# (b) without a command
for input in closed '1\n\n\n7\n' never; do
  n=$((n + 1)); go "b$n" "$input" run.azr
  [ "$(cat "$WORK/b$n/status")" = 2 ] && bounded "b$n" && grep -q "ERROR: No command" "$WORK/b$n/log" \
    && [ -z "$(find "$WORK/b$n/output" -type f)" ] \
    && ok "(b) no command, input '$input': refused at once, exit 2" \
    || bad "(b) no command, input '$input': status $(cat "$WORK/b$n/status"), $(head -2 "$WORK/b$n/log")"
done

# (c) each required option left out in turn
# leave NAME OPTION COMMAND ARGS...: COMMAND ARGS (without OPTION) is refused naming OPTION
leave() {
  local name="$1" option="$2"; shift 2
  go "$name" closed "$@" run.azr
  [ "$(cat "$WORK/$name/status")" = 2 ] && grep -q "needs .*$option" "$WORK/$name/log" \
    && [ -z "$(find "$WORK/$name/output" -type f)" ] \
    && ok "(c) $1 without $option: refused, exit 2, nothing calculated" \
    || bad "(c) $1 without $option: status $(cat "$WORK/$name/status"), $(grep -m1 ERROR "$WORK/$name/log")"
}
leave c1 --delta-chi2 minos
R="--rate-entrance=1 --rate-exit=-1"
T="--temperature-min=0.1 --temperature-max=1 --temperature-step=0.1"
leave c2 --rate-entrance rate --rate-exit=-1 $T
leave c3 --rate-exit rate --rate-entrance=1 $T
leave c4 --temperature-min rate $R --temperature-max=1 --temperature-step=0.1
leave c5 --temperature-max rate $R --temperature-min=0.1 --temperature-step=0.1
leave c6 --temperature-step rate $R --temperature-min=0.1 --temperature-max=1
if grep -q '^	mcmc:' <<< "$("$AZURE2_BIN" --help 2>&1)"; then
  leave c7 --walkers mcmc --steps=100 --sampler-threads=1 --basis=rwa
  leave c8 --steps mcmc --walkers=4 --sampler-threads=1 --basis=rwa
  leave c9 --sampler-threads mcmc --walkers=4 --steps=100 --basis=rwa
  leave c10 --basis mcmc --walkers=4 --steps=100 --sampler-threads=1
else
  echo "  (no sampler in this build: mcmc not checked)"
fi

[ "$fail" -eq 0 ] && echo "  PASS" || echo "  FAIL"
[ "$fail" -eq 0 ]
