#!/usr/bin/env bash
#
# End of input at a CLI prompt ends the run.
#
# Every prompt of the text interface re-asks until it gets a valid answer.  At
# end of file (a piped script one line short, a closed stdin) getline() fails
# at once on every pass, so such a loop printed its prompt forever: one MCMC
# run wrote 10 GB of "Number of Threads (1 or more):" and filled the disk.  A
# prompt that needs an answer now stops the run with an error naming it and
# exit status 1; a prompt whose blank answer is a documented default (the
# parameter and integral file names, the band question, the MCMC spreads and
# the overwrite question) takes that default at end of file, as before.
#
# Each menu path is cut after every answer it needs, with and without
# --no-readline (the file-name prompts read through readline when it is
# built).  Every run must end within seconds, write a bounded log and, where
# a needed answer is missing, exit non-zero with the message.
#
# Model: tests/identical_pp_res (p+p, one pair, no external capture), copied.
#
#   ./tests/cli_eof/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
# A run that loops is stopped here; a correct one ends in well under a second
# (the paths below stop before any calculation).
LIMIT="${EOF_TIME_LIMIT:-30}"
RUN="$(guard_command "$LIMIT")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
# Log cap: a looping prompt writes megabytes a second, so the cap is reached
# long before the time limit; a correct run writes a few kB.
CAP=2000000
MAXLOG=100000

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cli_eof.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
fail=0

SRC="$HERE/../identical_pp_res"
n=0
# eof NAME INPUT [FLAG]: run the menu with INPUT (printf %b) on stdin and
# nothing after it; NAME/status holds the exit status, NAME/log the output.
eof() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp -R "$SRC/data" "$d/data"
  cp "$SRC/identical_pp_res.azr" "$d/run.azr"
  local t0=$SECONDS
  (cd "$d" && printf '%b' "$2" | $RUN "$AZURE2_BIN" --no-gui ${3:-} run.azr 2>&1 | head -c "$CAP" > log
   echo "${PIPESTATUS[1]}" > status)
  echo $((SECONDS - t0)) > "$d/seconds"
}
# stops INPUT PROMPT [FLAG]: INPUT ends where PROMPT needs an answer.
stops() {
  n=$((n + 1))
  local name="r$n" input="$1" prompt="$2" flag="${3:-}"
  eof "$name" "$input" "$flag"
  local d="$WORK/$name" st size secs
  st="$(cat "$d/status")"; size="$(wc -c < "$d/log" | tr -d ' ')"; secs="$(cat "$d/seconds")"
  local label="'$(printf '%s' "$input" | sed 's/\\n/|/g')' $flag"
  if [ "$size" -gt "$MAXLOG" ]; then
    bad "$label: $size bytes of output (a prompt loops at end of file)"; return
  fi
  if [ "$st" = 124 ] || [ "$st" = 143 ] || [ "$secs" -ge "$LIMIT" ]; then
    bad "$label: still running after ${secs} s (status $st)"; return
  fi
  if [ "$st" = 0 ]; then bad "$label: exit status 0 at end of input"; return; fi
  if grep -qF "ERROR: the input ended (end of file) at the prompt \"$prompt" "$d/log"; then
    ok "$label: stops at \"$prompt\" (status $st, $size bytes, $secs s)"
  else
    bad "$label: no end-of-input error naming \"$prompt\" (status $st)"; tail -3 "$d/log" | sed 's/^/        /'
  fi
}
# defaults INPUT [FLAG]: INPUT ends where only defaulted prompts remain;
# the run must end on its own, without the end-of-input error.
defaults() {
  n=$((n + 1))
  local name="r$n"
  eof "$name" "$1" "${2:-}"
  local d="$WORK/$name" st size secs
  st="$(cat "$d/status")"; size="$(wc -c < "$d/log" | tr -d ' ')"; secs="$(cat "$d/seconds")"
  local label="'$(printf '%s' "$1" | sed 's/\\n/|/g')' ${2:-}"
  if [ "$size" -gt "$MAXLOG" ] || [ "$st" = 124 ] || [ "$st" = 143 ] || [ "$secs" -ge "$LIMIT" ]; then
    bad "$label: did not end (status $st, $size bytes, ${secs} s)"
  elif grep -q "ERROR: the input ended" "$d/log"; then
    bad "$label: blank defaults refused at end of file"
  else
    ok "$label: defaults taken at end of file, ends (status $st, ${secs} s)"
  fi
}

# The prompts come in this order: the menu, the MINOS variance (4), the
# parameter file (blank: from the .azr), then the rate (5) or sampler (6)
# questions.
F=--no-readline
echo "== $F"
stops ""                         "azure2: (menu option)" $F
stops "x\n9\n\n"                 "azure2: (menu option)" $F
stops "4\n"                      "Allowed Chi-Squared Variance" $F
stops "4\n-1\nabc\n"             "Allowed Chi-Squared Variance" $F
stops "5\n\n"                    "Reaction Rate Entrance Pair" $F
stops "5\n\n1\n"                 "Reaction Rate Exit Pair" $F
stops "5\n\n1\n1\n"              "Reaction Rate Entrance Pair" $F
stops "5\n\n1\n2\n"              "Use temperatures from file" $F
stops "5\n\n1\n2\nmaybe\n"       "Use temperatures from file" $F
stops "5\n\n1\n2\nyes\n"         "Temperature File Name" $F
stops "5\n\n1\n2\nyes\nnofile\n" "Temperature File Name" $F
stops "5\n\n1\n2\nno\n"          "Reaction Rate Min Temp" $F
stops "5\n\n1\n2\nno\n0.1\n"     "Reaction Rate Max Temp" $F
stops "5\n\n1\n2\nno\n0.1\n1\n"  "Reaction Rate Temp Step" $F
if grep -q "Perform MCMC" "$WORK/r1/log"; then
  stops "6\n\n"                   "Number of Walkers" $F
  stops "6\n\n1\n"                "Number of Walkers" $F
  stops "6\n\n10\n"               "Number of Steps" $F
  stops "6\n\n10\n100\n"          "Number of Threads" $F
  stops "6\n\n10\n100\n\n\n"      "Number of Threads" $F
  stops "6\n\n10\n100\n5\n1\n0\n" "Number of Threads" $F
  stops "6\n\n10\n100\n\n\n1\n"   "Use Reduced Width Amplitudes" $F
  stops "6\n\n10\n100\n\n\n1\nmaybe\n" "Use Reduced Width Amplitudes" $F
else
  echo "  (no MCMC in this build: menu 6 not checked)"
fi
n=$((n + 1)); eof "r$n" "7\n" $F
[ "$(cat "$WORK/r$n/status")" = 0 ] && ok "'7' $F: exit, status 0" || bad "'7' $F: status $(cat "$WORK/r$n/status")"
defaults "1\n" $F
defaults "3\n" $F

# Through readline where it is built (elsewhere --no-readline is the only
# reader and these repeat the above).  readline reads the descriptor while the
# menu reads the buffered std::cin, so only inputs that end at or before the
# first file-name prompt have a defined order.
echo "== readline (when built)"
stops ""      "azure2: (menu option)"
stops "4\n"   "Allowed Chi-Squared Variance"
stops "5\n"   "Reaction Rate Entrance Pair"
if grep -q "Perform MCMC" "$WORK/r1/log"; then
  stops "6\n" "Number of Walkers"
fi
defaults "1\n"

[ "$fail" -eq 0 ] && echo "  PASS" || echo "  FAIL"
[ "$fail" -eq 0 ]
