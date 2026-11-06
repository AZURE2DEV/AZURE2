# Sourced by tests/run_tests.sh and tests/*/check.sh: a time limit for one
# AZURE2 run, so that a hang fails the check instead of stalling the suite
# (an unparseable data line once spun ESegment::Fill forever at 100% CPU).
#
#   RUN="$(guard_command 600)"
#   printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr
#
# GNU timeout where there is one, gtimeout (Homebrew coreutils) on macOS,
# else run_guard below: plain bash, so the limit holds on every platform.

# run_guard SECONDS CMD [ARGS...]: CMD with this shell's stdin, sent SIGTERM
# after SECONDS; its exit status (143 when stopped).
run_guard() {
  local limit="$1" t=0
  shift
  "$@" <&0 &
  local pid=$!
  # The watchdog polls once a second and ends with the command; its output
  # goes nowhere, so a pipe after the command sees EOF when the command ends.
  (
    while kill -0 "$pid" 2>/dev/null; do
      if [ "$t" -ge "$limit" ]; then kill -TERM "$pid" 2>/dev/null; break; fi
      sleep 1
      t=$((t + 1))
    done
  ) </dev/null >/dev/null 2>&1 &
  wait "$pid"
}

# guard_command SECONDS: the prefix for $RUN.
guard_command() {
  if command -v timeout >/dev/null 2>&1; then
    echo "timeout $1"
  elif command -v gtimeout >/dev/null 2>&1; then
    echo "gtimeout $1"
  else
    echo "run_guard $1"
  fi
}
