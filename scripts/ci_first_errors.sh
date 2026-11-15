#!/usr/bin/env bash
#
# The first error lines of a configure or build log, for a run whose full log the reader cannot open: on GitHub
# Actions one ::error:: annotation (public with the run), elsewhere plain lines on stderr.
#
#   scripts/ci_first_errors.sh LOG LABEL
#
# From the first line that looks like an error (a compiler's or a linker's, CMake's, Ninja's FAILED:) 40 lines, each
# cut at 300 characters; the last 25 lines of the log when no line does.  Always exits 0: the caller fails the step.
# Used by .github/workflows/build.yml (Configure, Build) and scripts/sanitize_upgrade.sh.
set -uo pipefail
LOG=${1:?usage: ci_first_errors.sh LOG LABEL}
LABEL=${2:?usage: ci_first_errors.sh LOG LABEL}
[ -r "$LOG" ] || { echo "ci_first_errors: no $LOG" >&2; exit 0; }
awk -v label="$LABEL" -v gh="${GITHUB_ACTIONS:-}" '
  { gsub(/\r/, ""); line[NR] = substr($0, 1, 300) }
  !first && /FAILED:|[Ee]rror[: ]|undefined reference|Undefined symbols|ld: |fatal|Killed/ { first = NR }
  END {
    if (!NR) { line[1] = "(the log is empty)"; NR = 1 }
    if (!first) first = (NR > 25 ? NR - 24 : 1)
    last = (first + 39 < NR ? first + 39 : NR)
    head = sprintf("%s, first error lines (log lines %d-%d of %d):", label, first, last, NR)
    if (gh == "true") {
      out = ""
      for (i = first; i <= last; i++) { l = line[i]; gsub(/%/, "%25", l); out = out "%0A" l }
      printf "::error::%s%s\n", head, out
    } else {
      print head > "/dev/stderr"
      for (i = first; i <= last; i++) print "  " line[i] > "/dev/stderr"
    }
  }' "$LOG"
exit 0
