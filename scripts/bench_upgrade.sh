#!/usr/bin/env bash
#
# The T2 benchmark of the upgrade (docs/upgrade/ROADMAP.md, Phase 0; RULES.md §2.4): wall time and
# peak resident memory of the test project tests/upgrade/T2 at 1 and 8 threads, one line per run
# appended to a results file.
#
#   scripts/bench_upgrade.sh [path/to/AZURE2]
#
# Environment (all optional):
#   AZURE2_BINARY    the engine, if not given as the argument (default: build/src/AZURE2 of this repository)
#   BENCH_RESULTS    the results file (default: tests/upgrade/bench/T2.tsv of this repository)
#   BENCH_THREADS    thread counts to run (default: "1 8")
#   BENCH_COMMIT     the engine commit to record (default: git HEAD of this repository, "+dirty" if the
#                    work tree has uncommitted changes; set it when the binary was built from another tree)
#   BENCH_NOTE       free text for the last column (default: empty), e.g. the machine or "condor 8 cores"
#   BENCH_FIT_CALLS  function calls of the fit run (default 200)
#
# What is measured, per thread count:
#   calc  menu 1 (calculate with data) of T2.azr: all 13 segments, 340 points, every apparatus kernel once;
#   fit   menu 2 of T2.azr with --use-lm --max-calls=BENCH_FIT_CALLS from its start values: a fixed amount of fit
#         work (the THM segments have no analytic Jacobian, so this is MIGRAD with Minuit's numerical gradient, as
#         in tests/upgrade/T2/check_recovery.py; the full fit needs about 1100 calls, 585 s at 8 threads).  It
#         includes the sub-point grid re-placement and Coulomb functions at new sub-points that a calculation from
#         a fixed parameter set does not exercise.
# Both run with --no-coulomb-table (exact Coulomb functions), so a Coulomb table on one machine and not on
# another does not change the timing by its factor of about 4.  Peak RSS is that of the AZURE2 process,
# from GNU time (-v, Linux) or BSD time (-l, macOS); without either it is recorded as "na".
#
# Results file: tab-separated, '#' comment lines, one row per run with the columns
#   date_utc  commit  host  ncpu  threads  run  wall_s  peak_rss_MB  chi2  N  note
# run is "calc" or "fit<N>" (N = BENCH_FIT_CALLS).
# chi2 and N are the run's Total-Chi-Squared and Total-N (chiSquared.out): a benchmark whose chi2 moved is
# not timing the same calculation.  Compare rows of the same host and threads only.

set -uo pipefail
export LC_ALL=C
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${1:-${AZURE2_BINARY:-$REPO_ROOT/build/src/AZURE2}}"
[ -x "$BIN" ] || { echo "ERROR: AZURE2 binary '$BIN' not found (argument or AZURE2_BINARY)" >&2; exit 2; }
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
RESULTS="${BENCH_RESULTS:-$REPO_ROOT/tests/upgrade/bench/T2.tsv}"
FITCALLS="${BENCH_FIT_CALLS:-200}"
THREADS="${BENCH_THREADS:-1 8}"
NOTE="${BENCH_NOTE:-}"
T2="$REPO_ROOT/tests/upgrade/T2"
if [ -n "${BENCH_COMMIT:-}" ]; then
  COMMIT="$BENCH_COMMIT"
else
  COMMIT="$(git -C "$REPO_ROOT" rev-parse --short=8 HEAD 2>/dev/null || echo unknown)"
  git -C "$REPO_ROOT" diff --quiet HEAD 2>/dev/null || COMMIT="$COMMIT+dirty"
fi
NCPU="$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo na)"
HOST="$(hostname -s 2>/dev/null || hostname)"

# time_cmd OUTFILE CMD...: run CMD, write "wall_s rss_MB" to OUTFILE
if /usr/bin/time -v true >/dev/null 2>&1; then TIMEFLAG=-v
elif /usr/bin/time -l true >/dev/null 2>&1; then TIMEFLAG=-l
else TIMEFLAG=""; fi
time_cmd() {
  local out="$1"; shift
  local t0 t1 rss=na
  t0="$(date +%s.%N 2>/dev/null)"; case "$t0" in *N) t0="$(python3 -c 'import time; print(time.time())')";; esac
  if [ -n "$TIMEFLAG" ]; then /usr/bin/time "$TIMEFLAG" -o "$out.time" "$@"; else "$@"; fi
  t1="$(date +%s.%N 2>/dev/null)"; case "$t1" in *N) t1="$(python3 -c 'import time; print(time.time())')";; esac
  if [ "$TIMEFLAG" = -v ]; then
    rss="$(awk -F: '/Maximum resident set size/ { printf "%.1f", $2 / 1024 }' "$out.time")"
  elif [ "$TIMEFLAG" = -l ]; then
    rss="$(awk '/maximum resident set size/ { printf "%.1f", $1 / 1048576 }' "$out.time")"
  fi
  awk -v a="$t0" -v b="$t1" -v r="${rss:-na}" 'BEGIN { printf "%.2f %s\n", b - a, r }' > "$out"
}

WORK="$(mktemp -d "${TMPDIR:-/tmp}/bench_t2.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$(dirname "$RESULTS")"
if [ ! -s "$RESULTS" ]; then
  {
    echo "# T2 benchmark (scripts/bench_upgrade.sh): wall time and peak RSS of AZURE2 on tests/upgrade/T2"
    echo "# calc = menu 1, fit = menu 2 --use-lm --max-calls=N (numerical MIGRAD, N calls); --no-coulomb-table; compare same host/threads"
    printf 'date_utc\tcommit\thost\tncpu\tthreads\trun\twall_s\tpeak_rss_MB\tchi2\tN\tnote\n'
  } > "$RESULTS"
fi

status=0
for th in $THREADS; do
  for run in calc fit; do
    d="$WORK/$run-$th"
    mkdir -p "$d/output" "$d/checks"
    cp -r "$T2/T2.azr" "$T2/cov_diff.txt" "$T2/prior_cov.txt" "$T2/data" "$d/"
    menu='1\n\n\n7\n'; extra=""; [ "$run" = fit ] && { menu='2\nn\n\n\n7\n'; extra="--use-lm --max-calls=$FITCALLS"; }
    printf "$menu" > "$d/menu"
    ( cd "$d" && OMP_NUM_THREADS="$th" time_cmd "$d/meas" \
        "$BIN" --no-gui --no-readline --no-coulomb-table $extra T2.azr < "$d/menu" > "$d/log" 2>&1 )
    chi2="$(awk '/^Total-Chi-Squared:/ { print $2 }' "$d/output/chiSquared.out" 2>/dev/null)"
    n="$(awk '/^Total-Chi-Squared:/ { for (i = 1; i < NF; i++) if ($i == "Total-N:") print $(i + 1) }' \
          "$d/output/chiSquared.out" 2>/dev/null)"
    if [ -z "$chi2" ]; then
      echo "ERROR: $run at $th threads wrote no chiSquared.out; last lines:" >&2
      tail -5 "$d/log" >&2
      status=1; continue
    fi
    read -r wall rss < "$d/meas"
    label="$run"; [ "$run" = fit ] && label="fit$FITCALLS"
    row="$(date -u +%Y-%m-%dT%H:%M:%SZ)\t$COMMIT\t$HOST\t$NCPU\t$th\t$label\t$wall\t$rss\t$chi2\t$n\t$NOTE"
    printf "$row\n" >> "$RESULTS"
    printf "$row\n"
  done
done
exit $status
