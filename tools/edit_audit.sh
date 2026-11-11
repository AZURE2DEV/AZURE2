#!/usr/bin/env bash
# The audit of the checks' edits (docs/upgrade/format3_tree.md section 7): does every edit a tests/*/check.sh makes to
# a project, data or parameter file change that file, and does every variant run give a result of its own?
#
#   tools/edit_audit.sh <AZURE2> <out dir> [check name ...]     (default: every tests/*/check.sh, one after another)
#   python3 tools/edit_audit_report.py <out dir>                 the table: per check, edits that hit nothing and
#                                                                 variant runs whose outputs equal another run's
#
# Each check runs as tests/run_tests.sh runs it, with sed, awk, perl and python3 replaced on PATH by
# tools/edit_audit/shim.sh and the binary by tools/edit_audit/azure2_wrap.sh: both run the real program unchanged and
# log what it was given and what it changed to <out dir>/edits.tsv.  The check's own verdict is in <out dir>/status.tsv
# and its output in <out dir>/logs/<name>.log.  Threads: the checks' own (most take an exported OMP_NUM_THREADS, else
# 2); EDIT_AUDIT_CPU=<n> pins every check to that one CPU (taskset), whatever thread counts a check sets itself.

set -o pipefail
[ $# -ge 2 ] || { sed -n '2,13p' "$0"; exit 2; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BIN="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
OUTDIR="$2"
shift 2
mkdir -p "$OUTDIR" && OUTDIR="$(cd "$OUTDIR" && pwd)"
[ -x "$BIN" ] || { echo "no binary $BIN" >&2; exit 2; }

for t in sed awk perl python3; do
  r="$(command -v "$t")" || { echo "no $t on PATH" >&2; exit 2; }
  export "EDIT_AUDIT_REAL_$t=$r"
done
export EDIT_AUDIT_REAL_AZURE2="$BIN"
SHIMS="$OUTDIR/shims"
mkdir -p "$SHIMS" "$OUTDIR/azbin" "$OUTDIR/logs" "$OUTDIR/tmp" "$OUTDIR/shimtmp"
for t in sed awk perl python3; do ln -sf "$HERE/edit_audit/shim.sh" "$SHIMS/$t"; done
ln -sf "$HERE/edit_audit/azure2_wrap.sh" "$OUTDIR/azbin/AZURE2"
export EDIT_AUDIT_LOG="$OUTDIR/edits.tsv" EDIT_AUDIT_TMP="$OUTDIR/shimtmp"

if [ $# -eq 0 ]; then
  set -- $(cd "$ROOT/tests" && for c in */check.sh; do [ -x "$c" ] && dirname "$c"; done)
fi
for name in "$@"; do
  check="$ROOT/tests/$name/check.sh"
  [ -x "$check" ] || { printf '%s\tmissing\t0\n' "$name" >> "$OUTDIR/status.tsv"; continue; }
  scope="$OUTDIR/tmp/$name"
  rm -rf "$scope" && mkdir -p "$scope"
  t0=$(date +%s)
  EDIT_AUDIT_CHECK="$name" EDIT_AUDIT_SCOPE="$scope" EDIT_AUDIT_KEEP="$OUTDIR/keys/$name" TMPDIR="$scope" \
    PATH="$SHIMS:$PATH" ${EDIT_AUDIT_CPU:+taskset -c "$EDIT_AUDIT_CPU"} "$check" "$OUTDIR/azbin/AZURE2" > "$OUTDIR/logs/$name.log" 2>&1 < /dev/null
  st=$?
  printf '%s\t%s\t%s\n' "$name" "$st" "$(( $(date +%s) - t0 ))" >> "$OUTDIR/status.tsv"
  echo "$name: exit $st ($(( $(date +%s) - t0 )) s)"
  rm -rf "$scope"
done
