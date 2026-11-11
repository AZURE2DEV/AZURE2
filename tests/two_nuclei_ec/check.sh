#!/usr/bin/env bash
#
# External capture in the second compound nucleus of a project.  Two capture
# projects (models/A.azr: 12C(p,g0)13N, models/B.azr: 12C(n,g0)13C, both with an
# external-capture ground state) share one file (models/M.azr).  Each nucleus's
# capture extrapolation in the merged file must reproduce its standalone project
# to 1e-10: the nuclei do not couple.
#
# Before the fix the external-capture code looked the capture decay up by the
# exit pair's NUMBER in the entrance pair's decay list (EPoint::CalculateECAmplitudes,
# EPoint::GetECAmplitude, EData's intEC walk and writers).  That is the decay number
# only while every pair is a decay of the entrance pair, i.e. in a one-nucleus
# project; for the second nucleus (pair numbers 3 and 4, decay numbers 1 and 2)
# it indexed past the list and the run died with a segmentation fault.
#
#   ./tests/two_nuclei_ec/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/two_nuclei_ec.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {  # NAME MODEL FLAGS
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$2" "$dir/t.azr"
  (cd "$dir" && printf '3\nn\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
}
xs() { awk 'NF >= 4 { printf "%.12e %.12e\n", $1, $(NF - 1) }' "$1"; }
maxrel() {
  paste "$1" "$2" | awk '
    NF >= 4 { d = $2 - $4; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300;
              r = d / m; if (r > worst) worst = r; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t + 0) }'; }

M="$HERE/models"
for flags in "" "--use-brune" "--no-ec-table"; do
  tag="${flags:-default}"
  run "A_$tag" "$M/A.azr" "$flags" || true
  run "B_$tag" "$M/B.azr" "$flags" || true
  run "M_$tag" "$M/M.azr" "$flags" || true
  fa="$WORK/A_$tag/output/AZUREOut_aa=1_R=2.extrap"; fb="$WORK/B_$tag/output/AZUREOut_aa=1_R=2.extrap"
  m1="$WORK/M_$tag/output/AZUREOut_aa=1_R=2.extrap"; m2="$WORK/M_$tag/output/AZUREOut_aa=3_R=4.extrap"
  if [ ! -s "$fa" ] || [ ! -s "$fb" ]; then bad "$tag: a standalone extrapolation is missing"; continue; fi
  if [ ! -s "$m1" ] || [ ! -s "$m2" ]; then
    bad "$tag: the merged project wrote no extrapolation ($(tail -1 "$WORK/M_$tag/run.log" | cut -c1-80))"; continue
  fi
  xs "$fa" > "$WORK/a"; xs "$fb" > "$WORK/b"; xs "$m1" > "$WORK/m1"; xs "$m2" > "$WORK/m2"
  n="$(wc -l < "$WORK/b")"
  [ "$n" -gt 100 ] || bad "$tag: only $n points in B"
  # the two capture cross sections must differ, or the comparison could not tell the nuclei apart
  r="$(maxrel "$WORK/a" "$WORK/b")"
  if le "$r" 1e-3; then bad "$tag: A and B are the same curve (max rel $r)"; fi
  r="$(maxrel "$WORK/m1" "$WORK/a")"
  if le "$r" 1e-10; then ok "$tag: 12C(p,g)13N in the merged file == A.azr (max rel $r)"; else bad "$tag: merged 13N capture vs A.azr differ (max rel $r)"; fi
  r="$(maxrel "$WORK/m2" "$WORK/b")"
  if le "$r" 1e-10; then ok "$tag: 12C(n,g)13C in the merged file == B.azr (max rel $r)"; else bad "$tag: merged 13C capture vs B.azr differ (max rel $r)"; fi
done

if [ "$failures" -eq 0 ]; then echo "PASS: external capture in two compound nuclei"; else echo "FAIL: $failures check(s)"; exit 1; fi
