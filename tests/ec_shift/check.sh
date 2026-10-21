#!/usr/bin/env bash
#
# External capture at shifted energies (ECTable).  A capture segment with an
# energy shift S must give the same cross sections as the same data with their
# lab energies moved by S and no shift: both evaluate the external-capture
# amplitudes at the same c.m. energies, from the same energy tables.  Before the
# tables, a shifted point interpolated linearly between the amplitudes of the
# segment's own points and held them flat beyond its ends (still the behaviour
# of --no-ec-table, whose difference is printed for reference).
#
# Model: tests/13N (12C(p,g)), target integration off, the capture segments
# 1 (Artemov), 6 and 8 (Skowronski, LUNA HPGe and FELS) and 10 (Kettner, 55 deg)
# active, S = +30 and -30 keV lab, so the end points move past the segments' ends.
#
#   ./tests/ec_shift/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# The time limit of tests/lib/guard.sh: timeout, gtimeout or a shell watchdog.
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_shift.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

ACTIVE=" 1 6 8 10 "
# model SHIFT MOVED(0|1) OUT.azr: 13N with target integration off, only $ACTIVE on;
# MOVED=0 puts SHIFT (lab MeV) on those segments, MOVED=1 reads moved/<file> with no shift.
model() {
  awk -v S="$1" -v moved="$2" -v active="$ACTIVE" '
    /<targetInt>/ { ti = 1; print; next }
    /<\/targetInt>/ { ti = 0 }
    /<segmentsData>/ { sd = 1; print; next }
    /<\/segmentsData>/ { sd = 0 }
    ti && NF > 0 { sub(/^[ \t]*1/, "0") }
    sd && NF > 0 {
      n++
      if (index(active, " " n " ") == 0) { $1 = 0 }
      else if (moved == 1) { $12 = 0; f = $15; sub(/.*\//, "", f); $15 = "moved/" f }
      else { $12 = S }
    }
    { print }' "$SRC/13N.azr" > "$3"
}
# run NAME AZR FLAGS -> $WORK/NAME/output/AZUREOut_aa=1_R=2.out
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$SRC/data" "$dir/data"
  ln -s "$WORK/moved_$4" "$dir/moved" 2>/dev/null
  cp "$2" "$dir/t.azr"
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -s "$dir/output/AZUREOut_aa=1_R=2.out" ]; then
    bad "$1: no output written"; tail -15 "$dir/run.log"; return 1
  fi
}
# maxrel A B [any]: largest relative difference of the fitted cross sections (column 4),
# rows matched by c.m. energy and angle ("any": rows in order, energies may differ)
maxrel() {
  paste <(sort -g -k1,1 -k3,3 "$WORK/$1/output/AZUREOut_aa=1_R=2.out") \
        <(sort -g -k1,1 -k3,3 "$WORK/$2/output/AZUREOut_aa=1_R=2.out") | awk -v anyE="${3:-}" '
    NF >= 8 { n = NF / 2; de = $1 - $(n + 1); if (de < 0) de = -de; if (de > 1e-9 && anyE == "") bad = 1
              a = $4; b = $(n + 4); d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m < 1e-300) m = 1e-300
              r = d / m; if (r > worst) worst = r; rows++ }
    END { if (rows == 0 || bad) print "nan"; else printf "%.3e\n", worst }'
}
within() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r != "nan" && r + 0 <= t) }'; }

for S in 0.030 -0.030; do
  mkdir -p "$WORK/moved_$S"
  for f in "$SRC"/data/*.dat; do
    awk -v S="$S" 'NF >= 2 && $1 + 0 == $1 { $1 = sprintf("%.10e", $1 + S) } { print }' "$f" > "$WORK/moved_$S/$(basename "$f")"
  done
  model "$S" 0 "$WORK/shift_$S.azr"
  model "$S" 1 "$WORK/moved_$S.azr"
  run "tab_shift_$S" "$WORK/shift_$S.azr" "" "$S" || continue
  run "tab_moved_$S" "$WORK/moved_$S.azr" "" "$S" || continue
  r="$(maxrel "tab_shift_$S" "tab_moved_$S")"
  if within "$r" 1e-7; then ok "shift $S MeV: shifted segments = moved data (max rel diff $r)"
  else bad "shift $S MeV: shifted segments differ from moved data by $r (limit 1e-7)"; fi
  if run "old_shift_$S" "$WORK/shift_$S.azr" "--no-ec-table" "$S" && run "old_moved_$S" "$WORK/moved_$S.azr" "--no-ec-table" "$S"; then
    echo "        (--no-ec-table, interpolated shifted amplitudes: max rel diff $(maxrel "old_shift_$S" "old_moved_$S"))"
  fi
done

# The table against tightly converged direct integrals at every amplitude.
model 0 0 "$WORK/plain.azr"
if AZURE2_EC_TABLE_CHECK=1 run "check" "$WORK/plain.azr" "" "0.030"; then
  d="$(sed -n 's/.*largest relative difference from the direct integrals \([^ ]*\).*/\1/p' "$WORK/check/run.log" | head -1)"
  if within "${d:-nan}" 1e-7; then ok "table vs converged integrals: largest relative difference $d"
  else bad "table vs converged integrals: largest relative difference ${d:-missing} (limit 1e-7)"; fi
fi
# The shift must have moved the cross sections (else the comparisons above prove nothing).
for S in 0.030 -0.030; do
  r="$(maxrel "tab_shift_$S" check any)"
  if [ "$r" = nan ] || within "$r" 1e-3; then bad "shift $S MeV changed the cross sections by only $r"
  else ok "shift $S MeV moves the cross sections (max rel change $r)"; fi
done
# ... and against the direct integrals (relative tolerance 1e-4) of --no-ec-table.
if run "direct" "$WORK/plain.azr" "--no-ec-table" "0.030"; then
  r="$(maxrel check direct)"
  if within "$r" 1e-3; then ok "table vs --no-ec-table, no shift: max rel diff $r"
  else bad "table vs --no-ec-table, no shift: max rel diff $r (limit 1e-3)"; fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS ec_shift"; exit 0; fi
echo "FAIL ec_shift ($failures)"; exit 1
