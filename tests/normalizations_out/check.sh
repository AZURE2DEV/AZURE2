#!/usr/bin/env bash
#
# normalizations.out holds one COMPLETE row per segment whose normalization is fitted (audit item 18, 2026-10-27).
#
# EData::WriteOutputFiles wrote the first four fields (key, file, angle bounds) only for segments with a varied
# normalization but the remaining six (frames, energy bounds, norm, shift) for EVERY segment, so a project with fixed
# and varied norms got a six-field row without key, file or angles for each fixed one (11B+a wang_an_global:
# "Lab,2.0000e-01,2.0000e+00,Lab,1.0000e+00,0.0000e+00").  Now a segment is skipped unless its norm is varied or it
# belongs to a norm group, and every row has ten fields.
#
#   (a) tests/13N with the norm of segment 3 varied (5 %): the file is the header plus exactly one row, of ten fields,
#       starting "3,data/meyer_84.dat,8.4300e+01,8.4300e+01,Lab,";
#   (b) tests/norm_groups with the two normgroup D members' norms varied (the other segments fixed): the rows are
#       exactly those two, ten fields each.
# dev f18b36e: (a) 16 rows, 15 of them with six fields; (b) passes there too (both of that model's segments are varied).
#
#   ./tests/normalizations_out/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/normalizations_out.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
fields() { awk -F, 'NR > 1 { print NF }' "$1" | sort -n | uniq -c | awk '{ printf "%s rows x %s fields; ", $1, $2 }'; }

# ---- (a) 13N, one varied norm -----------------------------------------------------------------------------------
mkdir -p "$WORK/a/output" "$WORK/a/checks" && cp -r "$HERE/../13N/data" "$WORK/a/"
# data line: active ent exit minE maxE minA maxA isDiff norm varyNorm normErr ... ; the third data segment
awk 'BEGIN { n = 0 } /<segmentsData>/ { inb = 1 } /<\/segmentsData>/ { inb = 0 }
     { sub(/\r$/, ""); if (inb && $0 !~ /segmentsData/ && NF > 10) { n++; if (n == 3) { $10 = 1; $11 = 5; $1 = $1 } } print }' \
    "$HERE/../13N/13N.azr" > "$WORK/a/t.azr"
run a
f="$WORK/a/output/normalizations.out"
if [ -s "$f" ]; then
  nrows="$(awk 'NR > 1' "$f" | wc -l | tr -d ' ')"
  row="$(sed -n 2p "$f")"
  nf="$(echo "$row" | awk -F, '{ print NF }')"
  if [ "$nrows" -eq 1 ] && [ "$nf" -eq 10 ] && [[ "$row" == 3,data/meyer_84.dat,8.4300e+01,8.4300e+01,Lab,* ]]; then
    ok "(a) one varied norm: one complete row ($row)"
  else
    bad "(a) one varied norm: $nrows row(s), $(fields "$f")first row '$row'"
  fi
else
  bad "(a) no normalizations.out"; tail -3 "$WORK/a/run.log"
fi

# ---- (b) norm groups ----------------------------------------------------------------------------------------------
mkdir -p "$WORK/b/output" "$WORK/b/checks" && cp -r "$HERE/../norm_groups/data" "$WORK/b/"
# the two normgroup members with their norms varied (5 %); the other segments keep fixed norms
awk '/<segmentsData>/ { inb = 1 } /<\/segmentsData>/ { inb = 0 }
     { sub(/\r$/, ""); if (inb && $0 ~ /normgroup/) { $10 = 1; $11 = 5; $1 = $1 } print }' "$HERE/../norm_groups/norm_groups.azr" > "$WORK/b/t.azr"
run b
f="$WORK/b/output/normalizations.out"
want="$(awk '/<segmentsData>/ { inb = 1; next } /<\/segmentsData>/ { inb = 0 } inb && NF > 10 && $1 == 1 && $10 == 1 { n++ } END { print n + 0 }' "$WORK/b/t.azr")"
if [ -s "$f" ]; then
  nrows="$(awk 'NR > 1' "$f" | wc -l | tr -d ' ')"
  nbad="$(awk -F, 'NR > 1 && NF != 10' "$f" | wc -l | tr -d ' ')"
  if [ "$want" -ge 2 ] && [ "$nrows" -eq "$want" ] && [ "$nbad" -eq 0 ]; then
    ok "(b) norm groups: $nrows rows for the $want fitted norms, ten fields each"
  else
    bad "(b) norm groups: $nrows row(s) for $want fitted norms, $nbad malformed ($(fields "$f"))"
  fi
else
  bad "(b) no normalizations.out"; tail -3 "$WORK/b/run.log"
fi

if [ "$failures" -eq 0 ]; then echo "normalizations_out: PASS"; else echo "normalizations_out: $failures check(s) FAILED"; fi
[ "$failures" -eq 0 ]
