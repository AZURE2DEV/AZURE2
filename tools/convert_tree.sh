#!/usr/bin/env bash
#
# Convert a copy of the test tree to format 3 and summarise (docs/upgrade/project_format3.md section 8.5, the
# preparation of step 3b's gate).  Every .azr under tests/ and examples/ of SRC is converted with
#   AZURE2 --convert OLD OLD.f3.azr --to 3 --json-report OLD.report.json [CONVERT_OPTIONS]
# in a copy of the tree (the originals are never touched), then checked with AZURE2 --validate --json.  Every
# report and validation is kept beside its file.  Prints the per-file result and the table of section 8.1: per
# construct, the lines converted and the notes and questions they raised.
#
#   tools/convert_tree.sh <AZURE2 binary> <source tree (repository root)> <output dir>
#   JOBS=<n> (default 1; one thread per conversion)   CONVERT_OPTIONS="--exact-grids" (passed to every conversion)
#
# Exit 0 when every file converted and validated with zero errors, 1 otherwise.

set -o pipefail
export LC_ALL=C OMP_NUM_THREADS=1

BIN="${1:?usage: convert_tree.sh AZURE2 SRC OUT}"
SRC="${2:?usage: convert_tree.sh AZURE2 SRC OUT}"
OUT="${3:?usage: convert_tree.sh AZURE2 SRC OUT}"
JOBS="${JOBS:-1}"
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
if [ -e "$OUT" ]; then
  echo "$OUT exists; remove it first" >&2
  exit 2
fi
mkdir -p "$OUT"
for d in tests examples; do
  [ -d "$SRC/$d" ] && cp -r "$SRC/$d" "$OUT/$d"
done
# outputs of earlier runs are not projects
find "$OUT" -path '*/output/*' -name '*.azr' -delete

convert_one() {
  local f="$1"
  local new="${f%.azr}.f3.azr" report="${f%.azr}.report.json" val="${f%.azr}.validate.json"
  # shellcheck disable=SC2086
  if "$BIN" --convert "$f" "$new" --to 3 --json-report "$report" $CONVERT_OPTIONS > "${f%.azr}.convert.log" 2>&1; then
    if "$BIN" --validate "$new" --json > "$val" 2> "${f%.azr}.validate.log"; then
      echo "ok        $f"
    else
      echo "INVALID   $f"
    fi
  else
    echo "REFUSED   $f: $(grep -m1 '^ERROR' "${f%.azr}.convert.log")"
  fi
}
export -f convert_one
export BIN CONVERT_OPTIONS

find "$OUT" -name '*.azr' ! -name '*.f3.azr' | sort | xargs -P "$JOBS" -I{} bash -c 'convert_one "$@"' _ {} > "$OUT/results.txt"
sort -k2 "$OUT/results.txt" -o "$OUT/results.txt"
cat "$OUT/results.txt"

python3 - "$OUT" <<'EOF'
import json, os, sys
out = sys.argv[1]
table, files, notes, questions = {}, 0, 0, 0
for root, _, names in os.walk(out):
    for n in names:
        if not n.endswith(".report.json"):
            continue
        r = json.load(open(os.path.join(root, n)))
        files += 1
        notes += len(r["notes"])
        questions += len(r["questions"])
        for c, v in r["constructs"].items():
            t = table.setdefault(c, [0, 0, 0, 0])
            t[0] += v["converted"]; t[1] += v["notes"]; t[2] += v["questions"]; t[3] += 1
print()
print(f"{'construct (section 8.1)':40s} {'converted':>9s} {'notes':>6s} {'questions':>9s} {'files':>6s}")
for c in sorted(table):
    t = table[c]
    print(f"{c:40s} {t[0]:9d} {t[1]:6d} {t[2]:9d} {t[3]:6d}")
print(f"{'all':40s} {'':9s} {notes:6d} {questions:9d} {files:6d}")
EOF

total=$(wc -l < "$OUT/results.txt")
good=$(grep -c '^ok ' "$OUT/results.txt")
echo
echo "$good of $total files converted and validated with zero errors"
[ "$good" -eq "$total" ]
