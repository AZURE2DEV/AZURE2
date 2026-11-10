#!/usr/bin/env bash
#
# Convert a copy of the test tree to format 3 and summarise (docs/upgrade/project_format3.md section 8.5, the
# preparation of step 3b's gate).  Every .azr under tests/ and examples/ of SRC is converted with
#   AZURE2 --convert OLD OLD.f3.azr --to 3 --json-report OLD.report.json [CONVERT_OPTIONS]
# in a copy of the tree (the originals are never touched), then checked with AZURE2 --validate --json.  Every
# report and validation is kept beside its file.  Prints the per-file result and the table of section 8.1: per
# construct, the lines converted and the notes and questions they raised.
#
# Then every original is run by the old engine (format 2, menu 1, as tools/run_one.sh runs it, in a private copy of its
# top-level folder), and the files the old engine refuses (an ERROR line in its log) but the converter converts are
# listed with the old engine's message: the converter converts what is written and leaves the old engine's semantic
# refusals to the format-3 reader, which has to refuse each of them (step 3b; project_format3.md section 8.5).
#
#   tools/convert_tree.sh <AZURE2 binary> <source tree (repository root)> <output dir>
#   JOBS=<n> (default 1; one thread per run)   CONVERT_OPTIONS="--exact-grids" (passed to every conversion)
#   OLD_ENGINE=0 skips the old engine's runs.
#
# Exit 0 when every file converted and validated with zero errors, 1 otherwise, 2 on a usage error (or when OUT
# exists).

set -o pipefail
export LC_ALL=C OMP_NUM_THREADS=1

[ $# -eq 3 ] || { echo "usage: [JOBS=n] [CONVERT_OPTIONS=...] [OLD_ENGINE=0] convert_tree.sh AZURE2 SRC OUT" >&2; exit 2; }
BIN=$1; SRC=$2; OUT=$3
JOBS="${JOBS:-1}"
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
case "$OUT" in /*) ;; *) OUT="$PWD/$OUT" ;; esac
case "$SRC" in /*) ;; *) SRC="$PWD/$SRC" ;; esac
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
    # in the file's folder: data files are read relative to the working directory, as tools/run_one.sh runs the engine
    if (cd "$(dirname "$new")" && "$BIN" --validate "$(basename "$new")" --json) > "$val" 2> "${f%.azr}.validate.log"; then
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

# The old engine on every original: one job per top-level folder (its files share data and output folders).
old_engine_folder() {
  local top="$1" flat="${1//\//__}"
  local work="$OUT/.old-engine/$flat"
  mkdir -p "$work" && cp -r "$SRC/$top" "$work/" || return
  find "$work/$(basename "$top")" -name '*.azr' -not -path '*/output/*' | sort | while read -r f; do
    local rel="${f#$work/}" log
    log="$OUT/$(dirname "$top")/${rel%.azr}.oldengine.log"
    ( cd "$(dirname "$f")" && rm -rf output checks && mkdir -p output checks &&
      printf '1\n\n\n7\n' | OMP_NUM_THREADS=1 timeout 900 "$BIN" --no-gui --no-readline "$(basename "$f")" 2>&1 |
        head -c 2000000 > "$log" )
  done
  rm -rf "$work"
}
export -f old_engine_folder
export OUT SRC
if [ "${OLD_ENGINE:-1}" != 0 ]; then
  (cd "$SRC" && find tests examples -mindepth 1 -maxdepth 1 -type d) | sort |
    xargs -P "$JOBS" -I{} bash -c 'old_engine_folder "$@"' _ {}
  echo
  echo "converted although the old engine refuses them (the format-3 reader of step 3b must refuse each):"
  grep '^ok ' "$OUT/results.txt" | awk '{print $2}' | while read -r f; do
    log="${f%.azr}.oldengine.log"
    [ -f "$log" ] && grep -q '^ERROR' "$log" && echo "  ${f#$OUT/}: $(grep -m1 '^ERROR' "$log")"
  done
  echo "refused by the converter, and by the old engine:"
  grep -v '^ok ' "$OUT/results.txt" | awk '{print $2}' | tr -d ':' | while read -r f; do
    log="${f%.azr}.oldengine.log"
    [ -f "$log" ] && echo "  ${f#$OUT/}: $(grep -m1 '^ERROR' "$log" || echo '(no ERROR line)')"
  done
fi

total=$(wc -l < "$OUT/results.txt")
good=$(grep -c '^ok ' "$OUT/results.txt")
echo
echo "$good of $total files converted and validated with zero errors"
[ "$good" -eq "$total" ]
