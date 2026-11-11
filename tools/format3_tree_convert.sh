#!/usr/bin/env bash
#
# Replace old project files of the repository by their format-3 conversions, in place (docs/upgrade/format3_tree.md,
# docs/upgrade/project_format3.md section 8.5 item 6).  For each FILE (a path relative to the repository root):
#   AZURE2 --convert FILE <tmp> --to 3 --exact-grids --json-report $REPORTS/<FILE>.report.json      (one thread)
# in FILE's own folder (the converter reads the data files relative to it); then the original moves to
# tests/convert/fixtures/FILE (git mv, so the converter, its spec-equivalence test and the converted-tree gate keep
# reading it) and the conversion takes FILE's name.  Side files the converter writes (<stem>.f3<ext>) stay beside FILE
# and are added to git.  A file that is format 3 already, or that the converter refuses, is left alone and named.
#
#   tools/format3_tree_convert.sh ENGINE FILE...
#   REPORTS=<dir> (default $WORK/format3_tree_reports; WORK is required when REPORTS is not set)
#
# Run from anywhere; the repository is the one this script lies in.  Exit 1 when a file was refused, else 0.

[ $# -ge 2 ] || { echo "usage: [REPORTS=dir] format3_tree_convert.sh ENGINE FILE..." >&2; exit 2; }
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
shift
if [ -z "${REPORTS:-}" ]; then
  [ -n "${WORK:-}" ] || { echo "format3_tree_convert.sh: set REPORTS or WORK" >&2; exit 2; }
  REPORTS="$WORK/format3_tree_reports"
fi
export OMP_NUM_THREADS=1 LC_ALL=C
status=0
for rel in "$@"; do
  rel="${rel#./}"
  f="$ROOT/$rel"
  if [ ! -f "$f" ]; then echo "missing   $rel"; status=1; continue; fi
  if awk '/<config>/ { c = 1 } /<\/config>/ { c = 0 } c && /^format=3[[:space:]]*$/ { found = 1 } END { exit !found }' "$f"; then
    echo "format 3  $rel (left as it is)"
    continue
  fi
  dir="$(dirname "$f")"
  tmp="$dir/.format3_tree_convert.$$.azr"
  report="$REPORTS/$rel.report.json"
  mkdir -p "$(dirname "$report")"
  if (cd "$dir" && "$BIN" --convert "$(basename "$f")" "$(basename "$tmp")" --to 3 --exact-grids \
        --json-report "$report" > "$REPORTS/$rel.convert.log" 2>&1); then
    fixture="$ROOT/tests/convert/fixtures/$rel"
    mkdir -p "$(dirname "$fixture")"
    git -C "$ROOT" mv "$rel" "tests/convert/fixtures/$rel"
    mv "$tmp" "$f"
    git -C "$ROOT" add "$rel"
    # side files: <stem>.f3<ext> beside the project, written by this conversion
    for side in "$dir"/*.f3.*; do
      [ -f "$side" ] && [ "$side" -nt "$fixture" ] && git -C "$ROOT" add "${side#"$ROOT"/}" && echo "side file ${side#"$ROOT"/}"
    done
    echo "converted $rel"
  else
    rm -f "$tmp"
    echo "REFUSED   $rel: $(grep -m1 '^ERROR' "$REPORTS/$rel.convert.log")"
    status=1
  fi
done
exit $status
