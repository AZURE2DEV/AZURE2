#!/bin/bash
# The build-time assertion of the calculation plan on every project of tests/ (plan_test.cpp): the tests/ tree is
# copied to a scratch directory (set-up may write output/ and checks/), every .azr outside output/ is set up there.
#   plan_test.sh path/to/plan_test path/to/tests
set -u
export LC_ALL=C
[ $# -ge 2 ] || { echo "usage: plan_test.sh path/to/plan_test path/to/tests"; exit 2; }
BIN="$1"; SRC="$(cd "$2" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cp -r "$SRC" "$WORK/tests"
cd "$WORK/tests"
projects=()
while IFS= read -r f; do
  d="$(dirname "$f")"; mkdir -p "$d/output" "$d/checks"; projects+=("$WORK/tests/$f")
done < <(find . -name '*.azr' -not -path '*/output/*' | sed 's|^\./||' | sort)
echo "${#projects[@]} projects"
"$BIN" "${projects[@]}"
