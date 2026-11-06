#!/usr/bin/env bash
#
# Everything a change has to pass before it is merged into dev, in the order
# that fails fastest: repository consistency (a second), build, CTest without
# the slow label (minutes), the physics regression suite (tens of minutes) and
# the documentation with warnings as errors (if Sphinx is installed).
#
#   ./scripts/check.sh [--quick] [build-dir]
#
#   --quick     skip the physics regression suite and the docs
#   build-dir   default build/; must already be configured (see README.md)
#
# Exits non-zero at the first failing stage and says which one.

set -uo pipefail
export LC_ALL=C

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QUICK=0
BUILD_DIR=""
for arg in "$@"; do
  case "$arg" in
    --quick) QUICK=1 ;;
    -h|--help) sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) BUILD_DIR="$arg" ;;
  esac
done
BUILD_DIR="${BUILD_DIR:-$REPO_ROOT/build}"
case "$BUILD_DIR" in /*) ;; *) BUILD_DIR="$PWD/$BUILD_DIR" ;; esac

if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
  echo "ERROR: $BUILD_DIR is not a configured build directory." >&2
  echo "       cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -Dpybind11_DIR=\"\$(python3 -m pybind11 --cmakedir)\"" >&2
  exit 2
fi

stage() { echo; echo "==== $1"; }
fail() { echo; echo "CHECK FAILED at: $1"; exit 1; }

stage "repository consistency"
python3 "$REPO_ROOT/tests/consistency/repo_consistency_test.py" || fail "repository consistency"

stage "build ($BUILD_DIR)"
cmake --build "$BUILD_DIR" --parallel || fail "build"
if ! ls "$REPO_ROOT"/pyazr/_azure2* >/dev/null 2>&1; then
  echo "WARNING: pyazr/_azure2 was not built (pybind11 not found at configure time);"
  echo "         the pyazr tests that need the engine will skip their engine sections."
fi

stage "ctest -LE slow"
ctest --test-dir "$BUILD_DIR" -LE slow --output-on-failure -j "${CTEST_JOBS:-4}" || fail "ctest"

if [ "$QUICK" -eq 1 ]; then
  echo; echo "Quick check passed (physics regression suite and docs skipped)."
  exit 0
fi

stage "physics regression suite"
"$REPO_ROOT/tests/run_tests.sh" "$BUILD_DIR/src/AZURE2" || fail "physics regression suite"

stage "documentation"
if python3 -c 'import sphinx, sphinx_rtd_theme, sphinx_copybutton' 2>/dev/null; then
  python3 -m sphinx -W --keep-going -q -b html "$REPO_ROOT/docs/source" "$REPO_ROOT/docs/_build/html" \
    || fail "documentation (sphinx -W)"
else
  echo "Sphinx not installed (pip install -r docs/requirements.txt); docs not built. CI builds them."
fi

echo; echo "All checks passed."
