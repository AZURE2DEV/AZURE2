#!/usr/bin/env bash
#
# ThreadSanitizer or AddressSanitizer run of the upgrade test projects and of the tests that exercise point
# mapping, energy shifts and target integration, at OMP_NUM_THREADS=8 (docs/upgrade/ROADMAP.md, Phase 0;
# RULES.md §2.4).  Exits non-zero on any sanitizer report.  What runs, the toolchain it needs and the findings
# so far: docs/upgrade/sanitizer_findings.md.
#
#   scripts/sanitize_upgrade.sh thread|address [build-dir]
#
# The build directory (default build-san-<kind> in this repository) is configured with
# -DAZURE2_SANITIZE=<kind> (CMakeLists.txt) and built unless SAN_NO_BUILD=1; an existing directory configured
# for another sanitizer is refused.  The compiler is CMake's choice, or $CC/$CXX.
#
# Environment (all optional):
#   SAN_THREADS        OpenMP threads (default 8)
#   SAN_FIT_CALLS      function calls of the bounded T2 fit (default 30; the full fit needs about 1100)
#   SAN_SECTIONS       space-separated subset of: upgrade checks regress reference pyazr (default: all)
#   SAN_CMAKE_ARGS     extra configure arguments (word-split), e.g. a compiler resource directory
#   SAN_JOBS           parallel build jobs (default: CMake's)
#   SAN_PYTHON         interpreter for the pyazr tests and for pybind11 (default python3)
#   SAN_ARCHER         libarcher.so for TSan (default: next to the libomp the binary loads)
#   SAN_RUNTIME        sanitizer runtime to LD_PRELOAD into python (default: asked from the compiler)
#   SAN_REPORT_DIR     keep the reports and logs here (default: a temporary directory, removed at exit unless
#                      something failed; its path is printed)
#   SAN_LEAKS          1: LeakSanitizer on with ASan (default 0; see the findings document)
#   SAN_ALLOW_FAILURES 1: a test that fails without a sanitizer report does not fail the run (default 0)
#
# Exit status: 0 clean; 1 sanitizer reports; 2 test failures without reports; 3 setup error.
#
# Suppressions: scripts/sanitizers/<kind>.supp.  An entry is allowed only for a report traced to third-party
# code (Minuit2, GSL, the OpenMP runtime), with its reason above it.

set -uo pipefail
export LC_ALL=C

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
KIND="${1:-}"
case "$KIND" in
  thread) SANFLAG=thread; SHORT=tsan ;;
  address) SANFLAG=address; SHORT=asan ;;
  -h|--help|"") sed -n '2,32p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
  *) echo "usage: $0 thread|address [build-dir]" >&2; exit 3 ;;
esac
BUILD="${2:-$REPO_ROOT/build-san-$KIND}"
case "$BUILD" in /*) ;; *) BUILD="$PWD/$BUILD" ;; esac
THREADS="${SAN_THREADS:-8}"
FITCALLS="${SAN_FIT_CALLS:-30}"
SECTIONS="${SAN_SECTIONS:-upgrade checks regress reference pyazr}"
PY="${SAN_PYTHON:-python3}"
SUPP="$REPO_ROOT/scripts/sanitizers/$SHORT.supp"

die() { echo "sanitize_upgrade: $*" >&2; exit 3; }
want() { case " $SECTIONS " in *" $1 "*) return 0 ;; esac; return 1; }

# ------------------------------------------------------------------ build
if [ "${SAN_NO_BUILD:-0}" != 1 ]; then
  if [ -f "$BUILD/CMakeCache.txt" ]; then
    have="$(sed -n 's/^AZURE2_SANITIZE:STRING=//p' "$BUILD/CMakeCache.txt")"
    [ "$have" = "$SANFLAG" ] || die "$BUILD is configured with AZURE2_SANITIZE='$have', not '$SANFLAG'"
  fi
  PB=()
  if PBDIR="$("$PY" -m pybind11 --cmakedir 2>/dev/null)"; then
    PB=(-Dpybind11_DIR="$PBDIR" -DPYBIND11_FINDPYTHON=ON -DPython_EXECUTABLE="$(command -v "$PY")"
        -DPython3_EXECUTABLE="$(command -v "$PY")")
  else
    echo "WARNING: pybind11 not importable by $PY: no _azure2 module, the pyazr section will fail" >&2
  fi
  # shellcheck disable=SC2086
  cmake -S "$REPO_ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=RelWithDebInfo -DAZURE2_SANITIZE="$SANFLAG" \
    -DBUILD_GUI=OFF -DUSE_QWT=OFF -DUSE_API=ON -DPYAZR_OUTPUT_DIR="$BUILD/pyazr-module" \
    "${PB[@]}" ${SAN_CMAKE_ARGS:-} || die "configure failed"
  cmake --build "$BUILD" --parallel ${SAN_JOBS:-} || die "build failed"
fi
BIN="$BUILD/src/AZURE2"
[ -x "$BIN" ] || die "no $BIN"

# ------------------------------------------------------------------ runtime environment
WORK="$(mktemp -d "${TMPDIR:-/tmp}/sanitize_upgrade.XXXXXX")"
REP="${SAN_REPORT_DIR:-$WORK/reports}"
mkdir -p "$REP"
case "$REP" in /*) ;; *) REP="$PWD/$REP" ;; esac
keep=0
cleanup() { { [ "$keep" = 1 ] && [ -z "${SAN_REPORT_DIR:-}" ]; } || rm -rf "$WORK"; }   # reports outside WORK survive
trap cleanup EXIT
export OMP_NUM_THREADS="$THREADS"
export AZURE2_BIN="$BIN"

OMPLIB="$(ldd "$BIN" 2>/dev/null | awk '/libomp\.so|libgomp\.so|libiomp5\.so/ { print $3; exit }')"
echo "AZURE2: $BIN"
echo "sanitizer: $SANFLAG; threads: $THREADS; OpenMP runtime: ${OMPLIB:-none found}"
if [ "$KIND" = thread ]; then
  # libgomp's barriers are invisible to TSan (raw futexes): every parallel region would be reported.
  case "$OMPLIB" in
    *libgomp*) die "TSan needs LLVM's libomp with libarcher; this build loads libgomp ($OMPLIB). Build with clang." ;;
  esac
  ARCHER="${SAN_ARCHER:-}"
  if [ -z "$ARCHER" ] && [ -n "$OMPLIB" ]; then   # next to libomp, as loaded or behind its symlinks
    for d in "$(dirname "$OMPLIB")" "$(dirname "$(readlink -f "$OMPLIB")")"; do
      [ -f "$d/libarcher.so" ] && { ARCHER="$d/libarcher.so"; break; }
    done
  fi
  [ -n "$ARCHER" ] && [ -f "$ARCHER" ] || die "libarcher.so not found next to $OMPLIB; set SAN_ARCHER"
  export OMP_TOOL_LIBRARIES="$ARCHER"
  echo "archer: $ARCHER"
fi

# san_env STEP: the options of this sanitizer, each report to $REP/STEP/<kind>.<pid>
san_env() {
  local d="$REP/$1"
  mkdir -p "$d"
  if [ "$KIND" = thread ]; then
    export TSAN_OPTIONS="log_path=$d/tsan halt_on_error=0 exitcode=66 second_deadlock_stack=1 history_size=7 suppressions=$SUPP ${SAN_TSAN_EXTRA:-}"
  else
    export ASAN_OPTIONS="log_path=$d/asan detect_leaks=${SAN_LEAKS:-0} halt_on_error=1 abort_on_error=0 detect_stack_use_after_return=1 check_initialization_order=1 suppressions=$SUPP ${SAN_ASAN_EXTRA:-}"
    export LSAN_OPTIONS="suppressions=$REPO_ROOT/scripts/sanitizers/lsan.supp"
  fi
}

if [ "$KIND" = thread ]; then
  # The archer tool must really be loaded, or OpenMP barriers are invisible and reports are false.
  san_env archer-selftest
  mkdir -p "$WORK/selftest/output" "$WORK/selftest/checks"
  cp -R "$REPO_ROOT/tests/upgrade/T3/T3.azr" "$REPO_ROOT/tests/upgrade/T3/data" "$WORK/selftest/"
  (cd "$WORK/selftest" && ARCHER_OPTIONS=verbose=1 "$BIN" calculate --no-gui \
     --no-coulomb-table T3.azr > run.log 2>&1)
  grep -q "Archer detected OpenMP application" "$WORK/selftest/run.log" \
    || { tail -20 "$WORK/selftest/run.log" >&2; die "libarcher did not attach (OMP_TOOL_LIBRARIES=$ARCHER)"; }
  echo "archer: attached"
fi

results=()   # "STATUS name"
record() { results+=("$1 $2"); printf '  %-5s %s\n' "$1" "$2"; }

# ------------------------------------------------------------------ 1. T1, T2, T3
CALC=calculate; FIT=fit; EXTRAP=extrapolate
# run_project STEP PROJECT AZR MENU FLAGS...: one CLI run on a copy of tests/upgrade/PROJECT
run_project() {
  local step="$1" proj="$2" azr="$3" command="$4"; shift 4
  local d="$WORK/$step" src="$REPO_ROOT/tests/upgrade/$proj"
  mkdir -p "$d/output" "$d/checks"
  cp -R "$src"/*.azr "$src/data" "$d/"
  for f in temperatures.dat cov_diff.txt prior_cov.txt; do [ -f "$src/$f" ] && cp "$src/$f" "$d/"; done
  san_env "$step"
  (cd "$d" && "$BIN" $command --no-gui --no-coulomb-table "$@" "$azr" > run.log 2>&1)
  local rc=$?
  cp "$d/run.log" "$REP/$step/run.log"
  if [ "$rc" -eq 0 ] && [ -s "$d/output/chiSquared.out" -o "$command" != "$CALC" ]; then
    record ok "$step"
  else
    record FAIL "$step (exit $rc)"
  fi
}
if want upgrade; then
  echo "== T1, T2, T3"
  run_project T1-calc T1 T1.azr "$CALC"
  run_project T1-extrap T1 T1.azr "$EXTRAP"
  run_project T1-rate T1 T1.azr 'rate --rate-entrance=1 --rate-exit=2 --temperature-file=temperatures.dat'
  run_project T1-fit T1 T1_start.azr "$FIT" --use-gradient
  run_project T3-calc T3 T3.azr "$CALC"
  run_project T3-extrap T3 T3.azr "$EXTRAP"
  run_project T3-rate T3 T3.azr 'rate --rate-entrance=2 --rate-exit=1 --temperature-file=temperatures.dat'
  run_project T3-fit T3 T3_start.azr "$FIT" --use-gradient
  run_project T2-calc T2 T2.azr "$CALC"
  run_project "T2-fit$FITCALLS" T2 T2.azr "$FIT" --use-lm --max-calls="$FITCALLS"
fi

# ------------------------------------------------------------------ 2. behavioural checks
# Point mapping across segments and angles, energy shifts (groups, functions, sqrt, sub-points, penalties),
# target integration (thickness, density, windows, thick targets, straggling, profiles, convolution).
CHECKS="angle_mapping component_mapping component_subpoint_angle shared_point_tables shared_tables_race
shared_effect_pairs upos_pair_index
shift_fn shift_fn_groups shift_groups shift_penalty_zero shift_sqrt_subpoints shift_subpoint_angle energy_shift_sqrt
target_auto_windows target_density target_effect_shared target_thickness_exact target_window_fallback thick_target
convolution_unify convolution_window straggling_bohr natural_target natural_neon"
if want checks; then
  echo "== behavioural checks (tests/<name>/check.sh)"
  for t in $CHECKS; do
    san_env "check-$t"
    if "$REPO_ROOT/tests/$t/check.sh" "$BIN" > "$REP/check-$t/check.log" 2>&1; then record ok "check-$t"
    else record FAIL "check-$t"; fi
  done
fi

# ------------------------------------------------------------------ 3. regression projects without check.sh
REGRESS="beam_profile_kernel beam_spread_ay thick_target_resonance"
if want regress; then
  echo "== regression projects (calculate)"
  for t in $REGRESS; do
    d="$WORK/regress-$t"
    cp -R "$REPO_ROOT/tests/$t" "$d"; rm -rf "$d/output" "$d/checks"; mkdir -p "$d/output" "$d/checks"
    san_env "regress-$t"
    (cd "$d" && "$BIN" $CALC --no-gui "$t.azr" > run.log 2>&1)
    rc=$?
    cp "$d/run.log" "$REP/regress-$t/run.log"
    if [ "$rc" -eq 0 ] && [ -s "$d/output/chiSquared.out" ]; then record ok "regress-$t"
    else record FAIL "regress-$t (exit $rc)"; fi
  done
fi

# ------------------------------------------------------------------ 4. C++ reference tests of the integrators
if want reference; then
  echo "== reference tests (ctest)"
  for t in beam_profile_reference udr_reference adaptive_grid_reference; do
    san_env "reference-$t"
    if ctest --test-dir "$BUILD" -R "^$t\$" --output-on-failure > "$REP/reference-$t/ctest.log" 2>&1 \
       && grep -q "100% tests passed" "$REP/reference-$t/ctest.log"; then record ok "reference-$t"
    else record FAIL "reference-$t"; fi
  done
fi

# ------------------------------------------------------------------ 5. pyazr tests (in-process engine)
# The module is an instrumented shared object loaded by an uninstrumented python, so the sanitizer runtime is
# preloaded.  The tests import pyazr from the directory two above their own file, so they run from a mirror of
# the repository whose pyazr/ holds this build's module (the checkout's pyazr/ is left alone).
PYTESTS="shift_fn sqrt_shift shift_first_call jacobian_shift_columns target_density thickness natural_target
doppler data_covariance composite_terms parallel_points shared_state"
if want pyazr; then
  echo "== pyazr tests"
  MOD="$(ls "$BUILD"/pyazr-module/_azure2*.so 2>/dev/null | head -1)"
  RT="${SAN_RUNTIME:-}"
  if [ -z "$RT" ] && [ -n "$MOD" ]; then
    RT="$(ldd "$MOD" 2>/dev/null | awk -v k="$SHORT" '$1 ~ ("lib" k "\\.so") { print $3; exit }')"   # GCC: DT_NEEDED
  fi
  if [ -z "$RT" ]; then   # clang: the runtime is not linked into shared objects; ask the compiler
    CXXC="$(sed -n 's/^CMAKE_CXX_COMPILER:[A-Z]*=//p' "$BUILD/CMakeCache.txt")"
    for n in "libclang_rt.$SHORT.so" "libclang_rt.$SHORT-x86_64.so" "lib$SHORT.so"; do
      c="$("$CXXC" ${SAN_RT_FLAGS:-} -print-file-name="$n" 2>/dev/null)"
      [ -n "$c" ] && [ "$c" != "$n" ] && [ -f "$c" ] && { RT="$c"; break; }
    done
  fi
  if [ -z "$MOD" ] || [ -z "$RT" ]; then
    record FAIL "pyazr (module: ${MOD:-not built}; runtime: ${RT:-not found})"
  else
    echo "module: $MOD; preload: $RT"
    M="$WORK/mirror"
    mkdir -p "$M/pyazr"
    for f in "$REPO_ROOT"/* "$REPO_ROOT"/.[!.]*; do
      [ -e "$f" ] || continue
      case "$(basename "$f")" in pyazr|build*|.git) ;; *) ln -s "$f" "$M/" ;; esac
    done
    for f in "$REPO_ROOT"/pyazr/*; do
      case "$(basename "$f")" in _azure2*|__pycache__) ;; *) ln -s "$f" "$M/pyazr/" ;; esac
    done
    cp "$MOD" "$M/pyazr/"
    # A test that starts the CLI must not hand it the preloaded runtime: the executable carries its own
    # (clang links it statically), and two runtimes in one process abort.
    mkdir -p "$WORK/bin"
    printf '#!/bin/sh\nunset LD_PRELOAD\nexec "%s" "$@"\n' "$BIN" > "$WORK/bin/AZURE2"
    chmod +x "$WORK/bin/AZURE2"
    for t in $PYTESTS; do
      san_env "pyazr-$t"
      log="$REP/pyazr-$t/test.log"
      # python's own allocator and threads are not instrumented: plain malloc for ASan, and TSan ignores
      # accesses from uninstrumented modules.
      AZURE2_BIN="$WORK/bin/AZURE2" PYTHONMALLOC=malloc PYTHONDONTWRITEBYTECODE=1 LD_PRELOAD="$RT${LD_PRELOAD:+:$LD_PRELOAD}" \
        TSAN_OPTIONS="${TSAN_OPTIONS:-} ignore_noninstrumented_modules=1" \
        "$PY" "$M/tests/pyazr/${t}_test.py" > "$log" 2>&1
      rc=$?
      if [ "$rc" -eq 0 ] && ! grep -q "^skip:" "$log"; then record ok "pyazr-$t"
      elif grep -q "^skip:" "$log"; then record FAIL "pyazr-$t (skipped: $(grep -m1 '^skip:' "$log"))"
      else record FAIL "pyazr-$t (exit $rc)"; fi
    done
  fi
fi

# ------------------------------------------------------------------ summary
reports="$(find "$REP" -type f \( -name 'tsan.*' -o -name 'asan.*' \) -size +0 | sort)"
nrep=0; [ -n "$reports" ] && nrep="$(echo "$reports" | wc -l | tr -d ' ')"
nfail=0
for r in "${results[@]}"; do case "$r" in FAIL*) nfail=$((nfail + 1)) ;; esac; done
echo
echo "======================================"
if [ "$nrep" -gt 0 ]; then
  echo "SANITIZER REPORTS: $nrep file(s); distinct reports (count, summary line):"
  # shellcheck disable=SC2086
  grep -h '^SUMMARY' $reports | sort | uniq -c | sort -rn
  echo "first report of each file (at most 6 files; all are in $REP):"
  for f in $(echo "$reports" | head -6); do
    echo "---- ${f#"$REP"/}"
    sed -n '1,/^SUMMARY/p' "$f" | head -80
  done
fi
echo "steps: ${#results[@]}; failed: $nfail; sanitizer report files: $nrep"
for r in "${results[@]}"; do case "$r" in FAIL*) echo "  $r" ;; esac; done
if [ "$nrep" -gt 0 ] || [ "$nfail" -gt 0 ]; then keep=1; echo "reports and logs: $REP"; fi
[ "$nrep" -gt 0 ] && exit 1
[ "$nfail" -gt 0 ] && [ "${SAN_ALLOW_FAILURES:-0}" != 1 ] && exit 2
exit 0
