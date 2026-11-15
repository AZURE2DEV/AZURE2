#!/usr/bin/env bash
#
# Capture cascades: the secondary gamma ray of a capture final state after an UNOBSERVED primary capture photon
# (a UPOS segment on a capture pair).  The extension of Brune and deBoer, PRC 102, 024628 (2020) named in that paper's
# conclusions: an unobserved photon summed over direction and helicity leaves no interference between multipoles and
# orients the final state like a spin-0 particle with l' = L.  Models and expectations: make_models.py.
#
#   (a) resonant capture (13C + alpha): one level and one multipole; two primary multipoles of one level (incoherent)
#       with a mixed secondary transition; two levels of different J through one multipole, both signs;
#   (b) external capture alone (p + a 1/2- target): into an s state (isotropic only as a coherent sum over J) and a
#       d -> p single-particle transition;
#   (c) a resonance and external capture together, at three energies across the resonance;
#   (d) an angle-integrated UPOS segment on a capture pair = the capture cross section;
#   (e) analytic and Jacobian-path gradients = finite differences (AZURE_GRAD_CHECK), free level and ANC;
#   (f) "linewindow" on a cascade segment is refused;
#   (g) with a second binary: the non-UPOS segments of the same models are byte-identical between the two builds.
#
# Expectations are magnetic-substate sums that use neither the UPOS coefficient nor R_k; for (a) and (b) the PRIMARY
# distribution of the same amplitudes is checked too, which ties the conventions to code that existed before.
#
#   ./tests/capture_cascade/check.sh path/to/AZURE2 [path/to/reference/AZURE2]

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2 [reference AZURE2]}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
REF_BIN="${2:-}"
[ -n "$REF_BIN" ] && REF_BIN="$(cd "$(dirname "$REF_BIN")" && pwd)/$(basename "$REF_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PYTHON="$(command -v python3 || command -v python || true)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/capture_cascade.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run BIN DIR MODEL MODE [FLAGS] [ENV...] -> $WORK/DIR
run() {
  local bin="$1" dir="$WORK/$2" model="$3" mode="$4" flags="${5:-}"
  shift 5 2> /dev/null || shift $#
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$model.azr" "$dir/t.azr"
  cp -r "$HERE/models/data" "$dir/"
  # shellcheck disable=SC2086
  (cd "$dir" && for kv in "$@"; do export "$kv"; done
   printf "$mode\n\n\n7\n" | run_limited 300 "$bin" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}

MODELS="single_e1 multipoles twoj_plus twoj_minus ec_s ec_p ecres_both"
for m in $MODELS; do
  run "$AZURE2_BIN" "$m" "$m" 1 ""
  ls "$WORK/$m/output/"AZUREOut_* > /dev/null 2>&1 || { bad "$m: no output"; tail -15 "$WORK/$m/run.log"; }
done

# ---- (a)-(d) ---------------------------------------------------------------------------------------------------
if [ -n "$PYTHON" ] && "$PYTHON" -c "import numpy" > /dev/null 2>&1; then
  "$PYTHON" "$HERE/verify.py" "$WORK" || failures=$((failures + $?))
else
  echo "  skip  (a)-(d): no python3 with numpy"
fi

# ---- (e) gradients ---------------------------------------------------------------------------------------------
run "$AZURE2_BIN" grad grad 2 "--use-gradient" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_H=1e-7
awk '/^\[grad-check\] idx / { idx = $3; a = f = ""; j = "-"
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print idx, a, f, j }' "$WORK/grad/run.log" | sort -n > "$WORK/grad.rows"
nrows="$(wc -l < "$WORK/grad.rows" | tr -d ' ')"
if [ "$nrows" -lt 4 ]; then
  bad "(e) gradient check listed $nrows parameters"; grep -i "grad\|error" "$WORK/grad/run.log" | head -5
else
  worst="$(awk '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a; g[NR] = $2; f[NR] = $3; j[NR] = $4 }
             END { for (i = 1; i <= NR; i++) { d = g[i] - f[i]; if (d < 0) d = -d; e = j[i] - f[i]; if (e < 0) e = -e
                                               if (d > w) w = d; if (j[i] != "-" && e > w) w = e }
                   printf "%.2e", w / m }' "$WORK/grad.rows")"
  if awk -v w="$worst" 'BEGIN { exit !(w + 0 <= 1e-5) }'; then
    ok "(e) adjoint and Jacobian-path gradients = finite differences, $nrows parameters, worst $worst of the largest"
  else
    bad "(e) gradient vs finite differences, worst $worst of the largest (> 1e-5)"
  fi
fi

# ---- (f) refused combination -----------------------------------------------------------------------------------
run "$AZURE2_BIN" bad_linewindow bad_linewindow 1 ""
if grep -q "ERROR.*UPOS on a capture pair.*linewindow" "$WORK/bad_linewindow/run.log" &&
   ! ls "$WORK/bad_linewindow/output/"AZUREOut_* > /dev/null 2>&1; then
  ok "(f) a cascade segment with \"linewindow\" is refused"
else
  bad "(f) a cascade segment with \"linewindow\" was not refused"; tail -5 "$WORK/bad_linewindow/run.log"
fi

# ---- (g) non-UPOS segments byte-identical to a reference build ---------------------------------------------------
# rows 14-27 of a one-energy model: the angle-integrated cross section and the primary distribution
if [ -n "$REF_BIN" ]; then
  for m in single_e1 twoj_minus ec_p; do
    run "$REF_BIN" "ref_$m" "$m" 1 ""
    f="$(basename "$(ls "$WORK/$m/output/"AZUREOut_* | head -1)")"
    new="$(awk 'NF { n++; if (n >= 14 && n <= 27) print }' "$WORK/$m/output/$f" | tr -d '\r')"
    ref="$(awk 'NF { n++; if (n >= 14 && n <= 27) print }' "$WORK/ref_$m/output/$f" | tr -d '\r')"
    if [ -n "$new" ] && [ "$new" = "$ref" ]; then  # no cmp on the Windows runners
      ok "(g) $m: angle-integrated cross section and primary distribution byte-identical to the reference build"
    else
      bad "(g) $m: non-UPOS segments differ from the reference build"
    fi
  done
fi

echo
[ "$failures" -eq 0 ] && echo "capture_cascade: PASS" || echo "capture_cascade: $failures FAILURE(S)"
exit "$failures"
