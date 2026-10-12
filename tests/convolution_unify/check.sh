#!/usr/bin/env bash
#
# One Gaussian beam spread (2026-10): the constant-sigma field, an expression
# sigma_lab(E_lab) (keyword sigmalab) and a legacy c.m. expression must give the
# same result for the same physical width, alone and under target integration;
# plain target integration is unchanged; and the normalised kernel moves the
# constant-sigma result by exactly 1/erf(3/sqrt 2) relative to the pre-change
# binary (reference_old/, recorded 2026-10-11 with the last unnormalised build).
#
#   ./tests/convolution_unify/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# macOS has no timeout(1); run without a limit there
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/conv_unify.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {  # NAME -> $WORK/NAME/output/AZUREOut_aa=1_R=1.extrap
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$1.azr" "$dir/t.azr"
  (cd "$dir" && printf '3\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -s "$dir/output/AZUREOut_aa=1_R=1.extrap" ]; then bad "$1: no extrapolation written"; tail -8 "$dir/run.log"; return 1; fi
}
# maxrel A B -> largest relative difference of the cross-section column between two .extrap files
maxrel() {
  paste "$1" "$2" | awk 'NF >= 8 { n = NF / 2; a = $4; b = $(n + 4); d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m < 1e-300) m = 1e-300;
                                    r = d / m; if (r > worst) worst = r; rows++ }
                        END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
# ratiospread A B -> min and max of A/B over the rows, and the median
ratiostat() {
  paste "$1" "$2" | awk 'NF >= 8 { n = NF / 2; if ($(n + 4) > 0) print $4 / $(n + 4) }' | sort -g | awk '{ v[NR] = $1 } END { printf "%.7f %.7f %.7f\n", v[1], v[int((NR + 1) / 2)], v[NR] }'
}
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t) }'; }
within() { awk -v x="$1" -v c="$2" -v t="$3" 'BEGIN { d = x - c; if (d < 0) d = -d; exit !(d <= t) }'; }

for m in bare fixed legacy_cm_const legacy_cm_E lab_const lab_E tgt tgt_fixed tgt_legacy_const tgt_lab_const tgt_straggle; do run "$m" || exit 1; done
out() { echo "$WORK/$1/output/AZUREOut_aa=1_R=1.extrap"; }

# 1. the three ways of writing the same beam spread agree
r="$(maxrel "$(out fixed)" "$(out lab_const)")";            if le "$r" 1e-12; then ok "constant sigma == expression with sigmalab ($r)"; else bad "constant sigma vs sigmalab expression differ ($r)"; fi
r="$(maxrel "$(out fixed)" "$(out legacy_cm_const)")";      if le "$r" 1e-9;  then ok "constant sigma == legacy c.m. expression of the same width ($r)"; else bad "constant sigma vs legacy c.m. expression differ ($r)"; fi
r="$(maxrel "$(out lab_E)" "$(out legacy_cm_E)")";          if le "$r" 1e-9;  then ok "sigma_lab(E_lab) == its c.m. transcription ($r)"; else bad "lab and c.m. forms of sigma(E) differ ($r)"; fi
r="$(maxrel "$(out tgt_fixed)" "$(out tgt_lab_const)")";    if le "$r" 1e-12; then ok "under target integration: constant == sigmalab ($r)"; else bad "target + constant vs target + sigmalab differ ($r)"; fi
r="$(maxrel "$(out tgt_fixed)" "$(out tgt_legacy_const)")"; if le "$r" 1e-9;  then ok "under target integration: constant == legacy c.m. ($r)"; else bad "target + constant vs target + legacy c.m. differ ($r)"; fi

# 2. what did not change. The references were recorded on Linux; another platform's libm differs
#    from it at ~1e-11, so these are held to 1e-8 -- five orders below the smallest change in part 3.
r="$(maxrel "$(out bare)" "$HERE/reference_old/bare.extrap")"; if le "$r" 1e-8; then ok "no effect: unchanged ($r)"; else bad "bare calculation moved ($r)"; fi
r="$(maxrel "$(out tgt)"  "$HERE/reference_old/tgt.extrap")";  if le "$r" 1e-8; then ok "plain target integration: unchanged ($r)"; else bad "plain target integration moved ($r)"; fi

# 3. what changed, by how much, and why
#    Constant sigma: the old kernel was not renormalised over its +-3 sigma window, so every point was
#    low by erf(3/sqrt 2) = 0.9973002; the new result is the old one divided by that (as the 2-point
#    quadrature integrates the kernel, within 1e-5) and by the same factor at every energy.
read lo med hi <<< "$(ratiostat "$(out fixed)" "$HERE/reference_old/fixed.extrap")"
if within "$med" 1.0027069 1e-5 && within "$lo" "$med" 2e-6 && within "$hi" "$med" 2e-6; then ok "constant sigma: new/old = 1/erf(3/sqrt2) = 1.0027069 everywhere (min $lo, median $med, max $hi)"; else bad "constant sigma: new/old not the kernel normalisation (min $lo, median $med, max $hi; expected 1.0027069)"; fi
#    Legacy c.m. expression: already normalised before; the width is now taken at the point's energy
#    instead of at each sub-point, a change of order sigma'(E) * 3 sigma -- small for E^1.5.
r="$(maxrel "$(out legacy_cm_E)" "$HERE/reference_old/legacy_cm_E.extrap")"; if le "$r" 2e-2; then ok "legacy c.m. sigma(E): width now at the point energy, moved by at most $r"; else bad "legacy c.m. sigma(E) moved by $r (expected < 2e-2)"; fi
r="$(maxrel "$(out tgt_fixed)" "$HERE/reference_old/tgt_fixed.extrap")";       if le "$r" 2e-2; then ok "target + beam spread: inner kernel now normalised, moved by at most $r"; else bad "target + beam spread moved by $r (expected < 2e-2)"; fi
r="$(maxrel "$(out tgt_straggle)" "$HERE/reference_old/tgt_straggle.extrap")"; if le "$r" 2e-2; then ok "target + straggling: inner kernel now normalised, moved by at most $r"; else bad "target + straggling moved by $r (expected < 2e-2)"; fi

if [ "$failures" -eq 0 ]; then echo "PASS: one Gaussian beam spread"; else echo "FAIL: $failures check(s)"; exit 1; fi
