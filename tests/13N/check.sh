#!/usr/bin/env bash
#
# A fit of the 13N evaluation starts, and its analytic gradient is right.
# The working copies of a fit (EData::Clone) hold no external-capture
# amplitudes, and EPoint::GetECAmplitudeWithShift used to size its memo by the
# stored ones: every fit of a project with external capture crashed at its
# first evaluation (fixed in 2dcb846; tests/ec_shift runs a shifted fit).
# Here the full model is fitted with --use-gradient and AZURE_GRAD_CHECK=exit,
# which evaluates chi^2 and the gradient at the start (adjoint and
# residual-Jacobian paths) and compares them with central finite differences:
#
#   (a) the run exits normally and its chi^2 is the one of "calculate with data"
#       (expected/chiSquared.out);
#   (b) analytic and Jacobian-path gradients agree with finite differences to
#       1e-5 of the largest component, and with each other to 1e-10.
#
#   ./tests/13N/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# The time limit of tests/lib/guard.sh: timeout, gtimeout or a shell watchdog.
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/13N_fit.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

mkdir -p "$WORK/output" "$WORK/checks"
cp -r "$HERE/data" "$WORK/"
cp "$HERE/13N.azr" "$WORK/t.azr"
(cd "$WORK" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
   run_limited 1200 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
status=$?
if [ "$status" -ne 0 ]; then
  bad "the fit did not start (exit status $status)"; tail -5 "$WORK/run.log"
else
  ok "the fit starts (exit status 0)"
fi

chi="$(awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { sub(",", "", $3); print $3 }' "$WORK/run.log")"
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$HERE/expected/chiSquared.out")"
if [ -n "$chi" ] && awk -v a="$chi" -v b="$want" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(d <= 1e-5) }'; then
  ok "chi^2 at the start $chi = calculation $want"
else
  bad "chi^2 at the start '$chi', calculation $want"
fi

awk '{ sub(/\r$/, "") }
     /^\[grad-check\] idx / { a = f = j = ""
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print a, f, j }' "$WORK/run.log" > "$WORK/rows"
n="$(wc -l < "$WORK/rows" | tr -d ' ')"
if [ "$n" -lt 5 ]; then
  bad "gradient check listed $n parameters"
else
  read -r fd jac <<< "$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a
                                d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d
                                e = $1 - $3; if (e < 0) e = -e; if (e > v) v = e }
                              END { printf "%.2e %.2e", w / m, v / m }' "$WORK/rows")"
  awk -v w="$fd" 'BEGIN { exit !(w + 0 <= 1e-5) }' && ok "gradient = finite differences for $n parameters ($fd of the largest)" \
    || bad "gradient vs finite differences: $fd of the largest"
  awk -v w="$jac" 'BEGIN { exit !(w + 0 <= 1e-10) }' && ok "Jacobian-path gradient = adjoint gradient ($jac of the largest)" \
    || bad "Jacobian-path vs adjoint gradient: $jac of the largest"
fi

if [ "$failures" -eq 0 ]; then echo "13N: all checks passed"; else echo "13N: $failures check(s) FAILED"; fi
[ "$failures" -eq 0 ]
