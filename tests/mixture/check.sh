#!/usr/bin/env bash
#
# Mixture segments (natural-abundance / compound targets, 2026-10): a segment's
# model is w_self * (this project's value) + sum_i w_i * fold(T_i), the tables
# folded with the segment's own resolution function.  Checked on p+12C
# elastic with the models of make_models.py:
#
#   (a) a "mixture 1.0 0" block changes nothing (bit-identical)
#   (b) half the model plus half a table OF THE SAME MODEL (made here from its
#       bare extrapolation on a 0.5 keV grid) reproduces the plain folded
#       model, under a Gaussian beam spread and under a time-of-flight (UDR)
#       kernel: the table is folded on the same grid by the same kernel, so the
#       only difference is the table's linear interpolation
#   (c) linearity: two tables of 0.25, and a doubled table at 0.25, equal (b)
#   (d) the same as a DATA segment: chi2 against the plain folded values
#       vanishes, the .out fit column equals (b), the decomposition file adds up
#
#   ./tests/mixture/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/mixture.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MENU -> $WORK/NAME/output/AZUREOut_aa=1_R=1.{extrap|out}
run() {
  local dir="$WORK/$1" menu="${2:-3}"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$1.azr" "$dir/t.azr"
  for f in table.dat table2.dat data.dat; do [ -f "$WORK/$f" ] && cp "$WORK/$f" "$dir/"; done
  cp "$HERE/models/udr.txt" "$dir/"
  (cd "$dir" && printf '%s\n\n\n7\n' "$menu" | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  local want="$dir/output/AZUREOut_aa=1_R=1.extrap"; [ "$menu" = 1 ] && want="$dir/output/AZUREOut_aa=1_R=1.out"
  if [ ! -s "$want" ]; then bad "$1: no output written"; tail -8 "$dir/run.log"; return 1; fi
}
out() { echo "$WORK/$1/output/AZUREOut_aa=1_R=1.extrap"; }
# maxrel A B [colA colB] -> largest relative difference of the fit column (column 4 of either file by
# default; .out and .extrap files have different widths, so each column is extracted first)
maxrel() {
  local ca="${3:-4}" cb="${4:-4}"
  paste <(awk -v c="$ca" '$1 !~ /^#/ && NF >= c { print $c }' "$1") <(awk -v c="$cb" '$1 !~ /^#/ && NF >= c { print $c }' "$2") |
    awk 'NF == 2 { d = $1 - $2; if (d < 0) d = -d; m = $2 < 0 ? -$2 : $2; if (m < 1e-300) m = 1e-300;
                   r = d / m; if (r > worst) worst = r; rows++ }
         END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t) }'; }

# -- the tables: this model's own bare cross section, 0.2-1.8 MeV lab in 0.1 keV steps, E_lab = E_cm*(M1+M2)/M2
run bare_wide || exit 1
F="$(awk 'BEGIN { printf "%.12f", (1.00782 + 12.0) / 12.0 }')"
awk -v F="$F" 'NF >= 4 { printf "%.10e %.10e\n", $1 * F, $4 }' "$(out bare_wide)" > "$WORK/table.dat"
awk 'NF == 2 { printf "%s %.16e\n", $1, 2.0 * $2 }' "$WORK/table.dat" > "$WORK/table2.dat"
[ "$(wc -l < "$WORK/table.dat")" -gt 15000 ] && ok "table of the bare model written ($(wc -l < "$WORK/table.dat") rows)" || bad "table too short"

# gauss is the reference folded model (the convolution check's `fixed` effect on this grid)
for m in gauss identity mix_gauss mix_two mix_scaled plain_udr mix_udr; do run "$m" || exit 1; done

# (a)
# (cksum, not cmp: the Windows CI runner's MSYS2 has coreutils but no diffutils)
if [ "$(cksum < "$(out identity)")" = "$(cksum < "$(out gauss)")" ]; then ok "(a) mixture 1.0 0: bit-identical to the plain model"; else bad "(a) mixture 1.0 0 changed the result"; fi
grep -q "Segment #1: mixture, 1 x this model" "$WORK/identity/run.log" && ok "(a) the block is reported" || bad "(a) no mixture report in the log"

# (b) self-consistency: the table is folded on the same grid by the same kernel; the only difference is the
#     table's linear interpolation across the ~30 keV resonance on a 0.1 keV grid, ~(h/Gamma)^2 x curvature ~ 1e-4
r="$(maxrel "$(out mix_gauss)" "$(out gauss)")"; if le "$r" 3e-4; then ok "(b1) Gaussian: 0.5 model + 0.5 table == plain folded model ($r)"; else bad "(b1) Gaussian mixture differs from the folded model ($r)"; fi
r="$(maxrel "$(out mix_udr)" "$(out plain_udr)")";  if le "$r" 3e-4; then ok "(b2) UDR kernel: 0.5 model + 0.5 table == plain folded model ($r)"; else bad "(b2) UDR mixture differs from the folded model ($r)"; fi
#     and the kernels do something here: the 30 keV Gaussian and the ~20 keV (at 0.46 MeV) time-of-flight
#     kernel each move the resonance by percent and differ from each other by percent
r="$(maxrel "$(out plain_udr)" "$(out gauss)")"; le "$r" 1e-2 && bad "(b2) the UDR model is indistinguishable from the Gaussian one ($r): the kernel did not apply" || ok "(b2) the UDR kernel differs from the Gaussian by $r, so (b2) is an independent check"

# (c) linearity
r="$(maxrel "$(out mix_two)" "$(out mix_gauss)")";    if le "$r" 1e-12; then ok "(c) two tables of 0.25 == one of 0.5 ($r)"; else bad "(c) two tables of 0.25 differ from one of 0.5 ($r)"; fi
r="$(maxrel "$(out mix_scaled)" "$(out mix_gauss)")"; if le "$r" 1e-12; then ok "(c) a doubled table at 0.25 == the table at 0.5 ($r)"; else bad "(c) doubled table at 0.25 differs ($r)"; fi
dec="$WORK/mix_two/output/AZUREOut_mixture_seg=1.extrap"
if [ -s "$dec" ]; then
  r="$(awk '$1 !~ /^#/ && NF >= 5 { s = $3 + $4 + $5; d = $2 - s; if (d < 0) d = -d; if (d / $2 > w) w = d / $2 } END { printf "%.3e", w }' "$dec")"
  if le "$r" 1e-8; then ok "(c) decomposition file: total == w_self*R-matrix + table terms, to its 10-digit printing ($r)"; else bad "(c) decomposition columns do not add up ($r)"; fi
else bad "(c) no decomposition file written"; fi

# (d) as a data segment: data = the plain folded values (lab energy, c.m. angle, c.m. cross section, 1 %)
awk -v F="$F" 'NF >= 4 { printf "%.10e %.10e %.10e %.10e\n", $1 * F, $3, $4, 0.01 * $4 }' "$(out gauss)" > "$WORK/data.dat"
run mix_data 1 || exit 1
chi="$(awk '/Total-Chi-Squared:/ { for (i = 1; i <= NF; i++) if ($i == "Total-Chi-Squared:") print $(i + 1) }' "$WORK/mix_data/output/chiSquared.out" | head -1)"
# 426 points with 1 % errors and the table's <= 3e-4 interpolation error at the resonance: chi2 < N (3e-4/1e-2)^2 ~ 0.4
if le "$chi" 1; then ok "(d) data segment: chi2 of the mixture against the folded model = $chi (interpolation only)"; else bad "(d) data segment chi2 = $chi (expected < 1)"; fi
r="$(maxrel "$WORK/mix_data/output/AZUREOut_aa=1_R=1.out" "$(out mix_gauss)")"; if le "$r" 1e-8; then ok "(d) .out fit column == the extrapolation's ($r)"; else bad "(d) data-mode fit column differs from the extrapolation ($r)"; fi
[ -s "$WORK/mix_data/output/AZUREOut_mixture_seg=1.out" ] && ok "(d) decomposition file written in data mode" || bad "(d) no decomposition .out file"

# a table that does not cover the grid is refused with the required range
mkdir -p "$WORK/short/output" "$WORK/short/checks"; cp "$HERE/models/mix_gauss.azr" "$WORK/short/t.azr"
awk '$1 < 1.0' "$WORK/table.dat" > "$WORK/short/table.dat"
(cd "$WORK/short" && printf '3\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
if grep -q "ERROR: mixture: segment #1 needs lab energies" "$WORK/short/run.log" && [ ! -s "$WORK/short/output/AZUREOut_aa=1_R=1.extrap" ]; then ok "a table ending before the grid is refused, naming the range"; else bad "a short table was not refused"; tail -3 "$WORK/short/run.log"; fi

if [ "$failures" -eq 0 ]; then echo "PASS: mixture segments"; else echo "FAIL: $failures check(s)"; exit 1; fi
