#!/usr/bin/env bash
#
# Electron screening benchmark on real data: 6Li(p,a)3He from Cruz et al. 2005 (EXFOR O1282005-007) on Li metal,
# a Li2WO4 insulator and a PdLi alloy, against J. Skowronski's 7Be R-matrix evaluation held fixed (models and data
# from make_models.py; README.md has the source and the numbers).
#
#   1. screened / bare = exp(pi eta U_e / E) at every point of the three segments;
#   2. at the benchmark minimum the chi2 of each segment is the benchmark's;
#   3. MIGRAD and --use-gradient, with only U_e and the norms free, find that minimum from U_e = 500 eV:
#      U_e = 1211, 323, 3433 eV within 5 % of their uncertainties (published 1320 +- 70, 320 +- 110, 3760 +- 260).
#
#   ./tests/screening_6li/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# The time limit of tests/lib/guard.sh: timeout, gtimeout or a shell watchdog.
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/screening_6li.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME AZR MODE FLAGS -> $WORK/NAME (mode 1 calculate, 2 fit)
run() {
  local name="$1" azr="$2" mode="$3" flags="$4"
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$HERE/data" "$dir/data"
  cp "$azr" "$dir/t.azr"
  local input='1\n\n\n7\n'
  [ "$mode" = 2 ] && input='2\nn\n\n\n'
  (cd "$dir" && printf "$input" | run_limited 900 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}
out() { echo "$WORK/$1/output/AZUREOut_aa=2_R=1.out"; }
# Output files may end their lines with CR (Windows builds): strip it before reading fields.
chi() { awk -F, -v s="$2" '{ sub(/\r$/, "") } $1 == s { print $2; exit }' "$WORK/$1/output/chiSquared.out"; }
par() { awk -v n="$2" '{ sub(/\r$/, "") } $1 == n { print $2; exit }' "$WORK/$1/output/param.sav"; }
# |a/b - 1| <= tol
close() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a / b - 1; if (d < 0) d = -d; exit !(a != "" && d <= t) }'; }
# |a - b| <= t
near() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(a != "" && d <= t) }'; }

UE=(1211.1217 322.6758 3432.9061)          # make_models.py BEST_UE (metal, insulator, PdLi)
CHI=(12.491650 2.952962 8.261932)          # the benchmark's chi2 at that minimum
DUE=(65 98 231)                            # its Delta chi2 = 1 half-widths (eV); the fits must land within 5 % of them

# ---- 1. the factor, point by point -------------------------------------------------------
run bare "$HERE/bare.azr" 1 ""
run scr "$HERE/scr.azr" 1 ""
if [ -s "$(out bare)" ] && [ -s "$(out scr)" ]; then
  # Rows in segment order: metal (9), insulator (7), PdLi (7).  exp(pi eta Ue/E), eta of p + 6Li with AZURE2's
  # constants (Constants.h) and the pair's masses (1.00782, 6.0151).
  read -r worst nrow <<< "$(paste "$(out bare)" "$(out scr)" | awk -v u1="${UE[0]}" -v u2="${UE[1]}" -v u3="${UE[2]}" '
    BEGIN { pi = 3.141592653589793; uconv = 931.4940880; fs = 1 / 137.0359996790; mu = 1.00782 * 6.0151 / 7.02292 }
    { gsub(/\r/, "") }
    NF >= 8 { r++; m = NF / 2; E = $1; b = $4; s = $(m + 4)
      u = (r <= 9) ? u1 : (r <= 16) ? u2 : u3
      eta = sqrt(uconv / 2) * fs * 3 * sqrt(mu / E)
      d = s / b / exp(pi * eta * u * 1e-6 / E) - 1; if (d < 0) d = -d; if (d > w) w = d }
    END { printf "%.3e %d\n", w, r }')"
  if [ "$nrow" = 23 ] && awk -v w="$worst" 'BEGIN { exit !(w < 1e-8) }'; then
    ok "screened / bare = exp(pi eta U_e / E) at all 23 points (worst $worst)"
  else
    bad "screened / bare (worst $worst over $nrow rows)"
  fi
else
  bad "bare or screened run produced no output"; tail -5 "$WORK/bare/run.log" "$WORK/scr/run.log"
fi

# ---- 2. chi2 at the benchmark minimum -----------------------------------------------------
for s in 1 2 3; do
  c="$(chi scr "$s")"
  if close "$c" "${CHI[$((s - 1))]}" 1e-4; then ok "segment $s chi2 $c (benchmark ${CHI[$((s - 1))]})"
  else bad "segment $s chi2 '$c', benchmark ${CHI[$((s - 1))]}"; fi
done

# ---- 3. the fits find it ------------------------------------------------------------------
for flags in "" "--use-gradient"; do
  name="fit${flags:+_grad}"
  run "$name" "$HERE/fit.azr" 2 "$flags"
  got="$(for i in 1 2 3; do printf '%s ' "$(par "$name" "screening_$i")"; done)"
  read -r a b c <<< "$got"
  if near "$a" "${UE[0]}" "$(awk -v s="${DUE[0]}" 'BEGIN { print 0.05 * s }')" &&
     near "$b" "${UE[1]}" "$(awk -v s="${DUE[1]}" 'BEGIN { print 0.05 * s }')" &&
     near "$c" "${UE[2]}" "$(awk -v s="${DUE[2]}" 'BEGIN { print 0.05 * s }')"; then
    ok "fit ${flags:-(MIGRAD)}: U_e = $a, $b, $c eV (benchmark ${UE[*]})"
  else
    bad "fit ${flags:-(MIGRAD)}: U_e = '$got', benchmark ${UE[*]}"; tail -5 "$WORK/$name/run.log"
  fi
done

if [ "$failures" -eq 0 ]; then echo "PASS screening_6li"; exit 0; fi
echo "FAIL screening_6li ($failures)"; exit 1
