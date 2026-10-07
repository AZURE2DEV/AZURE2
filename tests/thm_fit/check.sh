#!/usr/bin/env bash
#
# A fit (menu 2, MIGRAD) through every THM experiment feature that changes
# the model: an experiment alone (profiled norm), distortion=coulomb (R(E)),
# vertexModel=dw, lineshape=on and a ps window.  The calculation of each is
# checked elsewhere (tests/thm_*); this checks that a fit runs through them.
#
# Model: tests/18O_p_a_thm without its folding (the <targetInt> line removed),
# the energy of the 8.6026 MeV level free and started 10 keV high.  For each
# variant:
#   - the fit ends and writes output/param.sav;
#   - it moves the level energy and does not end above the chi2 of the start;
#   - the chi2 it writes (output/chiSquared.out) is the chi2 of a
#     calculation from its param.sav (rel 1e-6): what the fit minimizes is
#     what a calculation of its result gives.
#
#   ./tests/thm_fit/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../18O_p_a_thm"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_fit.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-600}")"
KIN="beam=18O target=3He spectator=d Ebeam=115"
fail=0

# setup NAME BLOCK: a copy of the project, unfolded, the first level's energy
# free and 10 keV high, with BLOCK as its <thm> block.
setup() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  awk '/^<targetInt>/ { print; skip = 1; next } /^<\/targetInt>/ { skip = 0 } skip { next }
       $3 == "8.602600" { $3 = "8.612600"; $4 = 0 } { print }' "$SRC/18O_p_a_thm.azr" > "$d/run.azr"
  printf '\n<thm>\n%s\n</thm>\n' "$2" >> "$d/run.azr"
}
# total DIR: the total chi2 of DIR/output/chiSquared.out.
total() { tr -d '\r' < "$1/output/chiSquared.out" 2>/dev/null | grep -oE 'Total-Chi-Squared: [0-9.eE+-]+' | awk '{ print $2 }'; }
energy() { tr -d '\r' < "$1" | awk '$1 == "energy_1" { print $2 }'; }

for v in experiment distortion dw lineshape ps; do
  case $v in
    experiment) block="experiment[A] segments=1,2" ;;
    distortion) block="experiment[A] segments=1,2 $KIN distortion=coulomb" ;;
    dw) block="experiment[A] segments=1,2 $KIN distortion=coulomb vertexModel=dw" ;;
    lineshape) block="experiment[A] segments=1,2 $KIN lineshape=on" ;;
    ps) block="experiment[A] segments=1,2 $KIN ps=hulthen:0-40 psNodes=8" ;;
  esac
  echo "($v) $block"
  setup "$v.start" "$block"
  (cd "$WORK/$v.start" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
  c0="$(total "$WORK/$v.start")"
  setup "$v" "$block"
  (cd "$WORK/$v" && printf '2\nn\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
  if [ ! -f "$WORK/$v/output/param.sav" ] || [ -z "$c0" ]; then
    bad "$v: the fit did not finish (start chi2 '$c0')"; tail -5 "$WORK/$v/log" | sed 's/^/        /'
    continue
  fi
  c1="$(total "$WORK/$v")"
  e1="$(energy "$WORK/$v/output/param.sav")"
  awk -v a="$c0" -v b="$c1" -v e="$e1" 'BEGIN { d = e - 8.6126; if (d < 0) d = -d
      exit !(b != "" && b <= a * (1 + 1e-12) && d > 1e-6) }' \
    && ok "$v: E moved 8.6126 -> $e1 MeV, chi2 $c0 -> $c1" || bad "$v: chi2 $c0 -> $c1, E $e1"
  # A calculation from the fit's parameters gives the fit's chi2.
  setup "$v.check" "$block"
  cp "$WORK/$v/output/param.sav" "$WORK/$v.check/fit.par"
  (cd "$WORK/$v.check" && printf '1\nfit.par\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
  c2="$(total "$WORK/$v.check")"
  awk -v a="$c1" -v b="$c2" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(b != "" && d < 1e-6) }' \
    && ok "$v: a calculation from param.sav gives it ($c2)" || bad "$v: fit chi2 $c1, calculation from param.sav $c2"
done

echo
if [ "$fail" -eq 0 ]; then echo "PASS: fits through the THM experiment features"; else echo "FAIL"; exit 1; fi
