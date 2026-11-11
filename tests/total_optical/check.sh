#!/usr/bin/env bash
#
# Quantity sigma_total: the total cross section of a neutral entrance pair by the optical theorem,
# sigma_tot = (pi/k^2) g sum_J (2J+1) sum_c 2 Re T_cc (GenMatrixFunc::CalculateBareCrossSection), one T-matrix per
# point instead of one per exit pair of a sum; on a target with isotopes the abundance-weighted sum of its nuclei's
# totals, one elastic term per nucleus.  Built at run time from tests/natural_neon (format 3: n+20Ne, n+22Ne, n+21Ne;
# three compound nuclei, natural-abundance sums, beam spreads) with a free 1 keV width given to every open reaction
# channel (alpha+17O, n'+20Ne, n'+22Ne, n'+21Ne), which are all zero in the fitted model, and the extrapolations on
# 0.03-2.5 MeV in steps of 10 keV (as the old lines' grids: the same doubles):
#
#   sum    each total as the sum of every exit pair of its compound nucleus (quantity sigma, the suite's own records);
#   total  each total as quantity sigma_total (the natural sums: one elastic term per nucleus on the natural target);
#   elast  the elastic terms alone (quantity sigma).
#
#   (a) every cross section of `total` equals that of `sum` to 1e-9, data (mode 1) and extrapolation (mode 3), and
#       `elast` differs (the reaction channels matter);
#   (b) on data set 1 (Vaughn, 0-2.33 MeV, the natural sum of the three totals, inelastic channels open), the
#       reaction widths free: a fit with the analytic gradient starts from the chi2 of the calculation, and the
#       gradient agrees with finite differences and with the Jacobian path (AZURE_GRAD_CHECK);
#   (c) refused, with an ERROR, exit 1 and no output: sigma_total of a charged pair, of a reaction that is not elastic,
#       with an angle column (a differential data set), a sum on a target without isotopes, two terms on one nucleus,
#       a factor on a term.
#
#   ./tests/total_optical/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$HERE/../natural_neon"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# macOS has no timeout(1); run without a limit there
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/total_optical.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# 1 keV (alternating sign, free) on every reaction-channel row of a level above that channel's threshold
awk '/<levels>/ { f = 1; print; next } /<\/levels>/ { f = 0 }
     f && NF > 30 && ($6 == 1 || $6 == 3 || $6 == 5 || $6 == 7) && $3 > $23 + $26 + 0.05 { $11 = 0; $12 = (k++ % 2 ? -1000 : 1000) }
     { print }' "$REF/natural_neon.azr" > "$WORK/base.azr"
# models MODE: sum (the reactions as they are), total (sigma_total of the elastic terms), elast (sigma of them); every
# extrapolation's energies 0.03, 0.04, ... <= 2.5 MeV lab, accumulated as the old grid lines were
models() {
  python3 - "$WORK/base.azr" "$1" > "$WORK/$1.azr" <<'EOF'
import re, sys
mode = sys.argv[2]
energies = []
e = 0.03
while e <= 2.5:
    energies.append(repr(e))
    e += 0.01
out, kind = [], None
for line in open(sys.argv[1]).read().split("\n"):
    head = re.match(r"^(dataset|extrapolation) ", line)
    if head or not line.startswith("    "):
        kind = head.group(1) if head else (None if line.strip() == "" or line.startswith("<") else kind)
    words = line.split()
    if kind and len(words) > 1 and words[0] == "reaction" and mode != "sum":
        terms = [t.strip() for t in line.split(None, 1)[1].split(" + ")]
        elastic = [t for t in terms if re.fullmatch(r"(\w+)\(n,n\)\1", t)]
        line = "    reaction       " + " + ".join(elastic)
    elif kind and words[:2] == ["quantity", "sigma"] and mode == "total":
        line = "    quantity       sigma_total"
    elif kind == "extrapolation" and words[:1] == ["energies"]:
        line = "    energies       list " + " ".join(energies) + " MeV lab"
    out.append(line)
print("\n".join(out), end="")
EOF
}
for m in sum total elast; do models $m; done
[ "$(grep -c '^    quantity       sigma_total$' "$WORK/total.azr")" = 11 ] || { bad "models not built ($(grep -c 'sigma_total' "$WORK/total.azr") sigma_total records)"; exit 1; }
# The extrapolation copies: no apparatus on the extrapolations (the beams of x3, x4 would fold the grids, see the
# azure2-calculate skill)
for m in sum total elast; do
  awk '/^[^ ]/ { x = ($1 == "extrapolation") } x && /^    (beam|quadrature) / { next } { print }' \
    "$WORK/$m.azr" > "$WORK/x$m.azr"
done

# run NAME MODEL MODE -> $WORK/NAME
run() {
  local dir="$WORK/$1" input
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$REF/data" "$dir/data"
  cp "$WORK/$2" "$dir/t.azr"
  if [ "$3" = 3 ]; then input='3\n\n\n7\n'; else input='1\n\n\n7\n'; fi
  (cd "$dir" && printf "$input" | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}
# max relative difference of the calculated cross section (4th column) over every AZUREOut file of two runs
maxrel() {
  local a="$1" b="$2" f n=0 w=0 r
  for f in "$a"/output/AZUREOut_aa=*."$3"; do
    [ -e "$f" ] || continue
    [ -e "$b/output/$(basename "$f")" ] || { echo "missing"; return; }
    r="$(awk 'NR == FNR { if (NF >= 4) v[FNR] = $4; next }
              NF >= 4 { d = $4 - v[FNR]; if (d < 0) d = -d; s = v[FNR] < 0 ? -v[FNR] : v[FNR]; if (s > 0 && d / s > w) w = d / s; n++ }
              END { printf "%d %.3e", n, w }' "$f" "$b/output/$(basename "$f")")"
    n=$((n + ${r%% *})); w="$(awk -v a="$w" -v b="${r##* }" 'BEGIN { print (b > a ? b : a) }')"
  done
  echo "$n $w"
}

for m in sum total elast; do run "$m" "$m.azr" 1; run "x$m" "x$m.azr" 3; done
for kind in "data out ''" "extrapolation extrap x"; do
  set -- $kind; what="$1"; ext="$2"; p="${3//\'/}"
  read -r n w <<< "$(maxrel "$WORK/${p}sum" "$WORK/${p}total" "$ext")"
  if [ "$n" = missing ] || [ -z "$n" ] || [ "$n" -lt 100 ]; then bad "(a) $what: outputs missing ($n points)"
  elif awk -v w="$w" 'BEGIN { exit !(w <= 1e-9) }'; then ok "(a) $what: total = sum of every exit pair at $n points (max rel $w)"
  else bad "(a) $what: total vs sum of every exit pair: max rel $w at $n points"; fi
  read -r n w <<< "$(maxrel "$WORK/${p}sum" "$WORK/${p}elast" "$ext")"
  if awk -v w="$w" 'BEGIN { exit !(w > 1e-3) }'; then ok "(a) $what: elastic alone differs (max rel $w): the reaction channels matter"
  else bad "(a) $what: elastic alone is the sum (max rel $w): the test does not discriminate"; fi
done

# (b) the gradient check of tests/13N on the total model with only segment 1 active (Vaughn, 0-2.33 MeV: the natural
# sum of the three totals, inelastic channels open, no target integration), the reaction widths free.  (The
# target-integrated Sikkema segments disagree with finite differences by ~1e-2 on a few parameters with or without
# "total 1", from the integration grid.)
awk '/^dataset / { n++; print; if (n > 1) print "    active         no"; next } { print }' "$WORK/total.azr" > "$WORK/gtotal.azr"
run gcalc gtotal.azr 1
mkdir -p "$WORK/fit/output" "$WORK/fit/checks"; ln -s "$REF/data" "$WORK/fit/data"; cp "$WORK/gtotal.azr" "$WORK/fit/t.azr"
(cd "$WORK/fit" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
   run_limited 1200 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
chi="$(awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { sub(",", "", $3); print $3 }' "$WORK/fit/run.log")"
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/gcalc/output/chiSquared.out")"
if [ -n "$chi" ] && awk -v a="$chi" -v b="$want" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(d <= 1e-5) }'; then
  ok "(b) chi2 at the start of the fit $chi = calculation $want"
else bad "(b) chi2 at the start of the fit '$chi', calculation $want"; fi
awk '{ sub(/\r$/, "") }
     /^\[grad-check\] idx / { a = f = j = ""
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print a, f, j }' "$WORK/fit/run.log" > "$WORK/fit/rows"
nrow="$(awk 'END { print NR }' "$WORK/fit/rows")"
if [ "$nrow" -lt 5 ]; then bad "(b) gradient check listed $nrow parameters"
else
  read -r fd jac <<< "$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a
                               d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d
                               e = $1 - $3; if (e < 0) e = -e; if (e > v) v = e }
                             END { printf "%.2e %.2e", w / m, v / m }' "$WORK/fit/rows")"
  awk -v w="$fd" 'BEGIN { exit !(w + 0 <= 1e-5) }' && ok "(b) gradient = finite differences for $nrow parameters ($fd of the largest)" \
    || bad "(b) gradient vs finite differences: $fd of the largest"
  awk -v w="$jac" 'BEGIN { exit !(w + 0 <= 1e-10) }' && ok "(b) Jacobian-path gradient = adjoint gradient ($jac of the largest)" \
    || bad "(b) Jacobian-path vs adjoint gradient: $jac of the largest"
fi

# (c) refusals: data set Almeida1983_21Ne_ntot_lo (21Ne(n,n)21Ne, quantity sigma_total) or data set 1 (the natural
# sum) edited; a refused project is refused when it is read: its diagnostic, the refusal line, exit 1, no output
refuse() {
  local name="$1" set="$2" field="$3" value="$4" want="$5" dir="$WORK/$1"
  awk -v set="$set" -v field="$field" -v value="$value" '
    /^dataset / { on = ($2 == set) }
    on && $1 == field { printf "    %-14s %s\n", field, value; next } { print }' "$WORK/total.azr" > "$WORK/$name.azr"
  grep -qF -- "$value" "$WORK/$name.azr" || { bad "(c) $name: the edit did not apply"; return; }
  mkdir -p "$dir/output" "$dir/checks"; ln -s "$REF/data" "$dir/data"; cp "$WORK/$name.azr" "$dir/t.azr"
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); local rc=$?
  if grep -qF -- "$want" "$dir/run.log" && grep -q "^ERROR: t.azr (format 3) is refused" "$dir/run.log" && [ "$rc" = 1 ] \
    && [ ! -s "$dir/output/chiSquared.out" ]; then ok "(c) refused: $name"
  else bad "(c) not refused as expected: $name (exit $rc; $(grep -m1 ERROR "$dir/run.log"))"; fi
}
N21=Almeida1983_21Ne_ntot_lo
refuse charged $N21 reaction "17O(a,a)17O" "error E_REACTION_NOT_ALLOWED: dataset '$N21', field 'reaction': quantity 'sigma_total' with 17O(a,a)17O: the total cross section of charged particles is infinite"
refuse inelastic $N21 reaction "21Ne(n,n)21Ne[0.350725]" "quantity 'sigma_total' with 21Ne(n,n)21Ne[0.350725]: is defined for elastic scattering only"
refuse differential $N21 columns "E_lab[MeV] theta_lab[deg] value error" "quantity 'sigma_total' is angle-integrated"
refuse no_isotopes $N21 reaction "21Ne(n,n)21Ne + 20Ne(n,n)20Ne" "a sum of total cross sections is the total of a target with isotopes"
refuse twice Vaughn_ntot reaction "20Ne(n,n)20Ne + 22Ne(n,n)22Ne + 20Ne(n,n)20Ne" "two terms on pair 2: a target's total has one term per nucleus"
refuse factor Vaughn_ntot reaction "20Ne(n,n)20Ne + 2*22Ne(n,n)22Ne + 21Ne(n,n)21Ne" "a factor on a term of a target's total"

if [ "$failures" -eq 0 ]; then echo "PASS: total cross section (optical theorem)"; else echo "FAIL: $failures check(s)"; exit 1; fi
