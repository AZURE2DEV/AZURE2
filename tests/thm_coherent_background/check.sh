#!/usr/bin/env bash
#
# Coherent (interfering) THM background: `cbackground=` on an
# `experiment[<name>]` line of the <thm> block (ThmExperiment.h
# ThmCoherentBackground; docs/source/theory/thm_implementation.rst,
# "Coherent background").  A THM-only complex amplitude c(E) = c0 + c1 E,
# times the entrance vertex M_l of the bucket, is added to the resonant HOES
# amplitude of one (J^pi, entrance (s,l), exit (s',l')) combination before
# squaring; Re and Im of c0 (and c1) are ordinary fit parameters (cbkg_*).
#
# Models: toy/ (one 1/2+ level, neutral l = 0 entrance n + 16O and exit
# n + 12C channels, formal parameters: every factor is elementary) and
# tests/18O_p_a_thm (the 18O(p,a) 1/2+ doublet, two THM segments).
#
#   (a) without the key nothing moves: an experiment line without
#       cbackground gives the files of no line (tests/thm_experiment), and
#       the suite's pins are unchanged;
#   (b) a coherent background fixed at zero gives the model and chi2 of none
#       -- angle-integrated, with a theta window, with a ps window and a
#       linear incoherent background, and with the DW vertex (18O);
#   (c) the toy model against the closed form
#         m(E) = (2J+1) (k_f/mu_f) 2 P_f |M_0 (g_f g_c A + c(E))|^2,
#         A = 1 / (E_l - E - i (g_c^2 P_c + g_f^2 P_f)),  P = k a,
#         M_0 = (B - 1) j_0(pa) - pa j_0'(pa)  (B = 0; j_0' by the engine's
#         forward difference, step 1e-6),
#       for c = 0, a constant and a linear c; the ratio m(c)/m(0) =
#       |g_f g_c A + c|^2 / |g_f g_c A|^2 needs no vertex at all; a 0-180
#       window times 4 pi is the angle-integrated model with the background;
#   (d) synthetic 18O data made with c0 = 0.3 - 0.2 i are fitted (MIGRAD,
#       c0 starting at 0) back to c0 and chi2 ~ 0; param.sav and
#       thm_experiments.out carry the values;
#   (e) the Jacobian is checked in tests/pyazr/thm_coherent_background_test.py;
#   (f) refusals.
#
#   ./tests/thm_coherent_background/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
O18="$HERE/../18O_p_a_thm"
TOY="$HERE/toy"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_cbkg.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-300}"; else RUN=""; fi
OUT="AZUREOut_aa=1_R=2.out"
fail=0

# run NAME SRC AZR BLOCK [OPTION] -- a run of a copy of SRC/AZR with BLOCK
# (literal lines) as its <thm> block ("" = none); OPTION 1 (calculate with
# data, default) or 2 (fit).  Exit status in NAME/status.
run() {
  local d="$WORK/$1"
  rm -rf "$d"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$2/data" "$d/"
  cp "$2/$3" "$d/run.azr"
  [ -z "$4" ] || printf '\n<thm>\n%s\n</thm>\n' "$4" >> "$d/run.azr"
  (cd "$d" && printf '%s\n\n\n7\n' "${5:-1}" | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 | head -c 1000000 > log;
   echo "${PIPESTATUS[1]}" > status)
}
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; fail=1; }
same() { [ -f "$1" ] && [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
total() { tr -d '\r' < "$1/log" | awk '/Total Chi-Squared:/ { c = $3 } END { print c }'; }
# worst |a/b - 1| of column 4 (the model) of two output files, data lines only
worst4() { paste <(tr -d '\r' < "$1" | awk 'NF > 3 { print $4 }') <(tr -d '\r' < "$2" | awk 'NF > 3 { print $4 }') |
  awk -v f="${3:-1}" '{ d = ($1 * f - $2) / $2; if (d < 0) d = -d; if (d > w) w = d; n++ } END { printf "%.2e %d", w, n }'; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a != "" && a + 0 <= b + 0) }'; }
started() { [ "$(cat "$WORK/$1/status")" = 0 ] && [ -f "$WORK/$1/output/$OUT" ] || { bad "$1: the run failed"; tail -5 "$WORK/$1/log" | sed 's/^/        /'; }; }

# (a), (b) ------------------------------------------------------------------
echo "(a, b) no key, and a background fixed at zero, change nothing"
run none "$O18" 18O_p_a_thm.azr "experiment[E] segments=1,2"
run zero "$O18" 18O_p_a_thm.azr "experiment[E] segments=1,2 cbackground=1/2+:2=0f,0f"
started none; started zero
for f in chiSquared.out normalizations.out "$OUT"; do
  same "$WORK/none/output/$f" "$WORK/zero/output/$f" && ok "c = 0 (fixed): $f byte-identical" || bad "c = 0: $f differs"
done
grep -q "^cbkg cbkg_E_1/2+_2_1/2,0,1/2,1_re0 .* fixed" "$WORK/zero/output/thm_experiments.out" &&
  ok "thm_experiments.out lists the fixed parameters" || bad "thm_experiments.out has no cbkg lines"
grep -q "cbkg" "$WORK/none/output/thm_experiments.out" "$WORK/none/output/param.par" &&
  bad "without the key a cbkg line was written" || ok "without the key: no cbkg line in thm_experiments.out or param.par"
for w in "theta=50-70" "beam=18O target=d spectator=n Ebeam=54 ps=hulthen:0-40 background=linear" \
         "beam=18O target=d spectator=n Ebeam=54 distortion=coulomb vertexModel=dw"; do
  run w0 "$O18" 18O_p_a_thm.azr "experiment[E] segments=1,2 $w"
  run w1 "$O18" 18O_p_a_thm.azr "experiment[E] segments=1,2 $w cbackground=1/2+:2:linear=0f,0f,0f,0f"
  same "$WORK/w0/output/$OUT" "$WORK/w1/output/$OUT" && [ "$(total "$WORK/w0")" = "$(total "$WORK/w1")" ] &&
    ok "c = 0 with $w: model and chi2 identical" || bad "c = 0 with $w differs"
  run w2 "$O18" 18O_p_a_thm.azr "experiment[E] segments=1,2 $w cbackground=1/2+:2=0.3,-0.2"
  same "$WORK/w0/output/$OUT" "$WORK/w2/output/$OUT" && bad "c != 0 with $w: no effect" || ok "c != 0 with $w moves the model"
done

# (c) -----------------------------------------------------------------------
echo "(c) toy model: the closed form |M (g_f g_c A + c(E))|^2"
# formula CR0 CI0 CR1 CI1 FILE -> worst relative deviation of the model
# (column 4) from the closed form; RATIO=1: of m(c)/m(0) (FILE0 the c = 0 run)
formula() {
  tr -d '\r' < "$5" | awk -v cr0="$1" -v ci0="$2" -v cr1="$3" -v ci1="$4" -v ref="${6:-}" '
    BEGIN { hc = 197.32696310; u = 931.4940880; mn = 1.00866492
            mi = mn * 15.99052821 / (mn + 15.99052821) * u; mf = mn * 11.99670964 / (mn + 11.99670964) * u
            if (ref != "") while ((getline line < ref) > 0) { split(line, f, " "); if (f[4] != "") m0[++k] = f[4] } }
    NF > 3 { E = $1; n++
      Pc = sqrt(2 * mi * E) / hc * 5.0; kf = sqrt(2 * mf * (E + 2.143)) / hc; Pf = kf * 4.5
      x = sqrt(2 * mi * (E + 2.224566)) / hc * 5.0; h = 1e-6
      M = -sin(x) / x - x * (sin(x + h) / (x + h) - sin(x) / x) / h     # (B - 1) j0 - x j0", B = 0
      D = 0.5 - E; G = 0.09 * Pc + 0.04 * Pf; den = D * D + G * G       # E_l - E = 0.5 - E (c.m.)
      rr = 0.06 * D / den; ri = 0.06 * G / den                         # g_f g_c A
      ar = rr + cr0 + cr1 * E; ai = ri + ci0 + ci1 * E
      if (ref != "") { want = (ar * ar + ai * ai) / (rr * rr + ri * ri); got = $4 / m0[n] }
      else { want = 2 * (kf / mf) * 2 * Pf * M * M * (ar * ar + ai * ai); got = $4 }
      d = (got - want) / want; if (d < 0) d = -d; if (d > w) w = d }
    END { printf "%.2e", w }'
}
run t0 "$TOY" toy.azr "experiment[E] segments=1"
run t1 "$TOY" toy.azr "experiment[E] segments=1 cbackground=1/2+:2=0.05,-0.08"
run t2 "$TOY" toy.azr "experiment[E] segments=1 cbackground=1/2+:2:1/2,0,1/2,0:linear=0.05,-0.08,-0.03,0.1"
started t0; started t1; started t2
w="$(formula 0 0 0 0 "$WORK/t0/output/$OUT")"
lt "$w" 1e-7 && ok "c = 0: the closed form (worst rel $w; awk's own difference quotient limits it)" || bad "c = 0: $w"
w="$(formula 0.05 -0.08 0 0 "$WORK/t1/output/$OUT")"
lt "$w" 1e-7 && ok "c = 0.05 - 0.08 i: the closed form (worst rel $w)" || bad "constant c: $w"
w="$(formula 0.05 -0.08 -0.03 0.1 "$WORK/t2/output/$OUT")"
lt "$w" 1e-7 && ok "c = (0.05 - 0.08 i) + (-0.03 + 0.1 i) E: the closed form (worst rel $w)" || bad "linear c: $w"
w="$(formula 0.05 -0.08 0 0 "$WORK/t1/output/$OUT" "$WORK/t0/output/$OUT")"
lt "$w" 1e-9 && ok "m(c)/m(0) = |g_f g_c A + c|^2 / |g_f g_c A|^2 (worst rel $w)" || bad "ratio, constant c: $w"
w="$(formula 0.05 -0.08 -0.03 0.1 "$WORK/t2/output/$OUT" "$WORK/t0/output/$OUT")"
lt "$w" 1e-9 && ok "m(c)/m(0), linear c (worst rel $w)" || bad "ratio, linear c: $w"
run t3 "$TOY" toy.azr "experiment[E] segments=1 theta=0-180 cbackground=1/2+:2:linear=0.05,-0.08,-0.03,0.1"
started t3
read -r w n <<< "$(worst4 "$WORK/t3/output/$OUT" "$WORK/t2/output/$OUT" 12.566370614359172)"
lt "$w" 2e-10 && ok "theta=0-180 x 4 pi == angle-integrated, with the background ($n points, worst rel $w)" ||
  bad "0-180 window with the background: $w"

# (d) -----------------------------------------------------------------------
echo "(d) a coherent background is recovered by a fit"
run gen "$O18" 18O_p_a_thm.azr "experiment[E] segments=1 cbackground=1/2+:2=0.3,-0.2"
started gen
mkdir -p "$WORK/syn/data"
cp "$O18/18O_p_a_thm.azr" "$WORK/syn/"; cp "$O18/data/"* "$WORK/syn/data/"
F=data/lc723_thm_points.dat
n="$(tr -d '\r' < "$O18/$F" | grep -cv '^#')"
paste <(tr -d '\r' < "$O18/$F" | grep -v '^#') <(tr -d '\r' < "$WORK/gen/output/$OUT" | awk 'NF > 3' | head -n "$n" | awk '{ print $4 }') |
  awk '{ printf "%s %s %.12e %.12e\n", $1, $2, 250 * $5, 2.5 * $5 }' > "$WORK/syn/$F"
run fit "$WORK/syn" 18O_p_a_thm.azr "experiment[E] segments=1 cbackground=1/2+:2" 2
[ "$(cat "$WORK/fit/status")" = 0 ] && [ -f "$WORK/fit/output/param.sav" ] || { bad "fit run failed"; tail -5 "$WORK/fit/log"; }
re0="$(tr -d '\r' < "$WORK/fit/output/param.sav" | awk '$1 ~ /^cbkg_E_.*_re0$/ { print $2 }')"
im0="$(tr -d '\r' < "$WORK/fit/output/param.sav" | awk '$1 ~ /^cbkg_E_.*_im0$/ { print $2 }')"
near() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(a != "" && d <= t) }'; }
near "$re0" 0.3 1e-4 && near "$im0" -0.2 1e-4 && ok "c0 = $re0 + ($im0) i recovered (true 0.3 - 0.2 i)" || bad "c0 = $re0, $im0"
nrm="$(tr -d '\r' < "$WORK/fit/output/thm_experiments.out" | awk '$1 == "norm" { print $2; exit }')"
near "$nrm" 0.004 1e-4 && ok "norm 1/250 recovered: $nrm" || bad "norm $nrm (want 0.004)"
c1="$(tr -d '\r' < "$WORK/fit/output/chiSquared.out" | awk -F, 'NR == 2 { print $2 }')"
lt "$c1" 1e-3 && ok "segment 1 chi2 ~ 0 ($c1)" || bad "segment 1 chi2 $c1"
grep -q "^cbkg cbkg_E_1/2+_2_1/2,0,1/2,1_re0" "$WORK/fit/output/thm_experiments.out" && ok "thm_experiments.out carries the fitted values" ||
  bad "thm_experiments.out without the values"

# (f) -----------------------------------------------------------------------
echo "(f) refused"
refuse() {  # refuse NAME MESSAGE-FRAGMENT BLOCK
  run "$1" "$O18" 18O_p_a_thm.azr "$3"
  if [ "$(cat "$WORK/$1/status")" != 0 ] && grep -q "ERROR: <thm>" "$WORK/$1/log" && grep -q -- "$2" "$WORK/$1/log"; then
    ok "$1: $(grep -m1 'ERROR' "$WORK/$1/log" | tr -d '\r' | cut -c1-150)"
  else
    bad "$1: not refused with '$2' (status $(cat "$WORK/$1/status"))"; tail -3 "$WORK/$1/log" | sed 's/^/        /'
  fi
}
X="experiment[E] segments=1,2"
refuse no_parity "J^pi '1/2'" "$X cbackground=1/2:2"
refuse bad_spin "J^pi '1/3+'" "$X cbackground=1/3+:2"
refuse bad_exit "exit pair key 'x'" "$X cbackground=1/2+:x"
refuse bad_form "channels 'cubic'" "$X cbackground=1/2+:2:cubic"
refuse bad_channels "expected <s>,<l>,<s'>,<l'>" "$X cbackground=1/2+:2:1/2,0"
refuse few_values "expected 2 values" "$X cbackground=1/2+:2=1"
refuse linear_values "expected 4 values" "$X cbackground=1/2+:2:linear=1,2"
refuse bad_value "value 'a'" "$X cbackground=1/2+:2=1,a"
refuse empty_term "an empty term" "$X cbackground=1/2+:2;"
refuse no_group "no J^pi = 3/2+ group" "$X cbackground=3/2+:2"
refuse not_exit "no segment of the experiment has exit pair 1" "$X cbackground=1/2+:1"
refuse no_channel "no entrance channel (s,l) = (1/2,1)" "$X cbackground=1/2+:2:1/2,1,1/2,1"
refuse twice "is given twice" "$X cbackground=1/2+:2;1/2+:2:1/2,0,1/2,1"
refuse coherentL "entranceL=coherent" "entranceL=coherent
$X cbackground=1/2+:2"
refuse lineshape "lineshape=on and cbackground= cannot be combined" \
  "$X beam=18O target=d spectator=n Ebeam=54 lineshape=on cbackground=1/2+:2"

echo
if [ "$fail" -eq 0 ]; then echo "PASS: THM coherent background"; else echo "FAIL"; exit 1; fi
