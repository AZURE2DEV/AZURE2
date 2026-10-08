#!/usr/bin/env bash
#
# THM (HOES) segments under Park's parametrization (--use-park).
#
# Park's amplitudes are Brune's rescaled level by level, gamma_Park =
# gamma_Brune sqrt(J), J = 1 - sum_c gamma_Park^2 dS_c/dE, and Park's level
# matrix is Brune's with diag(sqrt J) on both sides.  The HOES amplitude
# sum_{l,m} gamma_{l f} A_{l m} v_m, v_m = sum_c gamma_{m c} M_l, is bilinear in
# the amplitudes, with level-diagonal factors only (a per-level vertex
# boundary, the line-shape factor N_C), so it is the same function of the
# physics in both parametrizations -- as long as everything that reads the
# amplitudes knows which they are: the line-shape width Gamma_lambda, the
# adaptive grid's resonance widths, and a channel entered as an amplitude
# (gammaIsRWA), which is Brune's in either mode.  Checked here:
#
#   1. the same .azr in both modes gives the same model at every point (rel
#      $TOL_MODEL), chi2, norms and THM experiment tables (rel $TOL_MODEL) and
#      physical widths (rel $TOL_WIDTH, the round trip of the input through
#      Brune's transformation; Park's is exact): every THM example, the
#      on-shell 6Li/7Li/15N examples, tests/18O_p_a_thm (amplitude input,
#      folding) and tests/17O, and 18O_p_a_thm with each THM feature that
#      reads the amplitudes or the levels -- line shape, R(E) and the DW vertex
#      (point Coulomb and optical), the spectator window, the fixed angle,
#      the coherent background, vertex=perlevel and onshell, the Coulomb
#      integral, entranceL=coherent with lambda32 kinematics and a spectator
#      energy;
#   2. parameter files convert on read: a Brune param.par read by a Park run
#      (and a Park one by a Brune run) reproduces the model (rel $TOL_FILE,
#      the 8 digits of param.par): 18O_p_a_thm, f19_pag_thm, tests/7Li_p_a
#      and tests/6Li_d with a fixed-angle window (identical particles) --
#      whose .azr widths lie beyond what Brune's transformation can reach
#      ("Denominator less than zero"), so only the file route compares them;
#   3. a fit (MIGRAD) in each mode from the same start: the same minimum (rel
#      1e-6) and physical widths within $TOL_FIT (the 8.6 MeV level's widths
#      are determined to ~30 %);
#   4. J <= 0 in a THM project: an observed width beyond the bound is reported
#      under Park and penalised (Total-Park-Chi-Squared); MCMC rejects every
#      such point (logP = -inf, acceptance 0) and samples a valid one.
#
# The analytic THM Jacobian under Park (chain rule through J) is checked by
# tests/pyazr/thm_park_test.py.  About ten minutes.
#
#   ./tests/thm_park/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_park.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-900}")"
TOL_MODEL="${TOL_MODEL:-1e-8}"
TOL_WIDTH="${TOL_WIDTH:-1e-5}"
TOL_FILE="${TOL_FILE:-1e-6}"
TOL_FIT="${TOL_FIT:-5e-3}"
KIN="beam=18O target=3He spectator=d Ebeam=115"
KIN18="beam=18O target=d spectator=n Ebeam=54"
fail=0

# stage NAME SRC [BLOCK [EDIT]]: a copy of project directory SRC (its first
# .azr, as run.azr) under $WORK/NAME; BLOCK replaces its <thm> block; EDIT is
# an awk program the .azr is passed through first.
stage() {
  local d="$WORK/$1" src="$2" azr
  rm -rf "$d"
  mkdir -p "$d/output" "$d/checks"
  cp -R "$src/data" "$d/data"
  azr="$(ls "$src"/*.azr | head -1)"
  if [ -n "${4:-}" ]; then tr -d '\r' < "$azr" | awk "$4" > "$d/run.azr"; else tr -d '\r' < "$azr" > "$d/run.azr"; fi
  if [ -n "${3:-}" ]; then
    awk '/^<thm>/ { T = 1 } !T { print } /^<\/thm>/ { T = 0 }' "$d/run.azr" > "$d/tmp.azr" && mv "$d/tmp.azr" "$d/run.azr"
    printf '\n<thm>\n%s\n</thm>\n' "$3" >> "$d/run.azr"
  fi
}
# run NAME MENU [FLAGS]: AZURE2 on $WORK/NAME/run.azr with the menu input
# MENU (printf format); status in NAME/status, output in NAME/log.
run() {
  local d="$WORK/$1"
  (cd "$d" && printf "$2" | $RUN "$AZURE2_BIN" --no-gui --no-readline ${3:-} run.azr 2>&1 | head -c 2000000 > log;
   echo "${PIPESTATUS[1]}" > status)
}
calc() { run "$1" '1\n\n\n7\n' "${2:-}"; }
calc_file() { run "$1" "1\n$2\n\n7\n" "${3:-}"; }
total() { tr -d '\r' < "$WORK/$1/output/chiSquared.out" 2>/dev/null | awk '/^Total-Chi-Squared:/ { print $2 }'; }
okrun() { [ "$(cat "$WORK/$1/status" 2>/dev/null)" = 0 ] && [ -f "$WORK/$1/output/chiSquared.out" ] || {
  bad "run $1 failed"; tail -5 "$WORK/$1/log" | sed 's/^/        /'; return 1; }; }

# models A B: the largest relative difference of the model (column 4) over
# every point of every AZUREOut file of run A and run B ("x" if the files or
# their energies differ).
models() {
  local f n=0 w=0 r
  for f in $(cd "$WORK/$1/output" && ls AZUREOut_* 2>/dev/null); do
    [ -f "$WORK/$2/output/$f" ] || { echo x; return; }
    r="$(paste <(tr -d '\r' < "$WORK/$1/output/$f" | awk 'NF > 4 { print $1, $4 }') \
               <(tr -d '\r' < "$WORK/$2/output/$f" | awk 'NF > 4 { print $1, $4 }') |
      awk '{ if ($1 != $3 || NF != 4) { print "x"; exit } d = $2 - $4; if (d < 0) d = -d; s = $2 < 0 ? -$2 : $2
             if (s > 0 && d / s > w) w = d / s; else if (s == 0 && d > 0) w = 1 } END { if (NR == 0) print "x"; else printf "%.3e\n", w }')"
    [ "$r" = x ] && { echo x; return; }
    w="$(awk -v a="$w" -v b="$r" 'BEGIN { print (b > a ? b : a) }')"
    n=$((n + 1))
  done
  [ "$n" -gt 0 ] && echo "$w" || echo x
}
# numclose FILE_A FILE_B TOL: the files agree token by token, numbers to
# relative TOL (" Total-Park-Chi-Squared: 0" ignored); prints the worst.
numclose() {
  [ -f "$1" ] && [ -f "$2" ] || { echo x; return 1; }
  awk -v tol="$3" '
    function num(s) { return s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?,?$/ }
    { sub(/ Total-Park-Chi-Squared: 0$/, "") }
    FNR == NR { a[FNR] = $0; n = FNR; next }
    { if (!(FNR in a)) { bad = 1; next }
      na = split(a[FNR], x, /[ \t,]+/); nb = split($0, y, /[ \t,]+/)
      if (na != nb) { bad = 1; next }
      for (i = 1; i <= na; i++) {
        if (num(x[i]) && num(y[i])) { d = x[i] - y[i]; if (d < 0) d = -d; s = x[i] < 0 ? -x[i] : x[i]
          r = s > 0 ? d / s : d; if (r > w) w = r }
        else if (x[i] != y[i]) bad = 1 } }
    END { if (FNR != n) bad = 1; printf "%s\n", bad ? "x" : sprintf("%.3e", w); exit !(!bad && w <= tol) }' \
    <(tr -d '\r' < "$1") <(tr -d '\r' < "$2")
}
# widths FILE: every G of parameters.out in eV.
widths() { tr -d '\r' < "$1" | awk '{ for (i = 1; i < NF - 1; i++) if ($i == "G" && $(i+1) == "=") { v = $(i+2); u = $(i+3)
  if (u == "keV") v *= 1e3; else if (u == "MeV") v *= 1e6; else if (u == "meV") v *= 1e-3; print v } }'; }
widthdiff() {
  paste <(widths "$WORK/$1/output/parameters.out") <(widths "$WORK/$2/output/parameters.out") |
    awk '{ if (NF != 2) { bad = 1; next } d = $1 - $2; if (d < 0) d = -d; s = $1 < 0 ? -$1 : $1
           if (s > 0 && d / s > w) w = d / s; n++ } END { if (bad || n == 0) print "x"; else printf "%.3e %d\n", w, n }'
}
below() { awk -v a="$1" -v t="$2" 'BEGIN { exit !(a != "x" && a + 0 <= t + 0) }'; }

# same_physics TAG SRC [BLOCK]: Brune and Park runs of the same .azr agree.
same_physics() {
  local tag="$1"
  stage "$tag.b" "$2" "${3:-}"; stage "$tag.p" "$2" "${3:-}"
  calc "$tag.b"; calc "$tag.p" --use-park
  okrun "$tag.b" && okrun "$tag.p" || return
  local m w f r
  m="$(models "$tag.b" "$tag.p")"
  read -r w _ <<< "$(widthdiff "$tag.b" "$tag.p")"
  below "$m" "$TOL_MODEL" && ok "$tag: model, Brune = Park at every point (max rel $m), chi2 $(total "$tag.b")" \
    || bad "$tag: model differs (max rel $m), chi2 Brune $(total "$tag.b"), Park $(total "$tag.p")"
  below "$w" "$TOL_WIDTH" && ok "$tag: physical widths agree (max rel $w)" || bad "$tag: widths differ (max rel $w)"
  for f in chiSquared.out normalizations.out thm_experiments.out; do
    [ -f "$WORK/$tag.b/output/$f" ] || continue
    r="$(numclose "$WORK/$tag.b/output/$f" "$WORK/$tag.p/output/$f" "$TOL_MODEL")" \
      && ok "$tag: $f agrees (max rel $r)" || bad "$tag: $f differs (max rel $r)"
  done
  grep -q "^ *#parametrization *2\." "$WORK/$tag.p/output/param.par" && ! grep -q "is not positive" "$WORK/$tag.p/log" \
    || bad "$tag: Park param.par not tagged, or a level with J <= 0"
}

# converted TAG SRC [BLOCK] [both]: Brune's param.par read in a Park run (and,
# with "both", Park's param.par in a Brune run) gives the Brune model.
converted() {
  local tag="$1" m
  stage "$tag.b" "$2" "${3:-}"; calc "$tag.b"; okrun "$tag.b" || return
  stage "$tag.bp" "$2" "${3:-}"; cp "$WORK/$tag.b/output/param.par" "$WORK/$tag.bp/brune.par"
  calc_file "$tag.bp" brune.par --use-park; okrun "$tag.bp" || return
  m="$(models "$tag.b" "$tag.bp")"
  grep -q "Converted the parameter file from Brune to Park" "$WORK/$tag.bp/log" && below "$m" "$TOL_FILE" \
    && ok "$tag: Brune param.par read under Park, converted, same model (max rel $m)" \
    || bad "$tag: Brune param.par under Park: model max rel $m (conversion message: $(grep -c 'Converted the parameter file' "$WORK/$tag.bp/log"))"
  [ "${4:-}" = both ] || return
  stage "$tag.p" "$2" "${3:-}"; calc "$tag.p" --use-park; okrun "$tag.p" || return
  stage "$tag.pb" "$2" "${3:-}"; cp "$WORK/$tag.p/output/param.par" "$WORK/$tag.pb/park.par"
  calc_file "$tag.pb" park.par; okrun "$tag.pb" || return
  m="$(models "$tag.b" "$tag.pb")"
  grep -q "Converted the parameter file from Park (observed) to Brune" "$WORK/$tag.pb/log" && below "$m" "$TOL_FILE" \
    && ok "$tag: Park param.par read under Brune, converted, same model (max rel $m)" \
    || bad "$tag: Park param.par under Brune: model max rel $m"
}

T18="$ROOT/tests/18O_p_a_thm"
echo "1. the same .azr, Brune and Park"
for e in o18_lacognata2010 o18_lacognata2008 f19_pag_thm c12c12_tumino2018 o17_guardo2017_fit \
         li6_pizzone2011 li7_tumino2006 n15_lacognata2007; do
  same_physics "$e" "$ROOT/examples/$e"
done
same_physics 18O_p_a_thm "$T18"
same_physics 17O "$ROOT/tests/17O"
same_physics 17O_perlevel "$ROOT/tests/17O" "vertex=perlevel"
same_physics 18O.lineshape "$T18" "experiment[A] segments=1,2 $KIN lineshape=on"
same_physics 18O.dw_coulomb_lineshape "$T18" "experiment[A] segments=1,2 $KIN distortion=coulomb vertexModel=dw lineshape=on"
same_physics 18O.dw_optical "$T18" "experiment[A] segments=1,2 $KIN18 distortion=optical opticalAA=ancai06 opticalSF=kd03:extrapolate vertexModel=dw background=linear"
same_physics 18O.distortion "$T18" "experiment[A] segments=1,2 $KIN distortion=coulomb"
same_physics 18O.window "$T18" "experiment[A] segments=1,2 $KIN ps=hulthen:0-40 psNodes=8"
same_physics 18O.angle_window "$T18" "experiment[A] segments=1,2 $KIN18 ps=hulthen:0-40 background=linear theta=50-70"
same_physics 18O.cbackground "$T18" "experiment[A] segments=1,2 background=linear cbackground=1/2+:2=0.3,-0.2"
same_physics 18O.cbackground_perlevel "$T18" "vertex=perlevel
experiment[A] segments=1,2 cbackground=1/2+:2=0.3,-0.2"
same_physics 18O.onshell "$T18" "vertex=onshell"
same_physics 18O.coulomb_integral "$T18" "coulombIntegral=1
kinematics=kf3body"
same_physics 18O.coherent_l "$T18" "entranceL=coherent
kinematics=lambda32
spectatorEnergy=0.3"

echo "2. parameter files converted on read"
converted 18O_p_a_thm.file "$T18" "" both
converted f19_pag_thm.file "$ROOT/examples/f19_pag_thm" "" both
converted 7Li_p_a.file "$ROOT/tests/7Li_p_a"
converted 6Li_d.angle.file "$ROOT/tests/6Li_d" "experiment[A] segments=1 theta=10-40"
converted 7Li_p_a.angle.file "$ROOT/tests/7Li_p_a" "experiment[A] segments=1 theta=30-60"
# Why the file route for these two: their widths are beyond Brune's reach.
grep -q "Denominator less than zero" "$WORK/7Li_p_a.file.b/log" && ok "7Li_p_a: the .azr widths are beyond Brune's transformation (warned)" \
  || bad "7Li_p_a: expected Brune's 'Denominator less than zero' warning"

echo "3. a fit in each mode"
# tests/18O_p_a_thm unfolded, the 8.6 MeV level 10 keV high and free, every
# amplitude free.
FREE='/^<targetInt>/ { print; skip = 1; next } /^<\/targetInt>/ { skip = 0 } skip { next }
      /^<levels>/ { L = 1 } /^<\/levels>/ { L = 0 }
      L && NF > 30 { if ($3 == "8.602600") { $3 = "8.612600"; $4 = 0 } $11 = 0 } { print }'
stage fit.b "$T18" "" "$FREE"; stage fit.p "$T18" "" "$FREE"
run fit.b '2\nn\n\n\n7\n'; run fit.p '2\nn\n\n\n7\n' --use-park
if okrun fit.b && okrun fit.p; then
  cb="$(tr -d '\r' < "$WORK/fit.b/log" | awk '/^Total Chi-Squared:/ { v = $3 } END { print v }')"
  cp="$(tr -d '\r' < "$WORK/fit.p/log" | awk '/^Total Chi-Squared:/ { v = $3 } END { print v }')"
  awk -v a="$cb" -v b="$cp" 'BEGIN { d = (a - b) / a; if (d < 0) d = -d; exit !(a != "" && b != "" && d < 1e-6 && a < 2228) }' \
    && ok "fit: the same minimum, Brune $cb, Park $cp (start 2228.9)" || bad "fit: minimum Brune '$cb', Park '$cp'"
  read -r w n <<< "$(widthdiff fit.b fit.p)"
  below "$w" "$TOL_FIT" && ok "fit: the $n physical widths agree (max rel $w)" || bad "fit: widths differ (max rel $w)"
  eb="$(tr -d '\r' < "$WORK/fit.b/output/param.sav" | awk '$1 == "energy_1" { print $2 }')"
  ep="$(tr -d '\r' < "$WORK/fit.p/output/param.sav" | awk '$1 == "energy_1" { print $2 }')"
  awk -v a="$eb" -v b="$ep" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(a != "" && d < 1e-5 && a != 8.6126) }' \
    && ok "fit: level energy moved to $eb (Brune), $ep (Park) MeV" || bad "fit: level energy Brune '$eb', Park '$ep'"
  grep -q "^ *#parametrization *2\." "$WORK/fit.p/output/param.sav" && ok "fit: param.sav tagged Park" || bad "fit: param.sav not tagged Park"
fi

echo "4. J <= 0"
# The 8.81 MeV level's p + 18O width given as an observed 3 MeV, far past the
# bound 2P/(dS/dE) at 5.1 fm.
WALL='/^<levels>/ { W = 1 } /^<\/levels>/ { W = 0 }
      W && NF > 30 && $3 == "8.805800" && $6 == 1 { $12 = 3000000; $33 = 0 } { print }'
stage_wall() {
  stage "$1" "$T18" "" "$FREE"
  awk "$WALL" "$WORK/$1/run.azr" > "$WORK/$1/tmp.azr" && mv "$WORK/$1/tmp.azr" "$WORK/$1/run.azr"
}
stage_wall wall; stage_wall wall.b
calc wall --use-park; calc wall.b
if okrun wall && okrun wall.b; then
  grep -q "Park norm J = .* is not positive" "$WORK/wall/log" && ok "wall: J <= 0 reported under Park" || bad "wall: no J <= 0 warning"
  pen="$(tr -d '\r' < "$WORK/wall/output/chiSquared.out" | awk '/^Total-Chi-Squared:/ { for (i = 1; i < NF; i++) if ($i == "Total-Park-Chi-Squared:") print $(i+1) }')"
  awk -v p="$pen" 'BEGIN { exit !(p != "" && p > 1e3) }' && ok "wall: penalised, Total-Park-Chi-Squared = $pen" \
    || bad "wall: Total-Park-Chi-Squared '$pen'"
  grep -q "Denominator less than zero" "$WORK/wall.b/log" && ok "wall: Brune warns about the same input" || bad "wall: no Brune warning"
fi
# MCMC (when built): 12 walkers, 100 steps, RWA, 1 % spread, 0.5 keV.
MCMC='6\n\n12\n100\n1\n0.5\n1\nyes\nyes\n'
if grep -q "Perform MCMC" "$WORK/wall/log"; then
  stage_wall wall.mcmc
  run wall.mcmc "$MCMC" --use-park
  acc="$(tr -d '\r' < "$WORK/wall.mcmc/log" | grep -o 'Acceptance fraction: [0-9.]*' | awk '{ print $3 }')"
  if [ "$acc" = "0.000" ] && tr -d '\r' < "$WORK/wall.mcmc/output/samples.mcmc" | awk -F, 'NR > 1 { n++; if ($3 != "-inf") b++ } END { exit !(n > 0 && b == 0) }'; then
    ok "MCMC: every point beyond J > 0 rejected (logP = -inf, acceptance $acc)"
  else
    bad "MCMC: points beyond J > 0 not rejected (acceptance '$acc')"
  fi
  stage mcmc "$T18" "" "$FREE"
  run mcmc "$MCMC" --use-park
  acc="$(tr -d '\r' < "$WORK/mcmc/log" | grep -o 'Acceptance fraction: [0-9.]*' | awk '{ print $3 }')"
  awk -v a="$acc" 'BEGIN { exit !(a != "" && a > 0.05) }' && grep -q "MCMC sampling completed successfully" "$WORK/mcmc/log" \
    && ok "MCMC: a THM project under Park samples (acceptance $acc)" || bad "MCMC: acceptance '$acc'"
else
  echo "  (MCMC not built: skipped)"
fi

echo
if [ "$fail" -eq 0 ]; then echo "PASS: THM under Park's parametrization"; else echo "FAIL"; exit 1; fi
