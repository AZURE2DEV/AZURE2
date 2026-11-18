#!/usr/bin/env bash
#
# Isospin-labelled levels across compound nuclei (the records of <isospin>; dev's test of 2026-10-17, re-written for
# format 3).  The models are dev's, converted by `AZURE2 --convert` (README.md): a toy A = 4 system in Hale's
# charge-independent scheme (p+3He, n+t, p+t, n+3He, d+d; a 2- T=1 level in 4Li, 4H and 4He, a 2- T=0 level in 4He),
# and dev's two-nucleus capture models (12C+p / 12C+n).
#
# (a) The factors printed at start-up are Hale's A=4 relations, and isospin.out obeys them: derived energies =
#     reference + threshold difference + shift, derived amplitudes = factor ratio x the free one, the T=1 level's d+d
#     amplitude 0, the T=0 level's n+3He amplitude = -(p+t).  parameters.out labels the members.
# (b) The 4H member equals a standalone 4H project built from its formal parameters (--no-transform) to 1e-10.
# (c) Refused: a cross-nucleus level without a record, unequal radii in one arrangement, a mix on an allowed channel,
#     T too small for a nucleus, Park's basis.
# (d) Pseudo-data from iso_mix.azr (shifts -0.86 / -0.64, mix 0.05): MIGRAD (numerical and analytic gradient) and LM
#     from -0.5 / -0.4 / 0 recover all three; param.sav holds the derived slots synced.
# (e) AZURE_GRAD_CHECK on iso_free.azr (every reference parameter free, priors): the analytic gradient and the
#     Jacobian path agree with finite differences through the ties.
# (g) Soft tie (iso_soft.azr): at the start every extrapolation equals the hard-tie model's; the gradient check (from
#     shifts -0.8 / -0.6) agrees for all 13 free parameters; a MIGRAD fit reports the soft penalty in
#     Total-Isospin-Chi-Squared and isospin.out; the least-squares fitter takes the soft tie as residual rows and
#     ends on the same objective.
# (h) relativistic=1: isospin.out keeps the energy relations, the 4H member equals the standalone 4H project run
#     relativistically, the key changes the tied cross section, the gradient check passes through the ties.
# (i) Brune basis (the default; (a)-(h) run the formal one): the same ties on observed energies and Brune amplitudes;
#     isospin.out prints no boundary values; the 4H member equals the standalone 4H project from its Brune parameters;
#     the reference nucleus has the formal run's p+3He cross section while the tied n+t one differs; the soft tie
#     starts on the hard tie; the gradient checks pass; a MIGRAD fit of Brune pseudo-data recovers shifts and mix.
# (j) Photon channels in a tied J-group stay outside the ties: every level of a one-nucleus capture model in a block
#     gives the model without a block in both bases (the formal block sets B = -l: capture to a final state that is
#     not at the boundary energy); 13N + 13C tied: the report lists the photon channels as not tied, the reference
#     nucleus's capture is the untied model's, the 13C capture follows the tie.
# (k) Excited particles are ranked: n + 12C*(7.654, 0+) alone in 13C is free (key 0+#2); with p + 12C*(7.654) in 13N
#     it is tied to it, and the ground-state pair stays tied.
#
# (f) The ensemble sampler (a build with it): a short chain in the fit's basis samples the 8 free parameters with the
#     isospin priors in the likelihood; a chain on observed parameters is refused.
#
#   ./tests/isospin_levels/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
M="$HERE/models"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }
command -v python3 >/dev/null 2>&1 || { echo "SKIP isospin_levels: python3 not found"; exit 0; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/isospin_levels.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

nok=0; nbad=0
ok()  { echo "  ok   $1"; nok=$((nok + 1)); }
bad() { echo "  FAIL $1"; nbad=$((nbad + 1)); }
# run NAME MODEL COMMAND FLAGS [ENV=VALUE...]: a fresh directory per run, with the pseudo-data written so far
run() {
  local dir="$WORK/$1" model="$2" cmd="$3" flags="$4"; shift 4
  mkdir -p "$dir/output" "$dir/checks"
  cp "$model" "$dir/t.azr"
  for f in nt.dat ptn.dat ddp.dat ddn.dat; do [ -f "$WORK/$f" ] && cp "$WORK/$f" "$dir/"; done
  # shellcheck disable=SC2086
  (cd "$dir" && env "$@" $(guard_command 900) "$AZURE2_BIN" $cmd $flags t.azr < /dev/null > run.log 2>&1)
}
xs() { awk 'NF >= 5 { printf "%.12e %.12e\n", $1, $4 }' "$1"; }
xsec() { awk 'NF >= 4 { printf "%.12e %.12e\n", $1, $(NF - 1) }' "$1"; }
maxrel() {
  paste "$1" "$2" | awk '
    NF >= 4 { d = $2 - $4; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300;
              r = d / m; if (r > worst) worst = r; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && r + 0 <= t + 0) }'; }
# Byte-identical files without cmp (the Windows MSYS2 runner has none); \r stripped so a CRLF-writing build compares alike.
same_file() { [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
# member FILE: the energy and the first amplitude of the 4H member of level 1 in an isospin.out
member() {
  awk '/^level ID 1/ { on = 1 } on && /Z=1 A=4/ { for (i = 1; i <= NF; i++) if ($i == "E") e = $(i + 2); getline;
       for (i = 1; i <= NF; i++) if ($i == "gamma") g = $(i + 2); print e, g; exit }' "$1"
}
# standalone E GAMMA: the standalone 4H model with these level parameters
standalone() { awk -v e="$1" -v g="$2" '/<levels>/ { print; getline; $3 = e; $12 = g; print; next } { print }' "$M/standalone_4H.azr"; }
rel() { awk '/<\/config>/ { print "relativistic=1" } { print }' "$1" > "$2"; }
# relations FILE LABEL WORD: the seven relations of (a) on an isospin.out whose energies are "<WORD> E ="
relations() {
  python3 - "$1" "$2" "$3" "${4:-7}" <<'PYEOF'
import sys, math
lv, cur, mem = {}, None, None
word = sys.argv[3] + " E ="
for raw in open(sys.argv[1]):
    t = raw.split()
    if raw.startswith("level ID"):
        cur = int(t[2]); lv[cur] = {}
    elif cur is not None and word in raw:
        mem = t[0] + t[1]; lv[cur][mem] = {"E": float(t[t.index("E") + 2]), "g": {}}
    elif cur is not None and raw.strip().startswith("pair"):
        lv[cur][mem]["g"][int(t[1])] = float(t[t.index("gamma") + 2])
li, h, he = "Z=3A=4", "Z=1A=4", "Z=2A=4"
r = math.sqrt(0.5)
try:
    checks = [
      ("E(4H)  = E(4Li) - 0.86", lv[1][h]["E"], lv[1][li]["E"] - 0.86),
      ("E(4He) = E(4Li) + 19.8136 - 0.64", lv[1][he]["E"], lv[1][li]["E"] + 19.8136 - 0.64),
      ("gamma(n+t) = gamma(p+3He)", lv[1][h]["g"][2], lv[1][li]["g"][1]),
      ("gamma(p+t, T=1) = gamma(p+3He)/sqrt2", lv[1][he]["g"][3], r * lv[1][li]["g"][1]),
      ("gamma(n+3He, T=1) = gamma(p+3He)/sqrt2", lv[1][he]["g"][4], r * lv[1][li]["g"][1]),
      ("gamma(d+d, T=1) = 0", lv[1][he]["g"][5], 0.0),
      ("gamma(n+3He, T=0) = -gamma(p+t, T=0)", lv[2][he]["g"][4], -lv[2][he]["g"][3]),
    ]
except KeyError as e:
    print("  FAIL %sisospin.out: member or channel missing (%s)" % (sys.argv[2], e)); sys.exit(7)
bad = 0
for name, a, b in checks[:int(sys.argv[4]) if len(sys.argv) > 4 else 7]:
    d = abs(a - b) / max(abs(b), 1.0)
    print(("  ok   " if d < 1e-12 else "  FAIL ") + sys.argv[2] + "isospin.out: " + name + " (rel %.1e)" % d)
    bad += d >= 1e-12
sys.exit(bad)
PYEOF
}
factors() { awk '/^Isospin levels:/ { on = 1; next } on && /factor/ { for (i = 1; i <= NF; i++) if ($i == "factor") printf "%s ", $(i + 1) } /^Filling Data/ { on = 0 }' "$1" | sed 's/ $//'; }
gradline() { grep -m1 "components with relative" "$1"; }
recovered() { awk -v v="$1" 'BEGIN { n = split(v, a, " "); exit !(n == 3 && (a[1]+0.86)^2 < 1e-8 && (a[2]+0.64)^2 < 1e-8 && (a[3]-0.05)^2 < 1e-8) }'; }
fitted() { awk '$1 == "shift_4H.shift" || $1 == "shift_4He.shift" || $1 == "mix_1.amplitude" { printf "%s ", $2 }' "$1" 2>/dev/null; }
mk() {  # extrap, lab Emin, step, file
  awk -v e0="$2" -v de="$3" 'NF >= 5 { e = e0 + (NR - 1) * de; printf "%.6f 0 %.12e %.12e\n", e, $4, 0.001 * $4 }' "$1" > "$WORK/$4"
}
pseudo_data() {  # the four data files from the extrapolations of a truth run
  mk "$1/AZUREOut_aa=2_R=2.extrap" 2.0 0.05 nt.dat
  mk "$1/AZUREOut_aa=3_R=4.extrap" 1.2 0.05 ptn.dat
  mk "$1/AZUREOut_aa=5_R=3.extrap" 0.1 0.02 ddp.dat
  mk "$1/AZUREOut_aa=5_R=4.extrap" 0.1 0.02 ddn.dat
}

echo "AZURE2: $AZURE2_BIN"

# ---- (a) factors and isospin.out ---------------------------------------------------------------
run iso "$M/iso.azr" extrapolate "--no-brune"
log="$WORK/iso/run.log"
want="1 1 0.707107 0.707107 0 0.707107 -0.707107 1"
got="$(factors "$log")"
if [ "$got" = "$want" ]; then ok "start-up factors are Hale's A=4 relations ($got)"; else bad "start-up factors: got '$got', want '$want'"; grep -i -m3 error "$log"; fi
iso="$WORK/iso/output/isospin.out"
if [ -s "$iso" ]; then
  if relations "$iso" "" formal; then nok=$((nok + 7)); else nbad=$((nbad + 1)); fi
else
  bad "iso.azr wrote no isospin.out"; tail -5 "$log"
fi
grep -q "ISOSPIN LEVEL 1 T=1 (derived)" "$WORK/iso/output/parameters.out" && ok "parameters.out labels the members" || bad "parameters.out has no isospin labels"

# ---- (b) the 4H member == a standalone 4H project with its formal parameters --------------------
if [ -s "$iso" ]; then
  read -r e4h g4h < <(member "$iso")
  standalone "$e4h" "$g4h" > "$WORK/standalone.azr"
  run alone "$WORK/standalone.azr" extrapolate "--no-brune --no-transform"
  a="$WORK/iso/output/AZUREOut_aa=2_R=2.extrap"; b="$WORK/alone/output/AZUREOut_aa=2_R=2.extrap"
  if [ -s "$a" ] && [ -s "$b" ]; then
    xs "$a" > "$WORK/a.xs"; xs "$b" > "$WORK/b.xs"
    r="$(maxrel "$WORK/a.xs" "$WORK/b.xs")"
    if le "$r" 1e-10; then ok "4H member == standalone 4H from its formal parameters (max rel $r, $(wc -l < "$WORK/a.xs") points)"; else bad "4H member vs standalone 4H: max rel $r"; fi
  else
    bad "standalone 4H comparison: an extrapolation is missing"; tail -5 "$WORK/alone/run.log"
  fi
fi

# ---- (c) refused models --------------------------------------------------------------------------
refused() {  # NAME MODEL COMMAND FLAGS PATTERN
  run "$1" "$2" "$3" "$4"
  if grep -q "$5" "$WORK/$1/run.log" && [ ! -s "$WORK/$1/output/chiSquared.out" ] && ! ls "$WORK/$1/output/"*.extrap >/dev/null 2>&1; then
    ok "refused: $1"
  else
    bad "not refused as expected: $1 (want '$5')"; grep -i -m3 "error" "$WORK/$1/run.log"
  fi
}
refused undeclared "$M/bad_undeclared.azr" extrapolate "--no-brune" "has channels of two compound nuclei"
refused radius "$M/bad_radius.azr" extrapolate "--no-brune" "need equal radii"
refused mix_allowed "$M/bad_mix.azr" extrapolate "--no-brune" "isospin-allowed"
refused T_too_small "$M/bad_T.azr" extrapolate "--no-brune" "has T=0 but lines in"
refused park "$M/iso.azr" extrapolate "--use-park" "not available in Park's basis"

# ---- (d) recovery fit ----------------------------------------------------------------------------
run truth "$M/iso_mix.azr" extrapolate "--no-brune"
T="$WORK/truth/output"
if [ -s "$T/AZUREOut_aa=2_R=2.extrap" ] && [ -s "$T/AZUREOut_aa=5_R=3.extrap" ]; then
  pseudo_data "$T"
  for fl in "" "--use-gradient" "--use-lm"; do
    tag="fit${fl//-/_}"
    run "$tag" "$M/iso_fit.azr" fit "--no-brune $fl"
    v="$(fitted "$WORK/$tag/output/param.sav")"
    if recovered "$v"; then
      ok "fit ${fl:-(MIGRAD)} recovers shift(4H) -0.86, shift(4He) -0.64, mix 0.05 from -0.5 / -0.4 / 0 ($v)"
    else
      bad "fit ${fl:-(MIGRAD)} gave '$v' (want -0.86 -0.64 0.05)"; tail -5 "$WORK/$tag/run.log"
    fi
  done
  sav="$WORK/fit/output/param.sav"
  if [ -s "$sav" ]; then
    # The derived slots are written synced to the fitted shift and mix (param.sav keeps 8 digits).
    if awk '{ v[$1] = $2 } END { d1 = v["energy_2"] - (v["energy_1"] + v["shift_4H.shift"]);
             d2 = v["energy_3"] - (v["energy_1"] + 19.8136 + v["shift_4He.shift"]); d3 = v["width_3_3"] - v["mix_1.amplitude"];
             exit !(("energy_2" in v) && ("mix_1.amplitude" in v) && d1 * d1 < 1e-12 && d2 * d2 < 1e-10 && d3 * d3 < 1e-14) }' "$sav"; then
      ok "param.sav: derived energies and the mixed d+d width are synced to the fitted shift and mix"
    else
      bad "param.sav: derived slots not synced"; grep -E "^ *(energy_[123]|width_3_3|shift_|mix_)" "$sav"
    fi
  else
    bad "recovery fit wrote no param.sav"; tail -8 "$WORK/fit/run.log"
  fi
else
  bad "truth model produced no extrapolations"; tail -5 "$WORK/truth/run.log"
fi

# ---- (e) gradient check through the ties ----------------------------------------------------------
HARD="energy 0/2 gamma 0/3 isospin_shift 0/2 isospin_mix 0/1; Jacobian path: energy 0/2 gamma 0/3 isospin_shift 0/2 isospin_mix 0/1"
SOFT="energy 0/4 gamma 0/7 isospin_shift 0/2; Jacobian path: energy 0/4 gamma 0/7 isospin_shift 0/2"
gradcheck() {  # NAME MODEL FLAGS PATTERN WORDS
  run "$1" "$2" fit "$3 --use-gradient" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50
  local line; line="$(gradline "$WORK/$1/run.log")"
  if grep -q "$4" <<< "$line"; then ok "$5"; else bad "$5: '$line'"; fi
}
# off the minimum: the soft model from shifts -0.8 / -0.6
sed 's/^\( *shift  *\)-0.86 MeV/\1-0.8 MeV/; s/^\( *shift  *\)-0.64 MeV/\1-0.6 MeV/' "$M/iso_soft.azr" > "$WORK/iso_soft_off.azr"
if [ -s "$WORK/nt.dat" ]; then
  gradcheck gc "$M/iso_free.azr" "--no-brune" "$HARD" "gradient check: reference energies, owner widths, shifts and mix agree with finite differences on both paths"
fi

# ---- (f) the ensemble sampler ---------------------------------------------------------------------
if [ -s "$WORK/nt.dat" ]; then
  help="$("$AZURE2_BIN" --help < /dev/null 2>&1 | tr -d '\r')"
  if grep -q "^[[:space:]]*mcmc:" <<< "$help"; then
    MC="--no-brune --walkers=8 --steps=100 --parameter-spread=1 --energy-spread=0.5 --sampler-threads=2 --overwrite-samples"
    run mc "$M/iso_free.azr" mcmc "$MC --basis=rwa"
    hdr="$(head -1 "$WORK/mc/output/samples.mcmc" 2>/dev/null | tr -d '\r')"
    if [ "$hdr" = "step,walker,logprob,loglikelihood,logprior,param0,param1,param2,param3,param4,param5,param6,param7" ] &&
       grep -a -q "3 numbers of gamma_decay or isospin records" "$WORK/mc/run.log"; then
      ok "sampler (fit's basis): 8 free parameters sampled, the isospin priors in the likelihood ($(($(wc -l < "$WORK/mc/output/samples.mcmc") - 1)) samples)"
    else
      bad "sampler (fit's basis): header '$hdr'"; grep -a -i -m3 "error\|isospin" "$WORK/mc/run.log"
    fi
    run mcp "$M/iso_free.azr" mcmc "$MC --basis=physical"
    if grep -a -q "switch the MCMC parameter type to RWA" "$WORK/mcp/run.log" && [ ! -s "$WORK/mcp/output/samples.mcmc" ]; then
      ok "refused: the sampler on observed parameters"
    else
      bad "the sampler on observed parameters was not refused"
    fi
  else
    echo "  (no sampler in this build: mcmc not checked)"
  fi
fi

# ---- (g) soft tie --------------------------------------------------------------------------------
run soft "$M/iso_soft.azr" extrapolate "--no-brune"
n=0; for f in "$WORK/iso/output/"*.extrap; do same_file "$f" "$WORK/soft/output/$(basename "$f")" || n=$((n + 1)); done
if [ "$n" = 0 ] && [ -n "$(ls "$WORK/soft/output/"*.extrap 2>/dev/null)" ]; then ok "soft tie at its start == hard tie (every extrapolation identical)"
else bad "soft tie at its start differs from the hard tie in $n file(s)"; fi
grep -q "SOFT (sigma E 0.05" "$WORK/soft/run.log" && grep -q "soft" "$WORK/soft/output/isospin.out" && ok "soft tie reported (start-up log, isospin.out)" || bad "soft tie not reported"
if [ -s "$WORK/nt.dat" ]; then
  gradcheck gcs "$WORK/iso_soft_off.azr" "--no-brune" "$SOFT" "soft-tie gradient check: 13 free parameters (soft slots included) agree with finite differences on both paths"
  run softfit "$M/iso_soft.azr" fit "--no-brune --use-gradient"
  if grep -q "Total-Isospin-Chi-Squared" "$WORK/softfit/output/chiSquared.out" 2>/dev/null && grep -q "soft-tie penalty" "$WORK/softfit/output/isospin.out"; then
    ok "soft-tie fit: $(grep -o 'Total-Isospin-Chi-Squared: [^ ]*' "$WORK/softfit/output/chiSquared.out"), $(grep 'soft-tie penalty' "$WORK/softfit/output/isospin.out")"
  else
    bad "soft-tie fit did not report the soft penalty"; tail -5 "$WORK/softfit/run.log"
  fi
  # The least-squares fitter takes the soft ties as residual rows (dev refused the combination): the same penalty.
  run softlm "$M/iso_soft.azr" fit "--no-brune --use-lm"
  p1="$(grep -o 'Total-Isospin-Chi-Squared: [^ ]*' "$WORK/softfit/output/chiSquared.out" 2>/dev/null | awk '{ print $2 }')"
  p2="$(grep -o 'Total-Isospin-Chi-Squared: [^ ]*' "$WORK/softlm/output/chiSquared.out" 2>/dev/null | awk '{ print $2 }')"
  if awk -v a="$p1" -v b="$p2" 'BEGIN { exit !(a ~ /^ *[-+]?\.?[0-9]/ && b ~ /^ *[-+]?\.?[0-9]/ && (a - b)^2 <= (1e-3 * a)^2) }'; then
    ok "soft tie with the least-squares fitter: the same soft penalty as MIGRAD ($p2)"
  else
    bad "soft tie with --use-lm: penalty '$p2', MIGRAD '$p1'"; tail -5 "$WORK/softlm/run.log"
  fi
fi

# ---- (h) relativistic kinematics (relativistic=1) through the ties ------------------------------------
# Level energies stay excitation energies under relativistic=1, so the ties are unchanged.
rel "$M/iso.azr" "$WORK/iso_rel.azr"
run isorel "$WORK/iso_rel.azr" extrapolate "--no-brune"
isor="$WORK/isorel/output/isospin.out"
if [ -s "$isor" ] && [ -s "$WORK/isorel/output/relativistic.out" ]; then
  if relations "$isor" "relativistic " formal 2; then nok=$((nok + 2)); else nbad=$((nbad + 1)); fi
  read -r e4h g4h < <(member "$isor")
  standalone "$e4h" "$g4h" > "$WORK/standalone_abs.azr"
  rel "$WORK/standalone_abs.azr" "$WORK/standalone_rel.azr"
  run alonerel "$WORK/standalone_rel.azr" extrapolate "--no-brune --no-transform"
  a="$WORK/isorel/output/AZUREOut_aa=2_R=2.extrap"; b="$WORK/alonerel/output/AZUREOut_aa=2_R=2.extrap"
  if [ -s "$a" ] && [ -s "$b" ] && [ -s "$WORK/a.xs" ]; then
    xs "$a" > "$WORK/ar.xs"; xs "$b" > "$WORK/br.xs"
    r="$(maxrel "$WORK/ar.xs" "$WORK/br.xs")"
    if le "$r" 1e-10; then ok "relativistic: 4H member == standalone relativistic 4H (max rel $r)"; else bad "relativistic: 4H member vs standalone: max rel $r"; fi
    r="$(maxrel "$WORK/ar.xs" "$WORK/a.xs")"
    if le "$r" 1e-7; then bad "relativistic=1 left the tied n+t cross section unchanged (max rel $r)"; else ok "relativistic=1 changes the tied n+t cross section (max rel $r against the non-relativistic run)"; fi
  else
    bad "relativistic: an extrapolation is missing"; tail -5 "$WORK/isorel/run.log" "$WORK/alonerel/run.log"
  fi
else
  bad "relativistic iso.azr wrote no isospin.out or relativistic.out"; tail -5 "$WORK/isorel/run.log"
fi
if [ -s "$WORK/nt.dat" ]; then
  rel "$M/iso_free.azr" "$WORK/iso_free_rel.azr"
  gradcheck gcr "$WORK/iso_free_rel.azr" "--no-brune" "$HARD" "relativistic gradient check: energies, widths, shifts and mix agree with finite differences on both paths"
fi

# ---- (i) Brune basis: ties on observed energies and Brune amplitudes -----------------------------------
run isoB "$M/iso.azr" extrapolate "--use-brune"
logB="$WORK/isoB/run.log"; isoB="$WORK/isoB/output/isospin.out"
got="$(factors "$logB")"
if [ "$got" = "$want" ] && grep -q "^Isospin levels:.*Brune basis" "$logB"; then ok "Brune: start-up factors are Hale's A=4 relations, the report names the Brune ties"
else bad "Brune start-up report: factors '$got'"; grep -i -m3 "error\|^Isospin levels" "$logB"; fi
if [ -s "$isoB" ]; then
  if relations "$isoB" "Brune " observed; then nok=$((nok + 7)); else nbad=$((nbad + 1)); fi
  if grep -q "Brune parameters" "$isoB" && ! grep -q "B=" "$isoB" && ! grep -q "formal E" "$isoB"; then ok "Brune isospin.out: Brune header, observed energies, no boundary values"
  else bad "Brune isospin.out: header or boundary values"; fi
  read -r e4h g4h < <(member "$isoB")
  standalone "$e4h" "$g4h" > "$WORK/standaloneB.azr"
  run aloneB "$WORK/standaloneB.azr" extrapolate "--use-brune --no-transform"
  a="$WORK/isoB/output/AZUREOut_aa=2_R=2.extrap"; b="$WORK/aloneB/output/AZUREOut_aa=2_R=2.extrap"
  if [ -s "$a" ] && [ -s "$b" ]; then
    xs "$a" > "$WORK/aB.xs"; xs "$b" > "$WORK/bB.xs"
    r="$(maxrel "$WORK/aB.xs" "$WORK/bB.xs")"
    if le "$r" 1e-10; then ok "Brune: 4H member == standalone 4H from its Brune parameters (max rel $r, $(wc -l < "$WORK/aB.xs") points)"; else bad "Brune: 4H member vs standalone 4H: max rel $r"; fi
    # The reference nucleus has one level in its group and the same physical input in both runs: same p+3He
    # cross section.  The tied 4H level sits elsewhere, because the shift now relates observed energies.
    if [ -s "$WORK/a.xs" ] && [ -s "$WORK/iso/output/AZUREOut_aa=1_R=1.extrap" ] && [ -s "$WORK/isoB/output/AZUREOut_aa=1_R=1.extrap" ]; then
      xs "$WORK/iso/output/AZUREOut_aa=1_R=1.extrap" > "$WORK/pF.xs"; xs "$WORK/isoB/output/AZUREOut_aa=1_R=1.extrap" > "$WORK/pB.xs"
      r="$(maxrel "$WORK/pB.xs" "$WORK/pF.xs")"
      if le "$r" 1e-6; then ok "Brune vs formal: same p+3He (reference nucleus, one level) cross section (max rel $r)"; else bad "Brune vs formal: p+3He differs (max rel $r)"; fi
      r="$(maxrel "$WORK/aB.xs" "$WORK/a.xs")"
      if le "$r" 1e-4; then bad "Brune vs formal: tied n+t cross section unchanged (max rel $r)"; else ok "Brune vs formal: the tied n+t cross section differs (max rel $r): the shift relates observed energies"; fi
    else
      bad "Brune vs formal: an extrapolation is missing"
    fi
  else
    bad "Brune standalone 4H comparison: an extrapolation is missing"; tail -5 "$WORK/aloneB/run.log"
  fi
else
  bad "iso.azr (Brune) wrote no isospin.out"; tail -5 "$logB"
fi
run softB "$M/iso_soft.azr" extrapolate "--use-brune"
n=0; for f in "$WORK/isoB/output/"*.extrap; do same_file "$f" "$WORK/softB/output/$(basename "$f")" || n=$((n + 1)); done
if [ "$n" = 0 ] && [ -n "$(ls "$WORK/softB/output/"*.extrap 2>/dev/null)" ]; then ok "Brune: soft tie at its start == hard tie (every extrapolation identical)"
else bad "Brune: soft tie at its start differs from the hard tie in $n file(s)"; fi
if [ -s "$WORK/nt.dat" ]; then
  # Gradient check on the formal pseudo-data (off the Brune model's minimum), hard and soft ties.
  gradcheck gcB "$M/iso_free.azr" "--use-brune" "$HARD" "Brune gradient check: energies, widths, shifts and mix agree with finite differences on both paths"
  gradcheck gcsB "$WORK/iso_soft_off.azr" "--use-brune" "$SOFT" "Brune soft-tie gradient check: 13 free parameters agree with finite differences on both paths"
fi
# Recovery fit on Brune pseudo-data (this replaces the data files, so it comes last).  The 4H shift starts at
# -0.9, not at iso_fit.azr's -0.5: the n+t pseudo-data pass through an interference zero (1e-6 of the peak, near
# 3.07 MeV c.m.) and, with 0.1 % relative errors, the points beside it give chi2 a second minimum at a shift of
# -0.8564, which MIGRAD reaches from above.  That is a property of the pseudo-data, not of the ties.
run truthB "$M/iso_mix.azr" extrapolate "--use-brune"
T="$WORK/truthB/output"
if [ -s "$T/AZUREOut_aa=2_R=2.extrap" ] && [ -s "$T/AZUREOut_aa=5_R=3.extrap" ]; then
  pseudo_data "$T"
  sed 's/^\( *shift  *\)-0.5 MeV/\1-0.9 MeV/' "$M/iso_fit.azr" > "$WORK/iso_fitB.azr"
  run fitB "$WORK/iso_fitB.azr" fit "--use-brune --use-gradient"
  v="$(fitted "$WORK/fitB/output/param.sav")"
  if recovered "$v"; then
    ok "Brune fit (MIGRAD, analytic gradient) recovers shift(4H) -0.86, shift(4He) -0.64, mix 0.05 from -0.9 / -0.4 / 0 ($v)"
  else
    bad "Brune fit gave '$v' (want -0.86 -0.64 0.05)"; tail -5 "$WORK/fitB/run.log"
  fi
else
  bad "Brune truth model produced no extrapolations"; tail -5 "$WORK/truthB/run.log"
fi

# ---- (j) photon channels in a tied J-group ------------------------------------------------------
for fl in "--use-brune" "--no-brune"; do
  run "ecA$fl" "$M/ecA.azr" extrapolate "$fl"; run "ecAiso$fl" "$M/ecA_iso.azr" extrapolate "$fl"
  a="$WORK/ecA$fl/output/AZUREOut_aa=1_R=2.extrap"; b="$WORK/ecAiso$fl/output/AZUREOut_aa=1_R=2.extrap"
  if [ -s "$a" ] && [ -s "$b" ]; then
    xsec "$a" > "$WORK/a.xs"; xsec "$b" > "$WORK/b.xs"; r="$(maxrel "$WORK/b.xs" "$WORK/a.xs")"
    if le "$r" 1e-10; then ok "capture, $fl: every level in a one-nucleus block == no block (max rel $r, $(wc -l < "$WORK/a.xs") points)"; else bad "capture, $fl: one-nucleus block changes 12C(p,g0) (max rel $r)"; fi
  else
    bad "capture, $fl: an extrapolation is missing"; grep -i -m3 "error" "$WORK/ecAiso$fl/run.log"
  fi
done
run ecM "$M/ecM.azr" extrapolate "--use-brune"; run ecMtied "$M/ecM_tied.azr" extrapolate "--use-brune"
a="$WORK/ecA--use-brune/output/AZUREOut_aa=1_R=2.extrap"
m1="$WORK/ecMtied/output/AZUREOut_aa=1_R=2.extrap"; m2="$WORK/ecMtied/output/AZUREOut_aa=3_R=4.extrap"; u2="$WORK/ecM/output/AZUREOut_aa=3_R=4.extrap"
if [ -s "$a" ] && [ -s "$m1" ] && [ -s "$m2" ] && [ -s "$u2" ]; then
  if grep -q "photon/beta channel(s), not tied" "$WORK/ecMtied/run.log"; then ok "tied 13N + 13C capture model: start-up report lists the photon channels as not tied"; else bad "tied capture model: no 'not tied' line in the start-up report"; fi
  xsec "$a" > "$WORK/a.xs"; xsec "$m1" > "$WORK/m1.xs"; xsec "$m2" > "$WORK/m2.xs"; xsec "$u2" > "$WORK/u2.xs"
  r="$(maxrel "$WORK/m1.xs" "$WORK/a.xs")"
  if le "$r" 1e-10; then ok "tied capture model: 12C(p,g0) of the reference nucleus == ecA.azr (max rel $r)"; else bad "tied capture model: reference-nucleus capture vs ecA.azr max rel $r"; fi
  r="$(maxrel "$WORK/m2.xs" "$WORK/u2.xs")"
  if le "$r" 1e-3; then bad "tied capture model: 12C(n,g0) is the untied curve (max rel $r); the tie did nothing"; else ok "tied capture model: 12C(n,g0) follows the tie (differs from the untied model, max rel $r)"; fi
else
  bad "tied capture model: an extrapolation is missing"; grep -i -m3 "error" "$WORK/ecMtied/run.log"
fi

# ---- (k) excited particles are not the ground-state arrangement ----------------------------------
lvl1() { awk '/^ *level ID 1 /{on = 1; next} /^ *level ID /{on = 0} on && /^ *pair /' "$1"; }
run ecExc1 "$M/ecM_exc1.azr" extrapolate "--use-brune"; run ecExc2 "$M/ecM_exc2.azr" extrapolate "--use-brune"
r1="$(lvl1 "$WORK/ecExc1/run.log")"; r2="$(lvl1 "$WORK/ecExc2/run.log")"
if [ -z "$r1" ] || [ -z "$r2" ]; then
  bad "excited-particle models: no start-up report for level 1"; grep -i -m3 "error" "$WORK/ecExc1/run.log" "$WORK/ecExc2/run.log"
else
  if grep -qE "^ *pair +5 .* free$" <<< "$r1"; then ok "n + 12C*(7.654, 0+) alone in 13C: free, not tied to the 13N ground-state channel"; else bad "n + 12C*(7.654) in 13C: $(echo "$r1" | grep -E '^ *pair +5')"; fi
  if grep -q "0+#2" "$WORK/ecExc1/output/isospin.out" 2>/dev/null; then ok "isospin.out key of the excited-12C channel carries the state rank (0+#2)"; else bad "isospin.out has no 0+#2 key"; fi
  if grep -qE "^ *pair +5 .*= 1 x free amplitude" <<< "$r2" && grep -qE "^ *pair +6 .* free$" <<< "$r2"; then ok "with p + 12C*(7.654) in 13N: n + 12C*(7.654) is tied to it"; else bad "excited-12C pair in both nuclei: $(echo "$r2" | grep -E '^ *pair +[56]')"; fi
  if grep -qE "^ *pair +3 .*= 1 x free amplitude" <<< "$r2"; then ok "the ground-state pair 3 is still tied to pair 1"; else bad "pair 3 lost its tie: $(echo "$r2" | grep -E '^ *pair +3')"; fi
fi

echo "isospin_levels: $nok ok, $nbad failed"
[ "$nbad" -eq 0 ]
