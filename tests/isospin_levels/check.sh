#!/usr/bin/env bash
#
# Isospin-labelled levels across compound nuclei (the <isospin> block).
#
#   tests/isospin_levels/check.sh [path/to/AZURE2]
#
# (a) The factors printed at start-up are Hale's A=4 relations, and isospin.out obeys them:
#     derived energies = reference + threshold difference + shift, derived amplitudes =
#     factor ratio x the free one, the T=1 level's d+d amplitude 0, the T=0 level's
#     n+3He amplitude = -(p+t).
# (b) The 4H member equals a standalone 4H project built from its formal parameters
#     (--no-transform) to 1e-10 in the n+t cross section.
# (c) Refused models: a cross-nucleus level without the block, unequal radii in one
#     arrangement, mix on an allowed channel, T too small for a nucleus, --use-brune.
# (d) Pseudo-data from iso_mix.azr (shifts -0.86 / -0.64, mix 0.05): MIGRAD (numerical and
#     analytic gradient), LM and GSL least squares from -0.5 / -0.4 / 0 recover all three.
# (e) AZURE_GRAD_CHECK on iso_free.azr (every reference parameter free, priors): the
#     analytic gradient and the Jacobian path agree with finite differences through the ties.
# (f) MCMC: a short RWA run samples the 8 free parameters; the physical-parameter run is refused.
# (g) Soft tie (iso_soft.azr): at the start (soft slots on their tie targets) every extrapolation
#     equals the hard-tie model's; the gradient check (from shifts -0.8 / -0.6, off the minimum)
#     agrees for all 13 free parameters, soft slots included; LM is refused; a MIGRAD fit reports
#     the soft penalty in Total-Isospin-Chi-Squared and isospin.out.
#
# The models are written by make_models.py (committed).

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"
AZURE2_BIN="${1:-}"
if [ -z "$AZURE2_BIN" ]; then
  for c in "$REPO_ROOT"/build*/src/AZURE2; do [ -x "$c" ] && { [ -z "$AZURE2_BIN" ] || [ "$c" -nt "$AZURE2_BIN" ]; } && AZURE2_BIN="$c"; done
fi
[ -x "$AZURE2_BIN" ] || { echo "ERROR: AZURE2 binary not found" >&2; exit 1; }
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
command -v python3 >/dev/null 2>&1 || { echo "SKIP isospin_levels: python3 not found"; exit 0; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/isospin_levels.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
M="$HERE/models"
nok=0; nbad=0
ok()  { echo "  ok   $1"; nok=$((nok + 1)); }
bad() { echo "  FAIL $1"; nbad=$((nbad + 1)); }
run_limited() { local t=$1; shift; if command -v timeout >/dev/null 2>&1; then timeout "$t" "$@"; else "$@"; fi; }
# run NAME MODEL MODE FLAGS: a fresh directory per run
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$2" "$dir/t.azr"
  for f in nt.dat ptn.dat ddp.dat ddn.dat; do [ -f "$WORK/$f" ] && cp "$WORK/$f" "$dir/"; done
  (cd "$dir" && printf "$3\n\n\n7\n" | run_limited 900 "$AZURE2_BIN" $4 --no-gui --no-readline t.azr > run.log 2>&1)
}
xs() { awk 'NF >= 5 { printf "%.12e %.12e\n", $1, $4 }' "$1"; }
maxrel() {
  paste "$1" "$2" | awk '
    NF >= 4 { d = $2 - $4; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300;
              r = d / m; if (r > worst) worst = r; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t + 0) }'; }

echo "AZURE2: $AZURE2_BIN"

# ---- (a) factors and isospin.out ---------------------------------------------------------------
run iso "$M/iso.azr" 3 "--no-brune"
log="$WORK/iso/run.log"
want="1 1 0.707107 0.707107 0 0.707107 -0.707107 1"
got="$(awk '/^Isospin levels:/ { on = 1; next } on && /factor/ { for (i = 1; i <= NF; i++) if ($i == "factor") printf "%s ", $(i + 1) } /^Filling Data/ { on = 0 }' "$log" | sed 's/ $//')"
if [ "$got" = "$want" ]; then ok "start-up factors are Hale's A=4 relations ($got)"; else bad "start-up factors: got '$got', want '$want'"; grep -i error "$log" | head -3; fi
iso="$WORK/iso/output/isospin.out"
if [ -s "$iso" ]; then
  if python3 - "$iso" <<'EOF'
import sys, math
lv, cur, mem = {}, None, None
for raw in open(sys.argv[1]):
    t = raw.split()
    if raw.startswith("level ID"):
        cur = int(t[2]); lv[cur] = {}
    elif cur is not None and "formal E =" in raw:
        mem = t[0] + t[1]; lv[cur][mem] = {"E": float(t[t.index("E") + 2]), "g": {}}
    elif cur is not None and raw.strip().startswith("pair"):
        lv[cur][mem]["g"][int(t[1])] = float(t[t.index("gamma") + 2])
li, h, he = "Z=3A=4", "Z=1A=4", "Z=2A=4"
r = math.sqrt(0.5)
checks = [
  ("E(4H)  = E(4Li) - 0.86", lv[1][h]["E"], lv[1][li]["E"] - 0.86),
  ("E(4He) = E(4Li) + 19.8136 - 0.64", lv[1][he]["E"], lv[1][li]["E"] + 19.8136 - 0.64),
  ("gamma(n+t) = gamma(p+3He)", lv[1][h]["g"][2], lv[1][li]["g"][1]),
  ("gamma(p+t, T=1) = gamma(p+3He)/sqrt2", lv[1][he]["g"][3], r * lv[1][li]["g"][1]),
  ("gamma(n+3He, T=1) = gamma(p+3He)/sqrt2", lv[1][he]["g"][4], r * lv[1][li]["g"][1]),
  ("gamma(d+d, T=1) = 0", lv[1][he]["g"][5], 0.0),
  ("gamma(n+3He, T=0) = -gamma(p+t, T=0)", lv[2][he]["g"][4], -lv[2][he]["g"][3]),
]
bad = 0
for name, a, b in checks:
    d = abs(a - b) / max(abs(b), 1.0)
    print(("  ok   " if d < 1e-12 else "  FAIL ") + "isospin.out: " + name + " (rel %.1e)" % d)
    bad += d >= 1e-12
sys.exit(bad)
EOF
  then nok=$((nok + 7)); else nbad=$((nbad + 1)); fi
else
  bad "iso.azr wrote no isospin.out"; tail -5 "$log"
fi
grep -q "ISOSPIN LEVEL 1 T=1 (derived)" "$WORK/iso/output/parameters.out" && ok "parameters.out labels the members" || bad "parameters.out has no isospin labels"

# ---- (b) the 4H member == a standalone 4H project with its formal parameters --------------------
if [ -s "$iso" ]; then
  read -r e4h g4h < <(awk '/^level ID 1/ { on = 1 } on && /Z=1 A=4/ { for (i = 1; i <= NF; i++) if ($i == "E") e = $(i + 2); getline; for (i = 1; i <= NF; i++) if ($i == "gamma") g = $(i + 2); print e, g; exit }' "$iso")
  awk -v e="$e4h" -v g="$g4h" '/<levels>/ { print; getline; $3 = e; $12 = g; print; next } { print }' "$M/standalone_4H.azr" > "$WORK/standalone.azr"
  run alone "$WORK/standalone.azr" 3 "--no-brune --no-transform"
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
refused() {  # NAME MODEL FLAGS MODE PATTERN
  run "$1" "$2" "$4" "$3"
  if grep -q "$5" "$WORK/$1/run.log" && [ ! -s "$WORK/$1/output/chiSquared.out" ] && ! ls "$WORK/$1/output/"*.extrap >/dev/null 2>&1; then
    ok "refused: $1"
  else
    bad "not refused as expected: $1 (want '$5')"; grep -i -m3 "error" "$WORK/$1/run.log"
  fi
}
refused undeclared "$M/bad_undeclared.azr" "--no-brune" 3 "has channels of two compound nuclei"
refused radius "$M/bad_radius.azr" "--no-brune" 3 "need equal radii"
refused mix_allowed "$M/bad_mix.azr" "--no-brune" 3 "isospin-allowed"
refused T_too_small "$M/bad_T.azr" "--no-brune" 3 "has T=0 but lines in"
refused brune "$M/iso.azr" "--use-brune" 3 "run with --no-brune"

# ---- (d) recovery fit ----------------------------------------------------------------------------
run truth "$M/iso_mix.azr" 3 "--no-brune"
T="$WORK/truth/output"
mk() {  # extrap, lab Emin, step, file
  awk -v e0="$2" -v de="$3" 'NF >= 5 { e = e0 + (NR - 1) * de; printf "%.6f 0 %.12e %.12e\n", e, $4, 0.001 * $4 }' "$1" > "$WORK/$4"
}
if [ -s "$T/AZUREOut_aa=2_R=2.extrap" ] && [ -s "$T/AZUREOut_aa=5_R=3.extrap" ]; then
  mk "$T/AZUREOut_aa=2_R=2.extrap" 2.0 0.05 nt.dat
  mk "$T/AZUREOut_aa=3_R=4.extrap" 1.2 0.05 ptn.dat
  mk "$T/AZUREOut_aa=5_R=3.extrap" 0.1 0.02 ddp.dat
  mk "$T/AZUREOut_aa=5_R=4.extrap" 0.1 0.02 ddn.dat
  for fl in "" "--use-gradient" "--use-lm" "--use-gsl-lm"; do
    tag="fit${fl//-/_}"
    run "$tag" "$M/iso_fit.azr" 2 "--no-brune $fl"
    s="$WORK/$tag/output/param.sav"
    v="$(awk '$1 == "isospin_shift_1_4" || $1 == "isospin_shift_2_4" || $1 == "isospin_mix_1" { printf "%s ", $2 }' "$s" 2>/dev/null)"
    if awk -v v="$v" 'BEGIN { n = split(v, a, " "); exit !(n == 3 && (a[1]+0.86)^2 < 1e-8 && (a[2]+0.64)^2 < 1e-8 && (a[3]-0.05)^2 < 1e-8) }'; then
      ok "fit ${fl:-(MIGRAD)} recovers shift(4H) -0.86, shift(4He) -0.64, mix 0.05 from -0.5 / -0.4 / 0 ($v)"
    else
      bad "fit ${fl:-(MIGRAD)} gave '$v' (want -0.86 -0.64 0.05)"; tail -5 "$WORK/$tag/run.log"
    fi
  done
  sav="$WORK/fit/output/param.sav"
  if [ -s "$sav" ]; then
    # The derived slots are written synced to the fitted shift and mix (param.sav keeps 8 digits).
    if awk '{ v[$1] = $2 } END { d1 = v["energy_2"] - (v["energy_1"] + v["isospin_shift_1_4"]);
             d2 = v["energy_3"] - (v["energy_1"] + 19.8136 + v["isospin_shift_2_4"]); d3 = v["width_3_3"] - v["isospin_mix_1"];
             exit !(d1 * d1 < 1e-12 && d2 * d2 < 1e-10 && d3 * d3 < 1e-14) }' "$sav"; then
      ok "param.sav: derived energies and the mixed d+d width are synced to the fitted shift and mix"
    else
      bad "param.sav: derived slots not synced"; grep -E "^ *(energy_[123]|width_3_3|isospin)" "$sav"
    fi
  else
    bad "recovery fit wrote no param.sav"; tail -8 "$WORK/fit/run.log"
  fi
else
  bad "truth model produced no extrapolations"; tail -5 "$WORK/truth/run.log"
fi

# ---- (e) gradient check through the ties ----------------------------------------------------------
if [ -s "$WORK/nt.dat" ]; then
  mkdir -p "$WORK/gc/output" "$WORK/gc/checks"; cp "$WORK"/*.dat "$WORK/gc/"; cp "$M/iso_free.azr" "$WORK/gc/t.azr"
  (cd "$WORK/gc" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50 run_limited 900 "$AZURE2_BIN" --no-brune --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
  line="$(grep "components with relative" "$WORK/gc/run.log" | head -1)"
  if echo "$line" | grep -q "energy 0/2 gamma 0/3 isospin 0/3; Jacobian path: energy 0/2 gamma 0/3 isospin 0/3"; then
    ok "gradient check: reference energies, owner widths, shifts and mix agree with finite differences on both paths"
  else
    bad "gradient check: '$line'"
  fi
fi

# ---- (f) MCMC ----------------------------------------------------------------------------------------
if [ -s "$WORK/nt.dat" ]; then
  mkdir -p "$WORK/mc/output" "$WORK/mc/checks"; cp "$WORK"/*.dat "$WORK/mc/"; cp "$M/iso_free.azr" "$WORK/mc/t.azr"
  (cd "$WORK/mc" && printf '6\n\n8\n100\n1\n0.5\n2\nyes\nyes\n' | run_limited 900 "$AZURE2_BIN" --no-brune --no-gui --no-readline t.azr > run.log 2>&1)
  hdr="$(head -1 "$WORK/mc/output/samples.mcmc" 2>/dev/null)"
  if [ "$hdr" = "step,walker,logprob,loglikelihood,logprior,param0,param1,param2,param3,param4,param5,param6,param7" ] &&
     grep -q "3 energy shifts taken from" "$WORK/mc/run.log"; then
    ok "MCMC (RWA): 8 free parameters sampled, isospin priors applied ($(($(wc -l < "$WORK/mc/output/samples.mcmc") - 1)) samples)"
  else
    bad "MCMC (RWA): header '$hdr'"; grep -a -i -m3 "error\|prior" "$WORK/mc/run.log"
  fi
  mkdir -p "$WORK/mcp/output" "$WORK/mcp/checks"; cp "$WORK"/*.dat "$WORK/mcp/"; cp "$M/iso_free.azr" "$WORK/mcp/t.azr"
  (cd "$WORK/mcp" && printf '6\n\n8\n100\n1\n0.5\n2\nno\nyes\n' | run_limited 300 "$AZURE2_BIN" --no-brune --no-gui --no-readline t.azr > run.log 2>&1)
  if grep -a -q "switch the MCMC parameter type to RWA" "$WORK/mcp/run.log" && [ ! -s "$WORK/mcp/output/samples.mcmc" ]; then
    ok "refused: MCMC on physical parameters"
  else
    bad "MCMC on physical parameters was not refused"
  fi
fi

# ---- (g) soft tie --------------------------------------------------------------------------------
run soft "$M/iso_soft.azr" 3 "--no-brune"
n=0; for f in "$WORK/iso/output/"*.extrap; do cmp -s "$f" "$WORK/soft/output/$(basename "$f")" || n=$((n + 1)); done
if [ "$n" = 0 ] && [ -n "$(ls "$WORK/soft/output/"*.extrap 2>/dev/null)" ]; then ok "soft tie at its start == hard tie (every extrapolation identical)"
else bad "soft tie at its start differs from the hard tie in $n file(s)"; fi
grep -q "SOFT (sigma E 0.05" "$WORK/soft/run.log" && grep -q "soft" "$WORK/soft/output/isospin.out" && ok "soft tie reported (start-up log, isospin.out)" || bad "soft tie not reported"
if [ -s "$WORK/nt.dat" ]; then
  mkdir -p "$WORK/gcs/output" "$WORK/gcs/checks"; cp "$WORK"/*.dat "$WORK/gcs/"
  sed 's/^shift 1 4 2 -0.86 1/shift 1 4 2 -0.8 1/; s/^shift 2 4 3 -0.64 1/shift 2 4 3 -0.6 1/' "$M/iso_soft.azr" > "$WORK/gcs/t.azr"
  (cd "$WORK/gcs" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50 run_limited 900 "$AZURE2_BIN" --no-brune --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
  line="$(grep "components with relative" "$WORK/gcs/run.log" | head -1)"
  if echo "$line" | grep -q "energy 0/4 gamma 0/7 isospin 0/2; Jacobian path: energy 0/4 gamma 0/7 isospin 0/2"; then
    ok "soft-tie gradient check: 13 free parameters (soft slots included) agree with finite differences on both paths"
  else
    bad "soft-tie gradient check: '$line'"
  fi
  refused soft_lm "$M/iso_soft.azr" "--no-brune --use-lm" 2 "a soft tie couples parameters"
  run softfit "$M/iso_soft.azr" 2 "--no-brune --use-gradient"
  if grep -q "Total-Isospin-Chi-Squared" "$WORK/softfit/output/chiSquared.out" 2>/dev/null && grep -q "soft-tie penalty" "$WORK/softfit/output/isospin.out"; then
    ok "soft-tie fit: $(grep -o 'Total-Isospin-Chi-Squared: [^ ]*' "$WORK/softfit/output/chiSquared.out"), $(grep 'soft-tie penalty' "$WORK/softfit/output/isospin.out")"
  else
    bad "soft-tie fit did not report the soft penalty"; tail -5 "$WORK/softfit/run.log"
  fi
fi

echo "isospin_levels: $nok ok, $nbad failed"
[ "$nbad" -eq 0 ]
