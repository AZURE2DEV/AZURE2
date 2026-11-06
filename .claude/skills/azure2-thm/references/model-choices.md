# Baseline model choices in detail (Stage 3)

Every key below is documented in `docs/source/theory/thm_implementation.rst`, cited here
as *impl*. Global keys go in `<thm>`; per-experiment keys go on an
`experiment[<name>]` line. Edit them in the GUI (THM Workspace) or with
`AzrModel.set_thm_option` / `set_thm_experiment`. Both apply the engine's rules and
raise `ValueError` on a refusal.

## Global options

**`entranceL=incoherent`** (default)
- Exact for a 4π-integrated, spin-summed observable: the l ≠ l′ cross terms vanish.
- In a restricted angular window they survive (7Li 50–70°: +77 % for s = 1, −9 % for s = 2).
  Use `theta=` there, which computes the Blatt–Biedenharn sum over all Jπ.
- `coherent` is neither the 4π nor the fixed-angle observable. It is kept to compare with older fits and is refused
  with `theta=`.
- `tests/7Li_p_a`: 2138 incoherent against 3197 coherent.

**`vertex=constant`** (default)
- B_c = S_c(E₁) of the lowest level of each Jπ group. M_l then factors out of the level
  sum, and the HOES amplitude is the same in every R-matrix representation (Brune and
  formal agree to 2e-11).
- `perlevel` puts S_c(E_λ) inside the sum. It is representation-dependent once two levels
  of one Jπ interfere (18O ApJ 723 band: 23 % against 5 % peak rms). Use it only to
  reproduce pre-Sep-2026 numbers.
- `onshell` (L = S + iP) is boundary-free and complex. It differs near and above the
  barrier. Use it as a variant.

**`kinematics=`**: the factor K(E) the data still contain.

| data were divided by | key | K(E) |
|---|---|---|
| nothing but \|φ\|² | `triple` | 1 |
| full three-body KF (Typel & Baur eq. 16) | `kf3body` | 1/(μ_f k_f) |
| λ₃/λ₂ with on-shell p (Pizzone 2011) | `lambda32` | 1/k_i |
| La Cognata working formula, or unknown | `lacognata` | k_f/μ_f |

- Only the energy dependence matters, since the norm is free. `lambda32` changes by ×7
  over 0.02–1 MeV.
- For several exit channels of different Q the convention moves their relative weight
  (12C+12C: ×5 between `lacognata` and `kf3body` for α against p).
- 17O ωγ(5⁻): only `lambda32` reproduces the published 0.046 eV; the other conventions
  give about twice that.

**`coulombIntegral=1`**
- Adds C_l, the x–A Coulomb interaction outside the channel radius (|C_l| = 20–45 % of
  the surface term below the barrier). It does nothing for neutrons and costs ×3–4
  run time. No published analysis includes it.
- Refused with `vertexModel=dw` and with R(E) on a distorted a+A wave (`coulomb`, or
  `optical` unless `opticalAA=plane`). Those already contain the force.

**`spectatorEnergy=` / `spectatorEnergy[k]=`**
- One mean T_s added to E + B. It is a crude stand-in for the `ps` window and excludes
  it for the same pair.
- A d spectator at 30–50 MeV/c cuts gives 0.5–1.4 MeV.

**`weight[k]=<file>`**
- w(E) multiplying the THM model of segment k before folding (log-linear
  interpolation).
- The hook for an external R(E) (finite-range DWBA, CDCC). The built-in R(E) is
  `distortion=`.

## Experiment-line options

**Kinematics** `beam= target= spectator= Ebeam=`: all four or none. Required by
everything below except `theta`, `background` and `cbackground`.

**`lineshape=on`**: the Coulomb line shape N_C of a charged spectator.
- Each level's exit amplitude is multiplied inside the coherent sum, with
  |N_C|² = exp[2ζ arctan(2(E_λ−E)/Γ_λ)] and ζ < 0, so peaks move up by −ζΓ/2.
- Brune or Park only (refused with the formal parametrization). Refused with
  `cbackground=`. N_C ≡ 1 for a neutron spectator.
- `session.thm_lineshape(name, E)` returns ζ, E_sF, η_sb (the approximation needs it
  ≪ 1) and |N_C|² per level.
- Sizes: 12C(14N,d) ζ ≈ −0.13 … −0.49 (large); 17O (p spectator) ζ ≈ −0.22, but
  η_sb ≈ 0.21 makes it marginal, with +10–17 % on strengths and Δχ² +3.6.

**`ps=`** (`delta` default, `hulthen:0-40`, `hulthen:a,b:0-40`, `gauss:FWHM:0-40`,
`table:<file>`) with **`psNodes=`**
- Averages the HOES cross section over the accepted spectator directions at fixed E.
- The weight is |φ(q)|² d cos θ_cm = |φ|² q dq over the reachable q only. It is not the
  isotropic q² dq measure that AZURE2 used before Oct 2026; refits made with the old
  measure are superseded.
- Matters where M_l has a node inside the window: the window fills the nodes.
- `session.thm_vertex(name, E, strict=True)` raises instead of silently returning the
  nearest point when E is out of reach.

**`spectatorAngles=`** (`[cm:]a-b`, `table:<file>`) with **`spectatorAngleNodes=`**
(default 8)
- The accepted spectator directions. With `ps` it forms one acceptance, over which the
  pw vertex, R(E) and the DW vertex are averaged.
- Needs `distortion=coulomb|optical` or a `ps` distribution. Excludes `spectatorAngle`
  and `psNodes`.

**`distortion=coulomb|optical|table:<file>`**
- Related keys: `opticalAA=`, `opticalSF=`, `spectatorAngle=qf|<lab deg>|cm:<deg>`,
  `distortionRef=`, `distortionRatio=dwpw|dw`, `boundState=whittaker|yukawa[:rmin]`.
- The zero-range DWBA factor R(E) = |M(E)|²/|M(E_ref)|² multiplies the pw model. It is
  the same as dividing a PWA-extracted S* by R.
- Use it for charged spectators below the s+F barrier.
- 12C+12C: it reproduces the Mukhamedzhanov & Pang 2019 curve with defaults; R(0.8)/R(2.7)
  ≈ 1/470. 18O(d,n): R = 0.93–1.08 (negligible).
- `session.thm_distortion(name, E)`.

**`vertexModel=dw`** (needs `distortion=coulomb|optical`)
- Replaces M_l by the surface term of the prior-form DWBA built from the a+A and s+F
  distorted waves and a finite-range s–x tail. R(E) is then not applied, because the DW
  vertex contains it.
- Mandatory as a variant when strengths of different l, or of levels far apart in ρ,
  come from peak areas.
- 19F |M1(213)/M0(324)|² at a_p = 4.1 / 5.1 / 6.1 fm:

  | vertex | ratio |
  |---|---|
  | pw, qf | 1.10 / 3.00 / 19.7 |
  | pw + window | 2.23 / 8.07 / 6.49 |
  | dw Coulomb + window | 0.36 / 0.52 / 0.81 |
  | dw optical + window | 0.89 / 0.94 / 1.06 |

  The radius dependence nearly disappears with dw; the potentials and the acceptance
  take its place.
- It does not always shrink: for 17O (neutron entrance), ωγ(5⁻) varies ×1.8 (pw) and
  ×2.7 (dw) over 3.64–5.64 fm.
- `session.thm_vertex` returns `model`, `M2`, `M2_qf` and `M2_pw`. Start-up takes 4–6 s.

**Optical potentials** `opticalAA=` / `opticalSF=`
- Global sets: `ancai06`, `daehnick80` (d), `kd03` (n, p), `bg71` (t, ³He), `liang09`
  (³He), `mcfadden66`, `avrigeanu94` (⁴He). Append `:extrapolate` outside the stated
  range (a WARNING). The other forms are `plane`, `coulomb`, or ten numbers
  `V,R,a,W,RW,aW,WD,RD,aD,RC`.
- There is no heavy-ion global set: 14N+12C takes ten numbers.
- The spread between standard sets is small (19F l ratios 1 %; d+24Mg ≤ 3 %), far below
  Coulomb → optical (×1.7–20).
- Light targets need `:extrapolate`. Treat the resulting DW variants as checks and say
  so.

**`theta=a-b`**
- The model becomes ⟨dσ/dΩ⟩ of the exit pair over the c.m. window, with all Jπ and l
  interfering.
- 7Li Tumino (θ_cm 50–70°): shape 14 % rms, 45 % max. 18O single 1/2⁺: isotropic, no
  change.

**`background=const|linear|quadratic`**
- A smooth incoherent term, added to the folded model and profiled linearly with the
  norm.
- It may come out negative: judge it on the plot.

**`cbackground=<Jπ>:<exit>[:s,l,s′,l′][:const|:linear][=Re,Im,...]`**
- A THM-only interfering amplitude c(E)·M_l. Its Re and Im are fit parameters (kind
  `cbkg`, Fitting tab > THM Background).
- Strongly correlated with the level parameters, with two phase solutions of equal χ².
- 18O doublet: it recovered 86 of 212 χ² units without reconciling THM and direct data.
  19F: it raised ωγ(828) only to ~40 % of the direct value.
- Use it as a diagnostic of a THM–direct tension. Physics present in the direct data
  belongs in a background pole: a broad level of the same Jπ on the Levels tab.

## Allowed combinations (engine-enforced)

| | `coulombIntegral=1` | `lineshape=on` | `distortion` (R) | `vertexModel=dw` |
|---|---|---|---|---|
| `coulombIntegral=1` | – | allowed | refused (`coulomb`; `optical` unless `opticalAA=plane`); `table` warned | refused |
| `lineshape=on` | allowed | – | allowed | allowed |
| `distortion` | see left | allowed | – | R not applied |
| `ps` window | allowed | allowed | `distortionRatio=dw` refused | allowed (per node) |

## Brune and Park

- AZURE2 inputs are Brune (observed) parameters. The CLI runs Brune by default
  (`--no-brune` for Lane–Thomas). pyazr: `use_brune=True` default; `use_park=True`.
- A formal published set (ApJ 723 Table 3): solve the Brune eigenproblem and enter
  amplitudes with `gammaIsRWA` (field 33). Check the on-shell σ against an independent
  formal calculation.
- Brune fails when γ²dS/dE is large ("Denominator less than zero"; 12C+12C θ² = 26 at
  6.41 fm). Use a larger radius or smaller seeds.
- Park (`--use-park`, GUI "Use Park parametrization", pyazr `use_park=True`):
  - It gives the same THM model as Brune (`tests/thm_park`: every example and option,
    ≤ 1e-8), and fits land on the same minimum.
  - The fit amplitude is the observed width (Γ = 2Pγ²), so fixing, bounding or putting
    a prior on one width acts on one parameter.
  - Keep J_λ > 0 (penalty and warning otherwise; MCMC rejects).
  - A `gammaIsRWA` value in the file is Brune's amplitude in both modes.
- A `.azr` whose widths Brune cannot reach is a different model in the two modes: fix
  the input first.

## Identical nuclei

- (1+δ₁₂) = 2 multiplies every cross section out of an identical entrance pair (since
  df8c2c8). Identical exit pairs get no factor.
- Two identical 0⁺ bosons fuse only through even-J natural-parity levels. THM can
  populate odd-J or unnatural-parity levels by transfer (12C+12C 0.877 MeV 1⁻).
  - AZURE2 keeps l = J with a Bose-symmetry warning and cannot average l = J ± 1.
  - Workaround: a second, identical entrance pair for the on-shell segments. Odd-J
    levels couple to the THM pair only; allowed levels couple to both with the same
    width. Cost: the second channel enters the level matrix twice (χ² 110.69 → 111.65
    near the Wigner limit).
  - Exact alternative: two sessions, with the forbidden amplitudes zeroed for the
    on-shell observables.

## Refit against estimate

- Re-profiling the THM scale at fixed parameters measures goodness of fit only (7Li
  `lambda32`: 70.7 re-profiled against 65.5 refit).
- The fixed-peak-area estimate scales γ_x² by [|M_l|²K]_old/[|M_l|²K]_new per level. It
  reproduced refits to ≤ 5 % for isolated 12C+12C levels below 1.75 MeV. It **fails
  near vertex nodes** (a 4⁺ level moved ×490) and for interfering levels.
- Mark estimates with `*`. Final numbers come from full refits only.
- A distortion weight is absorbed into the widths only by a full refit (12C+12C:
  restricted refit χ²_THM 147, full refit 56).

## Vertex nodes to know

- Plot |M_l(pa)|² over the data window at a and a ± 1 fm before choosing an anchor.
- p+18O, l = 0 node: 3.15 / 1.32 / 0.283 MeV at a_p = 4.1 / 5.1 / 6.1 fm.
- 12C+12C: l = 4 node at 1.98 MeV for 7.5 fm; l = 0 near 0.88–0.98 MeV for 6.5 fm (the 0.918 2⁺
  moves 35 keV there).
- Near-barrier Trojan horses (19F + d at 55 MeV): the entrance deceleration halves the
  local p at the surface, and plane-wave l ratios are off by an order of magnitude.
  Compare with `vertexModel=dw`.
