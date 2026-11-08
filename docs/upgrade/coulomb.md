# Coulomb functions and channel functions (`physics/coulomb/`, `channel_functions`)

Every particle channel of the R-matrix enters the collision matrix through four numbers at the channel energy: the
penetrability P_c, the shift function S_c, its energy derivative S′_c (Brune and Park bases) and the hard-sphere phase
φ_c, plus the Coulomb phase σ_l that relates the amplitudes of different l. All of them follow from the regular and
irregular Coulomb functions F_l, G_l and their derivatives at the channel radius, or, below threshold, from the
Whittaker function. This note fixes the interface of the part that computes them (`physics/coulomb/`, from Jakub
Skowroński's pycoulomb), the interface of `channel_functions` built on it, what in the engine they replace, the order
in which they replace it, and the tests. ARCHITECTURE.md §1 (exact energies, decision 11), §3, §4 and §6 set the
frame; this is step 3 of `cnuc_split.md`.

## 1. Coulomb functions

### 1.1 Definition

For a partition α of charge numbers Z₁, Z₂ (particles 1 and 2) and reduced mass μ (u), at channel energy E_c (MeV,
c.m., above the threshold of the pair in its states), the radial functions solve

    u″(ρ) + [1 − 2η/ρ − l(l+1)/ρ²] u(ρ) = 0,                                                   (1)

with ρ = k r (k in fm⁻¹, r in fm) and the Sommerfeld parameter η = Z₁Z₂ e² μ / (ħ² k). η > 0 is a repulsive pair,
η = 0 a neutral one (F, G are then ρ j_l(ρ) and −ρ y_l(ρ)). The normalisation is that of Abramowitz and Stegun ch. 14
and DLMF §33.2: F_l ~ sin θ_l, G_l ~ cos θ_l for ρ → ∞, θ_l = ρ − η ln 2ρ − lπ/2 + σ_l, with the Wronskian

    F′_l G_l − F_l G′_l = 1                                                                  (2)

(primes are d/dρ). The Coulomb phase is σ_l(η) = arg Γ(l + 1 + iη), the continuous branch, σ_0(0) = 0 and
σ_l = σ_{l−1} + arctan(η/l). The engine uses only e^{i(σ_l − σ_0)} (`conventions.rst`, Eq. conv-phases).

### 1.2 Interface

Only the dimensionless arguments enter; the kinematics (η and ρ from E_c, a_c and the masses) belong to the caller
(section 3.2). Namespace `azure2::physics::coulomb`, header `azure2/physics/coulomb/Coulomb.h`:

    enum class Method { coulsm, coul, coulfg, coulcc };
    enum class Status { ok, invalid_argument, failed };
    struct CoulombWaves { double F, Fp, G, Gp, sigma; Status status; };
    CoulombWaves coulomb(int l, double eta, double rho, Method method = Method::coulsm) noexcept;
    double coulomb_phase(int l, double eta) noexcept;
    const char *method_name(Method); bool method_from_name(const std::string &, Method &);

The rules are pycoulomb's (its `DESIGN.md` §2–3): a pure function of its arguments, no object, no static or global
mutable state, no cache, no printing, no exception leaving it; the same arguments give the same bits in any thread
and in any order of calls. Arguments are checked first (0 ≤ l ≤ 10⁴, η and ρ finite, ρ > 0); an algorithm that
throws, reports failure or returns a non-finite value gives `failed`. On any status but `ok` every number is NaN, so
a caller that ignores the status cannot use a wrong value without seeing NaN.

### 1.3 Algorithms, one file each behind the same signature

| `Method` | source (pycoulomb) | what it is |
|---|---|---|
| `coulsm` (default) | `libcoulomb/coulsm/` | N. Michel's `cwfcomp` (Comput. Phys. Commun. 176 (2007) 232) with the patch from the SAMMY side; AZURE2's `external/coul` since 2025-09-05. The default because it reproduces the engine's numbers bit for bit. |
| `coul` | `libcoulomb/coul/` | `cwfcomp` as in AZURE2 before 2025-09. |
| `coulfg` | `libcoulomb/coulfg/` | A. R. Barnett's COULFG (Comput. Phys. Commun. 27 (1982) 147), C translation; F′, G′ by the recurrence A&S 14.2.1 from l + 1. |
| `coulcc` | branch `coulcc-cpp`, `coulcc/` | COULCC (Thompson and Barnett, Comput. Phys. Commun. 36 (1985) 363), the C++17 translation bit-identical to the Fortran; called with real ρ, η and λ = l. |

Each algorithm is one adaptor, `Status <method>(int l, double eta, double rho, CoulombWaves &out)`, in
`physics/coulomb/src/method_<name>.cpp`; `Coulomb.cpp` checks the arguments, dispatches, sets σ_l and applies the
NaN rule. The third-party sources sit unchanged under `physics/coulomb/<method>/` apart from a namespace and the
`M_*` guards for MinGW/MSVC; `README.md` there records the pycoulomb commits and each local change. COULCC's complex
interface (complex l, η, ρ; Bessel functions) is not exposed: AZURE2 has no use for it yet.

GSL's `gsl_sf_coulomb_wave_FG_e` is not a method of `physics/coulomb/` (physics/ links nothing). It stays as an
independent cross-check in a test outside physics/ (section 5), and the engine's `--gsl-coul` path stays in the
engine (section 4, open question).

### 1.4 Whittaker functions (bound channels)

For E_c < 0 the decaying solution of (1) with κ = |k| and η_b = Z₁Z₂ e² μ / (ħ² κ) is the Whittaker function

    W_l(ρ) = W_{−η_b, l+½}(2ρ) = e^{−ρ} (2ρ)^{l+1} U(l + 1 + η_b, 2l + 2, 2ρ),   ρ = κ r,           (3)

and the bound-channel shift function is S_c = ρ W′/W at ρ = κ a_c. The interface is pycoulomb's: `whittaker(l, η_b,
ρ)` → W, W′, status, and `whittaker_scaled(l, η_b, ρ)` → mantissa, log₁₀ scale, ρW′/W, which stays finite just below a
charged threshold, where η_b grows like E_b^−½ and W underflows. pycoulomb computes U with GSL's `gsl_sf_hyperg_U`,
which physics/ cannot use.

**Design (2026-11-08, step 3b (i)).** Three ways were compared on 990 points (l = 0–8; η_b = 0, 0.01, 0.3, 1, 2, 3,
5, 10, 30, 100, 300; ρ = 10⁻³–50) against mpmath's `whitw` at 40 digits:

| method | worst relative error of W | points without a value |
|---|---|---|
| GSL `gsl_sf_hyperg_U` (the engine's `WhitFunc`) | 1.2 × 10⁻⁹ (l = 0, η_b = 100, ρ = 10⁻³) | 90 (GSL error, overflow or underflow) |
| COULCC at ρ = iκr, η = −iη_b (the earlier proposal) | — | 0/0 at η_b = l + 1, l + 2, … (below) |
| the integral representation below (Python prototype) | 4.6 × 10⁻¹³ (l = 5, η_b = 300, ρ = 10⁻³) | 0 |

COULCC was dropped: H⁺_l(−iη_b, iκr) = e^{−i(l+η_b)π/2} [Γ(l+1+η_b)/Γ(l+1−η_b)]^½ W_l(κr) (DLMF 33.2.7; checked
with mpmath to 30 digits), and the factor vanishes at η_b = l + 1, l + 2, …, so W = H⁺/factor loses digits there,
however accurate H⁺ is. Its logarithmic derivative is normalization-free, but the ANC conversions and N_f need W
itself.

**Method.** For a = l + 1 + η_b ≥ 1 (η_b ≥ 0: a neutral or repulsive pair; nuclear pairs are never attractive),
b = 2l + 2, z = 2ρ:

    U(a, b, z) = (1/Γ(a)) ∫₀^∞ e^{−zt} t^{a−1} (1+t)^{b−a−1} dt,    dU/dz = −a U(a+1, b+1, z)   (DLMF 13.4.4, 13.3.22)

so that ln W = −ρ + (l+1) ln 2ρ − ln Γ(a) + ln I(a, b, 2ρ) and ρW′/W = −ρ + (l+1) − 2aρ I(a+1, b+1, 2ρ)/I(a, b, 2ρ)
(primes d/dρ). The integrand is positive and has one maximum t₀ (a quadratic in t); with t = t₀ exp[(π/2) sinh s]
(double-exponential in both tails) the trapezoidal sum in s converges geometrically. The sum is formed in
logarithms (largest term factored out), so W is never formed before the end: `whittaker_scaled` returns mantissa
and power of ten, and `whittaker` returns W, W′ when they are inside the double range and `out_of_range` (NaN)
otherwise. The step is halved from 1/2 until two sums agree to 4ε; terms are summed while they exceed e⁻⁴⁰ of the
largest. No GSL, no table, no state: the same arguments give the same bits in any thread.

**Interface** (`Whittaker.h` in `physics/coulomb`, beside `Coulomb.h`):

    struct WhittakerWave { double W, Wp; Status status; };                 // W_l(ρ) = W_{−η_b, l+½}(2ρ), dW/dρ
    struct WhittakerScaled { double mantissa; long log10_scale; double rho_dlog; Status status; };
    WhittakerWave whittaker(int l, double eta_b, double rho) noexcept;    // out_of_range outside the double range
    WhittakerScaled whittaker_scaled(int l, double eta_b, double rho) noexcept;   // W = mantissa·10^scale, ρW′/W

`invalid_argument` for l < 0 or > 10000, η_b < 0 or not finite, ρ ≤ 0 or not finite. The binding exposes both with
array inputs.

**Tests.** `physics_whittaker` against a committed mpmath reference (`physics/tests/data/whittaker.csv`, made by
`make_whittaker_reference.py`, 40 digits: l = 0–10, η_b = 0–1000, ρ = 10⁻⁴–200, the scaled form where W leaves the
double range), with the tolerance set from the measured error; ρW′/W against mpmath's derivative; the same bits from
eight threads; `invalid_argument` cases. Outside physics/, a reference test compares with GSL (where GSL gives a
value) and times both. Speed decides the engine part: `ECIntegral` spends about 85 % of its time in GSL's U.

**Engine (step 3b (i), second part, after the timing).** `WhitFunc` becomes a thin adapter over `whittaker` /
`whittaker_scaled` (it keeps η_b and ρ from `Kinematics`); `ShftFunc`'s bound branch takes ρW′/W in place of
`CoulombTable::ExactBound` (GSL's U ratio, trusted for η_b ≤ 1000 and η_b ≤ 50 or z ≥ 0.05, linearly interpolated
to S(0) closer to threshold; the numerical derivative of W is only the opt-in `AZURE_BOUND_SHIFT_NUMERICAL` route).
The new function holds 10⁻¹² over that band and below it, so whether the interpolation band narrows (outputs near
threshold move) is decided with its own comparison. `CoulombTable`, `ThmNumerics` and `ThmFunc` call the same
function; GSL's U leaves the engine. In the data path `ECIntegral` evaluates W at every radius its adaptive rule
samples (13× GSL's time per call); the rate path already tabulates the final-state W once per pathway on a fixed
grid, and the per-integral table of step 3b (ii) does the same for the data path. Bound-state outputs (ANC ↔ γ conversions, N_f, external capture to bound states, sub-threshold shifts) move
in their last digits; the 179-project comparison classifies every change, and the T1–T3 truth checks must hold.

## 2. Channel functions

### 2.1 Definition

At channel energy E_c > 0 and channel radius a_c (fm), ρ = k a_c,

    P_c = ρ / (F_l² + G_l²),   S_c = ρ (F_l F′_l + G_l G′_l) / (F_l² + G_l²),
    e^{−iφ_c} = (G_l − i F_l) / (F_l² + G_l²)^½,   ω_c = σ_l − σ_0,                                (4)

all dimensionless, and S′_c = dS_c/dE_c (MeV⁻¹). The formulas and their order of operations are the engine's
(`CoulFunc.cpp`, `conventions.rst` Eqs. conv-ps, conv-phases), including the two guards for η ≳ 110 just above a
charged threshold, where F² + G² overflows: S and e^{−iφ} are then formed with F, G scaled by hypot(F, G), and only
then, so that every value the plain formula gives is unchanged. S′_c is the Richardson-extrapolated central
difference of S with step h = 10⁻³ E_c, (4 D(h/2) − D(h))/3, D(h) = [S(E + h) − S(E − h)]/2h; its error is about
10⁻⁸ relative, set by the 10⁻¹⁰ accuracy of S from `coulsm` (`cnuc_split.md` §5).

### 2.2 Interface

    struct ChannelKinematics { int z1z2; double mu_u, m1_MeV, m2_MeV; bool relativistic; };
    double eta(const ChannelKinematics &, double E_c);       double rho(…, double a_fm, double E_c);
    struct OpenChannelValues { double P, S, dS_dE; std::complex<double> exp_minus_i_phi; double sigma; Status status; };
    OpenChannelValues channel_functions(const ChannelKinematics &, int l, double a_fm, double E_c,
                                        coulomb::Method = coulsm, ChannelQuantities = all);

`eta` and `rho` are the engine's `Kinematics::Channel::Eta/Rho` (non-relativistic, or the relativistic forms of
Haeck et al., LA-UR-19-31084), moved into physics/ with the constants they use, so that the engine and physics/ share
one formula; the engine's `Kinematics` calls them. The same function is exposed in the binding with array inputs
(ARCHITECTURE.md §8) through the existing `coulomb_functions` call.

`channel_functions` covers the open channel, E_c > 0 with η ≤ 100. Below threshold, at threshold and in the
"Coulomb threshold" band η > 100 (P = 0 in double precision, S continued through E_c = 0), the engine's
`ChannelFunc` keeps its rules on the Whittaker branch until section 1.4 is done; it then becomes one function of
physics/ for all E_c.

### 2.3 The table

ARCHITECTURE.md §3–4 call for a read-only map of P, S, S′, φ, σ_l per (partition α, l, channel energy), built once per
evaluation-energy set and shared by every reaction, angle and J^π block (today's `EDependentTables::Slot(jGroup,
channel)` stores one (α, l) once per J^π group and channel spin). Its contract:

- keyed by (α, l, a_c) and the index of the evaluation energy E_n; values are those of `channel_functions` at
  E_c = (E_n − S_α) − e_α (MeV, c.m.; the order of `cnuc_split.md` §5);
- built in parallel, each thread writing only its own rows, so the values do not depend on the thread count;
- read-only once built; when a free energy shift ν moves the energies of a measurement, a new table for those
  energies is built single-threaded between evaluations and replaces the old one (decision 11); nothing writes into a
  table that an evaluation reads;
- no other cache: neither a process-wide memo nor a tolerance on the energy.

The table replaces the Coulomb rows of `EPoint::EDependentTables` (penetrability, shift, hard-sphere and Coulomb
phase), their copy-on-write sharing and `NodeTableCache` for those rows. That replacement changes who owns the
per-point energies (Phase 4, `EnergyGrid`), and is not part of this branch (section 6).

## 3. What this replaces

| today | role | becomes |
|---|---|---|
| `external/coul` (`cwfcomp`, `cwfcomp_accurate`) | the Coulomb routine | `physics/coulomb/coulsm/` (and `coul/`); deleted |
| direct `Coulomb_wave_functions_accurate` calls in `CoulFunc` (three: plain, hybrid boundary values, hybrid fallbacks), `CoulombTable::ExactOpen`, `ThmDistortion` | five copies of one call | `coulomb::coulomb(l, η, ρ, coulsm)` |
| `CoulFuncCache` (process-wide, per-thread memo; energies matched to 10⁻¹² MeV) | speed | deleted in step 3b, replaced by the channel table and the radial scratch (`channel_table.md`) |
| `CoulFunc::Penetrability`, `PEShift`, `HardSpherePhase`, `PEShift_dE` | Eq. (4) | thin calls of `channel_functions` (the Coulomb-table and hybrid paths stay in `CoulFunc`) |
| `Kinematics::Channel::Eta`, `Rho` | η, ρ | call physics/ |
| `WhitFunc`, `ShftFunc` | bound channels | stay (section 1.4) |
| `EDependentTables` Coulomb rows | per-point store | filled from the channel table (step 3b); the rows themselves stay until `EnergyGrid` (section 2.3) |
| `g_shiftDerivMemo` (exact-argument memo of S′, S″ at level energies) | speed | stays: level energies are fit parameters, not evaluation energies (`channel_table.md` §2) |

The memo's tolerance makes the engine's last digits depend on the order of the calls: an energy within 10⁻¹² MeV of
one computed earlier in the same thread gets the earlier one's functions. Deleting it lets every point get its own
values; outputs then move at the level allowed by S (10⁻¹⁰) and S′ (10⁻⁸) for an energy change of 10⁻¹² MeV, and the
pyazr/CLI identity checks of THM (HANDOFF, director's entry of 2026-11-07) pass. It is also most of the engine's speed
where energies recur, so it is deleted only with the table that replaces it (§7).

## 4. Engine and front ends

The engine always calls `Method::coulsm` in this branch: no new option. `--gsl-coul` keeps its present meaning (GSL's
routine, called from `CoulFunc`, which may use GSL), because removing or remapping a CLI flag is a user-visible change
to be decided (section 6). The memo stays in front of `CoulombWavesAt` as before, so `cache_stats()` is unchanged.

## 5. Tests

- `physics_coulomb_unit` (physics/ only): every method against mpmath at 40 digits on pycoulomb's reference set
  (`tests/data/coulomb.csv`, copied with its generator), with the tolerance each method meets there stated in the
  test; the Wronskian (2) to 10⁻¹⁰; σ_l against the arctan recurrence; the four methods against each other where
  each is accurate (identity test); argument checks (NaN and status); bit-identical results from 8 threads calling
  in different orders. It fails if a method is wired to the wrong source, if a sign or normalisation is changed, or
  if any state is shared between calls.
- `physics_channel_functions_unit`: P, S, φ against Eq. (4) evaluated with mpmath's F, G; the neutral closed forms
  (l = 0: P = ρ, S = 0, φ = ρ; l = 1: P = ρ³/(1 + ρ²), S = −1/(1 + ρ²), φ = ρ − arctan ρ); S′ against the derivative
  of the closed form; η and ρ against the non-relativistic formulas and the relativistic closed forms of
  `Kinematics.h`.
- `coulomb_gsl_crosscheck` (outside physics/, links GSL): `coulsm` against `gsl_sf_coulomb_wave_FG_e` where GSL is
  accurate.
- Engine: every `.azr` under `tests/` (179) in the Brune, formal (`--no-brune`) and Park (`--use-park`) bases, before
  and after, every output file compared number by number. The routing commits must be bit-identical (each changed file
  classified otherwise; anything above 10⁻⁹ relative on observables, unexplained, stops the step). `tests/upgrade/check.sh`, `ctest -LE slow`, `tests/run_tests.sh`, 1 vs 8 threads identical, T2 benchmark.

Code: `physics/coulomb/` (sources, adaptors, `Coulomb.h`), `physics/src/ChannelFunctions.cpp`,
`physics/include/azure2/physics/ChannelFunctions.h`, `physics/tests/coulomb_*`, `src/CoulFunc.cpp`, `include/Kinematics.h`.

## 6. Order and what is left

1. This note. 2. `physics/coulomb/` with its tests. 3. The engine calls it; `external/coul` deleted (bit-identical).
4. `channel_functions` and the kinematics in physics/; `CoulFunc` calls them (bit-identical). 5. `CoulFuncCache`
deleted (outputs classified) — moved to step 3b, after the table (§7). Left for step 3b (`cnuc_split.md` §3): the Whittaker functions without GSL and the
deletion of `WhitFunc`, `ShftFunc` and `ChannelFunc` (bound branch into `channel_functions`); the table of §2.3 in
place of the `EDependentTables` rows and `g_shiftDerivMemo`; the hybrid-potential channels ((α, l, j) entries of the
table); `CoulombTable` (an interpolation of these same functions, opt-in) and `--gsl-coul`, which are user-visible
options to keep, remap to a `Method`, or retire.

## 7. As built (branch `feature/upgrade-1-coulomb`, 2026-11-07)

Steps 1–4 of section 6 are done, one commit each. Step 5 (the memo's deletion) was built and measured, and moved to
the branch `feature/upgrade-1-coulomb-nomemo`, to be re-cut on top of the per-stage table (step 3b): without the
memo the engine is 3–17× slower and `convolution_unify` exceeds its time limit (below).

**Outputs.** Engine before (f549fa23) and after each step, GCC 11.5, Release, one thread, every `.azr` under
`tests/` (179) in three bases, every output file compared number by number (`cnuc_split.md` §6, step 3):

- Routing to `physics/coulomb` (step 3): identical except the projects with a hybrid-potential channel (`T3`,
  `T3_start`, `hybrid_potential`): observables ≤ 9.9 × 10⁻¹¹, one width in `parameters.out` 5.0 × 10⁻⁸ (formal
  basis, the formal → Brune iteration). Cause: `cwfcomp` returns F and G as complex numbers whose imaginary parts
  are not exactly 0 for real arguments (7 of 125 probe points, up to 1.0 × 10⁻¹⁵ relative); the hybrid model formed
  H⁺ = G + iF from the complex values and so carried them into the Numerov boundary values, `coulomb()` returns the
  real parts. Plain Coulomb channels never used the imaginary parts.
- `channel_functions` and the kinematics in physics/ (step 4): bit-identical to step 3 (2040 files).

**Speed of this branch** against f549fa23, one Condor host, both engines alternately, `--no-coulomb-table`, mean of
two: host d32cepyc086, mean of two, before (f549fa23) → this branch: T2 calc 1.56 → 1.36 s (1 thread; one
f549fa23 run at 1.89 s), 0.44 → 0.42 s (8); T2 fit200 80.5 → 81.9 s (1 thread), 36.5 → 37.0 s (8); menu-1
calculation, one thread, `13N` 10.05 → 10.16 s, `natural_neon` 12.35 → 12.45 s; `convolution_unify` `tgt_fixed`
(2 threads) 29.9 → 30.6 s. The routing costs at most 2 % (fit and `tgt_fixed`), within the run-to-run spread of the
calculations; peak RSS within 2 %.

**The memo's deletion, measured (branch `-nomemo`, c3dc22e5).** Observables ≤ 7.0 × 10⁻¹¹ relative, except the
rounding residue of an exact zero (E1 capture at θ = 0°, 1.5 × 10⁻²⁵ → 3.0 × 10⁻²⁵ beside 1.2 × 10⁻⁹ at 60°);
`intEC.dat` of 20 projects ≤ 4.0 × 10⁻¹¹; the THM checks `pyazr_thm_profiled` and `pyazr_thm_fixed_angle` pass (the
memo's call-order dependence was their failure). Speed (host qa-a10-001, same method): T2 calc 1.37 → 6.16 s
(1 thread), 0.46 → 1.40 s (8); fit200 93.1 → 114.4 s, 42.7 → 48.1 s; menu-1 calculation, one thread, `13N` 10.0 →
32.1 s, `beam_profile_kernel` 52.6 → 114.9 s, `natural_neon` 12.3 → 105.6 s, `18O_p_a_thm` 0.14 → 1.39 s;
`convolution_unify` `tgt_fixed` (2 threads) 84 → 1414 s, past its 600 s limit. The memo returned stored functions
wherever an energy recurred: shared point energies, target and resolution sub-points, external-capture radii, a
second pyazr session.

**Next (step 3b).** The table of section 2.3, one evaluation per (α, l, a_c, E_c) of the run's evaluation energies,
shared read-only, and a per-integral table of F, G at the quadrature radii for the external-capture integrals; then
the `-nomemo` commits re-cut on it, with the timings above repeated (they must come back near those of this branch).
