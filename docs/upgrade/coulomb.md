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
independent cross-check in a test outside physics/ (section 5). The engine's `--gsl-coul` path stayed in the
engine until 2026-11-16 (section 8 has what reached it and what it gave, section 9 its retirement).

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

**Part 1 verified (2026-11-09; build of bf8622e2, GCC 11.5 Release, pyazr on).** No warnings. `ctest -LE slow`
114/114, with `physics_whittaker`, `whittaker_gsl` and `pyazr_whittaker_wave` (the binding) passing in the full
build; 1 vs 8 threads byte-identical (T1–T3, 13N, three bases); the 179 projects × 3 bases against a18fe89e:
2040/2040 files identical (the engine does not call the new functions yet). `run_tests.sh` was not run: no engine
path changed.

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

**Engine part as built (2026-11-09, branch `feature/upgrade-1-whittaker-engine`).** Every call of GSL's U in the
engine goes through the two physics functions, with the engine's η_b and ρ = κr unchanged (non-relativistic or
`Kinematics` relativistic forms):

| engine call | was | becomes |
|---|---|---|
| `WhitFunc::operator()` | e^{−z/2} z^{l+1} `gsl_sf_hyperg_U` (z = 2κr) | `whittaker(l, η_b, κr).W`; outside the double range mantissa × 10^scale of `whittaker_scaled`, i.e. 0 on underflow and inf on overflow (IEEE), never a GSL error |
| `WhitFunc::Scaled` | `gsl_sf_hyperg_U_e10_e` | `whittaker_scaled` (mantissa, power of ten) |
| `CoulombTable::ExactBound` (ρW′/W; `ShftFunc`'s bound branch and the table's bound sheet) | ratio of two `gsl_sf_hyperg_U_e10_e` | `whittaker_scaled(...).rho_dlog` |
| `ThmBoundStateTail` | 2κ e^{−κr} `gsl_sf_hyperg_U_e`(1 + η_b, 2, 2κr) | `whittaker(0, η_b, κr).W / r` (same function: W_0/r) |

Invalid arguments (NaN η_b or ρ, ρ ≤ 0) keep the engine's present answer, 0 (the `WhitFunc` FIXME and `Scaled`'s
guard); they are not new failure modes. `ShftFunc`'s interpolation band near threshold (η_b ≤ 1000, and η_b ≤ 50 or
z ≥ 0.05) is kept unchanged here, so that this step changes only the function; narrowing it is a separate,
measured step. `ZeroEnergyLimit` (Bessel K) and the opt-in numerical route stay on GSL: they are not U.

The Coulomb table's bound sheet is filled from `ExactBound`, so its contents change: the table format becomes 2
(`kFormat`), and a format-1 file is refused with the existing warning ("table format 1, this build reads 2") and the
exact functions are used until a table is built. No file is overwritten.

**Comparison for this step.** Reference: the head of `feature/upgrade-1-whittaker` (engine identical to step 4).
Both builds run the 179 projects × 3 bases without a Coulomb table (`AZURE2_COULOMB_TABLE` pointing to a missing
file), so the difference is the function alone; the new build then also with a table of format 2 against itself
without it. Every changed file is classified; bound-state quantities (ANC ↔ γ, N_f, external capture to bound
states, sub-threshold S) are expected to move at ≲ 10⁻⁸ relative (GSL's error, measured above), anything larger is
explained. Timing: `ECIntegral`'s data path on a capture project (T1 or 13N), calc, one thread, before and after.

**Engine part verified (2026-11-09; build of 5ed85322, GCC 11.5 Release, pyazr on).** No warnings; no engine file
includes `gsl_sf_hyperg.h`. Without a table, against bf8622e2: Brune 595/689, formal 559/662, Park 605/689 files
identical; every changed observable within 2.2 × 10⁻⁹ of its column's largest value (worst: `relativistic/he5`);
`param.par`/`parameters.out` at most 9 × 10⁻⁸ relative in printed digits, larger relative changes only for values at
rounding (formal zero-input widths 10⁻¹⁸–10⁻³⁷; cross sections of 10⁻²² at 0° and 180° where they vanish), one
photon width of 0.99999999 eV printed as 1000 meV (`ec_threshold/ng`), and signs of −0.000000; all 20 `intEC.dat`
files change (every token, ~10⁻¹⁰), hence EC signature v6 with the stale-file case (s) of `tests/ec_m1_exact`
(fails without the bump). With a table built by each build (format 1 by bf8622e2, format 2 by this one), the
table-vs-exact differences are of the same size for both (worst `c12c12`, 5.7 and 7.6 × 10⁻⁵, interpolation), and
`c12c12` is identical between the two builds. `ctest -LE slow` 114/114 with a format-2 table (`pyazr_park_dsde`
compares a table session with an exact one and fails when the only table on the machine is refused);
`ec_m1_exact` PASS; `run_tests.sh` 133/133 (1421 s, 11 threads, format-2 table); 1 vs 8 threads byte-identical
(T1–T3, 13N, three bases). **Cost** (one thread, calc, mean of three, own table): 13N 7.10 → 9.14 s, T1 2.62 →
4.73 s; `AZURE_PROFILE` puts it all in set-up (the external-capture table: 1176 integrals, 1320 node values in T1,
recomputed every run, even with `intEC.dat` present), evaluation ~1 ms in both. This is the 13× cost per W call in
`ECIntegral`'s data path; the per-integral W table (step 3b (ii)) removes it and should come before merging.

**Step 3b (ii): the final-state W of the data path, memoized (design, 2026-11-09).** The cost above is the
final-state factor W_{l_f}(κ_f r) of `ECIntegral`: `fgEvaluate` (FW, GW) and the outgoing factor of `WWIntegrand`
evaluate it at every radius GSL's `qagiu` samples, for every node energy of every `ECTable` term, although it does
not depend on the entrance energy. The rate path's fixed Laguerre grid would remove it, but it changes the quadrature
of the data path (outputs move). Instead: `qagiu` maps [a, ∞) to t ∈ (0, 1] and bisects; a subinterval it creates is
the same bits whatever the integrand, and so are its 21 Gauss–Kronrod radii. The integrals of one final state at
different entrance energies therefore sample largely the same radii, and a memo keyed on the exact arguments returns
the bits a new evaluation would.

*Rule.* `WhitFunc::ScaledMemo(l, r, E_b, log10Scale)`: η_b and ρ formed as in `Scaled`, then a thread-local table
keyed on (l, η_b, ρ) — the arguments of `whittaker_scaled` themselves, so no key can be stale and no model or pair
identity enters — returns the stored (mantissa, scale) or computes and stores them. `whittaker_scaled` is a pure
function (same arguments, same bits, any thread), so results are bit-identical to `Scaled`. The table holds at most
2¹⁷ entries per thread (about 10 MB) and is cleared when full. It is used for the final-state factor in
`fgEvaluate` and in `WWIntegrand` (outgoing factor) only; the entrance factor of `WWIntegrand` depends on the
entrance energy and stays `Scaled`. `AZURE_PROFILE` reports hits and misses.

*Expected.* 179 × 3 byte-identical (with and without a table, against the head of this branch); threads identical;
the cost measured again (13N, T1, one thread, calc). If the hit rate is low (the trees of different energies share
few subintervals), the fallback is the fixed grid, as its own design with a comparison.

**Step 3b (ii) verified (2026-11-09; build of 82ca0985, GCC 11.5 Release, pyazr on).** No warnings. 179 × 3 against
5ed85322 (this branch before the memo): with a format-2 table 2040/2040 identical; without a table 1977/2040, the
63 others being `intEC.dat.sig` alone (the reference run without a table was made by 48fb454f, before the v6 bump;
the memo build's signatures equal the v6 build's, with and without a table). 1 vs 8 threads byte-identical (T1–T3,
13N, three bases); `ctest -LE slow` 114/114. Memo on T1 (calc): 394980 hits, 870 computed (99.8 %); 13N 394965 and
495. **Cost** (one thread, calc, mean of three, each build with its own table): T1 2.89 s (`-whittaker` head, GSL)
→ 5.14 s (engine switch) → 2.77 s (memo); 13N 7.60 → 9.71 → 7.08 s. The engine switch no longer costs time.

**The interpolation band near threshold (design, 2026-11-09).** `ShftFunc::operator()` evaluates ρW′/W exactly for
binding B ≥ B_min. Closer to threshold it interpolates linearly between S(−B_min) and S(0) (`ZeroEnergyLimit`, Bessel K).
B_min is the binding where η_b = 1000, or, if larger, the smaller of the bindings where η_b = 50 and where
z = 2κa = 0.05. The second condition described where GSL's U went wrong (off by 10⁻⁶ at η = 111, z = 1.2 × 10⁻³ for
p+p); `whittaker_scaled` has no such weakness. Measured on p+¹²C, p+p, α+¹²C and ¹²C+¹²C, l = 0 and 3, B = 10⁻¹⁴–1 MeV:
- status ok everywhere;
- S tends to S(0) linearly in B;
- below about 10⁻⁹ MeV (η_b ≳ 10⁴) it jitters at 10⁻⁹ relative, beyond the η_b ≤ 1000 its test covers to 10⁻¹²;
- each call takes 2–9 µs, but 16–230 µs at η_b ≳ 10⁷.

So the η_b = 1000 limit stays, and the GSL condition is dropped: B_min = B(η_b = 1000). For p+¹²C, B_min goes from
3.3 × 10⁻⁴ to 8.2 × 10⁻⁷ MeV. In what is left of the band, S − S(0) is linear to about 10⁻¹¹ (its curvature,
measured from the values above, is about 3 × 10⁻⁵ of the linear term at 10⁻³ MeV and shrinks with B).

*Expected.* S changes, without a table, only for bindings between the old and the new B_min. There the linear
interpolation was off by up to about 10⁻⁸ absolute, so outputs near threshold move at that level. dS/dE moves by up to
about 10⁻⁵ relative on a level within 3 × 10⁻⁴ MeV below threshold. The Coulomb table is filled from `ExactBound` (not
`operator()`), so its contents and format do not change. With a table, only the bindings the table does not cover
change. 179 × 3, with and without a table, against the head of this branch: every changed file classified. A change
beyond about 10⁻⁷ relative in an observable is examined before the step goes on. Then ctest, `run_tests.sh` and the
thread check.

**Band narrowed, verified (2026-11-09; build of 37e1c5ba).** 179 × 3 against the merged head db972c84:
- without a table, 2040/2040 identical: no project has a bound channel between the old and the new B_min;
- with a format-2 table, 2039/2040: `ec_threshold` `ng`, R = 3, one number, 1.2 × 10⁻¹¹ of its column (a point the
  table does not cover), class (b).

ctest 117/117, `run_tests.sh` 133/133, 1 vs 8 threads identical.

**Bound channel functions in physics/ (design, 2026-11-09).** After the band step, the bound branch of the channel
functions is the last engine code that does physics with GSL: `ShftFunc` (S(0) from `gsl_sf_bessel_Kn_scaled`, dS/dE
by `gsl_deriv_central`) and `ChannelFunc`/`CoulFunc` (the same derivative routine at a stencil across threshold and
above it). It moves to physics/:

- `physics/Derivatives.h`: `central_derivative(f, x, h)`, the algorithm of `gsl_deriv_central` (GSL 2.x
  `deriv/deriv.c`, GPL-3 like AZURE2; attribution in the source). It evaluates at x ± h and x ± h/2, then at most once
  more at the step that balances rounding and truncation, and gives the same bits as GSL for the same f.
  `second_derivative(f, x, h, side)` is the Richardson-extrapolated central second difference of `ShftFunc` and
  `CoulFunc`, with their one-sided stencil within 4h of threshold (side: below or above).
- `ChannelFunctions.h`:
  - `bessel_k_scaled(n, x)` = e^x K_n(x) for integer n ≥ 0, x > 0, by the trapezoid rule on
    ∫₀^∞ e^{−x(cosh t − 1)} cosh(nt) dt in long double (step 0.02, until a term is below 10⁻²² of the sum). Against
    mpmath (n = 0–21, x = 10⁻³–500) it agrees to 4.6 × 10⁻¹⁶;
  - `zero_energy_shift(c, a, l)` = −l − (x/2) K_{2l}(x)/K_{2l+1}(x), with x = (8ca)^½ (−l without Coulomb);
  - `neutral_bound_shift(l, ρ)`, the continued fraction of `ShftFunc::NeutralBoundShift`;
  - `bound_shift(kin, l, a, B)`: the order of `ShftFunc::operator()` without its table route. That is S(0) at B = 0,
    the neutral form, ρW′/W from `whittaker_scaled` for η_b ≤ 1000, and the linear interpolation to S(0) beyond.

The engine keeps the Coulomb table route, the stencil checks against the table (`forceExact_`), the opt-in numerical
W route (`AZURE_BOUND_SHIFT_NUMERICAL`, on GSL) and `--gsl-coul`. The functions are called with the engine's η_b, ρ and
c, so only S(0) changes: its K_n agree with GSL's to about 10⁻¹⁶.

*Verification.*
- physics tests: `bessel_k_scaled` against an mpmath table (`tests/data`, made by a script there);
  `zero_energy_shift` against the limit of `whittaker_scaled` at B = 10⁻⁸ MeV (linear extrapolation in B, about
  10⁻⁹); the derivatives on polynomials and exponentials, where they are exact or known.
- engine-side test: `central_derivative` against `gsl_deriv_central`, bit for bit, on 10⁴ cases (smooth, noisy and
  kinked functions; steps 10⁻⁶–10⁻²). This includes the second evaluation at the optimal step.
- 179 × 3 with and without a table against the band head. Expected: identical, except outputs that use S(0) (bound
  levels within B(η_b = 1000) of threshold, S at threshold), which move at 10⁻¹⁵. Classify, then ctest,
  `run_tests.sh` and the thread check.

**Bound channel functions verified (2026-11-09; build of 5c4fb3ed).** Physics test `physics_bound_functions` (K_n
against mpmath, 374 values, to 2 × 10⁻¹⁵; S(0) against the threshold limit; neutral closed forms; the band; exact
rules). Engine-side `physics_derivatives_gsl`: `central_derivative` against `gsl_deriv_central` bit for bit in
10⁴/10⁴ cases (5951 taking GSL's second, optimal step). 179 × 3 against the band head, with and without a table:
2040/2040 identical both ways (no project evaluates S(0) itself, and the derivative is a faithful port). ctest 119/119,
`run_tests.sh` 133/133, 1 vs 8 threads identical. `ShftFunc`, `ChannelFunc` and `CoulFunc` no longer include
`gsl_deriv.h` or `gsl_sf_bessel.h` (the opt-in `AZURE_BOUND_SHIFT_NUMERICAL` route keeps GSL in `ShftFunc`).

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
to be decided (section 6; decided on 2026-11-16: retired, section 9). The memo stays in front of `CoulombWavesAt` as
before, so `cache_stats()` is unchanged.

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
options to keep, remap to a `Method`, or retire (`--gsl-coul`: retired, section 9; the table: kept).

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

## 8. GSL's Coulomb functions in the engine: where they are still called (2026-11-16)

`FOR_THE_OWNER.md` section 10 item 6 said that GSL's `gsl_sf_coulomb_wave_FG_e` is still called in
`src/CoulFunc.cpp` and `gsl_sf_hyperg_U` in `src/ECIntegral.cpp`. This section records what reaches those calls,
what the two routines give at the arguments the test projects ask for, and what replacing the calls would change.
Nothing in the engine was changed for it. Engine: `feature/upgrade` at 3ac4a34a, GCC 11.5, Release, GSL 2.6.

### 8.1 The calls and the paths to them

- `gsl_sf_coulomb_wave_FG_e` in `CoulFunc::PlainWaves` (F, F′, G, G′ of a point-Coulomb or neutral pair): reached
  only with the option `use_gsl_coul`.
- `gsl_sf_coulomb_wave_FG_e` in `CoulFunc::computeHybrid` (H⁺ = G + iF at the matching radius of a pair with a
  potential): reached only with the option `use_gsl_coul`.
- `gsl_sf_coulomb_wave_FG_e` in `src/coulomb_table_main.cpp`: always in `azure2-coulomb-table validate`, at 5000
  random points, for one printed line ("references: Michel vs GSL Coulomb P …, S …") that is not part of its
  PASS/FAIL. `build` does not call it.
- `gsl_sf_hyperg_U`: not called by the engine since 2026-11-09 (section 1.4). `src/ECIntegral.cpp` only names it in
  a comment.

So with the default options the engine evaluates no GSL Coulomb function: `PlainWaves` and the hybrid boundary
values call `CoulombWavesAt`, the engine's one call of `physics/coulomb` (`coulsm`). The two sets of item 6 are one
set and one option.

**The option.** `use_gsl_coul` is a `<config>` key, the flag `--gsl-coul` and the pyazr keyword
`Session(use_gsl_coul=)`; default off (`io::ConfigKeys()`, bit `Config::USE_GSL_COULOMB_FUNC`). The bit is handed to
every `CoulFunc` and `ChannelFunc` the engine makes (22 places in `src/` and `api/src/`); the level doors of the
binding take it as the argument `use_gsl_coulomb` (`pyazr.endf` passes the session's option). With the option on,
GSL's routine gives:

- P, S, φ of the channel table and S, S′ at the level energies, **where the Coulomb table does not serve the
  point**: `CoulFunc` asks the table first, whatever the option, so with a table the option changes charged
  channels only outside the table, and without one it changes every open charged channel. What the option does to
  a result therefore depends on whether the machine has a table file;
- the same for neutral pairs and, through the boundary values, for pairs with a hybrid potential (never tabulated);
- F and G of the entrance channel at every radius of the external-capture integrals, and the hard-sphere phase and
  penetrability the capture amplitude uses at the channel radius;
- the Coulomb term of the THM vertex (`ThmCoulombTerm`) and the channel functions the binding returns.

The option does not reach `CoulombTable::ExactOpen` (the table's contents are always `physics/coulomb`'s), the
distorted waves of `physics/transfer`, or anything below threshold (Whittaker functions). Every store that depends
on it is keyed by it: the external-capture file's signature (`WorkingCopy::ECSignature`), the channel tables' pair
signature, `CoulRadialMemo::Key`, `ThmShiftMemo`. GSL's return status is not read at the two calls; an error
reaches the engine's handler (`GSLException::GSLErrorHandler`), which throws.

**Tests that use the option or GSL's routine.** `tests/config_options/check.sh` and
`tests/pyazr/config_options_test.py` (both values of the key), `tests/ec_integrals_signature/check.sh` (a file
computed with the default routine is refused with `--gsl-coul`), `experiment/tests/build_test.cpp` (the key is
read), `coulomb_gsl_crosscheck` (section 5: `coulsm` against GSL beyond the outer turning point). No test project
sets the option.

### 8.2 The two routines at the arguments the test projects ask for

*Method.* A copy of the engine with a probe in `CoulFunc` (not committed): at every evaluation of the default path
it also called GSL at the same (l, η, ρ), timed both calls, and kept the largest differences per band of η. Every
project file under `tests/` (208 files in 73 folders; `tests/convert` left out) was run with `calculate`, Brune
basis, one thread, once with `--no-coulomb-table` and once with the table. A sample of the arguments (every 499th
call, at most 400 per run) was compared with mpmath's `coulombf`/`coulombg` at 40 digits.

*Arguments.* Without the table 1,682,514 evaluations in 55 folders; with it 1,298,285 in 30 folders (capture
radii, neutral and hybrid pairs, points the table does not cover). l ≤ 8.

ΔP is the largest relative difference of P between GSL and `coulsm`, ΔS that of S per max(1, |S|); the times
are per call.

| η | calls (no table) | ρ | ΔP | ΔS | `coulsm`, µs | GSL, µs |
|---|---|---|---|---|---|---|
| 0 | 207 × 10³ | 7.5 × 10⁻⁴ – 5178 | 1.8 × 10⁻¹⁴ | 5.8 × 10⁻¹³ | 0.14 | 0.44 |
| 0 – 1 | 497 × 10³ | 0.081 – 47520 | 5.0 × 10⁻² | 6.6 × 10⁻² | 45 | 1.2 |
| 1 – 10 | 918 × 10³ | 0.041 – 30680 | 4.9 × 10⁻² | 4.9 × 10⁻² | 39 | 0.91 |
| 10 – 30 | 48 × 10³ | 0.023 – 563 | 2.0 × 10⁻² | 3.6 × 10⁻³ | 36 | 0.63 |
| 30 – 100 | 12 × 10³ | 0.022 – 315 | 2.0 × 10⁻² | 3.6 × 10⁻³ | 71 | 0.68 |
| > 100 | 66 | 0.006 – 0.067 | 2.3 × 10⁻² (7 points) | 5.6 × 10⁻³ | 153 | 0.67 |

For η > 100 (`ec_threshold`, capture just above a charged threshold) GSL returned an error at 27 of the 66 points
and F² + G² overflowed at 32. Summed over all runs the default routine took 60.9 s and GSL 1.6 s: for a charged
pair with η ≤ 100 GSL is 36 to 100 times faster per call, for a neutral pair three times slower.

*Which one is right.* Largest errors on 1500 distinct sampled arguments against mpmath (P relative, S per max(1, |S|);
ρ_tp = η + [η² + l(l+1)]^½ is the outer turning point):

| region | points | `coulsm`: P | `coulsm`: S | GSL: P | GSL: S |
|---|---|---|---|---|---|
| η = 0 | 192 | 9.9 × 10⁻¹⁶ | 2.5 × 10⁻¹⁴ | 4.6 × 10⁻¹⁵ | 2.0 × 10⁻¹⁴ |
| η > 0, ρ ≥ ρ_tp | 523 | 7.6 × 10⁻¹⁰ | 7.9 × 10⁻⁹ | 1.7 × 10⁻¹⁵ | 4.0 × 10⁻¹² |
| 0 < η ≤ 10, ρ < ρ_tp | 726 | 6.9 × 10⁻¹⁰ | 4.5 × 10⁻¹⁰ | 5.2 × 10⁻² | 6.2 × 10⁻² |
| 10 < η ≤ 100, ρ < ρ_tp | 58 | 8.7 × 10⁻¹⁰ | 3.0 × 10⁻¹³ | 1.6 × 10⁻² | 3.3 × 10⁻³ |

The differences of the first table are GSL's errors below the Coulomb barrier. Example: l = 3, η = 0.6449,
ρ = 1.2011 (`13N_capture_ay`): mpmath P = 2.35457878 × 10⁻³, `coulsm` the same to 3 × 10⁻¹¹, GSL 2.4771 × 10⁻³
(+5.2 %; its G is 2.5 % low, its own error estimate of G 1 × 10⁻³). Beyond the turning point GSL is the more accurate
of the two (10⁻¹⁵ against 10⁻⁹). Why GSL 2.6 loses accuracy inside the turning point was not examined. This agrees
with `docs/source/developer/coulomb_table.rst` ("not a usable reference") and with pycoulomb's measurements.

*A property of the default routine seen on the way.* For η = 0 and l ≥ 2 at small ρ, `coulsm`'s F alone is
inaccurate relative to itself (it is known: the comment of `coulomb_gsl_crosscheck`): against mpmath, at
ρ = 7.5 × 10⁻⁴ the relative error of F is 4 × 10⁻³ (l = 2), 3 × 10⁵ (l = 3), 3 × 10¹³ (l = 4); in the projects
the largest difference from GSL is 1.5 (η = 0, l ≤ 4, ρ ≥ 0.018: the neutral channels of `7Li_p_ay`, `7Li_p_a`,
`beam_spread_ay` and `data_file_comments`). The error is below 2.5 × 10⁻¹⁶ of
(F² + G²)^½, so P, S and φ do not see it (first row of the first table). F enters alone only in the
external-capture integrand: in `ec_threshold`, a capture project with a neutral pair (l ≤ 3, ρ ≥ 1.3), the
difference of F from GSL is 2.7 × 10⁻¹². A capture through a neutral d wave or higher at ρ ≪ 1 would read that F;
whether any project integrates it was not traced.

### 8.3 What the option changes in the outputs, and what it saves

Engine `work/b-trial95` (the same engine code as 3ac4a34a), `calculate`, Brune basis, one thread, default against
`--gsl-coul`, every output file compared with `tools/classify.py` (887 files of 208 project files):

| | files identical | project files that differ | observable files that differ | above 10⁻³ | above 10⁻² |
|---|---|---|---|---|---|
| with the table | 711 | 38 | 50 | 27 | 20 |
| `--no-coulomb-table` | 458 | 100 | 139 | 53 | 42 |

("above": the largest relative difference of a number in the observable file.)

χ² of the stored parameters, default → option: `13N` 336479.96 → 335217.54 (−0.38 %, with or without the table: its
change is the capture integrals); `T1` 482.2752 → 482.5713; `13N_capture_ay` 3.7395 → 3.7190; `T2` unchanged with
the table (10910.1338428, every channel is served by it) and 10910.1338 → 10897.1227 without. Largest changes of
an observable: with the table 3.2 × 10⁻² (`13N_capture_ay`, R = 2) apart from a cross section near a zero
(`ec_e1m1_phase` `s_1`, 0.67); without it 0.19 (`channel_space`), 0.17 (`beta_delayed_alpha` `16N`). All 23
external-capture files (`intEC.dat`) differ. No run failed with the option.

User-space instructions at one thread (`perf stat -e instructions:u`), default → `--gsl-coul`:

| run | default | `--gsl-coul` | change | χ² at the fit's end |
|---|---|---|---|---|
| `13N` calculation | 74.58 × 10⁹ | 16.75 × 10⁹ | −78 % | |
| `13N` calculation, no table | 76.50 × 10⁹ | 16.60 × 10⁹ | −78 % | |
| `T2` calculation, no table | 7.39 × 10⁹ | 1.23 × 10⁹ | −83 % | |
| `T2` gradient fit | 169.58 × 10⁹ | 169.12 × 10⁹ | −0.3 % | 297.860308551, unchanged |
| `T2` gradient fit, no table | 334.45 × 10⁹ | 305.38 × 10⁹ | −8.7 % | 297.8603 → 298.0371 |
| `13N` gradient fit | 183.85 × 10⁹ | 125.86 × 10⁹ | −32 % | 28100.17 → 28121.46 |

(Each a single run; two runs of the `T2` gradient fit with the same engine and options differed by 0.5 %, so its −0.3
% is no change. Wall times: `13N` calculation 7.6 → 3.2 s, `13N` fit 22.4 → 13.9 s.) What the option saves is the
default routine's time in the capture integrals and, without a table, in the channel functions; with a table and no
capture it saves nothing.

### 8.4 What calling `physics/coulomb` at the two places would change

The two calls would become `CoulombWavesAt`, and `useGSLFunctions_` with the constructor argument of `CoulFunc`
and `ChannelFunc`, the `gsl` members of the memo keys and the 22 hand-overs of the bit would go; `CoulFunc` would no
longer include a GSL header. Less than a day of code. Then:

- **No output of a run with the default options moves**: those runs do not reach GSL (expected: every project in
  three bases and the eight fits identical; no change of the external-capture signature's version).
- **Every run with the option gives the default's numbers**: it moves by the amounts of section 8.3 (percent in
  cross sections, 0.4 % in the χ² of `13N`), towards the values mpmath supports, and slows down by the factors
  there (a capture calculation 4.5 times in instructions).
- The option then means nothing. Removing a `<config>` key, a flag and a pyazr keyword, or keeping them without
  effect, changes what a user's project file or script does: `CLAUDE.md` leaves that to the owner. The places are
  the row of `io::ConfigKeys()`, `OptionStorages()`, `Session.__init__`, the tests of section 8.1, and the pages
  that name the flag (`command_line.rst`, `data_formats.rst`, `output_files.rst`, `pyazr.rst`,
  `coulomb_table.rst`, `README.md`, the skills' `physics-options.md` and `legacy-params.md`).

So the replacement was not made with this section: it is the retirement of a user-visible option, to be decided
(section 6 listed it). The PI chose the second of the choices below on 2026-11-16; section 9 has the design and
what was built. The choices, with what was measured for each:

1. *Keep the option as it is.* It is documented as fast and inaccurate; its numbers are wrong by up to 5 % in P
   below the barrier, and its meaning depends on the table file.
2. *Retire it* (the key refused with a message, or accepted without effect): the places listed above.
3. *Give it another routine of `physics/coulomb`.* Measured on 2387 of the sampled arguments (0 < η ≤ 100):
   `coulfg` (Barnett) takes 1.5 µs per call but has GSL's weakness (up to 1.1 × 10⁻¹ in P inside the turning
   point, 450 points against mpmath) and fails at 6 points (η ≥ 59); `coulcc` takes 15 µs against `coulsm`'s 48,
   with errors of 6.5 × 10⁻¹² in P and 7.2 × 10⁻¹² in S, but fails (NaN, status `failed`) at 11 points, all with
   η ≥ 48. None is both fast and safe as it stands; a faster default is a separate question, and it moves every
   output at about 10⁻⁹.

The cross-check line of `azure2-coulomb-table validate` compares over the whole table domain, where GSL is wrong
below the barrier, so it prints percent-size numbers that say nothing about the table; restricting it to
ρ ≥ ρ_tp (as `coulomb_gsl_crosscheck` does) would make it a check. Not changed here: it is the tool's output.

### 8.5 What was changed with this section

Two places of the engine, an include and a comment, no statement: `include/CoulFunc.h` no longer includes
`gsl/gsl_sf_coulomb.h` (the header declared nothing from it; `src/CoulFunc.cpp` includes it for the two calls), so
GSL's Coulomb functions are named in `src/CoulFunc.cpp` and `src/coulomb_table_main.cpp` only; and the comment of
`src/ECIntegral.cpp` that named `gsl_sf_hyperg_U` as the present cost says what it is now. The option, the two
calls and the tool's line are as they were.

## 9. The option `use_gsl_coul` retired; the tool's cross-check; F alone at η = 0 (design, 2026-11-16)

The PI decided on the three points of section 8 (2026-11-16, relayed by the director): the option is retired, the
cross-check line of `azure2-coulomb-table validate` is restricted to where GSL is right, and the F of `coulsm` at
η = 0 is measured, not changed.

### 9.1 Retiring `use_gsl_coul` / `--gsl-coul`

*What it is physically.* One routine for F_l, G_l is left in the engine, `physics/coulomb`'s `coulsm`. GSL's
routine is wrong by up to 5 % in P inside the outer turning point (section 8.2), which is where a charged channel
below its barrier lives; an option that selects it selects wrong numbers.

*The rule, for the three spellings.* The pattern is that of `useAdaptiveGrid` (`project_format3.md` section 3.1,
`E_CONFIG_KEY_REMOVED`): a removed key is refused with either value and the message says to delete it. A value of
0 is not accepted quietly: the line names an option that does not exist.

- One table in `io/`, `RemovedConfigKeys()` (`ConfigSection.h`): name, the flag it had, and why it is gone. Rows:
  `useAdaptiveGrid` (its text unchanged) and `use_gsl_coul` (`--gsl-coul`). One function gives the message in the
  spelling of the front end, `RemovedConfigKeyText(row, OptionSpelling)`: "use_gsl_coul is no longer a key of
  `<config>`: delete the line. …", "--gsl-coul is no longer an option of AZURE2: leave it out. …", "use_gsl_coul is
  no longer a keyword of pyazr.Session: leave it out. …", each followed by the same reason: GSL's Coulomb functions
  are wrong by up to 5 % in the penetrability inside the outer turning point; AZURE2's own functions
  (`physics/coulomb`) are the only ones, and they are what the default always used.
- *A format-3 file with the key* (either value): `E_CONFIG_KEY_REMOVED` at its line, from the reader
  (`experiment/src/ProjectReader.cpp`, which asks the table instead of comparing with one name), hence from
  `--validate` (exit 2), the run (exit 1, nothing computed), `pyazr.Project.validate` and a session. Before this
  change the key was a row of `io::ConfigKeys()`; a key that is no row gives `E_UNREAD_LINE` ("not a key of
  `<config>`"), which would not say why, so the removed keys are asked first.
- *The flag*: `parseOptions` throws "Unknown option --gsl-coul." followed by the message, as it does for
  `--no-readline` (exit 2, before the project is read). `--help` does not list it, as it does not list `--no-readline`;
  `command_line.rst` keeps its entry, as removed.
- *The keyword*: `pyazr.Session(..., use_gsl_coul=...)` raises `TypeError` with the message (a Python call with a
  keyword that does not exist is a `TypeError`); the signature loses the keyword. The binding gives the text
  (`_azure2.removed_option_text`), so Python holds no copy. The level doors of the binding lose their argument
  `use_gsl_coulomb` (private: `pyazr.endf` is their only caller).
- *The converter.* Old files (formats 1 and 2) have no such key: the GSL choice of the old setup program was a bit
  of `<lastRun>`, which the converter removes with the rest of that section (the PI, 2026-11-13, `qt_parity.md`):
  the converted project is computed with the program's own functions, as any project without the key. Nothing
  changes in `convert/`. This differs from what the old program did for a user who had ticked the box, by the
  amounts of section 8.3, and is row 19 of `FOR_THE_OWNER.md` section 3.

*What goes.* The two `gsl_sf_coulomb_wave_FG_e` calls and `useGSLFunctions_` in `CoulFunc`; the constructor
argument of `CoulFunc` and `ChannelFunc`; `Config::USE_GSL_COULOMB_FUNC` (bit 12 stays unused) and its 22
hand-overs; the `gsl` member of `CoulRadialMemo::Key`, `ShiftDerivativeMemo::Key` and `ThmShiftMemo`; the argument
of `EngineChannelFunctions`, `PhotonAmplitudes`, `FinalSqrtNF`; the argument `useGSL` of `ThmCoulombTerm`, which
had no effect; the rows of `io::ConfigKeys()` and `OptionStorages()` (no implication or contradiction named the
option), so `--describe`, `pyazr.kinds()["options"]` and the web page's form lose it with the row.

*Signatures (no default result may change).* `WorkingCopy::ECSignature` hashes the option's value as an integer
between the version text and the long-wavelength flag. A default run hashed 0 there; the line becomes the
constant 0 with a comment, so every signature of a default run is the one it was and the version stays v7. The
in-memory keys (the channel tables' pair signature, the external-capture table's term key) are compared only
within one run and lose the entry.

*What still links GSL.* The engine keeps GSL for the quadratures of the capture integrals and the folds
(`gsl_integration_*`), the fitter `FitGslTrustRegion`, Bessel and Gamma functions in `ThmFunc`,
`IntegratedFermiFunc` and others, and the opt-in numerical bound-shift route; `azure2-coulomb-table validate`
keeps its one call of GSL's Coulomb functions (section 9.2). The dependency stays.

*Tests.* A test of the option's two values becomes a test of the refusal:

- `experiment/tests/build_test.cpp`: `use_gsl_coul=0` and `=1` each give `E_CONFIG_KEY_REMOVED` alone, at the
  line, with "delete the line" and the reason; the every-option text loses the key.
- `tests/config_options/check.sh`: the two rows of the key-equals-flag table become three refusals: the key
  (`--validate` exit 2, the run exit 1 without output), and the flag (exit 2, the message, no output).
- `tests/pyazr/config_options_test.py`: the keyword with either value raises `TypeError` with the reason; a
  project with the key is refused by a session; the option is not in `pyazr.kinds()["options"]`.
- `tests/ec_integrals_signature` case 6 (the file of the other routine is refused) becomes: the run with
  `--gsl-coul` is refused, and the saved file is still taken by the next run without a warning.
- `tests/reference/thm_coulomb_term_test.cpp` ran every case with `useGSL` off and on to show it did not matter;
  with the argument gone each case runs once.
- The reason itself: `coulomb_gsl_crosscheck` gets stored mpmath values of P at points inside the turning point
  (among them l = 3, η = 0.6449, ρ = 1.2011): `coulsm` must agree to 10⁻⁹ and GSL must be off by more than 10⁻³
  (it is 5.2 × 10⁻² there). If a later GSL repairs its routine the test says so; the retirement would still hold.

*Gate.* 0 warnings; every project in three bases identical to `tab-trial94-1`, except the two `--describe`
stories under `tests/*/expected/` that list the options (one line less); the eight fits identical to
`work/b-trial95`; `ctest -LE slow`, `tests/run_tests.sh`; the affected tests also without a table file.

### 9.2 The cross-check line of `azure2-coulomb-table validate`

The 5000 random points are kept (same seed, same draws); a point enters the comparison only when
ρ ≥ ρ_tp = η + [η² + l(l+1)]^½, where GSL is good to 10⁻¹⁵ (section 8.2). The line says so and counts: "references
beyond the outer turning point (N of 5000 points; M more refused or scaled by GSL): Michel vs GSL Coulomb P …,
S …". The numbers become a check of the table's reference routine at the level of its own accuracy (10⁻⁹); they
stay printed and outside PASS/FAIL, whose limit is the interpolation error. No test compares the line.

### 9.3 F alone at η = 0 (measurement only)

Where the engine reads F by itself: the external-capture integrand (`ECIntegral`: `fgEvaluate` → FW, the rate's
fixed grid, the channel-capture integrand) and the Coulomb functions the binding returns to a caller
(`pyazr.physics.coulomb_functions`). Everything else forms P, S, φ or H⁺ = G + iF, where the error of F is below
2.5 × 10⁻¹⁶ of the amplitude. Measured with a scratch build over every test project: (i) every F of a neutral pair
read in an external-capture integrand, against GSL's at η = 0 (good to 2 × 10⁻¹⁵ against mpmath, section 8.2):
how many are off by more than 10⁻⁸ relative, in which project; (ii) the same build with those F replaced by GSL's,
every output file compared with the unreplaced run: the largest relative change of an integral (`intEC.dat`) and of
an observable. An observable moved by more than 10⁻¹⁰ stops the item and is reported with its project and size.

### 9.4 The retirement as built (2026-11-16)

Built as section 9.1 says, with these points settled on the way:

- The flag's message starts with the parser's words for a flag it does not take, "Unknown option --gsl-coul.",
  then the sentence of the table. `tests/cli_help` tells an accepted flag from a refused one by those words; the
  first build, without them, counted the flag as accepted and not printed.
- `Project.set_config(use_gsl_coul=...)` is refused by `io::ConfigEdits` with "'use_gsl_coul' is not an entry of
  `<config>` any more: it is not written." and the reason (it said "not an entry" with the list of entries before,
  as for a misspelt name; `useAdaptiveGrid` gets the same sentence now).
- A session given the name through the binding's `RuntimeOptions.options` gets the keyword's message
  (`NotAnOptionText`, `src/ProjectOptions.cpp`).
- The count of options is 10 where a test or a page said 11: `io_config_section`, `upgrade_project_options`,
  `tests/pyazr/config_options_test.py`, `describe_test.py`, `tests/config_options/check.sh` (the options block of
  `--describe`), `data_formats.rst`, `project_format3.md`. The two `--describe` stories under `tests/*/expected/`
  (`units_scales`, `ratio_to_rutherford`) lose the line `use_gsl_coul = 0 (default)` and nothing else.
- The description of `E_CONFIG_KEY_REMOVED` names both keys, so `tests/schema/schema_v3.json` and `format3.rst`
  were generated again (one line each).

**Verified** (the finished tree, with sections 9.5 and 9.6; `work/director/verify_coul_slice3.out`, build
`work/b-coul3`, GCC 11.5 Release, pyazr on). 0 warnings. Every project in three bases against `tab-trial94-1`: 887,
851 and 887 output files identical, none differs, none new (the `intEC.dat.sig` files among them: signature v7
unchanged). The eight fits identical to `work/b-trial95` file by file. `ctest -LE slow` 325 of 325;
`tests/run_tests.sh` 172 of 172. Without a Coulomb table file (build `work/b-coul2`, the same code): 42 of 42 tests
matching `pyazr_evaluation|coulomb|config_options|gsl|whittaker|ec_|project_options|io_config|experiment_|
pyazr_describe`, and the scripts `config_options` and `ec_integrals_signature`. Instructions at one thread against
`work/b-trial95`, no output differing: `T2` calculation 1,242,561,787 → 1,242,372,132; `13N` calculation
74,582,449,469 → 74,563,131,713; `13N` gradient fit 183,843,244,523 → 183,817,809,638; `T2` gradient fit
168,750,567,860 → 169,781,017,131 (+0.61 %) and, run twice more in another directory, 168,731,722,915 →
169,147,757,454 (+0.25 %, the same both times). The first full check of this step (`verify_coul_slice2.out`) had
two failures, both repaired before this one: `check_skills` (a date in two skill lines outside the skills' form
of a version note) and `tests/cli_help` (the flag's message without "Unknown option").

`src/CoulFunc.cpp` includes no GSL header any more. What calls GSL's Coulomb functions now: `azure2-coulomb-table
validate` (section 9.5) and the test `coulomb_gsl_crosscheck`, which prints "inside the outer turning point, 7
points against mpmath: coulsm within 1.52e-10 (P and S); GSL's P off by 4.93e-03 to 5.20e-02" (GSL 2.6).
